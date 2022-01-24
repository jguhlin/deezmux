use std::ffi::OsStr;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::sync::Arc;

use clap::{App, AppSettings, Parser, Subcommand};
use crossbeam::channel::{bounded, Receiver, Sender};
use crossbeam::thread;
use flate2::read::MultiGzDecoder;
use indicatif::ProgressStyle;
use indicatif::{MultiProgress, ProgressBar};
use wax::{Glob, Pattern};

mod fastq;
use fastq::*;

#[derive(Parser)]
#[clap(name = "deezmux")]
#[clap(author = "Joseph Guhlin <joseph.guhlin@gmail.com>")]
#[clap(version = "0.2.0")]
#[clap(about = "Fast demultiplexing of Illumina files where the barcode IDs are in the header", long_about = None)]
struct Cli {
    barcode_file: String,
    output_directory: String,
    read_prefix: String,
}

enum Mode {
    BarcodesInHeader,
    BarcodesInSeparateFile,
}

fn parse_barcode_file(barcode_file: &str) -> Vec<(String, String, String, String, String)> {
    let file = File::open(barcode_file).expect("Unable to open barcode file");
    let bufread = BufReader::new(file);

    let mut barcodes: Vec<(String, String, String, String, String)> = Vec::new();

    for line in bufread.lines().skip(1) {
        let line = line.unwrap();
        let j: Vec<&str> = line.split(",").collect();
        let z: Vec<&str> = j[1].split("+").collect();
        barcodes.push((
            z[0].to_string(), // Barcode 0
            z[1].to_string(), // Barcode 1
            j[0].to_string(), // Sample ID
            j[2].to_string(), // Barcode 0 sequence
            j[3].to_string(), // Barcode 1 sequence
        ));
    }

    barcodes
}

fn main() {
    let args = Cli::parse();

    /* let barcode_file = match matches.value_of("barcode_file") {
        Some(x) => x,
        None => panic!("No barcode file specified"),
    };

    let files = match matches.values_of("read_files") {
        Some(x) => x,
        None => panic!("No read files specified"),
    }; */
    let barcode_file = &args.barcode_file;

    println!("Parsing barcode file: {}", barcode_file);
    let barcodes = Arc::new(parse_barcode_file(barcode_file));

    let output_directory = &args.output_directory;
    let prefix = &args.read_prefix;

    // Figure out what type of files we are dealing with

    let prefix_directory = Path::new(prefix)
        .parent()
        .expect("Unable to determine parent directory of prefix path")
        .canonicalize()
        .expect("Unable to canonicalize prefix path");
    let prefix_files = Path::new(prefix).file_name().unwrap().to_str().unwrap();
    let glob = format!("{}*.gz", prefix_files).to_string();

    let glob = Glob::new(&glob).expect("Unable to create glob");
    let mut files = Vec::with_capacity(4);

    for i in glob.walk(prefix_directory, 1) {
        match i {
            Ok(x) => {
                files.push(x.into_path());
            }
            Err(e) => {
                println!("{}", e);
            }
        }
    }

    let mode = match files.len() {
        0 => panic!("No files found matching prefix provided"),
        1 => {
            panic!("Only one file found matching prefix. Expected paired reads")
        }
        2 => {
            println!("Found two files matching prefix. Assuming paired reads and barcode found in header");
            Mode::BarcodesInHeader
        }
        4 => {
            println!(
                "Found four files matching prefix. Assuming paired reads and barcode in I files"
            );
            Mode::BarcodesInSeparateFile
        }
        _ => {
            panic!("Found unusual number of files matching prefix. Expected paired reads")
        }
    };

    let index_files = match mode {
        Mode::BarcodesInHeader => None,
        Mode::BarcodesInSeparateFile => Some((
            files
                .iter()
                .find(|&x| x.file_name().unwrap().to_str().unwrap().contains("_I1"))
                .unwrap(),
            files
                .iter()
                .find(|&x| x.file_name().unwrap().to_str().unwrap().contains("_I2"))
                .unwrap(),
        )),
    };

    let files = [
        files
            .iter()
            .find(|&x| x.file_name().unwrap().to_str().unwrap().contains("_R1"))
            .unwrap(),
        files
            .iter()
            .find(|&x| x.file_name().unwrap().to_str().unwrap().contains("_R2"))
            .unwrap(),
    ];

    // let fqs = FastqSplitter::new().with_mm(2, 2);

    thread::scope(|s| {
        let m = MultiProgress::new();
        let mut t = Vec::new();

        for (i, file) in files.iter().enumerate() {
            let r = format!("r{}", i + 1);
            // println!("Processing {} as {}", file, r);

            /*
            let file_pb = File::open(file).expect("Unable to open file");
            let pb = ProgressBar::new(file_pb.metadata().unwrap().len());
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {eta}"),
            );
            let barcodes =
                fqs.match_barcodes(MultiGzDecoder::new(BufReader::new(pb.wrap_read(&file_pb))));

            */

            let file_pb = File::open(file).expect("Unable to open file");
            let pb = ProgressBar::new(file_pb.metadata().unwrap().len());
            pb.set_style(
                ProgressStyle::default_bar()
                    .template("{spinner:.blue}▕{bar:.green}▏{bytes:>4}/{total_bytes:4} {eta}")
                    // .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {eta}"),
                    .progress_chars("█▇▆▅▄▃▂▁  "),
            );

            let file_fh = pb.wrap_read(file_pb);

            let _pb = m.add(pb);

            let barcodes = Arc::clone(&barcodes);

            let handle = s.spawn(move |_| {
                split_by_barcodes(
                    MultiGzDecoder::new(BufReader::new(file_fh)),
                    r,
                    "output".to_string(),
                    barcodes,
                    index_files,
                );
            });
            t.push(handle);
        }

        m.join().unwrap();

        for i in t {
            i.join().unwrap();
        }
    })
    .expect("Unable to properly scope");
}
