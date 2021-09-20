use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use clap::{load_yaml, App};
use flate2::read::MultiGzDecoder;


mod fastq;
use fastq::FastqSplitter;

fn parse_barcode_file(barcode_file: &str) -> Vec<(String, String, String)> {
    let file = File::open(barcode_file).expect("Unable to open barcode file");
    let bufread = BufReader::new(file);

    let mut barcodes: Vec<(String, String, String)> = Vec::new();

    for line in bufread.lines().skip(1) {
        let line = line.unwrap();
        let j: Vec<&str> = line.split(",").collect();
        let z: Vec<&str> = j[1].split("+").collect();
        barcodes.push((z[0].to_string(), z[1].to_string(), j[0].to_string()));
    }

    barcodes
}

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from(yaml).get_matches();

    let barcode_file = match matches.value_of("barcode_file") {
        Some(x) => x,
        None => panic!("No barcode file specified"),
    };

    let files = match matches.values_of("read_files") {
        Some(x) => x,
        None => panic!("No read files specified"),
    };

    let barcodes = parse_barcode_file(barcode_file);

    let mut fqs = FastqSplitter::new().with_barcodes(barcodes).with_mm(2, 2);

    for file in files {    
        println!("Processing {}", file);
        let file = File::open(file).expect("Unable to open file");
        let pb = ProgressBar::new(file.metadata().unwrap().len());
        pb.set_style(ProgressStyle::default_bar().template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {eta}"));
        fqs.match_barcodes(MultiGzDecoder::new(pb.wrap_read(file)));
    }
}
