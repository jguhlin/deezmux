use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use indicatif::ProgressBar;

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
    let file = File::open("adapter_trimmed.collapsed.truncated").expect("Unable to open file");
    let pb = ProgressBar::new(file.metadata().unwrap().len());

    let barcodes = parse_barcode_file("j197_adapters.txt");

    let mut fqs = FastqSplitter::new().with_barcodes(barcodes).with_mm(2, 2);

    fqs.match_barcodes(pb.wrap_read(file));
}
