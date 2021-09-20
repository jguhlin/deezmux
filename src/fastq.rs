use bytelines::*;
use simdutf8::basic::from_utf8;
use triple_accel::*;

use std::collections::HashMap;
use std::io::{BufReader, Read};

pub struct FastqEntry<'fq> {
    pub id: &'fq [u8],
    pub scores: &'fq [u8],
    pub sequence: &'fq [u8],
}

pub struct FastqSplitter {
    files: Vec<String>,
    mm1: u32,
    mm2: u32,
    barcodes: Vec<(String, String, String)>,
    basename: String,
}

impl FastqSplitter {
    // Builder style

    pub fn new() -> FastqSplitter {
        FastqSplitter {
            files: Vec::new(),
            mm1: 2,
            mm2: 2,
            barcodes: Vec::new(),
            basename: "output".to_string(),
        }
    }

    pub fn with_files(mut self, files: Vec<String>) -> FastqSplitter {
        self.files = files;
        self
    }

    pub fn with_mm(mut self, mm1: u32, mm2: u32) -> FastqSplitter {
        self.mm1 = mm1;
        self.mm2 = mm2;
        self
    }

    pub fn with_barcodes(mut self, barcodes: Vec<(String, String, String)>) -> FastqSplitter {
        self.barcodes = barcodes;
        self
    }

    pub fn with_basename(mut self, basename: String) -> FastqSplitter {
        self.basename = basename;
        self
    }

    pub fn match_barcodes<R: Read>(&self, reader: R) {
        // -> Vec<(String, String, String)>{
        let mut lines = BufReader::with_capacity(2 * 1024 * 1024, reader).byte_lines();
        // let mut lines = BufReader::new(reader).byte_lines();
        let mut counts = HashMap::new();

        loop {
            let header;
            let id;

            match lines.next() {
                Some(Ok(line)) => {
                    header = from_utf8(line)
                        .expect("FASTQ Header line is not valid UTF-8")
                        .clone();
                    let n = header.len();

                    id = header[n - 17..n].to_string().clone();
                    let count = match counts.get_mut(&id) {
                        Some(count) => count,
                        None => {
                            counts.insert(id.clone(), 0);
                            counts.get_mut(&id).unwrap()
                        }
                    };

                    *count += 1;
                }
                _ => break,
            };

            let sequence = match lines.next() {
                Some(Ok(line)) => line.clone(),
                Some(Err(e)) => panic!("Invalid FASTQ: {}", e),
                _ => panic!("Out of data!")
            };

            match lines.next() {
                Some(Ok(_)) => (),
                _ => panic!("Invalid FASTQ"),
            };

            let scores = match lines.next() {
                Some(Ok(line)) => line.clone(),
                _ => panic!("Invalid FASTQ"),
            };
        }

        let mut hash_vec: Vec<(&String, &usize)> = counts.iter().collect();
        hash_vec.sort_by(|a, b| b.1.cmp(a.1));
        hash_vec.reverse();

        let mut assigned_reads: usize = 0;
        let mut unassigned_reads: usize = 0;
        let mut ambiguous_reads: usize = 0;

        let mut id_counts = HashMap::new();

        for (barcode, count) in hash_vec {
            let read_barcodes: Vec<&str> = barcode.split("+").collect();

            let mut scores = HashMap::new();

            for (bc0, bc1, id) in self.barcodes.iter() {
                let dist1 = levenshtein(read_barcodes[0].as_bytes(), bc0.as_bytes());
                let dist2 = levenshtein(read_barcodes[1].as_bytes(), bc1.as_bytes());

                scores.insert(id.clone(), dist1 + dist2);
            }

            let mut scores_vec: Vec<(&String, &u32)> = scores.iter().collect();
            scores_vec.sort_by(|a, b| b.1.cmp(a.1));
            scores_vec.reverse();

            let min = scores_vec[0].1;
            let min_id = scores_vec[0].0;

            // for i in scores_vec[1..].iter() {
            if min == scores_vec[1].1 && *min <= 4 {
                println!("Current Min: {} {} {}", min, min_id, barcode);
                println!("Ambiguous found: {} Affecting: {}", scores_vec[1].0, count);
                ambiguous_reads += count;
            } else if *min <= 4 {
                assigned_reads += count;

                let e = id_counts.entry(min_id.clone()).or_insert(0);
                *e += count;

            } else {
                unassigned_reads += count;
            }
            // }
        }

        println!("----Finished");
        println!();
        println!(
            "Assigned: {} Unassigned: {} Ambiguous Reads: {}",
            assigned_reads, unassigned_reads, ambiguous_reads
        );

        println!("{:#?}", id_counts);

        println!("----And done....");
        println!();

        //println!("{:#?}", hash_vec);

        //println!("{:#?}", counts.get("AAGCACTG+CGATGTTC"));
        //println!("{:#?}", counts.get("AACTGAGC+TCTTACGG"));
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn new_fastq() {
        let not_really_a_fastq = String::from(
            "@MT_E00516:746:HG3WYCCX2:6:1101:2229:1661 1:N:0:NGAGCTAG+NAGCCTGA\nNTG\n+\n#AA\n",
        );
        split_fastq_by_id(not_really_a_fastq.as_bytes(), "test");
    }
}
