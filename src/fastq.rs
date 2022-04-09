use bytelines::*;
use crossbeam::channel::{bounded, Receiver, Sender};
use crossbeam::thread;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use hashbrown::HashMap;
use simdutf8::basic::from_utf8;
use triple_accel::*;
use twox_hash::xxh3::RandomHashBuilder64;

// use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::sync::Arc;
// use std::thread;

/* pub struct FastqEntry<'fq> {
    pub id: &'fq [u8],
    pub scores: &'fq [u8],
    pub sequence: &'fq [u8],
} */

pub struct FastqSplitter {
    // files: Vec<String>,
    mm1: u32,
    mm2: u32,
    barcodes: Vec<(String, String, String, String, String)>,
    // basename: String,
}

impl FastqSplitter {
    // Builder style

    pub fn new() -> FastqSplitter {
        FastqSplitter {
            // files: Vec::new(),
            mm1: 2,
            mm2: 2,
            barcodes: Vec::new(),
            // basename: "output".to_string(),
        }
    }

    /*
    pub fn with_files(mut self, files: Vec<String>) -> FastqSplitter {
        self.files = files;
      self
    } */

    pub fn with_mm(mut self, mm1: u32, mm2: u32) -> FastqSplitter {
        self.mm1 = mm1;
        self.mm2 = mm2;
        self
    }

    pub fn with_barcodes(
        mut self,
        barcodes: Vec<(String, String, String, String, String)>,
    ) -> FastqSplitter {
        self.barcodes = barcodes;
        self
    }

    /*     pub fn with_basename(mut self, basename: String) -> FastqSplitter {
        self.basename = basename;
        self twox_hash::xxh3::RandomHashBuilder64
    } */

    pub fn match_barcodes<R: Read + Send + Sync>(
        &self,
        reader: R,
    ) -> HashMap<String, String, RandomHashBuilder64> {
        let (sender, receiver) = bounded(100000);

        let mut assigned_barcodes: HashMap<String, String, RandomHashBuilder64> =
            Default::default();
        let mut counts: HashMap<String, usize, RandomHashBuilder64> = Default::default();

        thread::scope(|s| {
            let mut lines = BufReader::new(reader).byte_lines();

            s.spawn(move |_| {
                let mut header;
                let mut id;

                loop {
                    match lines.next() {
                        Some(Ok(line)) => {
                            header = from_utf8(line).expect("FASTQ Header line is not valid UTF-8");
                            let n = header.len();
                            id = header[n - 17..n].to_string();
                            sender.send(Some(id)).expect("Error sending");
                        }
                        _ => {
                            sender.send(None).expect("Error sending complete");
                            break;
                        }
                    };

                    lines
                        .next()
                        .expect("Invalid FASTQ File")
                        .expect("Invalid FASTQ File");
                    lines
                        .next()
                        .expect("Invalid FASTQ File")
                        .expect("Invalid FASTQ File");
                    lines
                        .next()
                        .expect("Invalid FASTQ File")
                        .expect("Invalid FASTQ File");
                }
            });

            while let Ok(Some(id)) = receiver.recv() {
                let count = match counts.get_mut(&id) {
                    Some(count) => count,
                    None => {
                        counts.insert(id.clone(), 0);
                        counts.get_mut(&id).unwrap()
                    }
                };
                *count += 1;
            }
        })
        .unwrap();

        let mut hash_vec: Vec<(&String, &usize)> = counts.iter().collect();
        hash_vec.sort_by(|a, b| b.1.cmp(a.1));
        hash_vec.reverse();

        let mut assigned_reads: usize = 0;
        let mut unassigned_reads: usize = 0;
        let mut ambiguous_reads: usize = 0;

        let mut id_counts = HashMap::new();

        for (barcode, count) in hash_vec {
            let read_barcodes: Vec<&str> = barcode.split("+").collect();

            let mut scores: HashMap<String, u32, RandomHashBuilder64> = Default::default();

            for (bc0, bc1, id, _r1_file, _r2_file) in self.barcodes.iter() {
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
                // println!("Current Min: {} {} {}", min, min_id, barcode);
                // println!("Ambiguous found: {} Affecting: {}", scores_vec[1].0, count);
                ambiguous_reads += count;
                assigned_barcodes.insert(barcode.clone(), "AMBIGUOUS".to_string());
            } else if *min <= 4 {
                assigned_reads += count;

                let e = id_counts.entry(min_id.clone()).or_insert(0);
                *e += count;

                assigned_barcodes.insert(barcode.clone(), min_id.clone());
            } else {
                unassigned_reads += count;
                assigned_barcodes.insert(barcode.clone(), "UNASSIGNED".to_string());
            }
            // }
        }

        /*

        println!("----Finished");
        println!();
        println!(
            "Assigned: {} Unassigned: {} Ambiguous Reads: {}",
            assigned_reads, unassigned_reads, ambiguous_reads
        );

        println!("{:#?}", id_counts);

        println!("----And done....");
        println!();

        println!("{:#?}", assigned_barcodes); */

        //println!("{:#?}", hash_vec);

        //println!("{:#?}", counts.get("AAGCACTG+CGATGTTC"));
        //println!("{:#?}", counts.get("AACTGAGC+TCTTACGG"));

        /* Split files */

        assigned_barcodes
    }
}

pub fn split_by_barcodes<R: Read + Send + Sync>(
    reader: R,
    suffix: String,
    output_directory: String,
    barcodes: Arc<Vec<(String, String, String, String, String)>>,
    index_files: Option<(&PathBuf, &PathBuf)>,
) {
    let (sender, receiver) = bounded(8192);

    fs::create_dir_all(&output_directory).expect("Unable to create directory");

    let mut assigned_barcodes: HashMap<String, String, RandomHashBuilder64> = Default::default();

    thread::scope(|s| {
        let mut index_receiver1 = None;
        let mut index_receiver2 = None;

        if let Some((idx1, idx2)) = index_files {
            let (s1, ir1) = bounded(1);
            let (s2, ir2) = bounded(1);

            index_receiver1 = Some(ir1);
            index_receiver2 = Some(ir2);

            // Read the barcodes...
            s.spawn(move |_| {
                let fh1 = File::open(idx1).expect("Unable to open file");
                let fh1 = MultiGzDecoder::new(BufReader::new(fh1));
                let fh2 = File::open(idx2).expect("Unable to open file");
                let fh2 = MultiGzDecoder::new(BufReader::new(fh2));

                let mut lines1 = BufReader::new(fh1).byte_lines();
                let mut lines2 = BufReader::new(fh2).byte_lines();

                loop {
                    let header;

                    match lines1.next() {
                        Some(Ok(line)) => {
                            header = from_utf8(line)
                                .expect("FASTQ Header line is not valid UTF-8")
                                .clone();

                            let entry = (
                                header.to_string().clone(),
                                from_utf8(lines1.next().unwrap().unwrap())
                                    .unwrap()
                                    .to_string()
                                    .clone(),
                            );

                            s1.send(Some(entry)).expect("Error sending Fastq Barcode");
                            lines1.next();
                            lines1.next();
                        }
                        _ => {
                            s1.send(None).expect("Unable to send Done command");
                            break;
                        }
                    };

                    let header;

                    match lines2.next() {
                        Some(Ok(line)) => {
                            header = from_utf8(line)
                                .expect("FASTQ Header line is not valid UTF-8")
                                .clone();

                            let entry = (
                                header.to_string().clone(),
                                from_utf8(lines2.next().unwrap().unwrap())
                                    .unwrap()
                                    .to_string()
                                    .clone(),
                            );

                            s2.send(Some(entry)).expect("Error sending Fastq Barcode");
                            lines2.next();
                            lines2.next();
                        }
                        _ => {
                            s2.send(None).expect("Unable to send Done command");
                            break;
                        }
                    };
                }
            });
        }

        s.spawn(move |_| {
            let mut lines = BufReader::new(reader).byte_lines();

            loop {
                let header;
                let id;

                match lines.next() {
                    Some(Ok(line)) => {
                        header = from_utf8(line)
                            .expect("FASTQ Header line is not valid UTF-8")
                            .clone();

                        if index_receiver1.is_some() {
                            let i1 = index_receiver1.as_ref().unwrap().recv().expect("Error with barcode 1").expect("Error with barcode 1");
                            let i2 = index_receiver2.as_ref().unwrap().recv().expect("Error with barcode 2").expect("Error with barcode 2");

                            assert!(header.split_whitespace().next().unwrap() == i1.0.split_whitespace().next().unwrap());
                            assert!(header.split_whitespace().next().unwrap() == i2.0.split_whitespace().next().unwrap());

                            id = format!("{}+{}", i1.1, i2.1);

                        } else {
                            let n = header.len();
                            id = header[n - 17..n].to_string().clone();
                        }

                        let entry = (
                            id,
                            header.to_string().clone(),
                            from_utf8(lines.next().unwrap().unwrap())
                                .unwrap()
                                .to_string()
                                .clone(),
                            from_utf8(lines.next().unwrap().unwrap())
                                .unwrap()
                                .to_string()
                                .clone(),
                            from_utf8(lines.next().unwrap().unwrap())
                                .unwrap()
                                .to_string()
                                .clone(),
                        );

                        sender.send(Some(entry)).expect("Error sending Fastq entry");
                    }
                    _ => {
                        sender.send(None).expect("Unable to send Done command");
                        break;
                    }
                };
            }
        });

        // END OF THREAD

        let mut files: HashMap<String, Sender<Option<[String; 4]>>> = HashMap::new();

        let mut ids = (*barcodes)
            .iter()
            .map(|x| x.2.clone())
            .collect::<Vec<String>>();
        ids.push("AMBIGUOUS".to_string());
        ids.push("UNASSIGNED".to_string());

        for id in ids.into_iter() {
            let (send, r) = bounded(2048);
            files.insert(id.clone(), send);

            let output_directory = output_directory.clone();
            let suffix = suffix.clone();

            s.spawn(move |_| {
                let mut file = BufWriter::new(GzEncoder::new(
                    File::create(format!("{}/{}_{}.fq.gz", output_directory, id, suffix)).unwrap(),
                    Compression::fast(),
                ));

                while let Ok(Some(e1)) = r.recv() {
                    for e in e1.iter() {
                        file.write_all(e.as_bytes()).unwrap();
                        writeln!(file).unwrap();
                    }
                }
            });
        }

        while let Ok(Some((id, header, e1, e2, e3))) = receiver.recv() {
            let x = match assigned_barcodes.get(&id) {
                Some(x) => x,
                None => {
                    let read_barcodes: Vec<&str> = id.split("+").collect();
                    let mut scores: HashMap<String, u32, RandomHashBuilder64> = Default::default();
                    for (bc0, bc1, id, _r1_file, _r2_file) in barcodes.iter() {
                        let dist1 = levenshtein(read_barcodes[0].as_bytes(), bc0.as_bytes());
                        let dist2 = levenshtein(read_barcodes[1].as_bytes(), bc1.as_bytes());

                        scores.insert(id.clone(), dist1 + dist2);
                    }
                    let mut scores_vec: Vec<(&String, &u32)> = scores.iter().collect();
                    scores_vec.sort_by(|a, b| b.1.cmp(a.1));
                    scores_vec.reverse();

                    let min = scores_vec[0].1;
                    let min_id = scores_vec[0].0;

                    if min == scores_vec[1].1 && *min <= 4 {
                        assigned_barcodes.insert(id.clone(), "AMBIGUOUS".to_string());
                    } else if *min <= 4 {
                        assigned_barcodes.insert(id.clone(), min_id.clone());
                    } else {
                        assigned_barcodes.insert(id.clone(), "UNASSIGNED".to_string());
                    }

                    assigned_barcodes.get(&id).unwrap()
                }
            };

            let send = files.get_mut(x).unwrap();
            send.send(Some([header, e1, e2, e3]))
                .expect("Error sending Fastq entry");
        }

        for i in files.values_mut() {
            i.send(None).expect("Unable to send Finish command");
        }
    })
    .unwrap();
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn new_fastq() {
        let not_really_a_fastq = String::from(
            "@MT_E00516:746:HG3WYCCX2:6:1101:2229:1661 1:N:0:NGAGCTAG+NAGCCTGA\nNTG\n+\n#AA\n",
        );
        // split_fastq_by_id(not_really_a_fastq.as_bytes(), "test");
    }
}
