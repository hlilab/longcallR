use rust_htslib::{bam, bam::Record, bam::Format, bam::{Read, ext::BamRecordExtensions}};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use std::process;
use std::collections::{VecDeque, HashMap};
use threadpool::ThreadPool;
use crate::bam_reader::Region;
use crate::base_matrix::BaseMatrix;
use crate::base_matrix::{*};
use crate::isolated_region::find_isolated_regions;
use crate::profile::{*};

pub fn multithread_produce(bam_file: String, thread_size: usize, tx_low: mpsc::Sender<Region>, tx_high: mpsc::Sender<Region>) {
    let pool = ThreadPool::new(thread_size);
    let bam = bam::IndexedReader::from_path(bam_file.clone()).unwrap();
    let bam_header = bam.header().clone();
    let mut contig_names: VecDeque<String> = VecDeque::new();
    for ctg in bam_header.target_names() {
        contig_names.push_back(std::str::from_utf8(ctg).unwrap().to_string().clone());
    }
    while !contig_names.is_empty() {
        let ctg = contig_names.pop_front().unwrap();
        let tx_l = tx_low.clone();
        let tx_h = tx_high.clone();
        let bam_file_clone = bam_file.clone();
        pool.execute(move || {
            let (normal_depth_regions, high_depth_regions) = get_chrom_coverage_intervals(bam_file_clone.clone(), ctg.as_str(), 1000);
            for region in normal_depth_regions {
                tx_l.send(region).unwrap();
            }
            for region in high_depth_regions {
                tx_h.send(region).unwrap();
            }
        });
    }
    pool.join();
}

pub fn multithread_work(bam_file: String, ref_file: String, out_bam: String, thread_size: usize, rx_low: mpsc::Receiver<Region>, rx_high: mpsc::Receiver<Region>) {
    let pool = ThreadPool::new(thread_size);
    let bam_records_queue: Arc<Mutex<VecDeque<bam::Record>>> = Arc::new(Mutex::new(VecDeque::new()));
    let bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_file.clone()).unwrap();
    let header = bam::Header::from_template(bam.header());
    let bam_writer = bam::Writer::from_path(out_bam.clone(), &header, Format::Bam).unwrap();
    let mut bam_writer = Arc::new(Mutex::new(bam_writer));
    let ref_seqs = load_reference(ref_file.to_string());

    // multiplethreads for low coverage regions
    for normal_depth_region in rx_low {
        let bam_file_clone = bam_file.clone();
        let ref_seqs_clone = ref_seqs.clone();
        let bam_records_queue = Arc::clone(&bam_records_queue);
        let bam_writer_clone = Arc::clone(&bam_writer);
        pool.execute(move || {
            let mut base_matrix = BaseMatrix::new();
            base_matrix.load_data(bam_file_clone.clone().to_string(), normal_depth_region);
            base_matrix.load_ref_data(ref_seqs_clone.clone());
            base_matrix.expand_insertion();
            let (forward_donor_penalty,
                forward_acceptor_penalty,
                reverse_donor_penalty,
                reverse_acceptor_penalty) = get_donor_acceptor_penalty(&base_matrix.expanded_matrix, 30.0);
            let mut best_reduced_expanded_matrix: HashMap<String, Vec<u8>> = HashMap::new();
            let mut best_column_indexes: Vec<usize> = Vec::new();
            profile_realign(&base_matrix.expanded_matrix,
                            &forward_donor_penalty,
                            &forward_acceptor_penalty,
                            &reverse_donor_penalty,
                            &reverse_acceptor_penalty,
                            &mut best_reduced_expanded_matrix,
                            &mut best_column_indexes);
            update_expanded_matrix_from_realign(&mut base_matrix.expanded_matrix, &best_reduced_expanded_matrix, &best_column_indexes);
            update_bam_records_from_realign(&mut base_matrix.expanded_matrix, &mut base_matrix.bam_records, base_matrix.start_position, base_matrix.end_position);
            let mut queue = bam_records_queue.lock().unwrap();
            for (_, record) in base_matrix.bam_records {
                queue.push_back(record);
            }
            if queue.len() > 100000 {
                for record in queue.iter() {
                    let re = bam_writer_clone.lock().unwrap().write(&record);
                    if re != Ok(()) {
                        println!("write failed");
                        process::exit(1);
                    }
                }
                queue.clear();
            }
        });
    }
    pool.join();
    if bam_records_queue.lock().unwrap().len() > 0 {
        for record in bam_records_queue.lock().unwrap().iter() {
            let re = bam_writer.lock().unwrap().write(&record);
            if re != Ok(()) {
                println!("write failed");
                process::exit(1);
            }
        }
        bam_records_queue.lock().unwrap().clear();
    }

    let pool = ThreadPool::new(thread_size);
    for high_depth_region in rx_high {
        let bam_file_clone = bam_file.clone();
        let ref_seqs_clone = ref_seqs.clone();
        let bam_records_queue = Arc::clone(&bam_records_queue);
        let bam_writer_clone = Arc::clone(&bam_writer);
        pool.execute(move || {
            let mut base_matrix = BaseMatrix::new();
            base_matrix.load_data_without_extension(bam_file_clone.clone().to_string(), high_depth_region);
            base_matrix.load_ref_data(ref_seqs_clone.clone());
            base_matrix.expand_insertion();
            let (forward_donor_penalty,
                forward_acceptor_penalty,
                reverse_donor_penalty,
                reverse_acceptor_penalty) = get_donor_acceptor_penalty(&base_matrix.expanded_matrix, 30.0);
            let mut best_reduced_expanded_matrix: HashMap<String, Vec<u8>> = HashMap::new();
            let mut best_column_indexes: Vec<usize> = Vec::new();
            profile_realign(&base_matrix.expanded_matrix,
                            &forward_donor_penalty,
                            &forward_acceptor_penalty,
                            &reverse_donor_penalty,
                            &reverse_acceptor_penalty,
                            &mut best_reduced_expanded_matrix,
                            &mut best_column_indexes);
            update_expanded_matrix_from_realign(&mut base_matrix.expanded_matrix, &best_reduced_expanded_matrix, &best_column_indexes);
            update_bam_records_from_realign(&mut base_matrix.expanded_matrix, &mut base_matrix.bam_records, base_matrix.start_position, base_matrix.end_position);
            let mut queue = bam_records_queue.lock().unwrap();
            for (_, record) in base_matrix.bam_records {
                queue.push_back(record);
            }
            if queue.len() > 10000 {
                for record in queue.iter() {
                    let re = bam_writer_clone.lock().unwrap().write(&record);
                    if re != Ok(()) {
                        println!("write failed");
                        process::exit(1);
                    }
                }
                queue.clear();
            }
        });
    }
    pool.join();
    if bam_records_queue.lock().unwrap().len() > 0 {
        for record in bam_records_queue.lock().unwrap().iter() {
            let re = bam_writer.lock().unwrap().write(&record);
            if re != Ok(()) {
                println!("write failed");
                process::exit(1);
            }
        }
        bam_records_queue.lock().unwrap().clear();
    }
}


pub fn multithread_produce2(bam_file: String, thread_size: usize, tx_isolated_regions: mpsc::Sender<Region>) {
    let pool = ThreadPool::new(thread_size);
    let bam = bam::IndexedReader::from_path(bam_file.clone()).unwrap();
    let bam_header = bam.header().clone();
    let mut contig_names: VecDeque<String> = VecDeque::new();
    for ctg in bam_header.target_names() {
        contig_names.push_back(std::str::from_utf8(ctg).unwrap().to_string().clone());
    }
    while !contig_names.is_empty() {
        let ctg = contig_names.pop_front().unwrap();
        let tx_ir = tx_isolated_regions.clone();
        let bam_file_clone = bam_file.clone();
        pool.execute(move || {
            let isolated_regions = find_isolated_regions(bam_file_clone.clone().as_str(), 1, Some(ctg.as_str()));
            for region in isolated_regions {
                tx_ir.send(region).unwrap();
            }
        });
    }
    pool.join();
}

pub fn multithread_work2(bam_file: String, ref_file: String, out_bam: String, thread_size: usize, rx_isolated_regions: mpsc::Receiver<Region>) {
    let pool = ThreadPool::new(thread_size);
    let bam_records_queue: Arc<Mutex<VecDeque<bam::Record>>> = Arc::new(Mutex::new(VecDeque::new()));
    let bam: bam::IndexedReader = bam::IndexedReader::from_path(bam_file.clone()).unwrap();
    let header = bam::Header::from_template(bam.header());
    let bam_writer = bam::Writer::from_path(out_bam.clone(), &header, Format::Bam).unwrap();
    let mut bam_writer = Arc::new(Mutex::new(bam_writer));
    let ref_seqs = load_reference(ref_file.to_string());

    // multiplethreads for low coverage regions
    for reg in rx_isolated_regions {
        let bam_file_clone = bam_file.clone();
        let ref_seqs_clone = ref_seqs.clone();
        let bam_records_queue = Arc::clone(&bam_records_queue);
        let bam_writer_clone = Arc::clone(&bam_writer);
        pool.execute(move || {
            let mut profile = Profile::default();
            let mut readnames: Vec<String> = Vec::new();
            profile.init_with_pileup(bam_file_clone.clone().as_str(), &reg);
            profile.append_reference(&ref_seqs_clone);
            let mut parsed_reads = read_bam(bam_file_clone.clone().as_str(), &reg);
            for (rname, pr) in parsed_reads.iter_mut() {
                pr.init_parsed_seq(&profile);
                readnames.push(rname.clone());
            }
            profile.cal_intron_penalty();
            profile.cal_intron_intervals();
            realign(&mut profile, &mut parsed_reads, &readnames);
            let header = get_bam_header(bam_file_clone.clone().as_str());
            // write_bam(out_bam, &parsed_reads, &header);
            {
                let mut queue = bam_records_queue.lock().unwrap();
                for (_, pr) in parsed_reads.iter() {
                    queue.push_back(pr.bam_record.clone());
                }
                if queue.len() > 10000 {
                    for record in queue.iter() {
                        let re = bam_writer_clone.lock().unwrap().write(&record);
                        if re != Ok(()) {
                            println!("write failed");
                            process::exit(1);
                        }
                    }
                    queue.clear();
                }
            }
        });
    }
    pool.join();
    if bam_records_queue.lock().unwrap().len() > 0 {
        for record in bam_records_queue.lock().unwrap().iter() {
            let re = bam_writer.lock().unwrap().write(&record);
            if re != Ok(()) {
                println!("write failed");
                process::exit(1);
            }
        }
        bam_records_queue.lock().unwrap().clear();
    }
}