use rust_htslib::{bam, bam::{Read, ext::BamRecordExtensions}, bam::Record};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use std::collections::VecDeque;
use threadpool::ThreadPool;
use crate::base_matrix::get_chrom_coverage_intervals;


/*pub fn MultiThreadRun(thread_size: usize, func: fn(T) -> T) {
    let mut handles = vec![];
    for i in 0..thread_size {
        let handle = std::thread::spawn(move || {
            func(T);
        });
        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap();
    }
}*/

/*pub fn MultiThreadRun2(thread_size: usize, func: fn(T) -> T) -> Vec<T> {
    let mut handles = vec![];
    let mut results = vec![];
    let arc = Arc::new(Mutex::new(results));
    for i in 0..thread_size {
        let arc = arc.clone();
        let handle = std::thread::spawn(move || {
            let mut results = arc.lock().unwrap();
            results.push(func);
        });
        handles.push(handle);
    }

    for handle in handles {
        handle.join().unwrap();
    }

    let results = arc.lock().unwrap();
    return results.clone();
}*/

pub fn multithread_produce(bam_file: String, thread_size: usize, tx_low: mpsc::Sender<(String, i64, i64)>, tx_high: mpsc::Sender<(String, i64, i64)>) {
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