use std::collections::HashMap;

use mathru::algebra::abstr::cast::ToPrimitive;
use rust_htslib::{bam, bam::Read, bam::record::Record};

use crate::snp::AlleleClass;
use crate::snpfrags::SNPFrag;
use crate::util::Region;

pub fn calculate_prob_somatic(hap1_ref_baseqs: &Vec<u8>, hap1_alt_baseqs: &Vec<u8>, hap2_ref_baseqs: &Vec<u8>, hap2_alt_baseqs: &Vec<u8>, purity: f64) -> (AlleleClass, AlleleClass) {
    let mut hap1_allele_class;
    let mut hap2_allele_class;
    let som_rate = 5.0 / 1000000.0;    // each haplotype
    let het_rate = 1.0 / 2000.0;    // each haplotype
    let ref_rate = 1.0 - het_rate - som_rate;

    // for Hap1
    let mut prob_read_ref = 1.0;
    let mut prob_read_het = 1.0;
    let mut prob_read_som = 1.0;

    // P(read|ref), P(read|het), P(read|som)
    for q in hap1_ref_baseqs.iter() {
        let epsilon = 10.0_f64.powf(-(q.to_f64() / 10.0));  // error rate
        prob_read_ref *= 1.0 - epsilon; // ref->ref
        prob_read_het *= epsilon;   // alt->ref
        prob_read_som *= purity * epsilon + (1.0 - purity) * (1.0 - epsilon);   // purity: alt->ref, 1.0-purity: ref->ref
    }
    for q in hap1_alt_baseqs.iter() {
        let epsilon = 10.0_f64.powf(-(q.to_f64() / 10.0));
        prob_read_ref *= epsilon;   // ref->alt
        prob_read_het *= 1.0 - epsilon; // alt->alt
        prob_read_som *= purity * (1.0 - epsilon) + (1.0 - purity) * epsilon; // purity: alt->alt, 1.0-purity: ref->alt
    }

    // prob * prior
    let prob_read_ref_with_prior = prob_read_ref * ref_rate;
    let prob_read_het_with_prior = prob_read_het * het_rate;
    let prob_read_som_with_prior = prob_read_som * som_rate;
    let hap1_prob_ref = prob_read_ref_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    let hap1_prob_het = prob_read_het_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    let hap1_prob_som = prob_read_som_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    if hap1_prob_som > hap1_prob_ref && hap1_prob_som > hap1_prob_het {
        hap1_allele_class = AlleleClass { allcls: 2, prob: hap1_prob_som };
    } else if hap1_prob_het > hap1_prob_ref && hap1_prob_het > hap1_prob_som {
        hap1_allele_class = AlleleClass { allcls: 1, prob: hap1_prob_het };
    } else {
        hap1_allele_class = AlleleClass { allcls: 0, prob: hap1_prob_ref };
    }

    // for Hap2
    let mut prob_read_ref = 1.0;
    let mut prob_read_het = 1.0;
    let mut prob_read_som = 1.0;

    // P(read|ref), P(read|het), P(read|som)
    for q in hap2_ref_baseqs.iter() {
        let epsilon = 10.0_f64.powf(-(q.to_f64() / 10.0));  // error rate
        prob_read_ref *= 1.0 - epsilon; // ref->ref
        prob_read_het *= epsilon;   // alt->ref
        prob_read_som *= purity * epsilon + (1.0 - purity) * (1.0 - epsilon);   // purity: alt->ref, 1.0-purity: ref->ref
    }
    for q in hap2_alt_baseqs.iter() {
        let epsilon = 10.0_f64.powf(-(q.to_f64() / 10.0));
        prob_read_ref *= epsilon;   // ref->alt
        prob_read_het *= 1.0 - epsilon; // alt->alt
        prob_read_som *= purity * (1.0 - epsilon) + (1.0 - purity) * epsilon; // purity: alt->alt, 1.0-purity: ref->alt
    }

    // prob * prior
    let prob_read_ref_with_prior = prob_read_ref * ref_rate;
    let prob_read_het_with_prior = prob_read_het * het_rate;
    let prob_read_som_with_prior = prob_read_som * som_rate;
    let hap2_prob_ref = prob_read_ref_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    let hap2_prob_het = prob_read_het_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    let hap2_prob_som = prob_read_som_with_prior / (prob_read_ref_with_prior + prob_read_het_with_prior + prob_read_som_with_prior);
    if hap2_prob_som > hap2_prob_ref && hap2_prob_som > hap2_prob_het {
        hap2_allele_class = AlleleClass { allcls: 2, prob: hap2_prob_som };
    } else if hap2_prob_het > hap2_prob_ref && hap2_prob_het > hap2_prob_som {
        hap2_allele_class = AlleleClass { allcls: 1, prob: hap2_prob_het };
    } else {
        hap2_allele_class = AlleleClass { allcls: 0, prob: hap2_prob_ref };
    }
    return (hap1_allele_class, hap2_allele_class);
}


impl SNPFrag {
    pub fn get_somatic_haplotype_baseqs(&mut self, bam_path: &str, region: &Region, phased_fragments: &HashMap<String, i32>) {
        let mut bam_reader: bam::IndexedReader = bam::IndexedReader::from_path(bam_path).unwrap();
        bam_reader.fetch((region.chr.as_str(), region.start, region.end)).unwrap();
        let mut record = Record::new();
        // assert!(self.min_linkers >= 0, "Error: min_linkers <= 0");
        while let Some(result) = bam_reader.read(&mut record) {
            if self.somatic_snps.len() == 0 {
                continue;
            }
            if result.is_err() {
                panic!("BAM parsing failed...");
            }
            // TODO: filtering unmapped, secondary, supplementary reads?
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
                continue;
            }
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            if !phased_fragments.contains_key(&qname) {
                continue;
            }
            let assignment = phased_fragments[&qname];
            let pos = record.pos(); // 0-based
            if pos > self.candidate_snps[*self.somatic_snps.last().unwrap()].pos {
                continue;
            }
            let cigar = record.cigar();
            let seq = record.seq().as_bytes();
            let mut pos_on_ref = pos; // 0-based
            let mut pos_on_query = cigar.leading_softclips(); // 0-based
            let mut idx = 0; // index in self.somatic_snps
            let mut snp_pos = -1; // pre-computed position of candidate SNPs
            let mut alleles; // pre-computed alleles of candidate SNPs
            if pos <= self.candidate_snps[*self.somatic_snps.first().unwrap()].pos {
                snp_pos = self.candidate_snps[self.somatic_snps[idx]].pos;
                alleles = self.candidate_snps[self.somatic_snps[idx]].alleles.clone();
            } else {
                // find the first SNP in the read
                while idx < self.somatic_snps.len() {
                    if self.candidate_snps[self.somatic_snps[idx]].pos >= pos {
                        break;
                    }
                    idx += 1;
                }
                assert!(
                    idx < self.somatic_snps.len(),
                    "Error: idx < self.candidate_snps.len()"
                );
                snp_pos = self.candidate_snps[self.somatic_snps[idx]].pos;
                alleles = self.candidate_snps[self.somatic_snps[idx]].alleles.clone();
            }

            for cg in cigar.iter() {
                match cg.char() as u8 {
                    b'S' | b'H' => {
                        continue;
                    }
                    b'M' | b'X' | b'=' => {
                        for _ in 0..cg.len() {
                            if pos_on_ref == snp_pos {
                                let somatic_cand = &mut self.candidate_snps[self.somatic_snps[idx]];
                                let base = seq[pos_on_query as usize] as char;
                                let baseq = record.qual()[pos_on_query as usize];
                                let [allele1, allele2] = somatic_cand.alleles.clone();
                                let ref_allele = somatic_cand.reference;
                                if allele1 == ref_allele || allele2 == ref_allele {
                                    if base == allele1 || base == allele2 {
                                        if base == ref_allele {
                                            if assignment == 1 {
                                                somatic_cand.hap_quals.hap1_ref_baseqs.push(baseq);
                                            } else {
                                                somatic_cand.hap_quals.hap2_ref_baseqs.push(baseq);
                                            }
                                        } else {
                                            if assignment == 1 {
                                                somatic_cand.hap_quals.hap1_alt_baseqs.push(baseq);
                                            } else {
                                                somatic_cand.hap_quals.hap2_alt_baseqs.push(baseq);
                                            }
                                        }
                                    }
                                }
                                idx += 1;
                                if idx < self.somatic_snps.len() {
                                    snp_pos = self.candidate_snps[self.somatic_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.somatic_snps[idx]].alleles.clone();
                                }
                            }
                            pos_on_query += 1;
                            pos_on_ref += 1;
                        }
                    }
                    b'I' => {
                        pos_on_query += cg.len() as i64;
                    }
                    b'D' => {
                        for _ in 0..cg.len() {
                            if pos_on_ref == snp_pos {
                                idx += 1;
                                if idx < self.somatic_snps.len() {
                                    snp_pos = self.candidate_snps[self.somatic_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.somatic_snps[idx]].alleles.clone();
                                }
                            }
                            pos_on_ref += 1;
                        }
                    }
                    b'N' => {
                        for _ in 0..cg.len() {
                            if pos_on_ref == snp_pos {
                                idx += 1;
                                if idx < self.somatic_snps.len() {
                                    snp_pos = self.candidate_snps[self.somatic_snps[idx]].pos;
                                    alleles = self.candidate_snps[self.somatic_snps[idx]].alleles.clone();
                                }
                            }
                            pos_on_ref += 1;
                        }
                    }
                    _ => {
                        panic!("Error: unknown cigar operation: {}", cg.char());
                    }
                }
            }
        }
        // for idx in self.somatic_snps.iter() {
        //     println!("somatic snp:{:?}\n, hap1_ref_baseqs:{:?}\n, hap1_alt_baseqs:{:?}\n, hap2_ref_baseqs:{:?}\n, hap2_alt_baseqs:{:?}\n",
        //              self.candidate_snps[*idx].pos,
        //              self.candidate_snps[*idx].hap_quals.hap1_ref_baseqs,
        //              self.candidate_snps[*idx].hap_quals.hap1_alt_baseqs,
        //              self.candidate_snps[*idx].hap_quals.hap2_ref_baseqs,
        //              self.candidate_snps[*idx].hap_quals.hap2_alt_baseqs);
        // }
    }

    pub fn detect_somatic_by_het(&mut self, bam_path: &str, region: &Region) {
        if self.somatic_snps.len() == 0 {
            return;
        }
        // detect confident somatic mutation with phasing result of high allele fraction het snps.
        // 1. assign phased result for candidate somatic sites.
        let mut phased_fragments: HashMap<String, i32> = HashMap::new(); // read id, assignment
        for k in 0..self.fragments.len() {
            let frag = &self.fragments[k];
            if frag.assignment == 1 || frag.assignment == 2 {
                phased_fragments.insert(frag.read_id.clone(), frag.assignment);
            }
        }
        self.get_somatic_haplotype_baseqs(bam_path, region, &phased_fragments);
        // 2. find candidates meet the criteria of somatic mutation. haplotype-specific
        for i in 0..self.somatic_snps.len() {
            let som_cand = &mut self.candidate_snps[self.somatic_snps[i]];
            let (hap1_allele_class, hap2_allele_class) = calculate_prob_somatic(&som_cand.hap_quals.hap1_ref_baseqs, &som_cand.hap_quals.hap1_alt_baseqs, &som_cand.hap_quals.hap2_ref_baseqs, &som_cand.hap_quals.hap2_alt_baseqs, 0.3);
            if hap1_allele_class.allcls == 0 && hap2_allele_class.allcls == 2 {
                let somatic_score = -10.0_f64 * (1.0 - hap2_allele_class.prob).log10();
                som_cand.somatic = true;
                som_cand.variant_type = 1;
                som_cand.somatic_score = somatic_score;
                // println!("somatic snp:{}, score: {}", som_cand.pos, somatic_score);
                // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
                // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", som_cand.hap_quals.hap1_ref_baseqs, som_cand.hap_quals.hap1_alt_baseqs, som_cand.hap_quals.hap2_ref_baseqs, som_cand.hap_quals.hap2_alt_baseqs);
            } else if hap1_allele_class.allcls == 2 && hap2_allele_class.allcls == 0 {
                let somatic_score = -10.0_f64 * (1.0 - hap1_allele_class.prob).log10();
                som_cand.somatic = true;
                som_cand.variant_type = 1;
                som_cand.somatic_score = somatic_score;
                // println!("somatic snp:{}, score: {}", som_cand.pos, somatic_score);
                // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
                // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", som_cand.hap_quals.hap1_ref_baseqs, som_cand.hap_quals.hap1_alt_baseqs, som_cand.hap_quals.hap2_ref_baseqs, som_cand.hap_quals.hap2_alt_baseqs);
            }
        }
    }
}