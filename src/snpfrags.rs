use std::collections::{HashMap, HashSet};

use petgraph::algo::kosaraju_scc;
use petgraph::graphmap::GraphMap;
use petgraph::Undirected;
use rand::Rng;

use crate::snp::{CandidateSNP, Edge, Fragment};
use crate::somatic::calculate_prob_somatic;
use crate::util::Region;

#[derive(Debug, Clone, Default)]
pub struct SNPFrag {
    pub region: Region,
    pub candidate_snps: Vec<CandidateSNP>,
    // candidate SNPs
    pub homo_snps: Vec<usize>,
    // index of candidate homozygous SNPs
    pub edit_snps: Vec<usize>,
    // index of candidate rna editing SNPs
    pub low_frac_het_snps: Vec<usize>,
    // index of candidate low fraction het SNPs
    pub high_frac_het_snps: Vec<usize>,
    // index of candidate high fraction het SNPs
    pub somatic_snps: Vec<usize>,
    // index of candidate somatic mutation
    pub fragments: Vec<Fragment>,
    // multiple fragments
    pub phased: bool,
    // haplotype is phased or not
    pub edges: HashMap<[usize; 2], Edge>,
    // edges of the graph, key is [snp_idx of start_node, snp_idx of end_node]
    pub min_linkers: u32,
    // the number of links for snps can be phased
}

impl SNPFrag {
    pub unsafe fn init_haplotypes(&mut self) {
        // initialize haplotype of heterozygous snp
        let mut rng = rand::thread_rng();
        for i in self.high_frac_het_snps.iter() {
            let rg: f64 = rng.gen();
            if rg < 0.5 {
                self.candidate_snps[*i].haplotype = -1;
            } else {
                self.candidate_snps[*i].haplotype = 1;
            }
        }
    }

    pub unsafe fn init_assignment(&mut self) {
        let mut rng = rand::thread_rng();
        for k in 0..self.fragments.len() {
            if self.fragments[k].num_hete_links < self.min_linkers {
                continue;
            }
            let rg: f64 = rng.gen();
            if rg < 0.5 {
                self.fragments[k].haplotag = -1;
            } else {
                self.fragments[k].haplotag = 1;
            }
        }
    }

    pub fn cal_sigma_delta(sigma_k: i32, delta: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
        // calculate P(sigma_k | delta)
        // sigma_k: the assignment of read k, 1 or -1.
        // delta: the haplotypes of the SNPs covered by read k, each haplotype is 1 or -1.
        // ps: the allele of each base, 1,-1
        // probs: the probability of observing base at each SNP for read k, equals to 10^(-Q/10).
        let mut log_q1: f64 = 0.0;
        let mut log_q2: f64 = 0.0;
        let mut log_q3: f64 = 0.0;

        for i in 0..delta.len() {
            if sigma_k * delta[i] == ps[i] {
                log_q1 += (1.0 - probs[i]).log10();
            } else {
                log_q1 += probs[i].log10();
            }
        }

        for i in 0..delta.len() {
            if delta[i] == ps[i] {
                log_q2 += (1.0 - probs[i]).log10();
                log_q3 += probs[i].log10();
            } else {
                log_q2 += probs[i].log10();
                log_q3 += (1.0 - probs[i]).log10();
            }
        }
        let max_logq = log_q1.max(log_q2.max(log_q3));
        let q1 = 10.0_f64.powf(log_q1 - max_logq);
        let q2 = 10.0_f64.powf(log_q2 - max_logq);
        let q3 = 10.0_f64.powf(log_q3 - max_logq);
        // println!("sigma delta q1:{:?}, q2+q3:{:?}", q1, q2 + q3);
        return q1 / (q2 + q3);
    }

    pub fn cal_sigma_delta_log(
        sigma_k: i32,
        delta: &Vec<i32>,
        ps: &Vec<i32>,
        probs: &Vec<f64>,
    ) -> f64 {
        // same as call_sigma_delta, but return log10 value to avoid underflow
        let mut log_q1: f64 = 0.0;
        let mut log_q2: f64 = 0.0;
        let mut log_q3: f64 = 0.0;

        for i in 0..delta.len() {
            if sigma_k * delta[i] == ps[i] {
                log_q1 += (1.0 - probs[i]).log10();
            } else {
                log_q1 += probs[i].log10();
            }
        }

        for i in 0..delta.len() {
            if delta[i] == ps[i] {
                log_q2 += (1.0 - probs[i]).log10();
                log_q3 += probs[i].log10();
            } else {
                log_q2 += probs[i].log10();
                log_q3 += (1.0 - probs[i]).log10();
            }
        }

        // The exact formula: logP = log(\frac{A}{A+B})=logA-log(A+B)=logA-log(10^{logA}+10^{logB})
        // let log_p = log_q1 - f64::log10(10.0_f64.powf(log_q2) + 10.0_f64.powf(log_q3));
        // to avoid underflow, use approximate 1.0-log(A)/(log(A)+log(B)) as A/(A+B).
        // 0.99/(0.99+0.01) = 0.99, log(0.99)/(log(0.99)+log(0.01)) = 0.00217765
        // 0.01/(0.99+0.01) = 0.01, log(0.01)/(log(0.99)+log(0.01)) = 0.99782235
        let log_p = 1.0 - log_q1 / (log_q2 + log_q3);
        return log_p;
    }

    pub fn cal_delta_sigma(delta_i: i32, sigma: &Vec<i32>, ps: &Vec<i32>, probs: &Vec<f64>) -> f64 {
        // calculate P(delta_i | sigma)
        // delta_i: the haplotype of SNP i, 1 or -1.
        // sigma: the assignments of the reads cover SNP i, each haplotype is 1 or -1.
        // ps: the allele of each base, 1,-1
        // probs: the probability of observing base at SNP i for each read, equals to 10^(-Q/10).

        let mut log_q1: f64 = 0.0;
        let mut log_q2: f64 = 0.0;
        let mut log_q3: f64 = 0.0;

        for k in 0..sigma.len() {
            if delta_i * sigma[k] == ps[k] {
                log_q1 += (1.0 - probs[k]).log10();
            } else {
                log_q1 += probs[k].log10();
            }
        }

        for k in 0..sigma.len() {
            if sigma[k] == ps[k] {
                log_q2 += (1.0 - probs[k]).log10();
                log_q3 += probs[k].log10();
            } else {
                log_q2 += probs[k].log10();
                log_q3 += (1.0 - probs[k]).log10();
            }
        }
        let max_logq = log_q1.max(log_q2.max(log_q3));
        let q1 = 10.0_f64.powf(log_q1 - max_logq);
        let q2 = 10.0_f64.powf(log_q2 - max_logq);
        let q3 = 10.0_f64.powf(log_q3 - max_logq);
        // println!("delta sigma q1:{:?}, q2+q3:{:?}", q1, q2 + q3);
        return q1 / (q2 + q3);
    }

    pub fn cal_delta_sigma_log(
        delta_i: i32,
        sigma: &Vec<i32>,
        ps: &Vec<i32>,
        probs: &Vec<f64>,
    ) -> f64 {
        // same as call_delta_sigma, but return log10 value to avoid underflow
        let mut log_q1: f64 = 0.0;
        let mut log_q2: f64 = 0.0;
        let mut log_q3: f64 = 0.0;

        for k in 0..sigma.len() {
            if delta_i * sigma[k] == ps[k] {
                log_q1 += (1.0 - probs[k]).log10();
            } else {
                log_q1 += probs[k].log10();
            }
        }

        for k in 0..sigma.len() {
            if sigma[k] == ps[k] {
                log_q2 += (1.0 - probs[k]).log10();
                log_q3 += probs[k].log10();
            } else {
                log_q2 += probs[k].log10();
                log_q3 += (1.0 - probs[k]).log10();
            }
        }

        // logP = log(\frac{A}{A+B})=logA-log(A+B)=logA-log(10^{logA}+10^{logB})
        // let log_p = log_q1 - f64::log10(10.0_f64.powf(log_q2) + 10.0_f64.powf(log_q3));
        // to avoid underflow, use approximate 1.0-log(A)/(log(A)+log(B)) as A/(A+B).
        // 0.99/(0.99+0.01) = 0.99, log(0.99)/(log(0.99)+log(0.01)) = 0.00217765
        // 0.01/(0.99+0.01) = 0.01, log(0.01)/(log(0.99)+log(0.01)) = 0.99782235
        let log_p = 1.0 - log_q1 / (log_q2 + log_q3);
        return log_p;
    }

    pub fn cal_inconsistent_percentage(
        delta_i: i32,
        sigma: &Vec<i32>,
        ps: &Vec<i32>,
        probs: &Vec<f64>,
    ) -> f64 {
        let mut consisitent = 0;
        let mut inconsistent = 0;
        for k in 0..sigma.len() {
            if delta_i * sigma[k] == ps[k] {
                consisitent += 1;
            } else {
                inconsistent += 1;
            }
        }
        return 10e-6_f64.max((inconsistent as f64) / ((consisitent + inconsistent) as f64));
    }

    pub fn cal_overall_probability(snpfrag: &SNPFrag) -> f64 {
        // calculate the log10 probability of the current configuration of sigma and delta
        let mut logp = 0.0;
        for k in 0..snpfrag.fragments.len() {
            if snpfrag.fragments[k].haplotag == 0 {
                continue;
            }
            for fe in snpfrag.fragments[k].list.iter() {
                if fe.phase_site == false {
                    continue;
                }
                assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
                if snpfrag.fragments[k].haplotag * snpfrag.candidate_snps[fe.snp_idx].haplotype == fe.p {
                    logp += (1.0 - fe.prob).log10();
                } else {
                    logp += fe.prob.log10();
                }
            }
        }
        return logp;
    }

    pub fn cal_overall_probability_ase(snpfrag: &SNPFrag) -> f64 {
        // calculate the log10 probability of the current configuration of sigma and delta
        let mut logp = 0.0;
        for k in 0..snpfrag.fragments.len() {
            for fe in snpfrag.fragments[k].list.iter() {
                assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
                if snpfrag.fragments[k].haplotag * snpfrag.candidate_snps[fe.snp_idx].haplotype == fe.p {
                    logp += (1.0 - fe.prob).log10();
                } else {
                    logp += fe.prob.log10();
                }
            }
        }
        return logp;
    }

    pub fn check_new_haplotag(snpfrag: &SNPFrag, updated_haplotag: &HashMap<usize, i32>) -> i32 {
        // updated_haplotag: the index of the fragments will be updated
        let mut logp = 0.0;
        let mut pre_logp = 0.0;
        for (k, h) in updated_haplotag.iter() {
            let mut delta: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            if snpfrag.fragments[*k].haplotag == 0 {
                continue;
            }
            for fe in snpfrag.fragments[*k].list.iter() {
                if fe.phase_site == false {
                    continue;
                }
                assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
                ps.push(fe.p);
                probs.push(fe.prob);
                delta.push(snpfrag.candidate_snps[fe.snp_idx].haplotype);
            }
            if delta.len() == 0 {
                continue;
            }
            logp += SNPFrag::cal_sigma_delta_log(*h, &delta, &ps, &probs);
            pre_logp += SNPFrag::cal_sigma_delta_log(snpfrag.fragments[*k].haplotag, &delta, &ps, &probs);
        }

        let p = logp;
        let pre_p = pre_logp;
        // println!("haplotag p:{}, pre_p:{}", p, pre_p);
        let mut rv = 0;
        if p > pre_p {
            rv = 1;
        } else if p == pre_p {
            rv = 0;
        } else {
            rv = -1;
        }
        assert!(
            rv >= 0,
            "Error: update haplotag should not decrease the probability. {} -> {}",
            pre_logp,
            logp
        );
        return rv;
    }

    pub fn check_new_haplotag_ase(
        snpfrag: &SNPFrag,
        updated_haplotag: &HashMap<usize, i32>,
    ) -> i32 {
        // updated_haplotag: the index of the fragments will be updated
        let mut logp = 0.0;
        let mut pre_logp = 0.0;
        for (k, h) in updated_haplotag.iter() {
            let mut delta: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            for fe in snpfrag.fragments[*k].list.iter() {
                ps.push(fe.p);
                probs.push(fe.prob);
                delta.push(snpfrag.candidate_snps[fe.snp_idx].haplotype);
            }
            if delta.len() == 0 {
                continue;
            }
            logp += SNPFrag::cal_sigma_delta_log(*h, &delta, &ps, &probs);
            pre_logp += SNPFrag::cal_sigma_delta_log(snpfrag.fragments[*k].haplotag, &delta, &ps, &probs);
        }

        let p = logp;
        let pre_p = pre_logp;
        // println!("haplotag p:{}, pre_p:{}", p, pre_p);
        let mut rv = 0;
        if p > pre_p {
            rv = 1;
        } else if p == pre_p {
            rv = 0;
        } else {
            rv = -1;
        }
        assert!(
            rv >= 0,
            "Error: update haplotag should not decrease the probability. {} -> {}",
            pre_logp,
            logp
        );
        return rv;
    }

    pub fn check_new_haplotype(snpfrag: &SNPFrag, updated_haplotype: &HashMap<usize, i32>) -> i32 {
        let mut logp = 0.0;
        let mut pre_logp = 0.0;
        for (i, h) in updated_haplotype.iter() {
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            for k in snpfrag.candidate_snps[*i].snp_cover_fragments.iter() {
                if snpfrag.fragments[*k].haplotag == 0 {
                    continue;
                }
                for fe in snpfrag.fragments[*k].list.iter() {
                    if fe.snp_idx != *i {
                        continue;
                    }
                    if fe.phase_site == false {
                        continue;
                    }
                    assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
                    ps.push(fe.p);
                    probs.push(fe.prob);
                    sigma.push(snpfrag.fragments[*k].haplotag);
                }
            }
            if sigma.len() == 0 {
                // println!("SNP: {:?}", snpfrag.candidate_snps[*i]);
                // println!("SNP {} is not covered by any fragment.", snpfrag.candidate_snps[*i].pos);
                continue;
            }
            logp += SNPFrag::cal_delta_sigma_log(*h, &sigma, &ps, &probs);
            pre_logp += SNPFrag::cal_delta_sigma_log(
                snpfrag.candidate_snps[*i].haplotype,
                &sigma,
                &ps,
                &probs,
            );
        }
        let p = logp;
        let pre_p = pre_logp;
        // println!("haplotype p:{}, pre_p:{}", p, pre_p);
        let mut rv = 0;
        if p > pre_p {
            rv = 1;
        } else if p == pre_p {
            rv = 0;
        } else {
            rv = -1;
        }
        assert!(
            rv >= 0,
            "Error: update haplotype should not decrease the probability. {} -> {}",
            pre_logp,
            logp
        );
        return rv;
    }

    pub fn check_new_haplotype_ase(
        snpfrag: &SNPFrag,
        updated_haplotype: &HashMap<usize, i32>,
    ) -> i32 {
        let mut logp = 0.0;
        let mut pre_logp = 0.0;
        for (i, h) in updated_haplotype.iter() {
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            for k in snpfrag.candidate_snps[*i].snp_cover_fragments.iter() {
                for fe in snpfrag.fragments[*k].list.iter() {
                    if fe.snp_idx != *i {
                        continue;
                    }
                    assert_ne!(fe.p, 0, "Error: phasing with unexpected hete SNP.");
                    ps.push(fe.p);
                    probs.push(fe.prob);
                    sigma.push(snpfrag.fragments[*k].haplotag);
                }
            }
            if sigma.len() == 0 {
                // println!("SNP: {:?}", snpfrag.candidate_snps[*i]);
                // println!("SNP {} is not covered by any fragment.", snpfrag.candidate_snps[*i].pos);
                continue;
            }
            logp += SNPFrag::cal_delta_sigma_log(*h, &sigma, &ps, &probs);
            pre_logp += SNPFrag::cal_delta_sigma_log(
                snpfrag.candidate_snps[*i].haplotype,
                &sigma,
                &ps,
                &probs,
            );
        }
        let p = logp;
        let pre_p = pre_logp;
        // println!("haplotype p:{}, pre_p:{}", p, pre_p);
        let mut rv = 0;
        if p > pre_p {
            rv = 1;
        } else if p == pre_p {
            rv = 0;
        } else {
            rv = -1;
        }
        assert!(
            rv >= 0,
            "Error: update haplotype should not decrease the probability. {} -> {}",
            pre_logp,
            logp
        );
        return rv;
    }

    pub fn cross_optimize(&mut self) -> f64 {
        // Iteration:
        //     1. evaluate the assignment of each read based on the current SNP haplotype.
        //     2. evaluate the SNP haplotype based on the read assignment.
        // If P(sigma, delta) increase, repeat Iteration;
        // Else break;

        let mut phasing_increase: bool = true;
        let mut haplotag_increase: bool = true;
        let mut num_iters = 0;

        while phasing_increase | haplotag_increase {
            // optimize sigma
            let mut tmp_haplotag: HashMap<usize, i32> = HashMap::new();
            let mut processed_snps = HashSet::new(); // some snps in self.hete_snps may be filtered by previous steps, record the snps that covered by the fragments
            for k in 0..self.fragments.len() {
                let sigma_k = self.fragments[k].haplotag;
                let mut delta: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                if sigma_k == 0 {
                    continue;
                }
                for fe in self.fragments[k].list.iter() {
                    if fe.phase_site == false {
                        continue;
                    }
                    assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                    ps.push(fe.p);
                    probs.push(fe.prob);
                    delta.push(self.candidate_snps[fe.snp_idx].haplotype);
                    processed_snps.insert(fe.snp_idx);
                }

                let q = SNPFrag::cal_sigma_delta_log(sigma_k, &delta, &ps, &probs);
                let qn = SNPFrag::cal_sigma_delta_log(sigma_k * (-1), &delta, &ps, &probs);
                // println!("optimize sigma {} q:{}, qn:{}, sigma: {}", k, q, qn, sigma_k);

                if q < qn {
                    tmp_haplotag.insert(k, sigma_k * (-1));
                } else {
                    tmp_haplotag.insert(k, sigma_k);
                }
            }

            // assert!(SNPFrag::cal_overall_probability(&self, &processed_snps, &self.haplotype) >= SNPFrag::cal_overall_probability(&self, &self.haplotag, &self.haplotype));
            let check_val = SNPFrag::check_new_haplotag(&self, &tmp_haplotag);
            assert!(check_val >= 0, "ckeck val bug: {:?}", self.candidate_snps);
            for (k, h) in tmp_haplotag.iter() {
                // when prob is equal, we still perform the flip to avoid bug of underflow
                self.fragments[*k].haplotag = *h;
            }
            if check_val == 0 {
                haplotag_increase = false;
            } else {
                haplotag_increase = true;
                phasing_increase = true;
            }
            self.check_local_optimal_configuration(false, true);

            // optimize delta
            let mut tmp_haplotype: HashMap<usize, i32> = HashMap::new();
            for i in self.high_frac_het_snps.iter() {
                let delta_i = self.candidate_snps[*i].haplotype;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                    if self.fragments[*k].haplotag == 0 {
                        continue;
                    }
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == *i {
                            assert_ne!(fe.phase_site, false, "Error: phase for non-phase site.");
                            assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                            ps.push(fe.p);
                            probs.push(fe.prob);
                            sigma.push(self.fragments[*k].haplotag);
                        }
                    }
                }

                let q = SNPFrag::cal_delta_sigma_log(delta_i, &sigma, &ps, &probs);
                let qn = SNPFrag::cal_delta_sigma_log(delta_i * (-1), &sigma, &ps, &probs);
                // println!("optimize delta {} q:{:?}, qn:{:?}, delta: {}", i, q, qn, delta_i);
                if q < qn {
                    tmp_haplotype.insert(*i, delta_i * (-1));
                } else {
                    tmp_haplotype.insert(*i, delta_i);
                }
            }
            let check_val = SNPFrag::check_new_haplotype(&self, &tmp_haplotype);
            assert!(check_val >= 0, "ckeck val bug: {:?}", self.candidate_snps);
            for (i, h) in tmp_haplotype.iter() {
                // when prob is equal, we still perform the flip to avoid bug of underflow
                self.candidate_snps[*i].haplotype = *h;
            }
            if check_val == 0 {
                phasing_increase = false;
            } else {
                phasing_increase = true;
                haplotag_increase = true;
            }
            self.check_local_optimal_configuration(true, false);
            num_iters += 1;
            if num_iters > 20 {
                break;
            }
        }
        // sigma reaches the optimal solution first and then delta reaches the optimal solution. After this, equal probability flip of delta may destroy the optimum of sigma again.
        // self.check_local_optimal_configuration(true, true);
        let prob = SNPFrag::cal_overall_probability(&self);
        return prob;
    }

    fn save_best_configuration(
        &self,
        best_haplotype: &mut HashMap<usize, i32>,
        best_haplotag: &mut HashMap<usize, i32>,
    ) {
        best_haplotype.clear();
        best_haplotag.clear();
        for i in self.high_frac_het_snps.iter() {
            best_haplotype.insert(*i, self.candidate_snps[*i].haplotype);
        }
        for k in 0..self.fragments.len() {
            best_haplotag.insert(k, self.fragments[k].haplotag);
        }
    }

    fn load_best_configuration(
        &mut self,
        best_haplotype: &HashMap<usize, i32>,
        best_haplotag: &HashMap<usize, i32>,
    ) {
        for i in self.high_frac_het_snps.iter() {
            self.candidate_snps[*i].haplotype = best_haplotype.get(i).unwrap().clone();
        }
        for k in 0..self.fragments.len() {
            self.fragments[k].haplotag = best_haplotag.get(&k).unwrap().clone();
        }
    }

    pub fn phase(&mut self, max_enum_snps: usize, random_flip_fraction: f32, max_iters: i32) {
        let mut largest_prob = f64::NEG_INFINITY;
        let mut best_haplotype: HashMap<usize, i32> = HashMap::new();
        let mut best_haplotag: HashMap<usize, i32> = HashMap::new();

        if self.high_frac_het_snps.len() <= max_enum_snps {
            // enumerate the haplotype, then optimize the assignment
            let mut haplotype_enum: Vec<Vec<i32>> = Vec::new();
            let init_hap: Vec<i32> = vec![1; self.high_frac_het_snps.len()];
            haplotype_enum.push(init_hap.clone());
            for ti in 0..self.high_frac_het_snps.len() {
                for tj in 0..haplotype_enum.len() {
                    let mut tmp_hap = haplotype_enum[tj].clone();
                    tmp_hap[ti] = tmp_hap[ti] * (-1);
                    haplotype_enum.push(tmp_hap);
                }
            }
            assert!(haplotype_enum.len() == 2_usize.pow(self.high_frac_het_snps.len() as u32), "Error: Not all combinations included");
            for hap in haplotype_enum.iter() {
                for i in 0..self.high_frac_het_snps.len() {
                    self.candidate_snps[self.high_frac_het_snps[i]].haplotype = hap[i];
                }
                unsafe {
                    self.init_assignment();
                }
                let prob = self.cross_optimize();
                if prob > largest_prob {
                    largest_prob = prob;
                    self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                }
            }
            self.load_best_configuration(&best_haplotype, &best_haplotag);
        } else {
            // optimize haplotype and read assignment alternatively
            let mut max_iter: i32 = max_iters;
            while max_iter >= 0 {
                // random initialization of haplotype and haplotag at each iteration
                unsafe {
                    self.init_haplotypes();
                }
                unsafe {
                    self.init_assignment();
                }
                let prob = self.cross_optimize();
                if prob > largest_prob {
                    largest_prob = prob;
                    self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                }
                self.load_best_configuration(&best_haplotype, &best_haplotag);

                // when initial setting has reached to local optimal, flip all the snps after a specific position to jump out local optimization
                let mut unflipped_haplotype: Vec<i32> = Vec::new();
                for i in self.high_frac_het_snps.iter() {
                    unflipped_haplotype.push(self.candidate_snps[*i].haplotype);
                }
                for ti in 0..unflipped_haplotype.len() {
                    let mut tmp_hap: Vec<i32> = Vec::new();
                    for tj in 0..unflipped_haplotype.len() {
                        if tj < ti {
                            tmp_hap.push(unflipped_haplotype[tj]);
                        } else {
                            tmp_hap.push(unflipped_haplotype[tj] * (-1));
                        }
                    }
                    // block flip
                    {
                        assert_eq!(tmp_hap.len(), self.high_frac_het_snps.len());
                        for i in 0..self.high_frac_het_snps.len() {
                            self.candidate_snps[self.high_frac_het_snps[i]].haplotype = tmp_hap[i];
                        }
                    }
                    let prob = self.cross_optimize();
                    if prob > largest_prob {
                        largest_prob = prob;
                        self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                    }
                    self.load_best_configuration(&best_haplotype, &best_haplotag);

                    // when current block flip has reached to local optimal, flip a fraction of snps and reads to jump out local optimization
                    {
                        let mut rng = rand::thread_rng();
                        for ti in 0..self.high_frac_het_snps.len() {
                            let rg: f64 = rng.gen();
                            if rg < random_flip_fraction as f64 {
                                self.candidate_snps[self.high_frac_het_snps[ti]].haplotype = self.candidate_snps[self.high_frac_het_snps[ti]].haplotype * (-1);
                            }
                        }
                        for tk in 0..self.fragments.len() {
                            if self.fragments[tk].haplotag == 0 {
                                continue;
                            }
                            let rg: f64 = rng.gen();
                            if rg < random_flip_fraction as f64 {
                                self.fragments[tk].haplotag = self.fragments[tk].haplotag * (-1);
                            }
                        }
                    }
                    let prob = self.cross_optimize();
                    if prob > largest_prob {
                        largest_prob = prob;
                        self.save_best_configuration(&mut best_haplotype, &mut best_haplotag);
                    }
                    self.load_best_configuration(&best_haplotype, &best_haplotag);
                }
                self.load_best_configuration(&best_haplotype, &best_haplotag);
                max_iter -= 1;
            }
            self.load_best_configuration(&best_haplotype, &best_haplotag);
        }
    }

    pub fn assign_het_var_haplotype(
        &mut self,
        min_phase_score: f32,
        min_homozygous_freq: f32,
        somatic_allele_frac_cutoff: f32,
        somatic_allele_cnt_cutoff: u32,
    ) {
        // calculate phase score for each snp
        for ti in self.high_frac_het_snps.iter() {
            let snp = &mut self.candidate_snps[*ti];
            if !snp.for_phasing { continue; }
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
                continue;
            }
            let delta_i = snp.haplotype;
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            let mut assigns: Vec<i32> = Vec::new();
            let mut baseqs = Vec::new();
            let mut hap1_reads_num = 0;
            let mut hap2_reads_num = 0;
            for k in snp.snp_cover_fragments.iter() {
                if self.fragments[*k].assignment == 0 { continue; }
                if self.fragments[*k].num_hete_links < self.min_linkers { continue; }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == *ti {
                        if fe.base != '-' {
                            // ignore intron
                            if self.fragments[*k].assignment == 1 {
                                hap1_reads_num += 1;
                            } else if self.fragments[*k].assignment == 2 {
                                hap2_reads_num += 1;
                            }
                        }
                        assert_ne!(fe.phase_site, false, "Error: phase for non-phase site.");
                        assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                        ps.push(fe.p);
                        probs.push(fe.prob);
                        sigma.push(self.fragments[*k].haplotag);
                        assigns.push(self.fragments[*k].assignment);
                        baseqs.push(fe.baseq);
                    }
                }
            }

            let mut phase_score = 0.0;
            if sigma.len() > 0 || hap1_reads_num >= 2 || hap2_reads_num >= 2 {
                // each haplotype should have at least 2 reads
                phase_score = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(delta_i, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                if phase_score >= min_phase_score as f64 {
                    let mut haplotype_allele_expression: [u32; 4] = [0, 0, 0, 0];   // hap1_ref, hap1_alt, hap2_ref, hap2_alt
                    for k in 0..sigma.len() {
                        if sigma[k] == 1 {
                            // hap1
                            if ps[k] == 1 {
                                haplotype_allele_expression[0] += 1; // hap1 allele1, reference allele
                            } else if ps[k] == -1 {
                                haplotype_allele_expression[1] += 1; // hap1 allele2, alternative allele
                            }
                        } else if sigma[k] == -1 {
                            // hap2
                            if ps[k] == 1 {
                                haplotype_allele_expression[2] += 1; // hap2 allele1, reference allele
                            } else if ps[k] == -1 {
                                haplotype_allele_expression[3] += 1; // hap2 allele2, alternative allele
                            }
                        }
                    }
                    snp.germline = true;
                    snp.haplotype_expression = haplotype_allele_expression;
                    snp.phase_score = phase_score;
                } else {
                    snp.phase_score = phase_score;
                }
            }

            // TODO: het var with low phase score transfer to hom var
            if phase_score < min_phase_score as f64 {
                if snp.alleles[0] != snp.reference && snp.allele_freqs[0] >= min_homozygous_freq {
                    snp.variant_type = 2;
                    snp.hom_var = true;
                } else if snp.alleles[1] != snp.reference && snp.allele_freqs[1] >= min_homozygous_freq {
                    snp.variant_type = 2;
                    snp.hom_var = true;
                }
            }

            // TODO: het var with low phase score transfer to som var
            if phase_score < min_phase_score as f64 {
                // record HapQuals for somatic mutation detection
                for k in 0..assigns.len() {
                    if assigns[k] == 1 {
                        if ps[k] == 1 {
                            snp.hap_quals.hap1_ref_baseqs.push(baseqs[k]);
                        } else if ps[k] == -1 {
                            snp.hap_quals.hap1_alt_baseqs.push(baseqs[k]);
                        }
                    } else if assigns[k] == 2 {
                        if ps[k] == 1 {
                            snp.hap_quals.hap2_ref_baseqs.push(baseqs[k]);
                        } else if ps[k] == -1 {
                            snp.hap_quals.hap2_alt_baseqs.push(baseqs[k]);
                        }
                    }
                }
                let ref_allele_cnt = snp.hap_quals.hap1_ref_baseqs.len() + snp.hap_quals.hap2_ref_baseqs.len();
                let alt_allele_cnt = snp.hap_quals.hap1_alt_baseqs.len() + snp.hap_quals.hap2_alt_baseqs.len();
                if ref_allele_cnt * alt_allele_cnt > 0 && alt_allele_cnt as u32 >= somatic_allele_cnt_cutoff && alt_allele_cnt as f32 / (ref_allele_cnt + alt_allele_cnt) as f32 >= somatic_allele_frac_cutoff {
                    // calculate somatic mutation probability
                    let (hap1_allele_class, hap2_allele_class) = calculate_prob_somatic(&snp.hap_quals.hap1_ref_baseqs, &snp.hap_quals.hap1_alt_baseqs, &snp.hap_quals.hap2_ref_baseqs, &snp.hap_quals.hap2_alt_baseqs, 0.3);
                    if hap1_allele_class.allcls == 0 && hap2_allele_class.allcls == 2 {
                        let somatic_score = -10.0_f64 * (1.0 - hap2_allele_class.prob).log10();
                        snp.cand_somatic = true;
                        snp.somatic = true;
                        snp.variant_type = 1;
                        snp.somatic_score = somatic_score;
                        snp.phase_score = 0.0;
                        // println!("somatic snp:{}, score: {}", snp.pos, somatic_score);
                        // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
                        // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", snp.hap_quals.hap1_ref_baseqs, snp.hap_quals.hap1_alt_baseqs, snp.hap_quals.hap2_ref_baseqs, snp.hap_quals.hap2_alt_baseqs);
                    } else if hap1_allele_class.allcls == 2 && hap2_allele_class.allcls == 0 {
                        let somatic_score = -10.0_f64 * (1.0 - hap1_allele_class.prob).log10();
                        snp.cand_somatic = true;
                        snp.somatic = true;
                        snp.variant_type = 1;
                        snp.somatic_score = somatic_score;
                        snp.phase_score = 0.0;
                        // println!("somatic snp:{}, score: {}", snp.pos, somatic_score);
                        // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
                        // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", snp.hap_quals.hap1_ref_baseqs, snp.hap_quals.hap1_alt_baseqs, snp.hap_quals.hap2_ref_baseqs, snp.hap_quals.hap2_alt_baseqs);
                    }
                }
            }
        }
    }

    pub fn eval_low_frac_het_var_phase(&mut self, min_phase_score: f32,
                                       somatic_allele_frac_cutoff: f32,
                                       somatic_allele_cnt_cutoff: u32) {
        for ti in self.low_frac_het_snps.iter() {
            let snp = &mut self.candidate_snps[*ti];
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
                continue;
            }
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            let mut assigns: Vec<i32> = Vec::new();
            let mut baseqs = Vec::new();
            let mut hap1_reads_num = 0;
            let mut hap2_reads_num = 0;
            for k in snp.snp_cover_fragments.iter() {
                if self.fragments[*k].assignment == 0 { continue; }
                if self.fragments[*k].num_hete_links < self.min_linkers { continue; }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == *ti {
                        if fe.base != '-' {
                            // ignore intron
                            if self.fragments[*k].assignment == 1 {
                                hap1_reads_num += 1;
                            } else if self.fragments[*k].assignment == 2 {
                                hap2_reads_num += 1;
                            }
                        }
                        assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                        ps.push(fe.p);
                        probs.push(fe.prob);
                        sigma.push(self.fragments[*k].haplotag);
                        assigns.push(self.fragments[*k].assignment);
                        baseqs.push(fe.baseq);
                    }
                }
            }


            let mut phase_score1 = 0.0;
            let mut phase_score2 = 0.0;
            let mut phase_score = 0.0;
            if sigma.len() == 0 || hap1_reads_num < 2 || hap2_reads_num < 2 {
                // each haplotype should have at least 2 reads
                snp.single = true; // no surranding high_frac_het_snps
                continue;
            } else {
                phase_score1 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                phase_score2 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                snp.single = false;
            }

            if phase_score1.max(phase_score2) >= min_phase_score as f64 {
                phase_score = phase_score1.max(phase_score2);
                let mut haplotype_allele_expression: [u32; 4] = [0, 0, 0, 0];   // hap1_ref, hap1_alt, hap2_ref, hap2_alt
                for k in 0..sigma.len() {
                    if sigma[k] == 1 {
                        // hap1
                        if ps[k] == 1 {
                            haplotype_allele_expression[0] += 1; // hap1 allele1, reference allele
                        } else if ps[k] == -1 {
                            haplotype_allele_expression[1] += 1; // hap1 allele2, alternative allele
                        }
                    } else if sigma[k] == -1 {
                        // hap2
                        if ps[k] == 1 {
                            haplotype_allele_expression[2] += 1; // hap2 allele1, reference allele
                        } else if ps[k] == -1 {
                            haplotype_allele_expression[3] += 1; // hap2 allele2, alternative allele
                        }
                    }
                }
                snp.germline = true;
                snp.haplotype = if phase_score1 >= phase_score2 { 1 } else { -1 };
                snp.haplotype_expression = haplotype_allele_expression;
                snp.phase_score = phase_score;
            }

            // TODO: het var with low phase score transfer to som var
            if phase_score < min_phase_score as f64 {
                // record HapQuals for somatic mutation detection
                for k in 0..assigns.len() {
                    if assigns[k] == 1 {
                        if ps[k] == 1 {
                            snp.hap_quals.hap1_ref_baseqs.push(baseqs[k]);
                        } else if ps[k] == -1 {
                            snp.hap_quals.hap1_alt_baseqs.push(baseqs[k]);
                        }
                    } else if assigns[k] == 2 {
                        if ps[k] == 1 {
                            snp.hap_quals.hap2_ref_baseqs.push(baseqs[k]);
                        } else if ps[k] == -1 {
                            snp.hap_quals.hap2_alt_baseqs.push(baseqs[k]);
                        }
                    }
                }
                let ref_allele_cnt = snp.hap_quals.hap1_ref_baseqs.len() + snp.hap_quals.hap2_ref_baseqs.len();
                let alt_allele_cnt = snp.hap_quals.hap1_alt_baseqs.len() + snp.hap_quals.hap2_alt_baseqs.len();
                if ref_allele_cnt * alt_allele_cnt > 0 && alt_allele_cnt as u32 >= somatic_allele_cnt_cutoff && alt_allele_cnt as f32 / (ref_allele_cnt + alt_allele_cnt) as f32 >= somatic_allele_frac_cutoff {
                    // calculate somatic mutation probability
                    let (hap1_allele_class, hap2_allele_class) = calculate_prob_somatic(&snp.hap_quals.hap1_ref_baseqs, &snp.hap_quals.hap1_alt_baseqs, &snp.hap_quals.hap2_ref_baseqs, &snp.hap_quals.hap2_alt_baseqs, 0.3);
                    if hap1_allele_class.allcls == 0 && hap2_allele_class.allcls == 2 {
                        let somatic_score = -10.0_f64 * (1.0 - hap2_allele_class.prob).log10();
                        snp.cand_somatic = true;
                        snp.somatic = true;
                        snp.variant_type = 1;
                        snp.somatic_score = somatic_score;
                        snp.phase_score = 0.0;
                        // println!("somatic snp:{}, score: {}", snp.pos, somatic_score);
                        // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
                        // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", snp.hap_quals.hap1_ref_baseqs, snp.hap_quals.hap1_alt_baseqs, snp.hap_quals.hap2_ref_baseqs, snp.hap_quals.hap2_alt_baseqs);
                    } else if hap1_allele_class.allcls == 2 && hap2_allele_class.allcls == 0 {
                        let somatic_score = -10.0_f64 * (1.0 - hap1_allele_class.prob).log10();
                        snp.cand_somatic = true;
                        snp.somatic = true;
                        snp.variant_type = 1;
                        snp.somatic_score = somatic_score;
                        snp.phase_score = 0.0;
                        // println!("somatic snp:{}, score: {}", snp.pos, somatic_score);
                        // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
                        // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", snp.hap_quals.hap1_ref_baseqs, snp.hap_quals.hap1_alt_baseqs, snp.hap_quals.hap2_ref_baseqs, snp.hap_quals.hap2_alt_baseqs);
                    }
                }
            }
        }
    }
    pub fn eval_rna_edit_var_phase(&mut self, min_phase_score: f32) {
        for ti in self.edit_snps.iter() {
            let snp = &mut self.candidate_snps[*ti];
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
                continue;
            }
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            let mut assigns: Vec<i32> = Vec::new();
            let mut baseqs = Vec::new();
            let mut hap1_reads_num = 0;
            let mut hap2_reads_num = 0;
            for k in snp.snp_cover_fragments.iter() {
                if self.fragments[*k].assignment == 0 { continue; }
                if self.fragments[*k].num_hete_links < self.min_linkers { continue; }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == *ti {
                        if fe.base != '-' {
                            // ignore intron
                            if self.fragments[*k].assignment == 1 {
                                hap1_reads_num += 1;
                            } else if self.fragments[*k].assignment == 2 {
                                hap2_reads_num += 1;
                            }
                        }
                        assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                        ps.push(fe.p);
                        probs.push(fe.prob);
                        sigma.push(self.fragments[*k].haplotag);
                        assigns.push(self.fragments[*k].assignment);
                        baseqs.push(fe.baseq);
                    }
                }
            }

            let mut phase_score1 = 0.0;
            let mut phase_score2 = 0.0;
            let mut phase_score = 0.0;
            if sigma.len() == 0 || hap1_reads_num < 2 || hap2_reads_num < 2 {
                // each haplotype should have at least 2 reads
                snp.single = true; // no surranding high_frac_het_snps
                continue;
            } else {
                phase_score1 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                phase_score2 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                snp.single = false;
            }

            if phase_score1.max(phase_score2) >= min_phase_score as f64 {
                phase_score = phase_score1.max(phase_score2);
                let mut haplotype_allele_expression: [u32; 4] = [0, 0, 0, 0];   // hap1_ref, hap1_alt, hap2_ref, hap2_alt
                for k in 0..sigma.len() {
                    if sigma[k] == 1 {
                        // hap1
                        if ps[k] == 1 {
                            haplotype_allele_expression[0] += 1; // hap1 allele1, reference allele
                        } else if ps[k] == -1 {
                            haplotype_allele_expression[1] += 1; // hap1 allele2, alternative allele
                        }
                    } else if sigma[k] == -1 {
                        // hap2
                        if ps[k] == 1 {
                            haplotype_allele_expression[2] += 1; // hap2 allele1, reference allele
                        } else if ps[k] == -1 {
                            haplotype_allele_expression[3] += 1; // hap2 allele2, alternative allele
                        }
                    }
                }
                self.candidate_snps[*ti].germline = true;
                self.candidate_snps[*ti].rna_editing = false;
                if self.candidate_snps[*ti].alleles[0] != self.candidate_snps[*ti].reference && self.candidate_snps[*ti].alleles[1] != self.candidate_snps[*ti].reference {
                    // tri-allelic site
                    self.candidate_snps[*ti].variant_type = 3;
                    self.candidate_snps[*ti].hom_var = true;
                } else {
                    self.candidate_snps[*ti].variant_type = 1;
                }
                self.candidate_snps[*ti].haplotype = if phase_score1 >= phase_score2 { 1 } else { -1 };
                self.candidate_snps[*ti].haplotype_expression = haplotype_allele_expression;
                self.candidate_snps[*ti].phase_score = phase_score;
            }
        }
    }
    pub fn eval_som_var_phase(&mut self, somatic_allele_frac_cutoff: f32, somatic_allele_cnt_cutoff: u32) {
        for ti in self.somatic_snps.iter() {
            let snp = &mut self.candidate_snps[*ti];
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
                continue;
            }
            let mut ps: Vec<i32> = Vec::new();
            let mut assigns: Vec<i32> = Vec::new();
            let mut baseqs = Vec::new();
            for k in snp.snp_cover_fragments.iter() {
                if self.fragments[*k].assignment == 0 { continue; }
                if self.fragments[*k].num_hete_links < self.min_linkers { continue; }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == *ti {
                        assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                        ps.push(fe.p);
                        assigns.push(self.fragments[*k].assignment);
                        baseqs.push(fe.baseq);
                    }
                }
            }
            for k in 0..assigns.len() {
                if assigns[k] == 1 {
                    if ps[k] == 1 {
                        snp.hap_quals.hap1_ref_baseqs.push(baseqs[k]);
                    } else if ps[k] == -1 {
                        snp.hap_quals.hap1_alt_baseqs.push(baseqs[k]);
                    }
                } else if assigns[k] == 2 {
                    if ps[k] == 1 {
                        snp.hap_quals.hap2_ref_baseqs.push(baseqs[k]);
                    } else if ps[k] == -1 {
                        snp.hap_quals.hap2_alt_baseqs.push(baseqs[k]);
                    }
                }
            }
            let ref_allele_cnt = snp.hap_quals.hap1_ref_baseqs.len() + snp.hap_quals.hap2_ref_baseqs.len();
            let alt_allele_cnt = snp.hap_quals.hap1_alt_baseqs.len() + snp.hap_quals.hap2_alt_baseqs.len();
            if ref_allele_cnt * alt_allele_cnt > 0 && alt_allele_cnt as u32 >= somatic_allele_cnt_cutoff && alt_allele_cnt as f32 / (ref_allele_cnt + alt_allele_cnt) as f32 >= somatic_allele_frac_cutoff {
                // calculate somatic mutation probability
                let (hap1_allele_class, hap2_allele_class) = calculate_prob_somatic(&snp.hap_quals.hap1_ref_baseqs, &snp.hap_quals.hap1_alt_baseqs, &snp.hap_quals.hap2_ref_baseqs, &snp.hap_quals.hap2_alt_baseqs, 0.3);
                if hap1_allele_class.allcls == 0 && hap2_allele_class.allcls == 2 {
                    let somatic_score = -10.0_f64 * (1.0 - hap2_allele_class.prob).log10();
                    snp.cand_somatic = true;
                    snp.somatic = true;
                    snp.variant_type = 1;
                    snp.somatic_score = somatic_score;
                    snp.phase_score = 0.0;
                    // println!("somatic snp:{}, score: {}", snp.pos, somatic_score);
                    // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
                    // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", snp.hap_quals.hap1_ref_baseqs, snp.hap_quals.hap1_alt_baseqs, snp.hap_quals.hap2_ref_baseqs, snp.hap_quals.hap2_alt_baseqs);
                } else if hap1_allele_class.allcls == 2 && hap2_allele_class.allcls == 0 {
                    let somatic_score = -10.0_f64 * (1.0 - hap1_allele_class.prob).log10();
                    snp.cand_somatic = true;
                    snp.somatic = true;
                    snp.variant_type = 1;
                    snp.somatic_score = somatic_score;
                    snp.phase_score = 0.0;
                    // println!("somatic snp:{}, score: {}", snp.pos, somatic_score);
                    // println!("{:?},{:?}", hap1_allele_class, hap2_allele_class);
                    // println!("hap1_ref_baseqs:{:?}\nhap1_alt_baseqs:{:?}\nhap2_ref_baseqs:{:?}\nhap2_alt_baseqs:{:?}", snp.hap_quals.hap1_ref_baseqs, snp.hap_quals.hap1_alt_baseqs, snp.hap_quals.hap2_ref_baseqs, snp.hap_quals.hap2_alt_baseqs);
                }
            }
        }
    }
    pub fn eval_hom_var_phase(&mut self, min_phase_score: f32) {
        for ti in self.homo_snps.iter() {
            let snp = &mut self.candidate_snps[*ti];
            if snp.snp_cover_fragments.len() == 0 {
                // no surranding haplotype links
                snp.single = true;
                continue;
            }
            let mut sigma: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            let mut assigns: Vec<i32> = Vec::new();
            let mut baseqs = Vec::new();
            let mut hap1_reads_num = 0;
            let mut hap2_reads_num = 0;
            for k in snp.snp_cover_fragments.iter() {
                if self.fragments[*k].assignment == 0 { continue; }
                if self.fragments[*k].num_hete_links < self.min_linkers { continue; }
                for fe in self.fragments[*k].list.iter() {
                    if fe.snp_idx == *ti {
                        if fe.base != '-' {
                            // ignore intron
                            if self.fragments[*k].assignment == 1 {
                                hap1_reads_num += 1;
                            } else if self.fragments[*k].assignment == 2 {
                                hap2_reads_num += 1;
                            }
                        }
                        assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                        ps.push(fe.p);
                        probs.push(fe.prob);
                        sigma.push(self.fragments[*k].haplotag);
                        assigns.push(self.fragments[*k].assignment);
                        baseqs.push(fe.baseq);
                    }
                }
            }

            let mut phase_score1 = 0.0;
            let mut phase_score2 = 0.0;
            let mut phase_score = 0.0;
            if sigma.len() == 0 || hap1_reads_num < 2 || hap2_reads_num < 2 {
                // each haplotype should have at least 2 reads
                continue;
            } else {
                phase_score1 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
                phase_score2 = -10.0_f64 * (1.0 - SNPFrag::cal_delta_sigma_log(-1, &sigma, &ps, &probs)).log10(); // calaulate assignment score
            }

            if phase_score1.max(phase_score2) >= 3.0 * min_phase_score as f64 {
                // correct the genotype of this hom_var site, this site should be het_var
                phase_score = phase_score1.max(phase_score2);
                let mut haplotype_allele_expression: [u32; 4] = [0, 0, 0, 0];   // hap1_ref, hap1_alt, hap2_ref, hap2_alt
                for k in 0..sigma.len() {
                    if sigma[k] == 1 {
                        // hap1
                        if ps[k] == 1 {
                            haplotype_allele_expression[0] += 1; // hap1 allele1, reference allele
                        } else if ps[k] == -1 {
                            haplotype_allele_expression[1] += 1; // hap1 allele2, alternative allele
                        }
                    } else if sigma[k] == -1 {
                        // hap2
                        if ps[k] == 1 {
                            haplotype_allele_expression[2] += 1; // hap2 allele1, reference allele
                        } else if ps[k] == -1 {
                            haplotype_allele_expression[3] += 1; // hap2 allele2, alternative allele
                        }
                    }
                }
                self.candidate_snps[*ti].germline = true;
                if self.candidate_snps[*ti].variant_type != 3 {
                    // tri-allelic site will not be changed to het_var
                    self.candidate_snps[*ti].hom_var = false;
                    self.candidate_snps[*ti].variant_type = 1;
                }
                self.candidate_snps[*ti].haplotype = if phase_score1 >= phase_score2 { 1 } else { -1 };
                self.candidate_snps[*ti].haplotype_expression = haplotype_allele_expression;
                self.candidate_snps[*ti].phase_score = phase_score;
            }
        }
    }

    pub fn assign_reads_haplotype(&mut self, read_assignment_cutoff: f64) -> HashMap<String, i32> {
        let mut read_assignments: HashMap<String, i32> = HashMap::new();
        for k in 0..self.fragments.len() {
            let sigma_k = self.fragments[k].haplotag;
            let mut delta: Vec<i32> = Vec::new();
            let mut ps: Vec<i32> = Vec::new();
            let mut probs: Vec<f64> = Vec::new();
            for fe in self.fragments[k].list.iter() {
                if fe.phase_site == false { continue; }
                assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                ps.push(fe.p);
                probs.push(fe.prob);
                delta.push(self.candidate_snps[fe.snp_idx].haplotype);
            }
            if sigma_k == 0 {
                // unasigned haplotag, cluster the read into unknown group
                self.fragments[k].assignment = 0;
                self.fragments[k].assignment_score = 0.0;
                read_assignments.insert(self.fragments[k].read_id.clone(), 0);
            } else {
                let mut q = 0.0;
                let mut qn = 0.0;
                if delta.len() > 0 {
                    q = SNPFrag::cal_sigma_delta_log(sigma_k, &delta, &ps, &probs);
                    qn = SNPFrag::cal_sigma_delta_log(sigma_k * (-1), &delta, &ps, &probs);
                }

                if (q - qn).abs() >= read_assignment_cutoff {
                    if q >= qn {
                        if sigma_k == 1 {
                            self.fragments[k].assignment = 1;
                            self.fragments[k].assignment_score = q;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                        } else {
                            self.fragments[k].assignment = 2;
                            self.fragments[k].assignment_score = q;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                        }
                    } else {
                        if sigma_k == 1 {
                            self.fragments[k].assignment = 2;
                            self.fragments[k].assignment_score = qn;
                            self.fragments[k].haplotag = -1;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 2);
                        } else {
                            self.fragments[k].assignment = 1;
                            self.fragments[k].assignment_score = qn;
                            self.fragments[k].haplotag = 1;
                            read_assignments.insert(self.fragments[k].read_id.clone(), 1);
                        }
                    }
                } else {
                    // unknown which haplotype the read belongs to, cluster the read into unknown group
                    // panic!("Error: unexpected condition.");
                    self.fragments[k].assignment = 0;
                    self.fragments[k].assignment_score = 0.0;
                    read_assignments.insert(self.fragments[k].read_id.clone(), 0);
                }
            }
        }
        return read_assignments;
    }


    pub fn assign_phase_set(&mut self) -> HashMap<String, u32> {
        let mut phase_set: HashMap<String, u32> = HashMap::new();
        let mut graph: GraphMap<usize, Vec<usize>, Undirected> = GraphMap::new();  // node is index in candidate snp, edge is index in fragments
        // construct graph for hete snps
        for i in 0..self.candidate_snps.len() {
            let snp = &self.candidate_snps[i];
            if snp.variant_type == 0 || snp.variant_type == 2 || snp.variant_type == 3 {
                continue;
            }
            if snp.variant_type == 1 {
                // hete snps
                // let mut is_low_qual = false;
                // let mut is_dense = false;
                // let mut is_rna_edit = false;
                // let mut is_single_snp = false;
                // let mut is_unconfident_phased_snp = false;
                // let mut is_ase_snp = false;
                if snp.dense || snp.single || snp.rna_editing || snp.somatic || snp.phase_score == 0.0 {
                    continue;
                }
                // if snp.single == true {
                //     is_single_snp = true;
                // }
                // if snp.ase == true {
                //     is_ase_snp = true;
                // }
                // if !is_dense && !is_single_snp && snp.phase_score == 0.0 {
                //     is_unconfident_phased_snp = true;
                // }
                // if is_dense || is_single_snp || is_unconfident_phased_snp || snp.haplotype == 0 {
                //     continue;
                // }
                // // ase snp
                // if is_ase_snp && snp.phase_score < ase_ps_cutoff as f64 {
                //     continue;
                // }
                // // hete snp
                // if !is_ase_snp {
                //     if snp.variant_quality < min_qual_for_candidate as f64 {
                //         continue;
                //     }
                //     if snp.phase_score < min_phase_score as f64 {
                //         continue;
                //     }
                // }
                // ase snps > ase_ps_cutoff or hete snps > min_phase_score, construct graph
                graph.add_node(i);
            }
        }
        for k in 0..self.fragments.len() {
            let frag = &self.fragments[k];
            if frag.assignment == 0 { continue; }
            let mut node_snps = Vec::new();
            for fe in frag.list.iter() {
                if graph.contains_node(fe.snp_idx) {
                    node_snps.push(fe.snp_idx);
                }
            }
            if node_snps.len() >= 2 {
                for j in 0..node_snps.len() - 1 {
                    if !graph.contains_edge(node_snps[j], node_snps[j + 1]) {
                        graph.add_edge(node_snps[j], node_snps[j + 1], vec![k]);    // weight is a vector of fragment index, which is covered by the edge
                    } else {
                        graph.edge_weight_mut(node_snps[j], node_snps[j + 1]).unwrap().push(k);
                    }
                }
            }
        }
        let scc = kosaraju_scc(&graph);
        let region = self.region.clone().to_string();
        for component_nodes in scc.iter() {
            if component_nodes.len() <= 1 {
                continue;
            }
            let mut phase_id = 0;
            for node in component_nodes.iter() {
                if phase_id == 0 {
                    phase_id = (self.candidate_snps[*node].pos + 1) as u32;  // 1-based;
                }
                self.candidate_snps[*node].phase_set = phase_id;
                for edge in graph.edges(*node) {
                    let frag_idxes = edge.2;
                    for k in frag_idxes.iter() {
                        let fragment = &self.fragments[*k];
                        let read_id = fragment.read_id.clone();
                        if phase_set.contains_key(&read_id) {
                            continue;
                        }
                        phase_set.insert(read_id, phase_id);
                    }
                }
            }
        }
        return phase_set;
    }


    fn check_local_optimal_configuration(&self, used_for_haplotype: bool, used_for_haplotag: bool) {
        // check sigma
        if used_for_haplotag {
            for k in 0..self.fragments.len() {
                let sigma_k = self.fragments[k].haplotag;
                let mut delta: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                if sigma_k == 0 {
                    continue;
                }
                for fe in self.fragments[k].list.iter() {
                    if fe.phase_site == false { continue; }
                    assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                    ps.push(fe.p);
                    probs.push(fe.prob);
                    delta.push(self.candidate_snps[fe.snp_idx].haplotype);
                }
                if delta.len() == 0 {
                    continue;
                }
                let q = SNPFrag::cal_sigma_delta_log(sigma_k, &delta, &ps, &probs);
                let qn = SNPFrag::cal_sigma_delta_log(sigma_k * (-1), &delta, &ps, &probs);
                // println!("q:{}, qn:{}", q, qn);
                assert!(q >= qn, "{} Error: read assignment is not local optimal. {}->{}\n{:?}\ndelta:{:?}\nps:{:?}\nprobs:{:?}\nsigma:{}\n{:?}\n{:?}", k, q, qn, self.region, delta, ps, probs, sigma_k, used_for_haplotype, used_for_haplotag);
            }
        }

        // check delta
        if used_for_haplotype {
            for i in self.high_frac_het_snps.iter() {
                let delta_i = self.candidate_snps[*i].haplotype;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                    if self.fragments[*k].haplotag == 0 {
                        continue;
                    }
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == *i {
                            assert_ne!(fe.phase_site, false, "Error: phase for non-phase site.");
                            assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                            ps.push(fe.p);
                            probs.push(fe.prob);
                            sigma.push(self.fragments[*k].haplotag);
                        }
                    }
                }
                if sigma.len() == 0 {
                    continue;
                }
                let q = SNPFrag::cal_delta_sigma_log(delta_i, &sigma, &ps, &probs);
                let qn = SNPFrag::cal_delta_sigma_log(delta_i * (-1), &sigma, &ps, &probs);
                assert!(q >= qn, "{} Error: phase is not local optimal. {}->{}\n{:?}\nsigma:{:?}\nps:{:?}\nprobs:{:?}\ndelta:{}\n{:?}\n{:?}", i, q, qn, self.region, sigma, ps, probs, delta_i, used_for_haplotype, used_for_haplotag);
            }
        }
    }

    fn check_local_optimal_configuration_ase(
        &self,
        used_for_haplotype: bool,
        used_for_haplotag: bool,
    ) {
        // check sigma
        if used_for_haplotag {
            for k in 0..self.fragments.len() {
                let sigma_k = self.fragments[k].haplotag;
                let mut delta: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for fe in self.fragments[k].list.iter() {
                    ps.push(fe.p);
                    probs.push(fe.prob);
                    delta.push(self.candidate_snps[fe.snp_idx].haplotype);
                }
                if delta.len() == 0 {
                    continue;
                }
                let q = SNPFrag::cal_sigma_delta_log(sigma_k, &delta, &ps, &probs);
                let qn = SNPFrag::cal_sigma_delta_log(sigma_k * (-1), &delta, &ps, &probs);
                // println!("q:{}, qn:{}", q, qn);
                assert!(q >= qn, "{} Error: read assignment is not local optimal. {}->{}\n{:?}\ndelta:{:?}\nps:{:?}\nprobs:{:?}\nsigma:{}\n{:?}\n{:?}", k, q, qn, self.region, delta, ps, probs, sigma_k, used_for_haplotype, used_for_haplotag);
            }
        }

        // check delta
        if used_for_haplotype {
            for i in self.high_frac_het_snps.iter() {
                let delta_i = self.candidate_snps[*i].haplotype;
                let mut sigma: Vec<i32> = Vec::new();
                let mut ps: Vec<i32> = Vec::new();
                let mut probs: Vec<f64> = Vec::new();
                for k in self.candidate_snps[*i].snp_cover_fragments.iter() {
                    // k is fragment index
                    for fe in self.fragments[*k].list.iter() {
                        if fe.snp_idx == *i {
                            assert_ne!(fe.p, 0, "Error: phase for unexpected allele.");
                            ps.push(fe.p);
                            probs.push(fe.prob);
                            sigma.push(self.fragments[*k].haplotag);
                        }
                    }
                }
                if sigma.len() == 0 {
                    continue;
                }
                let q = SNPFrag::cal_delta_sigma_log(delta_i, &sigma, &ps, &probs);
                let qn = SNPFrag::cal_delta_sigma_log(delta_i * (-1), &sigma, &ps, &probs);
                assert!(q >= qn, "{} Error: phase is not local optimal. {}->{}\n{:?}\nsigma:{:?}\nps:{:?}\nprobs:{:?}\ndelta:{}\n{:?}\n{:?}", i, q, qn, self.region, sigma, ps, probs, delta_i, used_for_haplotype, used_for_haplotag);
            }
        }
    }
}