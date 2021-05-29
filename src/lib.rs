extern crate bio;

use std::collections::HashSet;
use bio::data_structures::bwt::{less, Occ, Less, BWT};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::alphabets;
use std::fs::File;
use std::io::Write;
use std::fmt::Pointer;
use std::ops::Index;

#[allow(dead_code)]
fn create_kmers(s: String, k: u32) -> Vec<String> {
    get_kmers(&mut vec![s], k)
}

#[allow(dead_code)]
fn get_kmers(reads: &mut Vec<String>, k: u32) -> Vec<String> {
    let mut v: Vec<String> = Vec::new();
    for s in reads {
        while s.len() >= k as usize {
            v.push(s[0..k as usize].to_string());
            *s = String::from(&s[1..]);
        }
    }
    v
}

#[allow(dead_code)]
fn create_kmers_unique(s: String, k: u32) -> Vec<String> {
    get_kmers_unique(&mut vec![s], k)
}

#[allow(dead_code)]
pub fn get_kmers_unique(reads: &mut Vec<String>, k: u32) -> Vec<String> {
    let mut v: Vec<String> = Vec::new();
    for s in reads {
        while s.len() >= k as usize {
            if !v.contains(&s[0..k as usize].to_string()) {
                v.push(s[0..k as usize].to_string());
            }
            *s = String::from(&s[1..]);
        }
    }
    v
}

#[allow(dead_code)]
fn create_kmers_succ(s: String, k: u32) -> Vec<String> {
    get_kmers_succ(&mut vec![s], k)
}

#[allow(dead_code)]
pub fn get_kmers_succ(reads: &mut Vec<String>, k: u32) -> Vec<String> {
    let mut v: Vec<String> = Vec::new();
    let mut dollars = String::new();
    for _i in 0..k - 1 {
        dollars += "$";
    }
    for s in &mut *reads {
        *s = format!("{}{}", dollars, s);
    }
    for s in reads {
        while s.len() >= k as usize {
            if !v.contains(&s[0..k as usize].to_string()) {
                v.push(s[0..k as usize].to_string());
            }
            *s = String::from(&s[1..]);
        }
    }
    v
}


pub struct SDbg {
    nodes: Vec<String>,
    last: Vec<usize>,
    out: Vec<char>,
    neg: Vec<bool>,
    neg_pos: Vec<usize>,
    nodes_rev: Vec<String>,
    fm: FMIndex<BWT, Less, Occ>,
}

impl SDbg {
    pub fn new(reads: &mut Vec<String>, k: u32) -> Self {
        let mut node_out: Vec<(String, char, String)> = Vec::new();
        let kmers = get_kmers_succ(reads, k);
        let mut set_check = HashSet::new();
        for kmer in &kmers {
            let kmer1 = kmer[0..3].to_string();
            node_out.push((kmer[0..3].to_string(),
                           kmer.chars().last().unwrap(),
                           kmer1.chars().rev().collect()));
            set_check.insert(kmer1);
        }
        for kmer in kmers {
            let kmer2 = kmer[1..].to_string();
            if !set_check.contains(&kmer2) {
                node_out.push((kmer[1..].to_string(),
                               '$',
                               kmer2.chars().rev().collect()));
            }
        }

        node_out.sort_by(|a, b| (a.2.cmp(&b.2)));
        let mut real_order = Vec::new();
        let mut tmp = Vec::new();
        for i in 0..node_out.len() {
            if i != node_out.len() - 1 && &node_out[i].0 == &node_out[i + 1].0 {
                tmp.push(node_out[i].clone());
            } else {
                tmp.push(node_out[i].clone());
                if tmp.len() != 1 {
                    tmp.sort_by(|a, b| (a.1.cmp(&b.1)));
                }
                for e in tmp.clone() {
                    real_order.push(e);
                }
                tmp.clear();
            }
        }
        let mut last = Vec::new();
        let mut nodes = Vec::new();
        let mut out = Vec::new();
        let mut nodes_rev = Vec::new();
        let mut neg = vec![false; real_order.len()];
        let mut neg_pos = Vec::new();
        for i in 0..real_order.len() - 1 {
            if &real_order[i].0 == &real_order[i + 1].0 {
                last.push(0);
            } else {
                last.push(1);
            }
            nodes.push(real_order[i].0.clone());
            out.push(real_order[i].1);
            nodes_rev.push(real_order[i].2.clone());
        }
        last.push(1);
        nodes.push(real_order[real_order.len() - 1].0.clone());
        out.push(real_order[real_order.len() - 1].1);
        nodes_rev.push(real_order[real_order.len() - 1].2.clone());
        let mut check_neg: HashSet<(String, char)> = HashSet::new();
        for i in 0..real_order.len() {
            if !check_neg.contains(&(nodes[i][1..].to_string(), out[i])) {
                check_neg.insert((nodes[i][1..].to_string(), out[i]));
            } else {
                neg[i] = true;
                neg_pos.push(i);
            }
        }
        let dna_alphabet = alphabets::Alphabet::new(b"$acgtACGT");
        let outb = out.iter().map(|c| *c as u8).collect::<Vec<_>>();
        let bwt: &[u8] = &outb[..];
        let c_fun = less(bwt, &dna_alphabet);
        let occ_fun = Occ::new(bwt, k as u32, &dna_alphabet);
        let fm = FMIndex::new(outb, c_fun, occ_fun);
        SDbg {
            nodes,
            last,
            out,
            nodes_rev,
            neg,
            fm,
            neg_pos,
        }
    }

    pub fn nodes(&self) -> &Vec<String> {
        &self.nodes
    }
    pub fn last(&self) -> &Vec<usize> {
        &self.last
    }
    pub fn out(&self) -> &Vec<char> {
        &self.out
    }
    pub fn nodes_rev(&self) -> &Vec<String> {
        &self.nodes_rev
    }
    pub fn fm(&self) -> &FMIndex<BWT, Less, Occ> {
        &self.fm
    }
    pub fn neg(&self) -> &Vec<bool> { &self.neg }
    pub fn neg_pos(&self) -> &Vec<usize> { &self.neg_pos }

    pub fn print(&self) {
        for i in 0..self.nodes.len() {
            println!("{}|{}|{} ({},{})", self.last[i], self.nodes[i], self.out[i],
                     self.nodes_rev[i], self.neg[i]);
        }
    }

    pub fn lf_function(&self, mut index: usize) -> isize {
        let symbol = self.fm().bwt()[index];
        let mut indexneg = index;
        let mut check = true;
        if self.neg()[indexneg] {
            while check {
                indexneg -= 1;
                if !self.neg()[indexneg] && self.fm().bwt()[indexneg] == symbol {
                    check = false;
                }
            }
        }
        println!("indexes: {}->{}", index, indexneg);
        let mut j = 0;
        let mut jump = 0;
        if index == 0 {
            j = (self.fm().less(self.fm().bwt()[index])) as isize;
        } else {
            j = (self.fm().less(self.fm().bwt()[indexneg]) +
                (self.fm().occ(indexneg, self.fm().bwt()[indexneg])) - 1) as isize
        }
        println!("pre jump: {}", j);
        if j > self.neg_pos()[0] as isize {
            for ind in 0..self.neg_pos().len() {
                if j <= self.neg_pos()[ind] as isize {
                    break;
                }
                jump += 1;
            }
        }
        println!("post jump: {}", j + jump);
        j + jump as isize
    }
    pub fn to_dot(&self, output: &str) {
        let mut fileout = File::create(output).expect("error");
        let mut nodes_tmp = self.nodes().clone();
        println!("{:?}", self.neg_pos());
        fileout
            .write("digraph sample{\n".as_bytes())
            .expect("error");
        let mut start = nodes_tmp[0].clone();
        let mut i = 0;
        let mut visited = vec![false; self.nodes.len()];
        //visited[0] = true;
        let mut not_visited = self.nodes.len();

        /*for i in 0..self.nodes.len(){
            print!("{}: ",i);
            for s in ['$', 'A', 'C', 'G', 'T']{
                print!("{} ",self.fm().occ(i, s as u8));
            }
            println!();
        }*/
        /*for s in ['$', 'A', 'C', 'G', 'T']{
            print!("{} {}\n",s, self.fm().less( s as u8));
        }*/
        while not_visited != 0 {
            println!("{} at {}", start, i);
            let mut j = self.lf_function(i);
            //println!("{},{},{},{}", start, i, j, &self.nodes()[(j) as usize]);

            if !visited[i as usize] && self.fm().bwt()[i] as char != '$' {
                fileout.write(
                    format!(
                        "\t\"{}\" -> \"{}\" [ label = \"{}\" ];\n",
                        start,
                        &self.nodes()[(j) as usize],
                        self.fm().bwt()[i] as char
                    )
                        .as_bytes(),
                )
                    .expect("error");
                visited[i as usize] = true;
                start = self.nodes()[(j) as usize].clone();
                not_visited -= 1;
                i = (j) as usize;
            } else {
                let mut o: isize = -1;
                for (ind, elem) in visited.clone().iter().enumerate() {
                    if !*elem && self.fm().bwt()[ind] as char == '$' {
                        visited[ind] = true;
                        not_visited -= 1;
                    } else if !*elem {
                        o = ind as isize;
                        break;
                    }
                }
                if o == -1 {
                    break;
                } else {
                    i = o as usize;
                }
                /*visited[i] = true;
                not_visited -= 1;*/
                start = self.nodes()[i as usize].clone();
            }
        }

        fileout.write("}".as_bytes()).expect("error");
    }
}

#[cfg(test)]
mod tests {
    use crate::SDbg;
    use bio::data_structures::fmindex::FMIndexable;

    #[test]
    fn test_sdbg() {
        let mut kmers = vec!["TACACT".to_string(),
                             "TACTCA".to_string(),
                             "GACTCG".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        //sdbg.print();
        assert_eq!(sdbg.nodes.len(), 16);
    }

    #[test]
    fn test_lf() {
        let mut kmers = vec!["TACACT".to_string(),
                             "TACTCA".to_string(),
                             "GACTCG".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        sdbg.to_dot("output/test.dot");
    }
}
