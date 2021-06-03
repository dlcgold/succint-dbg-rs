extern crate bio;

use std::collections::{HashSet, HashMap};
use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;
#[allow(unused_imports)]
use std::fs::File;
#[allow(unused_imports)]
use std::io::Write;
use math::round;

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
    kmersize: u32,
    n_nodes: usize,
    node_char: Vec<char>,
    last: Vec<usize>,
    out: Vec<char>,
    fvec: Vec<usize>,
    neg: Vec<bool>,
    rschar: HashMap<char, RankSelect>,
    rscharneg: HashMap<char, RankSelect>,
    rslast: RankSelect,
}


impl SDbg {
    pub fn new(reads: &mut Vec<String>, k: u32) -> Self {
        let mut node_out: Vec<(String, char, String)> = Vec::new();
        let kmers = get_kmers_succ(reads, k);
        let mut set_check = HashSet::new();
        for kmer in &kmers {
            let kmer1 = kmer[0..(k - 1) as usize].to_string();
            node_out.push((kmer[0..(k - 1) as usize].to_string(),
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
        let mut dollars = (0, 0);
        let mut dollbwt_pos = Vec::new();
        let mut doll_pos = Vec::new();
        let mut node_char = Vec::new();
        let mut fvec = vec![0; 256];
        for i in 0..real_order.len() {
            if i == 0 {
                fvec['$' as usize] = 0;
            } else if nodes[i - 1].chars().last().unwrap() != nodes[i].chars().last().unwrap() {
                fvec[nodes[i].chars().last().unwrap() as usize] = i;
            }

            if !check_neg.contains(&(nodes[i][1..].to_string(), out[i])) {
                check_neg.insert((nodes[i][1..].to_string(), out[i]));
            } else {
                neg[i] = true;
                neg_pos.push(i);
            }
            if nodes[i].chars().last().unwrap() == '$' {
                dollars.0 += 1;
                doll_pos.push(i);
            }
            if out[i] == '$' {
                dollars.1 += 1;
                dollbwt_pos.push(i);
            }
            node_char.push(nodes[i].chars().last().unwrap());
        }
        let mut bitvecs = HashMap::new();
        let mut bitvecsneg = HashMap::new();
        for s in "$acgtACGT".chars() {
            let mut bv = BitVec::new();
            let mut bvneg = BitVec::new();
            for (index, c) in out.iter().enumerate() {
                if c == &s && !neg[index] {
                    bv.push(true);
                } else {
                    bv.push(false);
                }
            }
            let krank = round::ceil(((bv.len() as f64).log2()).powf(2.) / (32 as f64), 0);
            let rank_select = RankSelect::new(bv.clone(), krank as usize);
            bitvecs.insert(s, rank_select);
            for (indexn, c) in out.iter().enumerate() {
                if c == &s && neg[indexn] {
                    bvneg.push(true);
                } else {
                    bvneg.push(false);
                }
            }
            let krankneg = round::ceil(((bvneg.len() as f64).log2()).powf(2.) / (32 as f64), 0);
            let rank_selectneg = RankSelect::new(bvneg.clone(), krankneg as usize);
            bitvecsneg.insert(s, rank_selectneg);
        }
        let mut bv = BitVec::new();

        for e in &last {
            match e {
                1 => bv.push(true),
                _ => bv.push(false),
            }
        }
        //println!("{:?}", bv);
        let krank = round::ceil(((bv.len() as f64).log2()).powf(2.) / (32 as f64), 0);
        let rslast = RankSelect::new(bv.clone(), krank as usize);
        let n_nodes = rslast.rank(out.len() as u64 - 1).unwrap() as usize;
        SDbg {
            kmersize: k,
            n_nodes,
            node_char,
            last,
            out,
            neg,
            rschar: bitvecs,
            rscharneg: bitvecsneg,
            fvec,
            rslast,
        }
    }

    pub fn kmersize(&self) -> u32 { self.kmersize }
    pub fn n_nodes(&self) -> usize { self.n_nodes }
    pub fn node_char(&self) -> &Vec<char> { &self.node_char }
    pub fn last(&self) -> &Vec<usize> { &self.last }
    pub fn out(&self) -> &Vec<char> { &self.out }
    pub fn neg(&self) -> &Vec<bool> { &self.neg }
    pub fn fvec(&self) -> &Vec<usize> { &self.fvec }
    pub fn rschar(&self) -> &HashMap<char, RankSelect> { &self.rschar }
    pub fn rscharneg(&self) -> &HashMap<char, RankSelect> { &self.rscharneg }
    pub fn rslast(&self) -> &RankSelect { &self.rslast }

    pub fn print(&self) {
        for i in 0..self.out().len() {
            println!("{}|{}|{} ({})", self.last[i], self.node_char[i], self.out[i], self.neg[i]);
        }
    }

    pub fn findsymbol(&self, i: isize) -> char {
        if i == 0 {
            return '$';
        }
        let mut symbol = '$';
        let mut lastcharnn = '$';
        let alphabet = "$ACGTacgt";
        for (index, elem) in alphabet.chars().enumerate() {
            if *&self.fvec[elem as usize] > i as usize {
                symbol = alphabet.chars().nth(index - 1).unwrap();
                break;
            }
            if *&self.fvec[elem as usize] != 0 {
                lastcharnn = elem;
            }
        }
        if symbol == '$' {
            symbol = lastcharnn;
        }
        symbol
    }

    pub fn first_edge(&self, i: isize) -> isize {
        let select;
        if i <= 0 {
            select = -1;
        } else {
            select = match &self.rslast.select(i as u64) {
                Some(t) => *t,
                None => 0,
            } as isize;
        }
        select as isize + 1
    }

    pub fn last_edge(&self, i: isize) -> isize {
        let select;
        if i + 1 <= 0 {
            select = -1;
        } else {
            select = match &self.rslast.select(i as u64 + 1) {
                Some(t) => *t,
                None => 0,
            } as isize;
        }
        select as isize
    }

    pub fn node_range(&self, i: isize) -> (isize, isize) {
        (self.first_edge(i), self.last_edge(i))
    }

    pub fn edge_to_node(&self, i: isize) -> isize {
        if i == 0 {
            return 0;
        }
        let rank = match self.rslast().rank(i as u64 - 1) {
            Some(t) => t,
            None => 0,
        };
        rank as isize
    }
    pub fn node_to_edge(&self, i: isize) -> isize {
        if i <= 0 {
            return 0;
        }
        let select = match self.rslast().select(i as u64) {
            Some(t) => t,
            None => 0,
        } + 1;
        select as isize
    }

    pub fn forward(&self, i: isize) -> isize {
        let symbol = self.out[i as usize];
        if symbol == '$' {
            return -1;
        }
        let relindex = match *&self.rschar()[&symbol].rank(i as u64) {
            Some(t) => t,
            None => 0,
        };
        let firstocc = self.fvec[symbol as usize];
        let ranklast = match *&self.rslast().rank(firstocc as u64 - 1) {
            Some(t) => t,
            None => 0,
        };
        let select;
        if ranklast + relindex <= 0 {
            select = -1;
        } else {
            select = match &self.rslast.select(ranklast + relindex) {
                Some(t) => *t,
                None => 0,
            } as isize;
        }
        select
    }

    pub fn backward(&self, i: isize) -> isize {
        let symbol = self.findsymbol(i);
        if symbol == '$' {
            return -1;
        }
        let firstocc = self.fvec[symbol as usize];
        let ranksymb = match *&self.rslast().rank(i as u64 - 1) {
            Some(t) => t,
            None => 0,
        } as isize;
        let ranklast = match *&self.rslast().rank(firstocc as u64 - 1) {
            Some(t) => t,
            None => 0,
        } as isize;
        let select;
        if ranksymb - ranklast + 1 <= 0 {
            select = -1;
        } else {
            select = match &self.rschar()[&symbol].select((ranksymb - ranklast + 1) as u64) {
                Some(t) => *t,
                None => 0,
            } as isize;
        }
        select as isize
    }

    pub fn outdegree(&self, i: isize) -> isize {
        if self.out()[i as usize] == '$' {
            return 0;
        }
        let range = self.node_range(i);
        if range.1 > range.0 {
            return range.1 - range.0 + 1;
        } else {
            return 1;
        }
        //range.1 - range.0 + 1
    }

    pub fn indegree(&self, i: isize) -> isize {
        let last = self.last_edge(i);
        let pred = self.backward(last);
        if pred == -1 {
            return 0;
        }
        let symbol = self.out()[pred as usize];
        let next;
        if pred + 1 <= 0 {
            next = -1
        } else {
            next = match &self.rschar()[&symbol].select(pred as u64 + 1) {
                Some(t) => *t,
                None => self.out().len() as u64 - 1,
            } as isize;
        }
        let negnextrank = match *&self.rscharneg[&symbol].rank(next as u64) {
            Some(t) => t,
            None => 0,
        } as isize;
        let negpredrank = match *&self.rscharneg[&symbol].rank(pred as u64) {
            Some(t) => t,
            None => 0,
        } as isize;
        negnextrank - negpredrank + 1
    }

    #[allow(unused_variables)]
    pub fn outgoing(&self, i: isize, symbol: char) -> isize {
        if symbol == '$' {
            return -1;
        }
        let range = self.node_range(i);
        let relindex = match *&self.rschar()[&symbol].rank(range.1 as u64) {
            Some(t) => t,
            None => 0,
        };
        let selectnode;
        if relindex <= 0 {
            selectnode = -1;
        } else {
            selectnode = match &self.rschar()[&symbol].select(relindex) {
                Some(t) => *t,
                None => 0,
            } as isize;
        }

        let relindexneg = match *&self.rscharneg()[&symbol].rank(range.1 as u64) {
            Some(t) => t,
            None => 0,
        };
        let selectnodeneg;
        if relindexneg <= 0 {
            selectnodeneg = -1;
        } else {
            selectnodeneg = match &self.rscharneg()[&symbol].select(relindexneg) {
                Some(t) => *t,
                None => 0,
            } as isize;
        }
        return if range.0 <= selectnode && selectnode <= range.1 {
            (match *&self.rslast.rank(self.forward(selectnode) as u64) {
                Some(t) => t,
                None => 0,
            } - 1) as isize
        } else if range.0 <= selectnodeneg && selectnodeneg <= range.1 {
            (match *&self.rslast.rank(self.forward(selectnodeneg) as u64) {
                Some(t) => t,
                None => 0,
            } - 1) as isize
        } else {
            -1
        };
    }

    pub fn successors(&self, i: isize) -> Vec<isize> {
        let mut succs = Vec::new();
        let range = self.node_range(i);
        for i in range.0..range.1 + 1 {
            succs.push(match self.rslast.rank(self.forward(i) as u64) {
                Some(t) => t,
                None => 0,
            } as isize - 1);
        }
        succs
    }

    pub fn label(&self, i: isize) -> String {
        if i == 0 {
            let mut countd = self.kmersize - 1;
            let mut dollars = "".to_string();
            while countd != 0 {
                dollars.push('$');
                countd -= 1;
            }
            return dollars;
        }
        let mut index = self.first_edge(i);
        //println!("{},{}", i, index);
        let mut symbol = self.findsymbol(index);
        let mut label = symbol.to_string();
        let mut count = self.kmersize - 2;

        while count != 0 {
            let new_index = self.backward(index);
            symbol = self.findsymbol(new_index);
            label.push(symbol);
            count -= 1;
            if symbol == '$' {
                while count != 0 {
                    label.push('$');
                    count -= 1;
                }
                break;
            }
            index = new_index;
        }
        let label = label.chars().rev().collect::<String>();
        label
    }

    pub fn to_dot(&self, output: &str) {
        let mut fileout = File::create(output).expect("error");
        fileout
            .write("digraph sample{\n".as_bytes())
            .expect("error");
        for i in 0..self.n_nodes() {
            for j in self.successors(i as isize) {
                if j != -1 {
                    let start = self.label(i as isize);
                    let end = self.label(j as isize);

                    fileout.write(
                        format!(
                            "\t\"{}\" -> \"{}\" [ label = \"{}\" ];\n",
                            start,
                            end,
                            end.chars().last().unwrap()
                        )
                            .as_bytes(),
                    )
                        .expect("error");
                }
            }
        }

        fileout.write("}".as_bytes()).expect("error");
    }
}


#[cfg(test)]
mod tests {
    use crate::SDbg;
    #[allow(unused_imports)]
    use std::path::Path;

    /* #[test]
     fn test_sdbg() {
         let mut kmers = vec!["TACACT".to_string(),
                              "TACTCA".to_string(),
                              "GACTCG".to_string()];
         let sdbg = SDbg::new(&mut kmers, 4);
         assert_eq!(sdbg.nodes.len(), 16);
     }*/

    #[test]
    fn test_sdbg_forward() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut forw = Vec::new();
        let forcheck = vec![10, 4, 5, 8, 11, 8, 9, 1, 12, 1, 2, -1, 6];
        for i in 0..sdbg.out().len() {
            forw.push(sdbg.forward(i as isize));
        }
        assert_eq!(forw, forcheck);
    }

    #[test]
    fn test_sdbg_backward() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut backw = Vec::new();
        let backcheck = vec![-1, 7, 10, 1, 1, 2, 12, 3, 3, 6, 0, 4, 8];
        for i in 0..sdbg.out().len() {
            backw.push(sdbg.backward(i as isize));
        }
        assert_eq!(backw, backcheck);
    }

    #[test]
    fn test_sdbg_outdegree() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut outd = Vec::new();
        let outcheck = vec![1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1];
        for i in 0..sdbg.n_nodes() {
            outd.push(sdbg.outdegree(i as isize));
        }
        assert_eq!(outd, outcheck);
    }

    #[test]
    fn test_sdbg_indegree() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut ind = Vec::new();
        let incheck = vec![0, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1];
        for i in 0..sdbg.n_nodes() {
            ind.push(sdbg.indegree(i as isize));
        }
        assert_eq!(ind, incheck);
    }

    #[test]
    fn test_sdbg_outgoing() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut outg = Vec::new();
        let outcheck = vec![-1, -1, -1, -1, 8, -1, -1, 3, -1, -1, -1, -1, 4, -1, -1, -1,
                            -1, -1, 6, 9, -1, -1, -1, 6, -1, -1, -1, -1, 7, -1, -1, 1, -1, -1, 10,
                            -1, 1, -1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 5,
                            -1, -1];
        for i in 0..sdbg.n_nodes() {
            for j in ['$', 'A', 'C', 'G', 'T'] {
                outg.push(sdbg.outgoing(i as isize, j));
            }
        }
        assert_eq!(outg, outcheck);
    }

    #[test]
    fn test_sdbg_label() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut label = Vec::new();
        let labelcheck = vec!["$$$", "CGA", "$TA", "GAC", "TAC", "GTC", "ACG", "TCG",
                              "$$T", "ACT", "CGT", "$$T", "CGA"];
        for i in 0..sdbg.out().len() {
            label.push(sdbg.label(i as isize));
        }
        assert_eq!(label, labelcheck);
    }

    #[test]
    fn test_sdbg_successors() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut outs = Vec::new();
        let succcheck = vec![vec![8], vec![3], vec![4], vec![6, 9], vec![6], vec![7],
                             vec![1, 10], vec![1], vec![2], vec![-1], vec![5]];
        for i in 0..sdbg.n_nodes() {
            outs.push(sdbg.successors(i as isize));
        }
        assert_eq!(outs, succcheck);
    }

    #[test]
    fn test_sdbg_to_dot() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        sdbg.to_dot("output/bowe.dot");
        assert!(Path::new("output/bowe.dot").exists());
    }

    #[test]
    fn test_dot_slide() {
        let mut kmers = vec!["TACACT".to_string(),
                             "TACTCA".to_string(),
                             "GACTCG".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        sdbg.to_dot("output/slide.dot");
        assert!(Path::new("output/slide.dot").exists());
    }

    #[test]
    fn test_dot_assign2() {
        let mut kmers = vec!["ATCTTGCATTACCGCCCCAATC".to_string(),
                             "ATCTTACATTACCGTCCCAACC".to_string()];
        let sdbg = SDbg::new(&mut kmers, 6);
        sdbg.to_dot("output/assign2.dot");
        assert!(Path::new("output/assign2.dot").exists());
    }
}
