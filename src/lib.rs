extern crate bio;

use std::collections::{HashSet, HashMap};
use bio::data_structures::bwt::{less, Occ, Less, BWT};
#[allow(unused_imports)]
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::data_structures::rank_select::RankSelect;
use bio::alphabets;
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
    nodes: Vec<String>,
    node_char: Vec<char>,
    last: Vec<usize>,
    out: Vec<char>,
    fvec: Vec<usize>,
    neg: Vec<bool>,
    neg_pos: Vec<usize>,
    nodes_rev: Vec<String>,
    dollars: (usize, usize),
    dollbwt_pos: Vec<usize>,
    doll_pos: Vec<usize>,
    rschar: HashMap<char, RankSelect>,
    rslast: RankSelect,
    fm: FMIndex<BWT, Less, Occ>,
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
        let dna_alphabet = alphabets::Alphabet::new(b"$acgtACGT");
        let outb = out.iter().map(|c| *c as u8).collect::<Vec<_>>();
        let bwt: &[u8] = &outb[..];
        let c_fun = less(bwt, &dna_alphabet);
        let occ_fun = Occ::new(bwt, k as u32, &dna_alphabet);
        let fm = FMIndex::new(outb, c_fun, occ_fun);
        let mut bitvecs = HashMap::new();
        for s in "$acgtACGT".chars() {
            let mut bv = BitVec::new();
            for c in &out {
                if c == &s {
                    bv.push(true);
                } else {
                    bv.push(false);
                }
            }
            let krank = round::ceil(((bv.len() as f64).log2()).powf(2.) / (32 as f64), 0);
            let rank_select = RankSelect::new(bv.clone(), krank as usize);
            bitvecs.insert(s, rank_select);
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

        SDbg {
            kmersize: k,
            nodes,
            node_char,
            last,
            out,
            nodes_rev,
            neg,
            fm,
            dollars,
            dollbwt_pos,
            doll_pos,
            rschar: bitvecs,
            neg_pos,
            fvec,
            rslast,
        }
    }

    pub fn kmersize(&self) -> u32 { self.kmersize }
    pub fn nodes(&self) -> &Vec<String> { &self.nodes }
    pub fn node_char(&self) -> &Vec<char> { &self.node_char }
    pub fn last(&self) -> &Vec<usize> { &self.last }
    pub fn out(&self) -> &Vec<char> { &self.out }
    pub fn nodes_rev(&self) -> &Vec<String> { &self.nodes_rev }
    pub fn fm(&self) -> &FMIndex<BWT, Less, Occ> { &self.fm }
    pub fn neg(&self) -> &Vec<bool> { &self.neg }
    pub fn dollars(&self) -> (usize, usize) { self.dollars }
    pub fn dollbwt_pos(&self) -> &Vec<usize> { &self.dollbwt_pos }
    pub fn doll_pos(&self) -> &Vec<usize> { &self.doll_pos }
    pub fn neg_pos(&self) -> &Vec<usize> { &self.neg_pos }
    pub fn fvec(&self) -> &Vec<usize> { &self.fvec }
    pub fn rschar(&self) -> &HashMap<char, RankSelect> { &self.rschar }
    pub fn rslast(&self) -> &RankSelect { &self.rslast }

    pub fn print(&self) {
        for i in 0..self.nodes.len() {
            println!("{}|{}|{} ({},{})", self.last[i], self.nodes[i], self.out[i],
                     self.nodes_rev[i], self.neg[i]);
        }
    }

    /*pub fn lf_function(&self, index: usize) -> isize {
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
        let mut jump = 0;
        println!("indexneg: {}", indexneg);
        let mut j;
        if index == 0 {
            j = (self.fm().less(self.fm().bwt()[indexneg])) as isize;
        } else {
            j = (self.fm().less(self.fm().bwt()[indexneg]) +
                (self.fm().occ(indexneg, self.fm().bwt()[indexneg])) - 1) as isize;
        }

        if j > self.neg_pos()[0] as isize {
            for ind in 0..self.neg_pos().len() {
                if j <= self.neg_pos()[ind] as isize {
                    break;
                }
                jump += 1;
            }
        }
        println!("jump: {}", jump);
        if index == 0 {
            j -= (self.dollars.1 - self.dollars.0) as isize;
        } else if self.dollars.1 != self.dollars.0 && j >= self.dollbwt_pos[0] as isize {
            let mut count_doll = 0;
            for ind in 0..self.dollbwt_pos().len() {
                if j < self.dollbwt_pos()[ind] as isize {
                    break;
                }
                count_doll += 1;
            }
            j -= count_doll;
        }
        j + jump as isize
    }*/

    /*pub fn to_dot(&self, output: &str) {
        let mut fileout = File::create(output).expect("error");
        let nodes_tmp = self.nodes().clone();
        fileout
            .write("digraph sample{\n".as_bytes())
            .expect("error");
        let mut start = nodes_tmp[0].clone();
        let mut i = 0;
        let mut visited = vec![false; self.nodes.len()];
        let mut not_visited = self.nodes.len();
        println!("{:?}", self.dollbwt_pos());
        for s in ['$', 'A', 'C', 'G', 'T'] {
            println!("C({}) = {}", s, self.fm().less(s as u8));
        }
        while not_visited != 0 {
            let j = self.lf_function(i);
            println!("{} -> {}", i, j);
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
                start = self.nodes()[i as usize].clone();
            }
        }

        fileout.write("}".as_bytes()).expect("error");
    }*/

    pub fn findsymbol(&self, i: usize) -> char {
        if i == 0 {
            return '$';
        }
        let mut symbol = '$';
        let mut lastcharnn = '$';
        let alphabet = "$ACGTacgt";
        for (index, elem) in alphabet.chars().enumerate() {
            if *&self.fvec[elem as usize] > i {
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

    pub fn first_edge(&self, i: usize) -> usize {
        let select = match &self.rslast.select(i as u64) {
            Some(t) => *t,
            None => 0 as u64,
        } + 1;
        select as usize
    }

    pub fn last_edge(&self, i: usize) -> usize {
        let select = match &self.rslast.select(i as u64 + 1) {
            Some(t) => *t,
            None => 0 as u64,
        };
        select as usize
    }

    pub fn node_range(&self, i: usize) -> (usize, usize) {
        (self.first_edge(i), self.last_edge(i))
    }

    pub fn edge_to_node(&self, i: usize) -> usize {
        if i == 0 {
            return 0;
        }
        let rank = self.rslast().rank(i as u64 + 1).unwrap();
        rank as usize
    }

    pub fn forward(&self, i: usize) -> usize {
        let symbol = self.out[i];
        if symbol == '$' {
            return self.out().len();
        }
        let relindex = *&self.rschar()[&symbol].rank(i as u64).unwrap();
        let firstocc = self.fvec[symbol as usize];
        let ranklast = *&self.rslast().rank(firstocc as u64 - 1).unwrap();
        let select = match &self.rslast.select(ranklast + relindex) {
            Some(t) => *t,
            None => 0 as u64,
        };
        select as usize
    }

    pub fn backward(&self, i: usize) -> usize {
        let symbol = self.findsymbol(i);
        if symbol == '$' {
            return self.out().len();
        }
        let firstocc = self.fvec[symbol as usize];
        let ranksymb = *&self.rslast().rank(i as u64 - 1).unwrap();
        let ranklast = *&self.rslast().rank(firstocc as u64 - 1).unwrap();
        let select = match &self.rschar()[&symbol].select(ranksymb - ranklast + 1) {
            Some(t) => *t,
            None => 0 as u64,
        };
        select as usize
    }

    pub fn outdegree(&self, i: usize) -> usize {
        let range = self.node_range(i);
        if range.1 > range.0 {
            return range.1 - range.0 + 1;
        } else {
            return 1;
        }
    }

    pub fn indegree(&self, i: usize) -> usize {
        let last = self.last_edge(i);
        let pred = self.backward(last);
        if pred == self.out().len() {
            return 0;
        }
        let symbol = self.node_char()[pred];
        let next = match &self.rschar()[&symbol].select(pred as u64 + 1) {
            Some(t) => *t,
            None => self.node_char().len() as u64,
        } as usize;
        let mut negpredrank = 0;
        let mut negnextrank = 0;
        for ind in (pred - 1)..0 {
            if self.out()[ind] == symbol && self.neg()[ind] == true {
                negpredrank += 1;
            }
        }
        for ind in (next - 1)..0 {
            if self.out()[ind] == symbol && self.neg()[ind] == true {
                negnextrank += 1;
            }
        }
        negnextrank - negpredrank + 1
    }

    #[allow(unused_variables)]
    pub fn outgoing(&self, i: usize, symbol: char) -> usize {
        if symbol == '$' {
            return self.out().len();
        }
        let range = self.node_range(i);
        let mut negrank = 0;
        let mut negselect = i;
        let mut check = true;
        for ind in (i - 1)..0 {
            if self.out()[ind] == symbol && self.neg()[ind] == true {
                if check {
                    negselect = ind;
                    check = false;
                }
                negrank += 1;
            }
        }
        let relindex = *&self.rschar()[&symbol].rank(range.1 as u64).unwrap();
        let selectnode = match &self.rschar()[&symbol].select(relindex) {
            Some(t) => *t,
            None => 0 as u64,
        } as usize;
        return if range.0 <= selectnode && selectnode <= range.1 {
            (*&self.rslast.rank(self.forward(selectnode) as u64).unwrap() - 1) as usize
        } else if range.0 <= negselect && negselect <= range.1 {
            (*&self.rslast.rank(self.forward(negselect) as u64).unwrap() - 1) as usize
        } else {
            self.out().len()
        };
    }

    pub fn successors(&self, i: usize) -> Vec<usize> {
        let mut succs = Vec::new();
        let range = self.node_range(i);
        for i in range.0..range.1 + 1 {
            succs.push(self.rslast.rank(self.forward(i) as u64).unwrap() as usize - 1);
        }
        succs
    }

    pub fn label(&self, i: usize) -> String {
        let nodepos = match &self.rslast().select(i as u64 + 1) {
            Some(t) => *t,
            None => 0 as u64,
        } as usize;
        let mut index = nodepos;
        let mut symbol = self.findsymbol(index);
        let mut label = symbol.to_string();
        let mut count = self.kmersize - 2;
        println!("{} at {}", symbol, index);
        if symbol == '$' {
            while count != 0 {
                label.push('$');
                count -= 1;
            }

            return label;
        }

        while count != 0 {
            let new_index = self.backward(index);
            symbol = self.findsymbol(new_index);
            label.push(symbol);
            println!("{}->", index);
            /*if self.last[new_index] == 0 {
                index = match &self.rslast().select(new_index as u64) {
                    Some(t) => *t,
                    None => 0 as u64,
                } as usize;
            }*/
            index = match &self.rslast().select(new_index as u64) {
                Some(t) => *t,
                None => 0 as u64,
            } as usize;
            println!("{}, {}, {}", index, new_index, symbol);
            count -= 1;
            if symbol == '$' {
                while count != 0 {
                    label.push('$');
                    count -= 1;
                }
                break;
            }
        }
        let label = label.chars().rev().collect::<String>();
        label
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
         sdbg.print()
         assert_eq!(sdbg.nodes.len(), 16);
     }*/

    #[test]
    fn test_sdbg_forward() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut forw = Vec::new();
        let forcheck = vec![10, 4, 5, 8, 11, 9, 10, 1, 12, 2, 4, 13, 6];
        for i in 0..sdbg.nodes.len() {
            forw.push(sdbg.forward(i));
        }
        assert_eq!(forw, forcheck);
    }

    #[test]
    fn test_sdbg_backward() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut backw = Vec::new();
        let backcheck = vec![13, 7, 9, 1, 1, 2, 12, 3, 3, 5, 0, 4, 8];
        for i in 0..sdbg.nodes.len() {
            backw.push(sdbg.backward(i));
        }
        assert_eq!(backw, backcheck);
    }

    #[test]
    fn test_sdbg_outdegree() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        println!("{:?}", sdbg.node_char());
        let mut outd = Vec::new();
        let outcheck = vec![1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1];
        for i in 1..sdbg.nodes().len() {
            outd.push(sdbg.outdegree(i));
        }
        assert_eq!(outd, outcheck);
    }

    /* #[test]
    fn test_sdbg_outgoing() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);

    }*/
    /*#[test]
    fn test_sdbg_label() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        for i in 0..sdbg.nodes.len() {
            let l = sdbg.label(i);
            println!("{} \n---------------", l);
        }
    }*/

    /*#[test]
    fn test_dot() {
        let mut kmers = vec!["TACACT".to_string(),
                             "TACTCA".to_string(),
                             "GACTCG".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        sdbg.to_dot("output/test.dot");
        assert!(Path::new("output/test.dot").exists());
    }*/


    /*#[test]
    fn test_dot_2() {
        let mut kmers = vec!["ATCTTGCATTACCGCCCCAATC".to_string(),
                             "ATCTTACATTACCGTCCCAACC".to_string()];
        let sdbg = SDbg::new(&mut kmers, 6);
        sdbg.print();

        sdbg.to_dot("output/test2.dot");
        //assert!(Path::new("output/test2.dot").exists());
    }*/
}
