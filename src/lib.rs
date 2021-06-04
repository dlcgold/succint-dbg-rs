//! Draft implementation of succint de Bruijn Graph of [Alex Bowe](https://alexbowe.com/succinct-debruijn-graphs/)
//! ([Python Code](https://github.com/alexbowe/debby))
//!
//! This implementation is slow and completely unoptimized
//!
//! TODO:
//! - [x] documentation (draft)
//! - [ ] optimizations
//! - [ ] refactoring
//! - [ ] more tests
//!
//! Examples will follow the graph in Alex Bowe's website
//!
//! The graph is derived from the string "TACGACGTCGACT" and has this rappresentation
//! (edge index, node index and label are explained only for simplicity regarding the examples,
//! they are not they are not stored):
//! ```rust,ignore
//! last edge negative edge_index node_index label
//! 1	 T	  false    0          0          $$$
//! 1	 C	  false    1          1          CGA
//! 1	 C	  false    2          2          $TA
//! 0	 G	  false    3          3          GAC
//! 1	 T	  false    4          3          GAC
//! 1	 G	  true     5          4          TAC
//! 1	 G	  false    6          5          GTC
//! 0	 A	  false    7          6          ACG
//! 1	 T	  false    8          6          ACG
//! 1	 A	  true     9          7          TCG
//! 1	 A	  false    10         8          $$T
//! 1	 $	  false    11         9          ACT
//! 1	 C	  false    12         10         CGT
//!
//! F array:
//! F($) = 0
//! F(A) = 1
//! F(C) = 3
//! F(G) = 7
//! F(T) = 10
//! ```
//! And the dot of the graph is:
//! <div>
//! <img src="../../../output/gh.png" />
//! </div>
//! <hr/>
extern crate bio;

use std::collections::{HashSet, HashMap};
use bio::data_structures::rank_select::RankSelect;
use bv::BitVec;
#[allow(unused_imports)]
use std::fs::File;
#[allow(unused_imports)]
use std::io::Write;
use math::round;

/// utility for create kmer-set (with correct number of $ at the begin) from a string.
#[allow(dead_code)]
fn create_kmers(s: String, k: u32) -> Vec<String> {
    get_kmers(&mut vec![s], k)
}

/// utility for create kmer-set (with correct number of $ at the begin) from Vec a string
#[allow(dead_code)]
fn get_kmers(reads: &mut Vec<String>, k: u32) -> Vec<String> {
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

/// struct to implement succint dbg with
/// - kmersize: the kmer size
/// - n_nodes: amount of nodes
/// - last: last array
/// - edge: edges array
/// - fvec: F array
/// - neg: array to recognize edges marked as negative
/// - rschar: map to obtain rank/select of edges
/// - rscharneg: map to obtain rank/select of edges marked as negative
/// - rslast: map to obtain rank/select of last array
pub struct SDbg {
    kmersize: u32,
    n_nodes: usize,
    last: Vec<usize>,
    edge: Vec<char>,
    fvec: Vec<usize>,
    neg: Vec<bool>,
    rschar: HashMap<char, RankSelect>,
    rscharneg: HashMap<char, RankSelect>,
    rslast: RankSelect,
}

impl SDbg {
    /// create new succint dbg from a set of strings
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// assert_eq!(sdbg.n_nodes(), 11)
    /// ```
    pub fn new(reads: &mut Vec<String>, k: u32) -> Self {
        let mut node_edge: Vec<(String, char, String)> = Vec::new();
        let kmers = get_kmers(reads, k);
        let mut set_check = HashSet::new();
        for kmer in &kmers {
            let kmer1 = kmer[0..(k - 1) as usize].to_string();
            node_edge.push((kmer[0..(k - 1) as usize].to_string(),
                            kmer.chars().last().unwrap(),
                            kmer1.chars().rev().collect()));
            set_check.insert(kmer1);
        }
        for kmer in kmers {
            let kmer2 = kmer[1..].to_string();
            if !set_check.contains(&kmer2) {
                node_edge.push((kmer[1..].to_string(),
                                '$',
                                kmer2.chars().rev().collect()));
            }
        }

        node_edge.sort_by(|a, b| (a.2.cmp(&b.2)));
        let mut real_order = Vec::new();
        let mut tmp = Vec::new();
        for i in 0..node_edge.len() {
            if i != node_edge.len() - 1 && &node_edge[i].0 == &node_edge[i + 1].0 {
                tmp.push(node_edge[i].clone());
            } else {
                tmp.push(node_edge[i].clone());
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
        let mut edge = Vec::new();
        let mut neg = vec![false; real_order.len()];
        let mut neg_pos = Vec::new();
        for i in 0..real_order.len() - 1 {
            if &real_order[i].0 == &real_order[i + 1].0 {
                last.push(0);
            } else {
                last.push(1);
            }
            nodes.push(real_order[i].0.clone());
            edge.push(real_order[i].1);
        }
        last.push(1);
        nodes.push(real_order[real_order.len() - 1].0.clone());
        edge.push(real_order[real_order.len() - 1].1);
        let mut check_neg: HashSet<(String, char)> = HashSet::new();
        let mut fvec = vec![0; 256];
        for i in 0..real_order.len() {
            if i == 0 {
                fvec['$' as usize] = 0;
            } else if nodes[i - 1].chars().last().unwrap() != nodes[i].chars().last().unwrap() {
                fvec[nodes[i].chars().last().unwrap() as usize] = i;
            }

            if !check_neg.contains(&(nodes[i][1..].to_string(), edge[i])) {
                check_neg.insert((nodes[i][1..].to_string(), edge[i]));
            } else {
                neg[i] = true;
                neg_pos.push(i);
            }
        }
        let mut bitvecs = HashMap::new();
        let mut bitvecsneg = HashMap::new();
        for s in "$acgtACGT".chars() {
            let mut bv = BitVec::new();
            let mut bvneg = BitVec::new();
            for (index, c) in edge.iter().enumerate() {
                if c == &s && !neg[index] {
                    bv.push(true);
                } else {
                    bv.push(false);
                }
                if c == &s && neg[index] {
                    bvneg.push(true);
                } else {
                    bvneg.push(false);
                }
            }
            let krank = round::ceil(((bv.len() as f64).log2()).powf(2.) / (32 as f64), 0);
            let krankneg = round::ceil(((bvneg.len() as f64).log2()).powf(2.) / (32 as f64), 0);
            let rank_select = RankSelect::new(bv.clone(), krank as usize);
            let rank_selectneg = RankSelect::new(bvneg.clone(), krankneg as usize);
            bitvecs.insert(s, rank_select);
            bitvecsneg.insert(s, rank_selectneg);
        }
        let mut bv = BitVec::new();
        for e in &last {
            match e {
                1 => bv.push(true),
                _ => bv.push(false),
            }
        }
        let krank = round::ceil(((bv.len() as f64).log2()).powf(2.) / (32 as f64), 0);
        let rslast = RankSelect::new(bv.clone(), krank as usize);
        let n_nodes = rslast.rank(edge.len() as u64 - 1).unwrap() as usize;
        SDbg {
            kmersize: k,
            n_nodes,
            last,
            edge,
            neg,
            rschar: bitvecs,
            rscharneg: bitvecsneg,
            fvec,
            rslast,
        }
    }

    /// create new succint dbg from a string
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let sdbg = SDbg::new_from_string(&"TACGACGTCGACT".to_string(), 4);
    /// assert_eq!(sdbg.n_nodes(), 11)
    /// ```
    pub fn new_from_string(read: &String, k: u32) -> Self {
        SDbg::new(&mut vec![read.clone()], k)
    }
    /// return the chosen k for kmers
    pub fn kmersize(&self) -> u32 { self.kmersize }
    /// return the amount of nodes in the succint dbg
    pub fn n_nodes(&self) -> usize { self.n_nodes }
    /// return the last array (as explained in Bowe's research)
    pub fn last(&self) -> &Vec<usize> { &self.last }
    /// return the vec of edges (as explained in Bowe's research)
    pub fn edge(&self) -> &Vec<char> { &self.edge }
    /// return the bool array that it specifies if edge in marked as negative
    pub fn neg(&self) -> &Vec<bool> { &self.neg }
    /// return the F array (as explained in Bowe's research)
    pub fn fvec(&self) -> &Vec<usize> { &self.fvec }
    /// return the map with the bitvec of a char, char used by key, with rank/select (as explained in Bowe's research)
    pub fn rschar(&self) -> &HashMap<char, RankSelect> { &self.rschar }
    /// return the map with the bitvec of a char marked as neg, char used by key, with rank/select (as explained in Bowe's research)
    pub fn rscharneg(&self) -> &HashMap<char, RankSelect> { &self.rscharneg }
    /// return the map with the bitvec of the last array with rank/select (as explained in Bowe's research)
    pub fn rslast(&self) -> &RankSelect { &self.rslast }

    /// print main structires of succint dbg
    pub fn print(&self) {
        println!("last edge negative");
        for i in 0..self.edge().len() {
            println!("{}\t{}\t{}", self.last[i], self.edge[i], self.neg[i]);
        }
        println!("F array:");
        for s in vec!['$', 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'] {
            if self.fvec[s as usize] != 0 || s == '$' {
                println!("F({}) = {}", s, self.fvec[s as usize]);
            }
        }
    }
    /// return the last char at index F array
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut symb = Vec::new();
    /// let symbcheck = vec!['$', 'A', 'A', 'C', 'C', 'C', 'C', 'G', 'G', 'G', 'T', 'T', 'T'];
    /// for i in 0..sdbg.edge().len() {
    ///     symb.push(sdbg.findsymbol(i as isize));
    /// }
    /// assert_eq!(symb, symbcheck);
    /// ```
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

    /// return first edge from a node
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut edged = Vec::new();
    /// let edgecheck = vec![0, 1, 2, 3, 5, 6, 7, 9, 10, 11, 12];
    /// for i in 0..sdbg.n_nodes() {
    ///     edged.push(sdbg.first_edge(i as isize));
    /// }
    /// assert_eq!(edged, edgecheck);
    ///```
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

    /// return last edge from a node
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut edged = Vec::new();
    /// let edgecheck = vec![0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12];
    /// for i in 0..sdbg.n_nodes() {
    /// edged.push(sdbg.last_edge(i as isize));
    /// }
    /// assert_eq!(edged, edgecheck);
    /// ```
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

    /// return range of edges from a node
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut ranged = Vec::new();
    /// let rangecheck = vec![(0, 0), (1, 1), (2, 2), (3, 4), (5, 5), (6, 6), (7, 8),
    ///                       (9, 9), (10, 10), (11, 11), (12, 12)];
    /// for i in 0..sdbg.n_nodes() {
    ///     ranged.push(sdbg.node_range(i as isize));
    /// }
    /// assert_eq!(ranged, rangecheck);
    /// ```
    pub fn node_range(&self, i: isize) -> (isize, isize) {
        (self.first_edge(i), self.last_edge(i))
    }

    /// change edge index in node index
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut edged = Vec::new();
    /// let edgecheck = vec![0, 1, 2, 3, 3, 4, 5, 6, 6, 7, 8, 9, 10];
    /// for i in 0..sdbg.edge().len() {
    ///     edged.push(sdbg.edge_to_node(i as isize));
    /// }
    /// assert_eq!(edged, edgecheck);
    /// ```
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

    /* change node index in edge index
    pub fn node_to_edge(&self, i: isize) -> isize {
        if i <= 0 {
            return 0;
        }
        let select = match self.rslast().select(i as u64) {
            Some(t) => t,
            None => 0,
        } - 1;
        select as isize
    }*/

    /// return index of the last edge of the node pointed to by edge as parameter (as explained in Bowe's research)
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut forw = Vec::new();
    /// let forcheck = vec![10, 4, 5, 8, 11, 8, 9, 1, 12, 1, 2, -1, 6];
    /// for i in 0..sdbg.edge().len() {
    ///     forw.push(sdbg.forward(i as isize));
    /// }
    /// assert_eq!(forw, forcheck);
    /// ```
    pub fn forward(&self, i: isize) -> isize {
        let symbol = self.edge[i as usize];
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

    /// return index of the first edge that points to the node that the edge at index exits (as explained in Bowe's research)
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut backw = Vec::new();
    /// let backcheck = vec![-1, 7, 10, 1, 1, 2, 12, 3, 3, 6, 0, 4, 8];
    /// for i in 0..sdbg.edge().len() {
    ///     backw.push(sdbg.backward(i as isize));
    /// }
    /// assert_eq!(backw, backcheck);
    /// ```
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

    /// return number of outgoing edges from a node (as explained in Bowe's research)
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut edged = Vec::new();
    /// let edgecheck = vec![1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1];
    /// for i in 0..sdbg.n_nodes() {
    ///     edged.push(sdbg.outdegree(i as isize));
    /// }
    /// assert_eq!(edged, edgecheck);
    /// ```
    pub fn outdegree(&self, i: isize) -> isize {
        if self.edge()[i as usize] == '$' {
            return 0;
        }
        let range = self.node_range(i);
        if range.1 > range.0 {
            return range.1 - range.0 + 1;
        } else {
            return 1;
        }
    }

    /// return number of incoming edges to a node (as explained in Bowe's research)
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut ind = Vec::new();
    /// let incheck = vec![0, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1];
    /// for i in 0..sdbg.n_nodes() {
    ///     ind.push(sdbg.indegree(i as isize));
    /// }
    /// assert_eq!(ind, incheck);
    /// ```
    pub fn indegree(&self, i: isize) -> isize {
        let last = self.last_edge(i);
        let pred = self.backward(last);
        if pred == -1 {
            return 0;
        }
        let symbol = self.edge()[pred as usize];
        let next;
        if pred + 1 <= 0 {
            next = -1
        } else {
            next = match &self.rschar()[&symbol].select(pred as u64 + 1) {
                Some(t) => *t,
                None => self.edge().len() as u64 - 1,
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

    /// return index obtained from a node, follow the edge labeled by a symbol (as explained in Bowe's research)
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut edgeg = Vec::new();
    /// let edgecheck = vec![-1, -1, -1, -1, 8, -1, -1, 3, -1, -1, -1, -1, 4, -1, -1, -1, -1, -1,
    ///                       6, 9, -1, -1, -1, 6, -1, -1, -1, -1, 7, -1, -1, 1, -1, -1, 10, -1, 1,
    ///                       -1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 5, -1, -1];
    /// for i in 0..sdbg.n_nodes() {
    ///     for j in ['$', 'A', 'C', 'G', 'T'] {
    ///         edgeg.push(sdbg.outgoing(i as isize, j));
    ///     }
    /// }
    /// assert_eq!(edgeg, edgecheck);
    /// ```
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
    /// return vec with index of successors from a node (as explained in Bowe's research)
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut edges = Vec::new();
    /// let succcheck = vec![vec![8], vec![3], vec![4], vec![6, 9], vec![6], vec![7],
    ///                      vec![1, 10], vec![1], vec![2], vec![-1], vec![5]];
    /// for i in 0..sdbg.n_nodes() {
    ///     edges.push(sdbg.successors(i as isize));
    /// }
    /// assert_eq!(edges, succcheck);
    /// ```
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

    /// return label of a node (as explained in Bowe's research)
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut label = Vec::new();
    /// let labelcheck = vec!["$$$", "CGA", "$TA", "GAC", "TAC", "GTC", "ACG", "TCG",
    ///                       "$$T", "ACT", "CGT", "$$T", "CGA"];
    /// for i in 0..sdbg.edge().len() {
    ///     label.push(sdbg.label(i as isize));
    /// }
    /// assert_eq!(label, labelcheck);
    /// ```
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

    fn firstchar(&self, i: isize) -> char {
        let node = self.edge_to_node(i);
        self.label(node).chars().next().unwrap()
    }

    fn selector(&self, i: isize, symbol: char, flag: isize, pred: isize) -> isize {
        if i > 0 {
            let select;
            if flag + i <= 0 {
                select = -1
            } else {
                select = match &self.rscharneg()[&symbol].select((flag + i) as u64) {
                    Some(t) => *t,
                    None => 0,
                } as isize;
            }
            return select;
        } else {
            return pred;
        }
    }

    fn accessor(&self, i: isize, symbol: char, flag: isize, pred: isize) -> char {
        self.firstchar(self.selector(i, symbol, flag, pred))
    }

    /// return predecessor node starting with a symbol, that has an edge to a node (as explained in Bowe's research)
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// let mut ing = Vec::new();
    /// let incheck = vec![-1, -1, -1, -1, -1, -1, 6, -1, -1, 7, 8, -1, -1, -1, -1, -1, -1,
    ///                    1, -1, -1, 2, -1, -1, -1, -1, -1, -1, 10, -1, -1, -1, -1, -1, 3, 4, -1,
    ///                    -1, -1, 5, -1, 0, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, 6, -1, -1, -1];
    /// for i in 0..sdbg.n_nodes() {
    ///     for j in ['$', 'A', 'C', 'G', 'T'] {
    ///         ing.push(sdbg.incoming(i as isize, j));
    ///     }
    /// }
    /// assert_eq!(ing, incheck);
    /// ```
    pub fn incoming(&self, i: isize, symbolf: char) -> isize {
        let last = self.last_edge(i);
        let pred = self.backward(last);
        if pred == -1 {
            return -1;
        }
        let symbol = self.edge()[pred as usize];
        let next;
        if pred + 1 <= 0 {
            next = -1
        } else {
            next = match &self.rschar()[&symbol].select(pred as u64 + 1) {
                Some(t) => *t,
                None => self.edge().len() as u64 - 1,
            } as isize;
        }
        let flag = match self.rscharneg[&symbol].rank(pred as u64) {
            Some(t) => t,
            None => 0
        } as isize;
        let indeg = match self.rscharneg[&symbol].rank(next as u64) {
            Some(t) => t,
            None => 0
        } as isize - flag + 1;

        let mut subind: isize = -1;
        let mut tmp = Vec::new();
        for ind in 0..indeg {
            tmp.push(self.accessor(ind, symbol, flag, pred));
        }
        tmp.sort();
        for (ind, elem) in tmp.iter().enumerate() {
            if ind != tmp.len() && elem == &symbolf {
                subind = ind as isize;
            }
        }
        if subind == -1 {
            return -1;
        }
        self.edge_to_node(self.selector(subind, symbol, flag, pred))
    }

    /// print to file of the .dot of tne graph
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// use std::path::Path;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// sdbg.to_dot("output/bowe.dot");
    /// assert!(Path::new("output/bowe.dot").exists());
    /// ```
    pub fn to_dot(&self, edgeput: &str) {
        let mut fileedge = File::create(edgeput).expect("error");
        fileedge
            .write("digraph sample{\n".as_bytes())
            .expect("error");
        for i in 0..self.n_nodes() {
            for j in self.successors(i as isize) {
                if j != -1 {
                    let start = self.label(i as isize);
                    let end = self.label(j as isize);

                    fileedge.write(
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
        fileedge.write("}".as_bytes()).expect("error");
    }

    /// print to file of the .dot of tne graph without $ nodes
    /// # Example
    /// ```
    /// use succint_dbg_rs::SDbg;
    /// use std::path::Path;
    /// let mut kmers = vec!["TACGACGTCGACT".to_string()];
    /// let sdbg = SDbg::new(&mut kmers, 4);
    /// sdbg.to_dot_no_dollar("output/bowend.dot");
    /// assert!(Path::new("output/bowend.dot").exists());
    /// ```
    pub fn to_dot_no_dollar(&self, edgeput: &str) {
        let mut fileedge = File::create(edgeput).expect("error");
        fileedge
            .write("digraph sample{\n".as_bytes())
            .expect("error");
        for i in 0..self.n_nodes() {
            if self.label(i as isize).chars().next().unwrap() != '$' {
                for j in self.successors(i as isize) {
                    if j != -1 {
                        let start = self.label(i as isize);
                        let end = self.label(j as isize);

                        fileedge.write(
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
        }
        fileedge.write("}".as_bytes()).expect("error");
    }
}


#[cfg(test)]
mod tests {
    use crate::SDbg;
    #[allow(unused_imports)]
    use std::path::Path;

    #[test]
    fn test_sdbg() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let lastcheck: Vec<usize> = vec![1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1];
        let edgecheck = vec!['T', 'C', 'C', 'G', 'T', 'G', 'G', 'A', 'T', 'A', 'A',
                             '$', 'C'];
        let negcheck = vec![false, false, false, false, false, true, false, false,
                            false, true, false, false, false];
        let fcheck = vec![0, 1, 3, 7, 10];
        let mut fvec = Vec::new();
        for elem in vec!['$', 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'] {
            if sdbg.fvec[elem as usize] != 0 || elem == '$' {
                fvec.push(sdbg.fvec[elem as usize]);
            }
        }
        assert_eq!(sdbg.last(), &lastcheck);
        assert_eq!(sdbg.edge(), &edgecheck);
        assert_eq!(sdbg.neg(), &negcheck);
        assert_eq!(&fvec, &fcheck);
    }

    #[test]
    fn test_sdbg_findsymbol() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut symb = Vec::new();
        let symbcheck = vec!['$', 'A', 'A', 'C', 'C', 'C', 'C', 'G', 'G', 'G', 'T', 'T', 'T'];
        for i in 0..sdbg.edge().len() {
            symb.push(sdbg.findsymbol(i as isize));
        }
        assert_eq!(symb, symbcheck);
    }

    #[test]
    fn test_sdbg_forward() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut forw = Vec::new();
        let forcheck = vec![10, 4, 5, 8, 11, 8, 9, 1, 12, 1, 2, -1, 6];
        for i in 0..sdbg.edge().len() {
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
        for i in 0..sdbg.edge().len() {
            backw.push(sdbg.backward(i as isize));
        }
        assert_eq!(backw, backcheck);
    }

    #[test]
    fn test_sdbg_outdegree() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut edged = Vec::new();
        let edgecheck = vec![1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1];
        for i in 0..sdbg.n_nodes() {
            edged.push(sdbg.outdegree(i as isize));
        }
        assert_eq!(edged, edgecheck);
    }

    #[test]
    fn test_sdbg_firstedge() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut edged = Vec::new();
        let edgecheck = vec![0, 1, 2, 3, 5, 6, 7, 9, 10, 11, 12];
        for i in 0..sdbg.n_nodes() {
            edged.push(sdbg.first_edge(i as isize));
        }
        assert_eq!(edged, edgecheck);
    }

    #[test]
    fn test_sdbg_lastedge() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut edged = Vec::new();
        let edgecheck = vec![0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12];
        for i in 0..sdbg.n_nodes() {
            edged.push(sdbg.last_edge(i as isize));
        }
        assert_eq!(edged, edgecheck);
    }

    #[test]
    fn test_sdbg_noderange() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut ranged = Vec::new();
        let rangecheck = vec![(0, 0), (1, 1), (2, 2), (3, 4), (5, 5), (6, 6), (7, 8), (9, 9), (10, 10), (11, 11), (12, 12)];
        for i in 0..sdbg.n_nodes() {
            ranged.push(sdbg.node_range(i as isize));
        }
        assert_eq!(ranged, rangecheck);
    }

    #[test]
    fn test_sdbg_edgenode() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut edged = Vec::new();
        let edgecheck = vec![0, 1, 2, 3, 3, 4, 5, 6, 6, 7, 8, 9, 10];
        for i in 0..sdbg.edge().len() {
            edged.push(sdbg.edge_to_node(i as isize));
        }
        assert_eq!(edged, edgecheck);
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
        let mut edgeg = Vec::new();
        let edgecheck = vec![-1, -1, -1, -1, 8, -1, -1, 3, -1, -1, -1, -1, 4, -1, -1, -1,
                             -1, -1, 6, 9, -1, -1, -1, 6, -1, -1, -1, -1, 7, -1, -1, 1, -1, -1, 10,
                             -1, 1, -1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 5,
                             -1, -1];
        for i in 0..sdbg.n_nodes() {
            for j in ['$', 'A', 'C', 'G', 'T'] {
                edgeg.push(sdbg.outgoing(i as isize, j));
            }
        }
        assert_eq!(edgeg, edgecheck);
    }

    #[test]
    fn test_sdbg_label() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut label = Vec::new();
        let labelcheck = vec!["$$$", "CGA", "$TA", "GAC", "TAC", "GTC", "ACG", "TCG",
                              "$$T", "ACT", "CGT", "$$T", "CGA"];
        for i in 0..sdbg.edge().len() {
            label.push(sdbg.label(i as isize));
        }
        assert_eq!(label, labelcheck);
    }

    #[test]
    fn test_sdbg_successors() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut edges = Vec::new();
        let succcheck = vec![vec![8], vec![3], vec![4], vec![6, 9], vec![6], vec![7],
                             vec![1, 10], vec![1], vec![2], vec![-1], vec![5]];
        for i in 0..sdbg.n_nodes() {
            edges.push(sdbg.successors(i as isize));
        }
        assert_eq!(edges, succcheck);
    }

    #[test]
    fn test_sdbg_incoming() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        let mut ing = Vec::new();
        let incheck = vec![-1, -1, -1, -1, -1, -1, 6, -1, -1, 7, 8, -1, -1, -1, -1, -1, -1,
                           1, -1, -1, 2, -1, -1, -1, -1, -1, -1, 10, -1, -1, -1, -1, -1, 3, 4, -1,
                           -1, -1, 5, -1, 0, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, 6, -1, -1, -1];
        for i in 0..sdbg.n_nodes() {
            for j in ['$', 'A', 'C', 'G', 'T'] {
                ing.push(sdbg.incoming(i as isize, j));
            }
        }
        assert_eq!(ing, incheck);
    }

    #[test]
    fn test_sdbg_to_dot() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        sdbg.to_dot("output/bowe.dot");
        assert!(Path::new("output/bowe.dot").exists());
    }

    #[test]
    fn test_sdbg_to_dot_nodol() {
        let mut kmers = vec!["TACGACGTCGACT".to_string()];
        let sdbg = SDbg::new(&mut kmers, 4);
        sdbg.to_dot_no_dollar("output/bowend.dot");
        assert!(Path::new("output/bowend.dot").exists());
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
