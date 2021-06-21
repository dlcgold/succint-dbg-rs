# succint-dbg-rs

Porting in Rust of succint De Bruijn graph from [Alex Bowe](https://alexbowe.com/succinct-debruijn-graphs/)
([Python Code](https://github.com/alexbowe/debby))

Succint rappresentation support:
- forward(i)
- backward(i)
- outdegree(n)
- indegree(n)
- outgoing(n, s)
- incoming(n, s)
- successors(n)
- label(n)
- creation from a String, from a String Vector, from a FASTQ file and from a FASTA file
- print of the vectors that represent the dbg
- dot creation

Examples will follow the graph in Alex Bowe's website.

```Rust
    let mut reads = vec!["TACGACGTCGACT".to_string()];
    let sdbg = SDbg::new(&mut reads, 4);
```

 The graph is derived from the string "TACGACGTCGACT" and has this rappresentation
 (edge index, node index and label are explained only for simplicity regarding the examples,
 they are not they are not stored):
 ```rust
 last   edge    negative edge_index node_index label
 1      'T'     false    0          0          "$$$"
 1      'C'     false    1          1          "CGA"
 1      'C'     false    2          2          "$TA"
 0      'G'     false    3          3          "GAC"
 1      'T'     false    4          3          "GAC"
 1      'G'     true     5          4          "TAC"
 1      'G'     false    6          5          "GTC"
 0      'A'     false    7          6          "ACG"
 1      'T'     false    8          6          "ACG"
 1      'A'     true     9          7          "TCG"
 1      'A'     false    10         8          "$$T"
 1      '$'     false    11         9          "ACT"
 1      'C'     false    12         10         "CGT"

 F array:
 F('$') = 0
 F('A') = 1
 F('C') = 3
 F('G') = 7
 F('T') = 10
 ```
![](output/gh.png)

More examples in Doc:
```bash
cargo doc --open
```

## Todo
- [x] documentation (draft)
- [ ] optimizations
- [ ] refactoring
- [ ] more tests
