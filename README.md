# Rust-bwa 

Rust-bwa is a Rust wrapper of the BWA api from 10x Genomics. Pass read-pair information in, and get Rust-htslib BAM records back.
Get started quickly with default settings & a reasonable paired-end model. See docs for more details on customing
parameters or the paired-end model.

```
extern crate bwa;
extern crate seq_io;

use bwa::BwaAligner;
use seqio::fastq

fn main() {
    let bwa = BwaAligner::from_path(&"tests/test_ref.fa").unwrap();

    let r1 = b"GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGG";
    let q1 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";
    let r2 = b"TGCTGCGTAGCAGATCGACCCAGGCATTCCCTAGCGTGCTCATGCTCTGGCTGGTAAACGCACGGATGAGGGCAAAAAT";
    let q2 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";

    let (r1_alns, _r2_alns) = bwa.align_read_pair(b"read_name", r1, q1, r2, q2).unwrap();
    println!("r1 mapping -- tid: {}, pos: {}", r1_alns[0].tid(), r1_alns[0].pos());

    let fastq = b"@id1
    GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGG
    +
    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    @id2
    TGCTGCGTAGCAGATCGACCCAGGCATTCCCTAGCGTGCTCATGCTCTGGCTGGTAAACGCACGGATGAGGGCAAAAAT
    +
    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII";
    
    // collect vector of seqio::fastq records
    let mut reader = fastq::Reader::new(&fastq[..]);
    let seqs = reader.records().collect::<Result<Vec<_>, _>>().unwrap();
    
    let alns = bwa.align_fastq_records(&seqs, /* paired */ false).unwrap();
    
    println!("mapping id1 -- tid: {}, pos: {}", alns[0].tid(), alns[0].pos());
    println!("mapping id2 -- tid: {}, pos: {}", alns[1].tid(), alns[1].pos());
}
```

Pre-built rust bindings were generated using `bindgen` for linux using the command:

```
cd bwa-sys
bindgen \
  --no-doc-comments \
  --allowlist-function mem_align1_core \
  --allowlist-function mem_sam_pe \
  --allowlist-function mem_opt_init \
  --allowlist-function bwa_idx_build \
  --allowlist-function bwa_idx_load \
  --allowlist-function bwa_idx_destroy \
  --allowlist-function mem_process_seq_pe \
  --allowlist-function mem_process_seqs \
  --allowlist-function bwa_fill_scmat \
  --allowlist-var "BWA_IDX_.*" \
  --allowlist-var "MEM_F_.*" \
  wrapper.h \
  -o src/lib.rs
```

`bindgen` can be installed using `cargo install bindgen-cli`. See the documentation [here](https://rust-lang.github.io/rust-bindgen/command-line-usage.html).

