// make -C bwa-sys/bwa/ -n libbwa.a | grep -o -E "[A-Za-z0-9_]+\.c"
const FILES: &[&str] = &[
    "bwa/utils.c",
    "bwa/kthread.c",
    "bwa/kstring.c",
    "bwa/ksw.c",
    "bwa/bwt.c",
    "bwa/bntseq.c",
    "bwa/bwa.c",
    "bwa/bwamem.c",
    "bwa/bwamem_pair.c",
    "bwa/bwamem_extra.c",
    "bwa/malloc_wrap.c",
    "bwa/QSufSort.c",
    "bwa/bwt_gen.c",
    "bwa/rope.c",
    "bwa/rle.c",
    "bwa/is.c",
    "bwa/bwtindex.c"
];

// make -C bwa-sys/bwa/ -nd libbwa.a | grep -o -E "[A-Za-z0-9_]+\.h" | sort | uniq
const HEADERS: &[&str] = &[
    "bwa/QSufSort.h",
    "bwa/bntseq.h",
    "bwa/bwa.h",
    "bwa/bwamem.h",
    "bwa/bwt.h",
    "bwa/kbtree.h",
    "bwa/khash.h",
    "bwa/kseq.h",
    "bwa/ksort.h",
    "bwa/kstring.h",
    "bwa/ksw.h",
    "bwa/kvec.h",
    "bwa/malloc_wrap.h",
    "bwa/neon_sse.h",
    "bwa/rle.h",
    "bwa/rope.h",
    "bwa/scalar_sse.h",
    "bwa/utils.h"
];

fn main() {
    for file in FILES {
        println!("cargo:rerun-if-changed={}", file);
    }
    for file in HEADERS {
        println!("cargo:rerun-if-changed={}", file);
    }
    cc::Build::new()
        .define("COMPILATION_TIME_PLACE", "\"build.rs\"")
        .warnings(false)
        .extra_warnings(false)
        .files(FILES)
        .flag("-fPIC")
        .compile("bwa");
}
