import gzip
import os
import tempfile
import time
from prody import parseMSA

STOCKHOLM = "Stockholm"


def parse_stockholm(path):
    if path.endswith(".gz"):
        # ProDy's parser can choke on blank lines in compressed inputs; filter them.
        with gzip.open(path, "rt") as handle, tempfile.NamedTemporaryFile(
            mode="wb", suffix=".sth", delete=False
        ) as tmp:
            for line in handle:
                if line.strip() == "":
                    continue
                tmp.write(line.encode())
        try:
            return parseMSA(tmp.name, format=STOCKHOLM)
        finally:
            os.remove(tmp.name)
    return parseMSA(path, format=STOCKHOLM)


def bench_read(path, label, **kwargs):
    start = time.perf_counter()
    if kwargs.get("format") == STOCKHOLM:
        parse_stockholm(path)
    else:
        parseMSA(path, **kwargs)
    elapsed = time.perf_counter() - start
    print(f"[BENCH] ProDy Read {label}: {elapsed}")


bench_read("../data/PF00089.stockholm.gz", "compressed Stockholm PF00089", compressed=True, format=STOCKHOLM)
bench_read("../data/PF00089.sth", "uncompressed Stockholm PF00089", format=STOCKHOLM)
bench_read("../data/PF00089.fasta.gz", "compressed FASTA PF00089", compressed=True)
bench_read("../data/PF00089.fasta", "uncompressed FASTA PF00089")

bench_read("../data/PF16957.stockholm.gz", "compressed Stockholm PF16957", compressed=True, format=STOCKHOLM)
bench_read("../data/PF16957.sth", "uncompressed Stockholm PF16957", format=STOCKHOLM)
bench_read("../data/PF16957.fasta.gz", "compressed FASTA PF16957", compressed=True)
bench_read("../data/PF16957.fasta", "uncompressed FASTA PF16957")
