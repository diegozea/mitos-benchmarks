import gzip
import time
from Bio import AlignIO


def read_stockholm(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as handle:
            return list(AlignIO.parse(handle, "stockholm"))
    with open(path) as handle:
        return list(AlignIO.parse(handle, "stockholm"))


if __name__ == "__main__":
    start = time.perf_counter()
    alignments = read_stockholm("../data/PF08171.sth")
    elapsed = time.perf_counter() - start
    print("[BENCH] Biopython Stockholm PF08171 (parse+annotations):", elapsed)

    if not alignments:
        raise RuntimeError("No alignments parsed from Stockholm file")

    alignment = alignments[0]
    ann = alignment.annotations
    col_ann = alignment.column_annotations
    seqs = list(alignment)

    ss_cons_key = "consensus secondary structure"
    if ss_cons_key not in col_ann and "secondary_structure" in col_ann:
        ss_cons_key = "secondary_structure"
    residue_ss_key = "secondary structure"
    if residue_ss_key not in seqs[0].letter_annotations and "secondary_structure" in seqs[0].letter_annotations:
        residue_ss_key = "secondary_structure"

    n_seq = len(seqs)
    n_seq_with_ac = sum(1 for record in seqs if "accession" in record.annotations)
    has_ss_cons = ss_cons_key in col_ann
    n_with_residue_ss = sum(1 for record in seqs if residue_ss_key in record.letter_annotations)

    start = time.perf_counter()
    _ = ann.get("accession")
    _ = ann.get("identifier")
    _ = col_ann.get(ss_cons_key)
    _ = [record.annotations.get("accession") for record in seqs]
    _ = [record.letter_annotations.get(residue_ss_key) for record in seqs]
    access_elapsed = time.perf_counter() - start
    print("[BENCH] Biopython Stockholm PF08171 annotation access:", access_elapsed)

    print("Sequences:", n_seq)
    print("Sequences with accession:", n_seq_with_ac)
    print("Has consensus secondary structure:", has_ss_cons)
    print("Sequences with residue secondary structure:", n_with_residue_ss)
