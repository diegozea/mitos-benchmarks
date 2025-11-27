import sys
import unittest
from pathlib import Path

import numpy as np
from Bio import AlignIO
from prody import buildSeqidMatrix, parseMSA

ROOT = Path(__file__).resolve().parents[1]
DATA_FASTA = ROOT / "tests" / "data" / "test_alignment.fasta"
DATA_STO = ROOT / "tests" / "data" / "test_alignment.sth"

# Allow importing helper functions from the benchmark scripts
sys.path.insert(0, str(ROOT / "BioPython"))
from Entropy_MI import alignment_array, mutual_information, shannon_entropy  # noqa: E402
from PID import pid_matrix  # noqa: E402

EXPECTED_PID = np.array(
    [
        [100.0, 75.0, 75.0],
        [75.0, 100.0, 50.0],
        [75.0, 50.0, 100.0],
    ]
)
EXPECTED_ENTROPY = np.array([0.0, 0.0, 0.9182958340544896, 0.9182958340544896])


def normalize_pid(mat: np.ndarray) -> np.ndarray:
    """Scale PID matrix to percentage using its diagonal as reference."""
    diag = np.diag(mat)
    scale = float(np.median(diag))
    if scale == 0:
        raise ValueError("PID matrix diagonal is zero; cannot normalize.")
    return mat / scale * 100.0


class TestPythonBenchmarks(unittest.TestCase):
    def test_biopython_pid_matches_expected(self):
        aln = AlignIO.read(DATA_FASTA, "fasta")
        observed = pid_matrix(aln)
        self.assertTrue(
            np.allclose(observed, EXPECTED_PID),
            "Biopython PID does not match expected percentages",
        )

    def test_prody_pid_matches_expected_after_scaling(self):
        msa = parseMSA(str(DATA_FASTA))
        observed = np.array(buildSeqidMatrix(msa), dtype=float)
        normalized = normalize_pid(observed)
        self.assertTrue(
            np.allclose(normalized, EXPECTED_PID, atol=1e-6),
            "ProDy PID (scaled) does not match expected percentages",
        )

    def test_biopython_entropy_matches_expected(self):
        aln = AlignIO.read(DATA_FASTA, "fasta")
        ent = shannon_entropy(alignment_array(aln))
        self.assertTrue(
            np.allclose(ent, EXPECTED_ENTROPY, atol=1e-6),
            "Biopython Shannon entropy differs from expected values",
        )

    def test_biopython_mi_matrix_invariants(self):
        arr = alignment_array(AlignIO.read(DATA_STO, "stockholm"))
        mi = mutual_information(arr)
        self.assertEqual(mi.shape, (arr.shape[1], arr.shape[1]))
        self.assertTrue(np.allclose(mi, mi.T, atol=1e-12), "MI matrix not symmetric")
        self.assertTrue(np.allclose(np.diag(mi), 0.0, atol=1e-12), "MI diagonal not zero")
        self.assertFalse(np.isnan(mi).any())
        self.assertGreater(mi[2, 3], 0.0, "Expected non-zero MI for varying columns")

    def test_stockholm_annotations_are_parsed(self):
        aln = AlignIO.read(DATA_STO, "stockholm")
        self.assertEqual(len(aln), 3)
        self.assertEqual(aln.get_alignment_length(), 4)
        # Column-level annotations may be stored under either Stockholm-style or Biopython keys
        ss_cons = (
            aln.column_annotations.get("consensus secondary structure")
            or aln.column_annotations.get("secondary_structure")
        )
        self.assertIsNotNone(ss_cons)
        self.assertEqual("".join(ss_cons), "HHHH")
        first = aln[0]
        self.assertEqual(first.annotations.get("accession"), "ACC1")
        ss_key = "secondary structure" if "secondary structure" in first.letter_annotations else "secondary_structure"
        self.assertIn(ss_key, first.letter_annotations)
        self.assertEqual("".join(first.letter_annotations[ss_key]), "HHHH")


if __name__ == "__main__":
    unittest.main()
