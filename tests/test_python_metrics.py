import sys
import unittest
from pathlib import Path

import numpy as np
from Bio import AlignIO
from Bio.Align import stockholm
from prody import buildSeqidMatrix, parseMSA

ROOT = Path(__file__).resolve().parents[1]
DATA_FASTA = ROOT / "tests" / "data" / "test_alignment.fasta"
DATA_STO = ROOT / "tests" / "data" / "test_alignment.sth"

# Allow importing helper functions from the benchmark scripts
sys.path.insert(0, str(ROOT / "BioPython"))
from PID import pid_matrix  # noqa: E402

EXPECTED_PID = np.array(
    [
        [100.0, 75.0, 75.0],
        [75.0, 100.0, 50.0],
        [75.0, 50.0, 100.0],
    ]
)


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

    def test_stockholm_annotations_are_parsed(self):
        alignments = list(stockholm.AlignmentIterator(DATA_STO))
        self.assertGreaterEqual(len(alignments), 1, "Stockholm file not parsed")
        aln = alignments[0]
        self.assertEqual(aln.annotations.get("identifier"), "TESTALIGN")
        self.assertEqual(aln.annotations.get("accession"), "PFTEST")
        self.assertEqual(
            aln.column_annotations.get("consensus secondary structure"),
            list("HHHH"),
        )
        first = aln.sequences[0]
        self.assertEqual(first.annotations.get("accession"), "ACC1")
        self.assertIn("secondary structure", first.letter_annotations)
        self.assertEqual(len(first.letter_annotations["secondary structure"]), 4)


if __name__ == "__main__":
    unittest.main()
