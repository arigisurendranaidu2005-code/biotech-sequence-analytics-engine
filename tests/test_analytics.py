import unittest
from app.models.sequence import SequenceRecord, SequenceType
from app.services.bio_analyzer import BioAnalyzerService

class TestBioAnalyzerService(unittest.TestCase):

    def test_gc_content_calculation(self):
        # 5 G/C out of 10 bases = 50.0%
        dna = "ATGCGATCGT"
        gc = BioAnalyzerService.calculate_gc_content(dna)
        self.assertEqual(gc, 50.0)

    def test_transcription_dna_to_rna(self):
        dna = "ATGCGATCGT"
        rna = BioAnalyzerService.transcribe(dna)
        self.assertEqual(rna, "AUGCGAUCGU")

    def test_translation_to_protein(self):
        # AUG = M (Methionine), GUC = V (Valine), UAA = Stop (*)
        mrna = "AUGGUCUAA"
        protein = BioAnalyzerService.translate(mrna)
        self.assertEqual(protein, "MV")

    def test_motif_search(self):
        dna = "CCCTATAAAGGGCTATAAACCC"
        positions = BioAnalyzerService.find_motifs(dna, "TATAAA")
        self.assertEqual(positions, [4, 14])

    def test_single_nucleotide_variant_detection(self):
        ref = "ATGC"
        query = "ACGT"  # Pos 2: T->C (Transition pyrimidine to pyrimidine), Pos 3: G->G, Pos 4: C->T (Transition)
        variants = BioAnalyzerService.detect_single_nucleotide_variants(ref, query)
        self.assertEqual(len(variants), 2)
        self.assertEqual(variants[0].position, 2)
        self.assertEqual(variants[0].reference_base, "T")
        self.assertEqual(variants[0].mutant_base, "C")
        self.assertEqual(variants[0].mutation_type, "TRANSITION")

    def test_pairwise_alignment_score(self):
        s1 = "ATGC"
        s2 = "ATGC"
        result = BioAnalyzerService.align_pairwise_simple(s1, s2)
        self.assertEqual(result.identity_percentage, 100.0)
        self.assertEqual(result.matches, 4)
        self.assertEqual(result.mismatches, 0)

    def test_full_record_analysis(self):
        record = SequenceRecord(
            id="SEQ-GENE-101",
            description="Synthetic Test Gene",
            sequence="ATGTATAAAGGG",
            seq_type=SequenceType.DNA
        )
        res = BioAnalyzerService.analyze_record(record)
        self.assertEqual(res.sequence_id, "SEQ-GENE-101")
        self.assertEqual(res.length, 12)
        self.assertGreater(res.gc_content_pct, 0)
        self.assertEqual(len(res.motifs_found), 1)

if __name__ == '__main__':
    unittest.main()
