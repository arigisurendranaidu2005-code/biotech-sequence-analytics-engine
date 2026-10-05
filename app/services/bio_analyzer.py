import re
from typing import Dict, List, Tuple, Optional
from app.models.sequence import (
    SequenceRecord,
    SequenceType,
    MutationVariant,
    GenomicAnalysisResult,
    AlignmentResult
)

class BioAnalyzerService:
    """
    Bioinformatics sequence analysis algorithms and pipeline routines.
    Implements standard NCBI codon tables, GC calculation, SNP discovery,
    and pairwise sequence alignment scoring.
    """

    # Standard genetic code dictionary (RNA Codons -> Amino Acids)
    CODON_TABLE: Dict[str, str] = {
        'AUG': 'M',  # Start codon (Methionine)
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',  # Stop codons
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        """Calculates percentage of G and C nucleotides in DNA/RNA."""
        cleaned = sequence.upper().strip()
        if not cleaned:
            return 0.0
        gc_count = cleaned.count('G') + cleaned.count('C')
        return round((gc_count / len(cleaned)) * 100, 2)

    @staticmethod
    def transcribe(dna_seq: str) -> str:
        """Transcribes DNA to mRNA by replacing Thymine (T) with Uracil (U)."""
        return dna_seq.upper().replace('T', 'U')

    @classmethod
    def translate(cls, mrna_seq: str) -> str:
        """Translates mRNA sequence into protein peptide string using standard codon table."""
        mrna = mrna_seq.upper().replace('T', 'U')
        protein = []
        for i in range(0, len(mrna) - 2, 3):
            codon = mrna[i:i+3]
            amino_acid = cls.CODON_TABLE.get(codon, 'X')
            if amino_acid == '*':
                break  # Stop codon encountered
            protein.append(amino_acid)
        return "".join(protein)

    @staticmethod
    def find_motifs(sequence: str, motif: str) -> List[int]:
        """Finds 1-based start positions of a target sequence motif."""
        positions = []
        seq = sequence.upper()
        target = motif.upper()
        start = 0
        while True:
            idx = seq.find(target, start)
            if idx == -1:
                break
            positions.append(idx + 1)
            start = idx + 1
        return positions

    @staticmethod
    def detect_single_nucleotide_variants(ref_seq: str, query_seq: str) -> List[MutationVariant]:
        """Identifies point mutations between reference and query sequences."""
        variants = []
        min_len = min(len(ref_seq), len(query_seq))
        purines = {'A', 'G'}
        pyrimidines = {'C', 'T'}

        for i in range(min_len):
            r_base = ref_seq[i].upper()
            q_base = query_seq[i].upper()
            if r_base != q_base:
                if (r_base in purines and q_base in purines) or (r_base in pyrimidines and q_base in pyrimidines):
                    m_type = "TRANSITION"
                else:
                    m_type = "TRANSVERSION"

                variants.append(MutationVariant(
                    position=i + 1,
                    reference_base=r_base,
                    mutant_base=q_base,
                    mutation_type=m_type
                ))
        return variants

    @staticmethod
    def align_pairwise_simple(seq1: str, seq2: str) -> AlignmentResult:
        """Calculates global match/mismatch alignment score between two sequences."""
        s1 = seq1.upper()
        s2 = seq2.upper()
        min_len = min(len(s1), len(s2))
        matches = 0
        mismatches = 0

        for i in range(min_len):
            if s1[i] == s2[i]:
                matches += 1
            else:
                mismatches += 1

        length_diff = abs(len(s1) - len(s2))
        score = (matches * 2) - (mismatches * 1) - (length_diff * 2)
        total_eval = max(len(s1), len(s2))
        identity_pct = round((matches / total_eval) * 100, 2) if total_eval > 0 else 0.0

        return AlignmentResult(
            seq1_id="SeqA",
            seq2_id="SeqB",
            alignment_score=score,
            matches=matches,
            mismatches=mismatches + length_diff,
            identity_percentage=identity_pct
        )

    @classmethod
    def analyze_record(cls, record: SequenceRecord) -> GenomicAnalysisResult:
        seq = record.sequence
        gc = cls.calculate_gc_content(seq)
        composition = {b: seq.count(b) for b in set(seq)}

        rna = cls.transcribe(seq) if record.seq_type == SequenceType.DNA else seq
        protein = cls.translate(rna)

        # Common promoter TATA box motif search
        tata_positions = cls.find_motifs(seq, "TATAAA")
        motifs = [{"motif": "TATAAA (TATA-box)", "position": p} for p in tata_positions]

        return GenomicAnalysisResult(
            sequence_id=record.id,
            seq_type=record.seq_type.value,
            length=record.length,
            gc_content_pct=gc,
            base_composition=composition,
            transcribed_rna=rna,
            translated_protein=protein,
            open_reading_frames_count=1 if "M" in protein else 0,
            motifs_found=motifs
        )
