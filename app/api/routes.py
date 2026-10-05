from fastapi import APIRouter, HTTPException
from typing import Dict, List, Any
from app.models.sequence import SequenceRecord, SequenceType
from app.services.bio_analyzer import BioAnalyzerService

router = APIRouter(prefix="/api/v1/bio", tags=["Genomic Sequence Analytics"])

@router.post("/analyze", summary="Analyze DNA/RNA Sequence (GC-Content, Transcription, Translation)")
async def analyze_sequence(payload: Dict[str, Any]):
    """
    Analyzes raw genetic sequences for:
    - GC content percentage
    - mRNA transcription
    - Codon translation to amino acid protein sequence
    - Promoter motif identification
    """
    seq_str = payload.get("sequence", "")
    seq_id = payload.get("id", "SEQ-001")
    seq_type_str = payload.get("type", "DNA").upper()

    if not seq_str:
        raise HTTPException(status_code=400, detail="Sequence string cannot be empty")

    seq_type = SequenceType.DNA if seq_type_str == "DNA" else SequenceType.RNA
    record = SequenceRecord(id=seq_id, description="Client Query Sequence", sequence=seq_str, seq_type=seq_type)
    result = BioAnalyzerService.analyze_record(record)
    return result.__dict__

@router.post("/variants", summary="Scan Point Mutations & SNPs Against Reference")
async def scan_variants(payload: Dict[str, str]):
    """Identifies nucleotide transition and transversion variants."""
    ref = payload.get("reference", "")
    query = payload.get("query", "")
    if not ref or not query:
        raise HTTPException(status_code=400, detail="Both reference and query sequences are required")

    variants = BioAnalyzerService.detect_single_nucleotide_variants(ref, query)
    return {
        "variant_count": len(variants),
        "variants": [v.__dict__ for v in variants]
    }

@router.post("/align", summary="Pairwise Sequence Alignment Scoring")
async def align_sequences(payload: Dict[str, str]):
    """Computes pairwise sequence alignment score and percentage identity."""
    s1 = payload.get("seq1", "")
    s2 = payload.get("seq2", "")
    if not s1 or not s2:
        raise HTTPException(status_code=400, detail="Both sequences are required for alignment")

    result = BioAnalyzerService.align_pairwise_simple(s1, s2)
    return result.__dict__
