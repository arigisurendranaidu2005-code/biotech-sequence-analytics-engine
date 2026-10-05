from dataclasses import dataclass, field
from typing import List, Dict, Optional
from enum import Enum

class SequenceType(str, Enum):
    DNA = "DNA"
    RNA = "RNA"
    PROTEIN = "PROTEIN"

@dataclass
class SequenceRecord:
    id: str
    description: str
    sequence: str
    seq_type: SequenceType = SequenceType.DNA
    length: int = field(init=False)

    def __post_init__(self):
        self.sequence = self.sequence.upper().strip()
        self.length = len(self.sequence)

@dataclass
class MutationVariant:
    position: int
    reference_base: str
    mutant_base: str
    mutation_type: str  # TRANSITION, TRANSVERSION, INSERTION, DELETION

@dataclass
class GenomicAnalysisResult:
    sequence_id: str
    seq_type: str
    length: int
    gc_content_pct: float
    base_composition: Dict[str, int]
    transcribed_rna: Optional[str] = None
    translated_protein: Optional[str] = None
    open_reading_frames_count: int = 0
    motifs_found: List[Dict[str, int]] = field(default_factory=list)

@dataclass
class AlignmentResult:
    seq1_id: str
    seq2_id: str
    alignment_score: int
    matches: int
    mismatches: int
    identity_percentage: float
