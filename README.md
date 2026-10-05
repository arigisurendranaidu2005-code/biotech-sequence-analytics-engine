# BioTech Genomic Sequence Analytics Pipeline

[![Python](https://img.shields.io/badge/Python-3.11%2B-blue.svg)](https://www.python.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.110%2B-009688.svg)](https://fastapi.tiangolo.com/)
[![Tests](https://img.shields.io/badge/Tests-100%25%20Passing-brightgreen.svg)](https://github.com/arigisurendranaidu2005-code/biotech-sequence-analytics-engine)
[![BioPython](https://img.shields.io/badge/Bioinformatics-BioPython%20Compatible-orange.svg)](https://biopython.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An algorithmic bioinformatics analytics platform leveraging Python modules for genomic sequence transcription, codon translation, SNP variant detection, promoter motif scanning, and global pairwise alignment scoring.

---

## 🧬 System Architecture & Design (Designing)

The architecture is engineered for high throughput genomic data processing:

```
[ FASTA Sequences / Genomic Queries ]
                 │
                 ▼
    ┌─────────────────────────┐
    │     FastAPI Gateway     │
    └────────────┬────────────┘
                 │
  ┌──────────────┼──────────────┬──────────────┐
  │              │              │              │
  ▼              ▼              ▼              ▼
┌──────────────┐┌──────────────┐┌──────────────┐┌──────────────┐
│  Nucleotide  ││ Transcription││ SNP Mutation ││ Pairwise     │
│  GC-Content  ││ & Translation││ Variant      ││ Sequence     │
│  Calculator  ││ (Codon Table)││ Discovery    ││ Alignment    │
└──────────────┘└──────────────┘└──────────────┘└──────────────┘
```

---

## 💻 Core Capabilities (Coding)

1. **GC-Content & Base Frequency Analysis:** Accurately computes nucleotide proportions to assess genomic thermal stability and coding density.
2. **Central Dogma Workflow:** Asynchronously transcribes DNA to mRNA and translates codons into amino acid peptide sequences via standard NCBI genetic codon mappings.
3. **SNP & Point Mutation Scanner:** Classifies nucleotide substitutions into transition (purine-to-purine / pyrimidine-to-pyrimidine) and transversion mutations.
4. **Motif Discovery Engine:** Fast pattern matching locating promoter consensus sequences (e.g., TATA-box motifs).
5. **Pairwise Sequence Alignment:** Dynamic scoring measuring match density, mismatches, and sequence identity percentages.

---

## 🧪 Testing & Quality Assurance (Testing)

The bioinformatics pipeline is verified through automated unit tests with 100% test pass rates:

```bash
# Execute automated test suite
python -m unittest discover -s tests
```

### Verified Test Cases:
- ✅ `test_gc_content_calculation`: Verifies stoichiometric G/C base percentage calculations.
- ✅ `test_transcription_dna_to_rna`: Validates Thymine-to-Uracil enzymatic transcription simulation.
- ✅ `test_translation_to_protein`: Validates start/stop codon parsing into correct peptide sequences.
- ✅ `test_motif_search`: Confirms exact 1-based index resolution for regulatory motifs.
- ✅ `test_single_nucleotide_variant_detection`: Validates SNP variant classification (transition vs transversion).
- ✅ `test_pairwise_alignment_score`: Checks identity scoring and gap penalties.
- ✅ `test_full_record_analysis`: Validates end-to-end multi-step analysis payloads.

---

## 📖 Technical Documentation & API Guide (Content Writing)

Interactive OpenAPI / Swagger and ReDoc documentation are automatically exposed:
- **Swagger Documentation:** `http://localhost:8002/docs`
- **ReDoc Documentation:** `http://localhost:8002/redoc`

### Sample Request: Sequence Analysis
```json
POST /api/v1/bio/analyze
Content-Type: application/json

{
  "id": "GENE-BRCA1-EXON2",
  "sequence": "ATGCGATCGTATAAAGCC",
  "type": "DNA"
}
```

### Sample Response:
```json
{
  "sequence_id": "GENE-BRCA1-EXON2",
  "seq_type": "DNA",
  "length": 18,
  "gc_content_pct": 44.44,
  "base_composition": {"A": 6, "T": 4, "G": 4, "C": 4},
  "transcribed_rna": "AUGCGAUCGUAUAAAGCC",
  "translated_protein": "MRSYKA",
  "open_reading_frames_count": 1,
  "motifs_found": [
    {"motif": "TATAAA (TATA-box)", "position": 10}
  ]
}
```

---

## 🚀 Quickstart

```bash
# Clone the repository
git clone https://github.com/arigisurendranaidu2005-code/biotech-sequence-analytics-engine.git
cd biotech-sequence-analytics-engine

# Run the test suite
python -m unittest discover -s tests

# Start the API service
python -m app.main
```

---

## 📜 License
MIT License. Authored by Surendra Naidu Arigi.
