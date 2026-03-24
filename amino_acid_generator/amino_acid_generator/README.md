# 🧬 Amino Acid Generator

Real-time protein synthesis visualization that receives translation data from JCVI Genome Analyzer.

## Overview

This application connects to the JCVI Genome Analyzer via ZeroMQ and displays:
1. **Primary Structure** - Amino acid sequences as genes are being translated
2. **Secondary Structure** - Predicted protein folding (α-helix, β-sheet, turns, coils)

```
JCVI Genome Analyzer          Amino Acid Generator
┌─────────────────┐           ┌─────────────────┐
│  Translation    │  ZeroMQ   │  Receives data  │
│  Animation      │ ────────► │  Shows proteins │
│  (port 5555)    │           │  in real-time   │
└─────────────────┘           └─────────────────┘
```

## Features

- **Real-time connection** via ZeroMQ (port 5555)
- **Live amino acid translation** - CDS → Protein sequence
- **Progressive display** - Shows partial protein based on translation %
- **Three-letter format** - Met-Lys-Arg-Phe...
- **Two tabs:**
  - 🔤 **Primary Structure** - Amino acid sequence
  - 🌀 **Secondary Structure** - Chou-Fasman prediction
- **Two sections per tab:**
  - 🔄 Active Translation - Genes currently being translated
  - ✅ Completed Proteins - Fully translated proteins

## Secondary Structure Prediction

Uses the **Chou-Fasman method** (~65% accuracy) to predict:

| Symbol | Structure | Description |
|--------|-----------|-------------|
| **H** | α-Helix | Right-handed coiled structure |
| **E** | β-Sheet | Extended strand structure |
| **T** | Turn | Reverse direction |
| **C** | Coil | Random/unstructured |

### How It Works

Each amino acid has propensity values for forming different structures:
- High helix formers: Ala, Glu, Leu, Met
- High sheet formers: Val, Ile, Tyr, Phe
- Turn formers: Gly, Asn, Pro, Ser

## Requirements

```bash
pip install PyQt5 pyzmq
```

## Usage

1. **Start JCVI Genome Analyzer** and load your genome
2. **Start Amino Acid Generator:**
   ```bash
   python amino_acid_app.py
   ```
3. **Start Translation Animation** in JCVI Analyzer
4. Watch proteins appear in real-time! 🧬

## How It Works

1. JCVI Genome Analyzer broadcasts translation data via ZeroMQ
2. This app subscribes to the data stream
3. For each update:
   - Receives: gene_id, protein_level (%), cds_sequence
   - Calculates: How many amino acids to show based on %
   - Displays: Partial or complete protein sequence

## Example

When a gene is at 30% translation:
```
JCVISYN3_0294 [30%]: Met-Lys-Arg-Phe-Leu-Ser...
```

When translation completes (100%):
```
✅ JCVISYN3_0294: Met-Lys-Arg-Phe-Leu-Ser-Thr-Pro-Val-Ala... (88 aa)
```

## Data Format (JSON received)

```json
{
    "timestamp": 0.3,
    "gene_id": "JCVISYN3_0294",
    "status": "translating",
    "protein_level": 30.5,
    "strand": "+",
    "start": 176396,
    "end": 176663,
    "cds_sequence": "ATGAAACGA..."
}
```

## Genetic Code

The app uses the standard genetic code to translate codons:
- ATG → Met (Start)
- TAA, TAG, TGA → Stop
- All 64 codons mapped to 20 amino acids

## Troubleshooting

**"Not Connected" status:**
- Make sure JCVI Genome Analyzer is running
- Check that ZeroMQ publisher is active on port 5555
- Try clicking "Connect" button

**No data appearing:**
- Start the Translation Animation in JCVI Analyzer
- Check that genes have CDS sequences loaded

**PyZMQ not installed:**
```bash
pip install pyzmq
```
