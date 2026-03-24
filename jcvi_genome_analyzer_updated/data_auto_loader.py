#!/usr/bin/env python3
"""
Data Auto-Discovery Loader
===========================
Automatically discovers and maps required genomic data files from a
data directory, eliminating the need for manual file selection.

Usage:
    from data_auto_loader import auto_discover_data, DataManifest
    manifest = auto_discover_data("/path/to/Data")
"""

import os
import sys
import glob
from dataclasses import dataclass, field
from typing import Optional, List, Dict


@dataclass
class DataManifest:
    """Complete mapping of discovered data files."""

    # Required files (will raise error if missing)
    fasta_path: Optional[str] = None
    gff3_path: Optional[str] = None

    # BED files (optional but enriching)
    transcription_bed_path: Optional[str] = None
    translation_bed_path: Optional[str] = None
    likely_promoters_path: Optional[str] = None
    internal_operons_path: Optional[str] = None
    regulation_bed_path: Optional[str] = None
    complex_cases_path: Optional[str] = None

    # Enrichment datasets
    enrichment_files: Dict[str, str] = field(default_factory=dict)

    # Auto-generated model files
    embeddings_cache: Optional[str] = None
    perturbation_dataset: Optional[str] = None
    ai_model_weights: Optional[str] = None

    # Discovery log
    discovered: List[str] = field(default_factory=list)
    missing_required: List[str] = field(default_factory=list)
    missing_optional: List[str] = field(default_factory=list)


# JCVI file patterns by version
VERSIONED_PATTERNS = {
    "1.0": {
        "fasta":           ("JCVI_1.0.fasta",                         True ),
        "gff3":            ("JCVI_1.0.gff3",                          True ),
        "transcription":   ("JCVI_Transcription_Dynamics_1.0#.bed",   False),
        "translation":     ("JCVI_Translation_Dynamics_1.0#.bed",     False),
        "promoters":       ("Likely_Promoters_1.0#.bed",              False),
        "operons":         ("Internal_Operons_1.0#.bed",              False),
        "regulation":      ("JCVI_Gene_Regulation_1.0#.bed",          False),
        "complex_cases":   ("Complex_Cases_1.0#.bed",                 False),
    },
    "3.0": {
        "fasta":           ("JCVI_3.0.fasta",                         True ),
        "gff3":            ("JCVI_3.0.gff3",                          True ),
        "transcription":   ("JCVI_Transcription_Dynamics_3.0#.bed",   False),
        "translation":     ("JCVI_Translation_Dynamics_3.0#.bed",     False),
        "promoters":       ("Likely_Promoters_3.0#.bed",              False),
        "operons":         ("Internal_Operons_3.0#.bed",              False),
        "regulation":      ("JCVI_Gene_Regulation_3.0#.bed",          False),
        "complex_cases":   ("Complex_Cases_3.0#.bed",                 False),
    }
}

ENRICHMENT_FILES = {
    "kinetic_params":      "Kinetic_Parameters.tsv",
    "mrna_counts":         "mRNA_counts.csv",
    "proteomics":          "proteomics_data.xlsx",
    "full_genome":         "JCVI_syn3A_full_genome_dataset.csv",
    "degradation_rates":   "syn3A_degradation_rates.csv",
    "orthologs_rbs":       "syn3A_gene_orthologs_rbs.csv",
    "global_params":       "global_parameters.csv",
    "transport_fluxes":    "transport_fluxes.tsv",
    "morphology":          "syn3A_morphology_params.csv",
    "lipid_data":          "syn3A_lipid_data.csv",
    "kinetic_summary":     "syn3A_kinetic_parameters.csv",
    "stress_response":     "syn3A_stress_response_data.csv",
    "division_genes":      "syn3A_division_genes_data.csv",
    "gene_dataset":        "syn3A_gene_dataset.csv",
    "metabolic_model":     "metabolic_model_iMB155.xml",
    "metabolic_recon":     "metabolic_reconstruction.xlsx",
}

AUTO_GENERATED = {
    "embeddings_cache":     "embeddings_cache.h5",
    "perturbation_dataset": "perturbation_dataset_50k.h5",
}


def auto_discover_data(data_dir: str = None, genome_version: str = "3.0") -> DataManifest:
    """
    Auto-discover all required and optional data files.

    Args:
        data_dir: Path to the data directory.  If None, searches common
                  locations relative to the script directory.
        genome_version: "3.0" or "1.0" to select the right file set.

    Returns:
        DataManifest with all discovered paths.
    """
    if data_dir is None:
        data_dir = _find_data_dir()

    manifest = DataManifest()

    print("┌─────────────────────────────────────────────────────────┐")
    print(f"│      DATA AUTO-DISCOVERY ENGINE v1.1 - [JCVI {genome_version}]      │")
    print("└─────────────────────────────────────────────────────────┘")
    print(f"  Scanning: {os.path.abspath(data_dir)}")
    print()

    if not os.path.isdir(data_dir):
        manifest.missing_required.append(f"Data directory not found: {data_dir}")
        print(f"  ❌ Data directory not found: {data_dir}")
        return manifest

    # Get patterns for this version
    patterns = VERSIONED_PATTERNS.get(genome_version, VERSIONED_PATTERNS["3.0"])

    # ── Core Genomic Files ────────────────────────────────────────
    print("  ── Core Genomic Files ──")
    _discover_core(data_dir, manifest, patterns)

    # ── BED Annotation Files ──────────────────────────────────────
    print("\n  ── BED Annotation Files ──")
    _discover_beds(data_dir, manifest, patterns)

    # ── Enrichment Datasets ───────────────────────────────────────
    # Enrichment files are currently shared or specific to 3.0
    print("\n  ── Enrichment Datasets ──")
    _discover_enrichment(data_dir, manifest)

    # ── Auto-Generated Artifacts ──────────────────────────────────
    print("\n  ── Auto-Generated Artifacts ──")
    _discover_autogen(data_dir, manifest)

    # ── Summary ───────────────────────────────────────────────────
    print("\n  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
    print(f"  Total discovered: {len(manifest.discovered)}")
    if manifest.missing_required:
        print(f"  ❌ MISSING REQUIRED ({len(manifest.missing_required)}):")
        for m in manifest.missing_required:
            print(f"       • {m}")
    if manifest.missing_optional:
        print(f"  ⚠️  Missing optional: {len(manifest.missing_optional)}")
    print()

    return manifest


def _find_data_dir() -> str:
    """Search common locations for the Data directory."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, "..", "Data"),
        os.path.join(script_dir, "Data"),
        os.path.join(script_dir, "data"),
        os.path.join(os.getcwd(), "Data"),
        os.path.join(os.getcwd(), "..", "Data"),
    ]
    for c in candidates:
        if os.path.isdir(c):
            return os.path.abspath(c)
    return os.path.join(script_dir, "..", "Data")


def _discover_core(data_dir: str, manifest: DataManifest, patterns: dict):
    fasta = os.path.join(data_dir, patterns["fasta"][0])
    gff3 = os.path.join(data_dir, patterns["gff3"][0])

    if os.path.exists(fasta):
        manifest.fasta_path = fasta
        manifest.discovered.append(fasta)
        print(f"    ✅ FASTA: {os.path.basename(fasta)}")
    else:
        manifest.missing_required.append(patterns["fasta"][0])
        print(f"    ❌ FASTA: {patterns['fasta'][0]} NOT FOUND")

    if os.path.exists(gff3):
        manifest.gff3_path = gff3
        manifest.discovered.append(gff3)
        print(f"    ✅ GFF3:  {os.path.basename(gff3)}")
    else:
        manifest.missing_required.append(patterns["gff3"][0])
        print(f"    ❌ GFF3:  {patterns['gff3'][0]} NOT FOUND")


def _discover_beds(data_dir: str, manifest: DataManifest, patterns: dict):
    bed_map = {
        "transcription":  "transcription_bed_path",
        "translation":    "translation_bed_path",
        "promoters":      "likely_promoters_path",
        "operons":        "internal_operons_path",
        "regulation":     "regulation_bed_path",
        "complex_cases":  "complex_cases_path",
    }
    for key, attr in bed_map.items():
        filename = patterns[key][0]
        path = os.path.join(data_dir, filename)
        if os.path.exists(path):
            setattr(manifest, attr, path)
            manifest.discovered.append(path)
            print(f"    ✅ {key:20s}: {filename}")
        else:
            manifest.missing_optional.append(filename)
            print(f"    ⚠️  {key:20s}: {filename} (not found)")


def _discover_enrichment(data_dir: str, manifest: DataManifest):
    for key, filename in ENRICHMENT_FILES.items():
        path = os.path.join(data_dir, filename)
        if os.path.exists(path):
            manifest.enrichment_files[key] = path
            manifest.discovered.append(path)
            print(f"    ✅ {key:22s}: {filename}")
        else:
            manifest.missing_optional.append(filename)
            print(f"    ⚠️  {key:22s}: {filename} (not found)")


def _discover_autogen(data_dir: str, manifest: DataManifest):
    for key, filename in AUTO_GENERATED.items():
        path = os.path.join(data_dir, filename)
        if os.path.exists(path):
            setattr(manifest, key, path)
            manifest.discovered.append(path)
            size_mb = os.path.getsize(path) / (1024 * 1024)
            print(f"    ✅ {key:22s}: {filename} ({size_mb:.1f} MB)")
        else:
            print(f"    ⚠️  {key:22s}: {filename} (not generated yet)")

    # AI model weights live in the code directory, not Data
    script_dir = os.path.dirname(os.path.abspath(__file__))
    model_path = os.path.join(script_dir, "ai_virtual_instrument_best.pt")
    if os.path.exists(model_path):
        manifest.ai_model_weights = model_path
        manifest.discovered.append(model_path)
        size_mb = os.path.getsize(model_path) / (1024 * 1024)
        print(f"    ✅ {'ai_model_weights':22s}: ai_virtual_instrument_best.pt ({size_mb:.1f} MB)")
    else:
        print(f"    ⚠️  {'ai_model_weights':22s}: not trained yet")


# ── Standalone test ───────────────────────────────────────────────
if __name__ == "__main__":
    manifest = auto_discover_data()
    if manifest.missing_required:
        print("\n⛔ Cannot proceed: required files are missing.")
        sys.exit(1)
    else:
        print("\n✅ All required files discovered. Ready to launch.")
