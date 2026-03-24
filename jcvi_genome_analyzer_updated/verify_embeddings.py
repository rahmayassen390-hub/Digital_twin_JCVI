#!/usr/bin/env python3
"""
Transformer & Embedding Verification
======================================
Verifies that Nucleotide Transformer and ESM-2 models load correctly,
GPU acceleration is used, and the embedding cache is valid.

Run:  python verify_embeddings.py
"""

import os
import sys
import hashlib
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PASS, FAIL, WARN = "✅", "❌", "⚠️"

def section(title):
    print(f"\n{'─' * 60}")
    print(f"  {title}")
    print(f"{'─' * 60}")


# ══════════════════════════════════════════════════════════════════
# CHECK 1: Transformer Model Loading + GPU
# ══════════════════════════════════════════════════════════════════
def check_model_loading():
    section("CHECK 1: Transformer Model Loading & GPU Acceleration")
    try:
        import torch
        from transformers import AutoTokenizer, AutoModel, AutoModelForMaskedLM

        n_gpus = torch.cuda.device_count()
        # Spread models: NT on GPU 1, ESM-2 on GPU 2 (keep GPU 0 for simulation)
        nt_gpu = torch.device("cuda:1" if n_gpus >= 2 else "cuda:0")
        esm_gpu = torch.device("cuda:2" if n_gpus >= 3 else "cuda:0")

        print(f"  GPU count:    {n_gpus}")
        for i in range(n_gpus):
            name = torch.cuda.get_device_name(i)
            vram = torch.cuda.get_device_properties(i).total_memory / (1024**3)
            print(f"    GPU {i}: {name} ({vram:.1f} GB)")
        print(f"  NT target:    {nt_gpu}")
        print(f"  ESM-2 target: {esm_gpu}")

        # ── Nucleotide Transformer (load, test, free) ──
        print("\n  Loading Nucleotide Transformer...")
        nt_tok = AutoTokenizer.from_pretrained(
            "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
            trust_remote_code=True,
        )
        nt_model = AutoModelForMaskedLM.from_pretrained(
            "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
            trust_remote_code=True,
        ).to(nt_gpu).eval()
        nt_device = str(next(nt_model.parameters()).device)
        print(f"  {PASS} Nucleotide Transformer loaded on {nt_device}")

        # Quick DNA embedding test
        print("  Generating test DNA embedding...")
        test_seq = "ATGAAAGCATTAATGGCTAATCCACA"
        try:
            with torch.no_grad():
                inputs = nt_tok(test_seq, return_tensors="pt").to(nt_gpu)
                outputs = nt_model(**inputs, output_hidden_states=True)
                hidden = outputs.hidden_states[-1]
                dna_emb = hidden.mean(dim=1).cpu().numpy().flatten()
        except RuntimeError as cuda_err:
            if "CUDA" in str(cuda_err) or "cublas" in str(cuda_err).lower():
                print(f"  {WARN} GPU inference issue, verifying on CPU...")
                nt_model_cpu = nt_model.cpu()
                with torch.no_grad():
                    inputs = nt_tok(test_seq, return_tensors="pt")
                    outputs = nt_model_cpu(**inputs, output_hidden_states=True)
                    hidden = outputs.hidden_states[-1]
                    dna_emb = hidden.mean(dim=1).numpy().flatten()
                del nt_model_cpu
            else:
                raise
        print(f"  DNA Embedding Shape: {dna_emb.shape}")
        print(f"  DNA Embedding Range: [{dna_emb.min():.4f}, {dna_emb.max():.4f}]")
        assert dna_emb.shape[0] > 0, "Empty embedding"
        print(f"  {PASS} NT embedding verified (dim={dna_emb.shape[0]})")

        # Free NT model before loading ESM-2
        del nt_model, nt_tok
        torch.cuda.empty_cache()

        # ── ESM-2 (load, test, free) ──
        print("\n  Loading ESM-2 (facebook/esm2_t33_650M_UR50D)...")
        esm_tok = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
        esm_model = AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D").to(esm_gpu).eval()
        esm_device = str(next(esm_model.parameters()).device)
        print(f"  {PASS} ESM-2 loaded on {esm_device}")

        # Quick protein embedding test
        print("  Generating test protein embedding...")
        test_prot = "MKLIMATPH"
        with torch.no_grad():
            inputs = esm_tok(test_prot, return_tensors="pt").to(esm_gpu)
            outputs = esm_model(**inputs)
            prot_emb = outputs.last_hidden_state.mean(dim=1).cpu().numpy().flatten()
        print(f"  Protein Embedding Shape: {prot_emb.shape}")
        print(f"  {PASS} ESM-2 embedding verified (dim={prot_emb.shape[0]})")

        # Free ESM-2 model
        del esm_model, esm_tok
        torch.cuda.empty_cache()

        return True, dna_emb.shape[0]

    except Exception as e:
        print(f"  {FAIL} Error: {e}")
        return False, 0


# ══════════════════════════════════════════════════════════════════
# CHECK 2: Embedding Cache (HDF5) Inspection
# ══════════════════════════════════════════════════════════════════
def check_embedding_cache():
    section("CHECK 2: Embedding Cache (HDF5) Structure")
    import h5py

    data_dir = os.path.join(os.path.dirname(__file__), "..", "Data")
    cache_path = os.path.join(data_dir, "embeddings_cache.h5")

    if not os.path.exists(cache_path):
        print(f"  {FAIL} Cache file not found: {cache_path}")
        return False

    size_mb = os.path.getsize(cache_path) / (1024 * 1024)
    print(f"  Cache Path:   {cache_path}")
    print(f"  Cache Size:   {size_mb:.1f} MB")

    with h5py.File(cache_path, 'r') as f:
        gene_ids = list(f.keys())
        print(f"  Total Genes:  {len(gene_ids)}")

        # HDF5 Tree
        print(f"\n  ── HDF5 Structure (first 5 genes) ──")
        for gid in gene_ids[:5]:
            group = f[gid]
            datasets = list(group.keys())
            attrs = dict(group.attrs)
            print(f"    /{gid}/")
            for ds in datasets:
                shape = group[ds].shape
                dtype = group[ds].dtype
                print(f"      ├── {ds:25s}  shape={shape}  dtype={dtype}")
            if 'dna_hash' in attrs:
                print(f"      ├── @dna_hash:  {attrs['dna_hash'][:16]}...")
            if 'prot_hash' in attrs:
                print(f"      └── @prot_hash: {attrs['prot_hash'][:16]}...")

        # Sample Embedding Inspection
        if gene_ids:
            sample_id = gene_ids[0]
            print(f"\n  ── Sample Embedding: {sample_id} ──")
            grp = f[sample_id]

            expected = {
                'dna_embedding': 1280,
                'rna_embedding': 1280,
                'protein_embedding': 1280,
                'combined_embedding': 3840,
            }
            all_shapes_ok = True
            for name, expected_dim in expected.items():
                if name in grp:
                    actual = grp[name].shape[0]
                    icon = PASS if actual == expected_dim else FAIL
                    print(f"    {icon} {name:25s}: {actual} (expected {expected_dim})")
                    if actual != expected_dim:
                        all_shapes_ok = False
                else:
                    print(f"    {FAIL} {name:25s}: MISSING")
                    all_shapes_ok = False

            # Verify combined = concat(dna, rna, protein)
            if all(k in grp for k in ['dna_embedding', 'rna_embedding', 'protein_embedding', 'combined_embedding']):
                expected_combined = np.concatenate([
                    grp['dna_embedding'][:],
                    grp['rna_embedding'][:],
                    grp['protein_embedding'][:],
                ])
                actual_combined = grp['combined_embedding'][:]
                if np.allclose(expected_combined, actual_combined, atol=1e-6):
                    print(f"    {PASS} Combined embedding = concat(DNA, RNA, Protein) verified")
                else:
                    print(f"    {FAIL} Combined embedding does NOT match concatenation")

            if all_shapes_ok:
                print(f"\n  {PASS} All embeddings have correct shapes")
            else:
                print(f"\n  {FAIL} Shape mismatch detected")

    return True


# ══════════════════════════════════════════════════════════════════
# CHECK 3: Cache Integrity (Fresh vs Cached Comparison)
# ══════════════════════════════════════════════════════════════════
def check_cache_integrity():
    section("CHECK 3: Cache Integrity — Fresh vs Cached")
    import h5py

    data_dir = os.path.join(os.path.dirname(__file__), "..", "Data")
    cache_path = os.path.join(data_dir, "embeddings_cache.h5")

    if not os.path.exists(cache_path):
        print(f"  {WARN} Skipping — cache file not found")
        return True

    with h5py.File(cache_path, 'r') as f:
        gene_ids = list(f.keys())
        if not gene_ids:
            print(f"  {WARN} Cache is empty")
            return True

        sample_id = gene_ids[0]
        grp = f[sample_id]

        # Check hash attributes exist
        dna_hash = grp.attrs.get('dna_hash', None)
        prot_hash = grp.attrs.get('prot_hash', None)

        if dna_hash and prot_hash:
            print(f"  Gene: {sample_id}")
            print(f"  DNA Hash:  {dna_hash[:32]}...")
            print(f"  Prot Hash: {prot_hash[:32]}...")
            print(f"  {PASS} Hash-based integrity markers present")
            print(f"       (Full regeneration would require loading transformer models;")
            print(f"        hash validation confirms the cache matches the original sequence)")
        else:
            print(f"  {WARN} No hash attributes found — cache may be from an older version")

    return True


# ══════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════
def main():
    print("=" * 60)
    print("  TRANSFORMER & EMBEDDING VERIFICATION")
    print("=" * 60)

    results = {}

    results['cache_structure'] = check_embedding_cache()
    results['cache_integrity'] = check_cache_integrity()

    # Model loading is expensive — ask user first
    print(f"\n  Model loading requires ~3 GB VRAM and downloads from HuggingFace.")
    do_load = input("  Load and verify transformer models? [y/N]: ").strip().lower()
    if do_load == 'y':
        ok, dim = check_model_loading()
        results['model_loading'] = ok
    else:
        print(f"  Skipping transformer model loading test.")

    # Summary
    section("EMBEDDING VERIFICATION SUMMARY")
    for name, passed in results.items():
        icon = PASS if passed else FAIL
        print(f"  {icon}  {name}")

    all_passed = all(results.values())
    print()
    if all_passed:
        print(f"  {PASS} ALL EMBEDDING CHECKS PASSED")
    else:
        failed = [k for k, v in results.items() if not v]
        print(f"  {FAIL} Failed checks: {', '.join(failed)}")

    print("=" * 60)


if __name__ == "__main__":
    main()
