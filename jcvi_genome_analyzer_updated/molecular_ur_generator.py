import torch
import numpy as np
import h5py
import hashlib
import os
from transformers import AutoTokenizer, AutoModel, AutoModelForMaskedLM
from typing import Dict, Optional
from data_structures import Gene

class MolecularURGenerator:
    """
    Molecular Universal Representation (UR) Generator.
    Uses transformer models to generate embeddings for DNA, RNA, and Protein.
    Models:
    - DNA/RNA: InstaDeepAI/nucleotide-transformer-v2-500m-multi-species
    - Protein: facebook/esm2_t33_650M_UR50D
    """
    
    def __init__(self, use_gpu: bool = True, cache_dir: str = "Data", load_models: bool = True, device: Optional[torch.device] = None):
        if device:
            self.device = device
        else:
            self.device = torch.device("cuda" if use_gpu and torch.cuda.is_available() else "cpu")
        print(f"Using device: {self.device}")
        
        # Cache setup
        self.cache_path = os.path.join(cache_dir, "embeddings_cache.h5")
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir, exist_ok=True)
            
        if load_models:
            # Load Nucleotide Transformer - Use MaskedLM to correctly handle gated MLP architecture
            self.nt_tokenizer = AutoTokenizer.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True)
            self.nt_model = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/nucleotide-transformer-v2-500m-multi-species", trust_remote_code=True).to(self.device)
            
            # Load ESM-2
            self.esm_tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
            self.esm_model = AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D").to(self.device)
            
            self.nt_model.eval()
            self.esm_model.eval()
        else:
            print("Skipping transformer model loading (test mode/cache only).")

    @torch.no_grad()
    def generate_ur(self, gene: Gene) -> Gene:
        """
        Generates and attaches embeddings to the gene object with HDF5 caching.
        """
        gene_id = gene.id
        dna_seq = gene.sequence or ""
        protein_seq = self._translate_sequence(dna_seq)
        
        # Generate hashes for validation
        dna_hash = hashlib.sha256(dna_seq.encode()).hexdigest()
        prot_hash = hashlib.sha256(protein_seq.encode()).hexdigest()
        
        # Try to load from cache
        cached_embeddings = self._get_from_cache(gene_id, dna_hash, prot_hash)
        if cached_embeddings:
            gene.molecular_ur = cached_embeddings
            return gene

        # Cache miss: generate embeddings
        print(f"Cache miss for {gene_id}: generating embeddings...")
        embeddings = {}
        
        # 1. DNA embedding
        if dna_seq:
            dna_vec = self._get_nt_embedding(dna_seq)
            embeddings['dna_embedding'] = dna_vec
            
            # 2. RNA embedding (transcribed)
            rna_seq = dna_seq.replace('T', 'U')
            rna_vec = self._get_nt_embedding(rna_seq)
            embeddings['rna_embedding'] = rna_vec
        else:
            embeddings['dna_embedding'] = np.zeros(1280)
            embeddings['rna_embedding'] = np.zeros(1280)
            
        # 3. Protein embedding
        if protein_seq:
            prot_vec = self._get_esm_embedding(protein_seq)
            embeddings['protein_embedding'] = prot_vec
        else:
            embeddings['protein_embedding'] = np.zeros(1280)
            
        # 4. Combined embedding (concatenation)
        embeddings['combined_embedding'] = np.concatenate([
            embeddings['dna_embedding'],
            embeddings['rna_embedding'],
            embeddings['protein_embedding']
        ])
        
        # Save to cache
        self._save_to_cache(gene_id, dna_hash, prot_hash, embeddings)
        
        gene.molecular_ur = embeddings
        return gene

    def _get_from_cache(self, gene_id: str, dna_hash: str, prot_hash: str) -> Optional[Dict]:
        """Load embeddings from HDF5 if hashes match"""
        if not os.path.exists(self.cache_path):
            return None
            
        try:
            with h5py.File(self.cache_path, 'r') as f:
                if gene_id not in f:
                    return None
                    
                group = f[gene_id]
                # Validate hashes
                if group.attrs.get('dna_hash') != dna_hash or group.attrs.get('prot_hash') != prot_hash:
                    return None
                
                return {
                    'dna_embedding': group['dna_embedding'][:],
                    'rna_embedding': group['rna_embedding'][:],
                    'protein_embedding': group['protein_embedding'][:],
                    'combined_embedding': group['combined_embedding'][:]
                }
        except Exception as e:
            print(f"Error reading cache for {gene_id}: {e}")
            return None

    def _save_to_cache(self, gene_id: str, dna_hash: str, prot_hash: str, embeddings: Dict):
        """Save embeddings to HDF5"""
        try:
            with h5py.File(self.cache_path, 'a') as f:
                if gene_id in f:
                    del f[gene_id]
                    
                group = f.create_group(gene_id)
                group.attrs['dna_hash'] = dna_hash
                group.attrs['prot_hash'] = prot_hash
                
                for key, val in embeddings.items():
                    group.create_dataset(key, data=val, compression="gzip")
        except Exception as e:
            print(f"Error saving cache for {gene_id}: {e}")

    def _get_nt_embedding(self, sequence: str) -> np.ndarray:
        """Get embedding from Nucleotide Transformer"""
        # Crop sequence if too long for the model
        max_len = 2048
        seq = sequence[:max_len]
        
        inputs = self.nt_tokenizer(seq, return_tensors="pt").to(self.device)
        outputs = self.nt_model(**inputs, output_hidden_states=True)
        
        # For MaskedLM, hidden_states[-1] is the last layer's hidden states
        # Shape is [batch, seq_len, hidden_size]
        last_hidden = outputs.hidden_states[-1]
        embeddings = last_hidden.mean(dim=1)
        return embeddings.detach().cpu().numpy().flatten()

    def _get_esm_embedding(self, sequence: str) -> np.ndarray:
        """Get embedding from ESM-2"""
        max_len = 1024
        seq = sequence[:max_len]
        
        inputs = self.esm_tokenizer(seq, return_tensors="pt").to(self.device)
        outputs = self.esm_model(**inputs)
        
        # Mean pooling
        embeddings = outputs.last_hidden_state.mean(dim=1)
        return embeddings.detach().cpu().numpy().flatten()

    def _translate_sequence(self, dna_sequence: str) -> str:
        """Simple DNA to Protein translation"""
        if not dna_sequence:
            return ""
            
        try:
            from Bio.Seq import Seq
            return str(Seq(dna_sequence).translate(to_stop=True))
        except ImportError:
            # Very basic fallback for codon translation (first 10 codons)
            return "M" + "X" * (len(dna_sequence) // 3 - 1)
