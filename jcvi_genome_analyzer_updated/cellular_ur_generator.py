"""
Cellular Universal Representation (UR) Generator
=================================================

Aggregates gene-level Molecular URs into a single cell-state vector
that represents the complete state of the virtual cell.

Components:
    1. Gene Aggregation Layer - attention-weighted mean of 498 gene embeddings
    2. Temporal Dynamics - ODE engine state (mRNA/protein concentrations)
    3. Environmental Context - from global_parameters.csv
    4. Resource State - ribosome/polymerase/ATP saturation

Output: CellStateVector (1888-dim)
"""

import os
import numpy as np
import pandas as pd
from typing import Dict, Optional
from data_structures import Gene, CellStateVector


class CellularURGenerator:
    """
    Generates the Cellular Universal Representation by integrating:
    - All gene Molecular URs (from MolecularURGenerator)
    - Current ODE simulation state (from CellSimulationEngine)
    - Environmental parameters (from global_parameters.csv)
    - Cellular resource state (from CellularResources)
    """
    
    # Target dimensions for each component
    GENE_EMBED_DIM = 1280
    DYNAMICS_DIM = 512
    ENV_DIM = 64
    RESOURCE_DIM = 32
    COMBINED_DIM = GENE_EMBED_DIM + DYNAMICS_DIM + ENV_DIM + RESOURCE_DIM  # 1888
    
    def __init__(self, data_dir: str = "Data"):
        self.data_dir = data_dir
        self.env_params = self._load_environment_params()
    
    def _load_environment_params(self) -> dict:
        """Load global environmental parameters"""
        path = os.path.join(self.data_dir, "global_parameters.csv")
        params = {
            'temperature': 37.0,    # °C (Mycoplasma optimal)
            'ph': 7.2,
            'osmolarity': 300.0,    # mOsm
            'doubling_time': 180.0, # minutes (~3hr for syn3A)
            'growth_rate': 0.0039,  # 1/min (ln2/180)
            'atp_concentration': 5.0,  # mM
            'glucose_available': 1.0,  # normalized
            'oxygen_available': 0.0,   # anaerobic (Mycoplasma)
        }
        
        if os.path.exists(path):
            try:
                df = pd.read_csv(path)
                for _, row in df.iterrows():
                    param_name = str(row.get('parameter', '')).lower().replace(' ', '_')
                    if param_name in params:
                        params[param_name] = float(row.get('value', params[param_name]))
                print(f"Loaded {len(df)} global parameters from {path}")
            except Exception as e:
                print(f"Warning: Could not load global_parameters.csv: {e}")
        
        return params
    
    def generate_cellular_ur(
        self,
        enriched_genes: Dict[str, Gene],
        simulation_engine=None,
        cellular_resources=None,
        simulation_time: float = 0.0,
    ) -> CellStateVector:
        """
        Generate the complete Cellular UR from all components.
        
        Args:
            enriched_genes: Dict[gene_id, Gene] with molecular_ur populated
            simulation_engine: CellSimulationEngine (optional, for dynamics)
            cellular_resources: CellularResources (optional, for resource state)
            simulation_time: Current simulation time in seconds
            
        Returns:
            CellStateVector with all components populated
        """
        # 1. Gene Aggregation Layer
        gene_embed = self._aggregate_gene_embeddings(enriched_genes)
        
        # 2. Temporal Dynamics Integration
        dynamics = self._extract_dynamics(enriched_genes, simulation_engine)
        dynamics = np.nan_to_num(dynamics.astype(np.float32))
        
        # 3. Environmental Context Vector
        env_vec = self._encode_environment()
        
        # 4. Resource State Vector
        if cellular_resources is None and simulation_engine is not None:
             cellular_resources = getattr(simulation_engine, 'cellular_resources', None)
        resource_vec = self._encode_resources(cellular_resources)
        resource_vec = np.nan_to_num(resource_vec.astype(np.float32))
        
        # 5. Concatenate into combined Cellular UR
        combined = np.concatenate([gene_embed, dynamics, env_vec, resource_vec])
        combined = np.nan_to_num(combined.astype(np.float32))
        
        # Count active genes
        active_count = 0
        if simulation_engine:
            for state in simulation_engine.gene_states.values():
                if not state.is_transcription_paused and state.mrna_count > 0:
                    active_count += 1
        else:
            active_count = len(enriched_genes)
        
        # Determine growth phase
        growth_phase = self._determine_growth_phase(simulation_engine, simulation_time)
        
        return CellStateVector(
            gene_embedding_aggregate=gene_embed,
            dynamics_vector=dynamics,
            environment_vector=env_vec,
            resource_vector=resource_vec,
            combined_cellular_ur=combined,
            timestamp=simulation_time,
            num_genes=len(enriched_genes),
            num_active_genes=active_count,
            growth_phase=growth_phase,
        )
    
    def _aggregate_gene_embeddings(self, genes: Dict[str, Gene]) -> np.ndarray:
        """
        Attention-weighted mean pooling of gene Molecular URs.
        
        Weights are derived from expression level:
            w_i = softmax(log(1 + mrna_expression_i))
        
        This ensures highly expressed genes contribute more to the
        cell-state representation.
        """
        embeddings = []
        weights = []
        
        for gene in genes.values():
            if gene.molecular_ur and 'combined_embedding' in gene.molecular_ur:
                emb = gene.molecular_ur['combined_embedding']
                # Use only the first 1280 dims (DNA embedding) for aggregation
                # to keep the cell-state compact
                if 'dna_embedding' in gene.molecular_ur:
                    emb = gene.molecular_ur['dna_embedding']
                
                if isinstance(emb, np.ndarray) and emb.shape[0] > 0:
                    embeddings.append(emb)
                    # Expression-weighted attention
                    w = np.log1p(max(gene.mrna_expression, 1.0))
                    weights.append(w)
        
        if not embeddings:
            return np.zeros(self.GENE_EMBED_DIM)
        
        embeddings = np.array(embeddings)
        weights = np.array(weights)
        
        # Softmax normalization
        weights = np.exp(weights - np.max(weights))
        weights = weights / (weights.sum() + 1e-8)
        
        # Weighted mean
        aggregated = np.average(embeddings, axis=0, weights=weights)
        
        # Ensure correct dimensionality
        if len(aggregated) < self.GENE_EMBED_DIM:
            aggregated = np.pad(aggregated, (0, self.GENE_EMBED_DIM - len(aggregated)))
        elif len(aggregated) > self.GENE_EMBED_DIM:
            aggregated = aggregated[:self.GENE_EMBED_DIM]
        
        return aggregated.astype(np.float32)
    
    def _extract_dynamics(self, genes: Dict[str, Gene], engine=None) -> np.ndarray:
        """
        Extract temporal dynamics vector from ODE engine state.
        
        Captures the instantaneous state of molecular concentrations:
        - Per-gene: mRNA count, protein count, status encoding
        - Aggregated: mean/std/min/max of mRNA and protein distributions
        - Trends: rising/declining fraction, steady-state fraction
        """
        dynamics = np.zeros(self.DYNAMICS_DIM)
        
        if engine is None:
            return dynamics
        
        states = list(engine.gene_states.values())
        if not states:
            return dynamics
        
        n = len(states)
        
        # === Distribution statistics (indices 0-63) ===
        mrna_counts = np.array([s.mrna_count for s in states])
        protein_counts = np.array([s.protein_count for s in states])
        mrna_pct = np.array([s.mrna_percentage for s in states])
        protein_pct = np.array([s.protein_percentage for s in states])
        
        # Debug trace
        if np.isnan(mrna_pct).any():
            nan_indices = np.where(np.isnan(mrna_pct))[0]
            print(f"DEBUG: NaN in mrna_pct for genes: {[states[i].gene_id for i in nan_indices[:5]]}")
        
        # mRNA statistics
        dynamics[0] = np.nanmean(mrna_counts)
        dynamics[1] = np.std(mrna_counts)
        dynamics[2] = np.median(mrna_counts)
        dynamics[3] = np.min(mrna_counts) if len(mrna_counts) > 0 else 0
        dynamics[4] = np.max(mrna_counts) if len(mrna_counts) > 0 else 0
        dynamics[5] = np.sum(mrna_counts)  # Total mRNA pool
        
        # Protein statistics
        dynamics[10] = np.mean(protein_counts)
        dynamics[11] = np.std(protein_counts)
        dynamics[12] = np.median(protein_counts)
        dynamics[13] = np.min(protein_counts) if len(protein_counts) > 0 else 0
        dynamics[14] = np.max(protein_counts) if len(protein_counts) > 0 else 0
        dynamics[15] = np.sum(protein_counts)  # Total protein pool
        
        # Percentage statistics
        dynamics[20] = np.nanmean(mrna_pct)
        dynamics[21] = np.nanmean(protein_pct)
        
        if np.isnan(dynamics[20:22]).any():
            print(f"DEBUG: NaN in Dynamics indices 20-21. mrna_pct has NaNs: {np.isnan(mrna_pct).any()}")

        # === Status fractions (indices 64-95) ===
        status_map = {
            'initializing': 0, 'rising': 1, 'steady_state': 2,
            'declining': 3, 'declining_fast': 4, 'depleted': 5,
            'paused': 6, 'affected_by_knockout': 7
        }
        for s in states:
            idx = status_map.get(s.status, 0)
            dynamics[64 + idx] += 1.0 / n  # Normalize to fractions
        
        # === Mass balance trends (indices 96-127) ===
        productions = np.nan_to_num(np.array([s.last_mrna_production for s in states]))
        degradations = np.nan_to_num(np.array([s.last_mrna_degradation for s in states]))
        net_balance = productions - degradations
        
        dynamics[96] = np.nanmean(net_balance)
        dynamics[97] = np.nanstd(net_balance)
        dynamics[98] = np.sum(net_balance > 0) / max(n, 1)  # Fraction growing
        dynamics[99] = np.sum(net_balance < 0) / max(n, 1)  # Fraction declining
        
        if np.isnan(dynamics[96:100]).any():
            print(f"DEBUG: NaN in Dynamics indices 96-99 (Mass balance)")

        # === Top-N gene concentrations (indices 128-255) ===
        # Encode top 64 most abundant proteins
        valid_prot = np.nan_to_num(protein_counts)
        top_protein_indices = np.argsort(valid_prot)[-64:]
        for i, idx in enumerate(top_protein_indices):
            dynamics[128 + i * 2] = protein_counts[idx]
            dynamics[129 + i * 2] = mrna_counts[idx]
            
        if np.isnan(dynamics[128:256]).any():
             print(f"DEBUG: NaN in Dynamics indices 128-255 (Top genes)")
        
        # === Kinetic parameter distribution (indices 256-319) ===
        half_lives_mrna = np.array([s.mrna_half_life for s in states])
        half_lives_prot = np.array([s.protein_half_life for s in states])
        tx_rates = np.array([s.transcription_rate for s in states])
        tl_rates = np.array([s.translation_rate for s in states])
        
        dynamics[256] = np.mean(half_lives_mrna)
        dynamics[257] = np.std(half_lives_mrna)
        dynamics[260] = np.mean(half_lives_prot)
        dynamics[261] = np.std(half_lives_prot)
        dynamics[264] = np.mean(tx_rates)
        dynamics[265] = np.std(tx_rates)
        dynamics[268] = np.mean(tl_rates)
        dynamics[269] = np.std(tl_rates)
        
        # === Simulation metadata (indices 320-511) ===
        dynamics[320] = engine.simulation_time
        dynamics[321] = float(n)
        dynamics[322] = engine.noise_level
        dynamics[323] = 1.0 if engine.enable_noise else 0.0
        dynamics[324] = 1.0 if engine.enable_resource_limits else 0.0
        
        return dynamics.astype(np.float32)
    
    def _encode_environment(self) -> np.ndarray:
        """
        Encode environmental parameters into a fixed-size vector.
        Uses normalized values for numerical stability.
        """
        env = np.zeros(self.ENV_DIM)
        
        # Core parameters (normalized)
        env[0] = self.env_params.get('temperature', 37.0) / 50.0
        env[1] = self.env_params.get('ph', 7.2) / 14.0
        env[2] = self.env_params.get('osmolarity', 300.0) / 1000.0
        env[3] = self.env_params.get('doubling_time', 180.0) / 360.0
        env[4] = self.env_params.get('growth_rate', 0.0039) * 100.0
        env[5] = self.env_params.get('atp_concentration', 5.0) / 10.0
        env[6] = self.env_params.get('glucose_available', 1.0)
        env[7] = self.env_params.get('oxygen_available', 0.0)
        
        # One-hot encoded organism type (Mycoplasma = minimal)
        env[10] = 1.0  # is_minimal_cell
        env[11] = 0.0  # is_ecoli_like
        env[12] = 0.0  # is_eukaryotic
        
        # Growth condition encoding
        env[20] = 1.0  # in_vitro
        env[21] = 0.0  # in_vivo
        
        return env.astype(np.float32)
    
    def _encode_resources(self, resources=None) -> np.ndarray:
        """
        Encode cellular resource state into a fixed-size vector.
        """
        res = np.zeros(self.RESOURCE_DIM)
        
        if resources is None:
            return res.astype(np.float32)
        
        # Ribosome state
        res[0] = resources.total_ribosomes / 5000.0  # Normalized
        res[1] = resources.ribosomes_in_use / max(resources.total_ribosomes, 1)
        res[2] = resources.ribosome_saturation
        
        # Polymerase state
        res[5] = resources.total_rna_polymerases / 500.0
        res[6] = resources.polymerases_in_use / max(resources.total_rna_polymerases, 1)
        res[7] = resources.polymerase_saturation
        
        # Energy state
        res[10] = min(1.0, resources.total_atp / 1e6)
        res[11] = min(1.0, resources.total_amino_acids / 1e7)
        res[12] = resources.atp_regeneration_rate / 1e4
        
        # Genome derivation metadata
        res[15] = resources.ribosomal_gene_count / 60.0
        res[16] = resources.polymerase_gene_count / 10.0
        res[17] = resources.atp_synthase_gene_count / 10.0
        res[18] = 1.0 if resources.resource_source == "genome" else 0.0
        
        return res.astype(np.float32)
    
    def _determine_growth_phase(self, engine, simulation_time: float) -> str:
        """Determine current growth phase from simulation state."""
        if engine is None:
            return "unknown"
        
        stats = engine.get_statistics()
        if not stats:
            return "unknown"
        
        avg_protein_pct = stats.get('avg_protein_percentage', 0)
        steady_fraction = stats.get('steady_state_genes', 0) / max(stats.get('total_genes', 1), 1)
        
        if simulation_time < 60:
            return "lag"
        elif avg_protein_pct < 30:
            return "lag"
        elif steady_fraction > 0.7:
            return "stationary"
        elif avg_protein_pct > 50:
            return "exponential"
        else:
            return "lag"
    
    def update_cellular_ur(
        self,
        previous_state: CellStateVector,
        enriched_genes: Dict[str, Gene],
        simulation_engine=None,
        cellular_resources=None,
        simulation_time: float = 0.0,
    ) -> CellStateVector:
        """
        Efficiently update the Cellular UR (only recompute dynamics + resources,
        keep gene embeddings cached from previous state).
        """
        # Gene embeddings are static — reuse from previous state
        gene_embed = previous_state.gene_embedding_aggregate
        
        # Recompute dynamic components
        dynamics = self._extract_dynamics(enriched_genes, simulation_engine)
        env_vec = self._encode_environment()
        resource_vec = self._encode_resources(cellular_resources)
        
        combined = np.concatenate([gene_embed, dynamics, env_vec, resource_vec])
        
        active_count = 0
        if simulation_engine:
            for state in simulation_engine.gene_states.values():
                if not state.is_transcription_paused and state.mrna_count > 0:
                    active_count += 1
        else:
            active_count = previous_state.num_active_genes
        
        growth_phase = self._determine_growth_phase(simulation_engine, simulation_time)
        
        return CellStateVector(
            gene_embedding_aggregate=gene_embed,
            dynamics_vector=dynamics,
            environment_vector=env_vec,
            resource_vector=resource_vec,
            combined_cellular_ur=combined,
            timestamp=simulation_time,
            num_genes=previous_state.num_genes,
            num_active_genes=active_count,
            growth_phase=growth_phase,
        )
