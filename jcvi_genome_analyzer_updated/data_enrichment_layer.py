import os
import pandas as pd
import numpy as np
import torch
from typing import Dict, List, Optional
from data_structures import Gene

class DataEnrichmentLayer:
    """
    Data Enrichment Layer for AI Virtual Cell.
    Loads supplementary data and enriches Gene objects.
    Fills gaps using AI models.
    """
    
    def __init__(self, data_dir: str = "Data"):
        self.data_dir = data_dir
        self.datasets = {}
        self._load_all_datasets()
        
    def _load_all_datasets(self):
        """Load all biological datasets into memory"""
        files_to_load = {
            'kinetic_params': 'Kinetic_Parameters.tsv',
            'mrna_counts': 'mRNA_counts.csv',
            'proteomics': 'proteomics_data.xlsx',
            'full_genome': 'JCVI_syn3A_full_genome_dataset.csv',
            'degradation_rates': 'syn3A_degradation_rates.csv',
            'orthologs_rbs': 'syn3A_gene_orthologs_rbs.csv',
            'global_params': 'global_parameters.csv',
            'transport_fluxes': 'transport_fluxes.tsv',
            'morphology': 'syn3A_morphology_params.csv',
            'lipid_data': 'syn3A_lipid_data.csv',
            'summary_kinetics': 'syn3A_kinetic_parameters.csv'
        }
        
        for key, filename in files_to_load.items():
            path = os.path.join(self.data_dir, filename)
            if not os.path.exists(path):
                print(f"Warning: Data file not found: {path}")
                continue
                
            try:
                if filename.endswith('.csv'):
                    self.datasets[key] = pd.read_csv(path)
                elif filename.endswith('.tsv'):
                    self.datasets[key] = pd.read_csv(path, sep='\t')
                elif filename.endswith('.xlsx'):
                    self.datasets[key] = pd.read_excel(path)
                print(f"Loaded {filename}")
            except Exception as e:
                print(f"Error loading {filename}: {e}")

    def enrich_gene(self, gene: Gene) -> Gene:
        """
        Enriched a single gene object with multi-source data.
        Maps data using locus_tag (gene.id).
        
        PRIORITY SYSTEM:
        1. Specialized Files (mRNA counts, degradation rates, RBS orthologs)
        2. General File (JCVI_syn3A_full_genome_dataset.csv) as fallback
        """
        locus_tag = gene.id
        norm_tag = locus_tag.replace('JCVISYN3A', 'JCVISYN3')
        
        # --- PHASE 1: LOAD BASELINE FROM GENERAL FILE ---
        if 'full_genome' in self.datasets:
            df = self.datasets['full_genome']
            gene_data = df[df['locus_tag'].str.replace('JCVISYN3A', 'JCVISYN3') == norm_tag]
            
            if not gene_data.empty:
                row = gene_data.iloc[0]
                gene.essentiality_status = row.get('essentiality_status', gene.essentiality_status)
                gene.function_class = row.get('function_class', gene.function_class)
                gene.kegg_id = row.get('kegg_id', gene.kegg_id)
                gene.protein_abundance = float(row.get('protein_abundance', gene.protein_abundance))
                gene.rbs_strength = float(row.get('rbs_strength_index', gene.rbs_strength))
                gene.mrna_half_life = float(row.get('mrna_half_life_min', gene.mrna_half_life))
                gene.protein_half_life = float(row.get('protein_half_life_hr', gene.protein_half_life))
                gene.ortholog_id = row.get('ortholog_syn1', gene.ortholog_id)

        # --- PHASE 2: OVERWRITE WITH SPECIALIZED DATA (HIGHER ACCURACY) ---
        
        # 2.1 mRNA Expression (Specialized mRNA counts)
        if 'mrna_counts' in self.datasets:
            df = self.datasets['mrna_counts']
            if 'LocusTag' in df.columns:
                gene_data = df[df['LocusTag'].str.replace('JCVISYN3A', 'JCVISYN3') == norm_tag]
                if not gene_data.empty:
                    # OVERWRITE: Specialized count is more accurate
                    gene.mrna_expression = float(gene_data.iloc[0].get('Count', gene.mrna_expression))

        # 2.2 Degradation Rates (Specialized stability data)
        if 'degradation_rates' in self.datasets:
            df = self.datasets['degradation_rates']
            if 'locus_tag' in df.columns:
                gene_data = df[df['locus_tag'].str.replace('JCVISYN3A', 'JCVISYN3') == norm_tag]
                if not gene_data.empty:
                    # OVERWRITE: Measured half-life from specialized table
                    gene.mrna_half_life = float(gene_data.iloc[0].get('mrna_half_life_min', gene.mrna_half_life))

        # 2.3 RBS Strength / Orthologs (Specialized supplementary)
        if 'orthologs_rbs' in self.datasets:
            df = self.datasets['orthologs_rbs']
            col_name = 'locus_tag_syn3A' if 'locus_tag_syn3A' in df.columns else 'locus_tag'
            if col_name in df.columns:
                gene_data = df[df[col_name].str.replace('JCVISYN3A', 'JCVISYN3') == norm_tag]
                if not gene_data.empty:
                    # OVERWRITE: Specialized RBS strength index
                    gene.rbs_strength = float(gene_data.iloc[0].get('rbs_strength_relative', gene.rbs_strength))
                    gene.ortholog_id = gene_data.iloc[0].get('locus_tag_syn1.0', 
                                       gene_data.iloc[0].get('ortholog_syn1', gene.ortholog_id))

        # --- PHASE 3: AI GAP FILLING (LAST RESORT) ---
        if gene.mrna_half_life <= 0.0 or gene.mrna_half_life == 5.0: # 5.0 is current default
             gene.mrna_half_life = self.predict_mrna_half_life(gene)
        if gene.protein_half_life <= 0.0 or gene.protein_half_life == 25.0: # 25.0 is current default
             gene.protein_half_life = self.predict_protein_half_life(gene)

        return gene

    def predict_rbs_strength(self, gene_sequence: str) -> float:
        """
        AI GAP FILLING: Predict RBS strength using Shine-Dalgarno analysis.
        Standard anti-SD: CCUCCU (complement: AGGAGG)
        """
        if not gene_sequence:
            return 1.0
            
        # Simplified: check for variants in the start of the sequence
        # (Assuming the sequence provided starts near the RBS)
        upstream = gene_sequence[:30]
        
        # Scoring based on consensus severity
        sd_scores = {
            'AGGAGG': 10.0,
            'AAGGAG': 8.0,
            'AGGAG': 6.0,
            'GGAGG': 6.0,
            'GAGG': 4.0,
            'AGGA': 4.0
        }
        
        max_score = 1.0
        for variant, score in sd_scores.items():
            if variant in upstream:
                max_score = max(max_score, score)
        
        return max_score

    def predict_mrna_half_life(self, gene: Gene) -> float:
        """
        AI GAP FILLING: Predict mRNA half-life using Linear Regression.
        Linear Regression is used instead of Random Forest to prevent overfitting 
        on the extremely small available dataset (syn3A_degradation_rates.csv).
        """
        if 'degradation_rates' not in self.datasets:
            return 2.4 # Global mean for JCVI syn3A
            
        try:
            from sklearn.linear_model import LinearRegression
            
            df = self.datasets['degradation_rates'].copy()
            df = df[df['locus_tag'] != 'default_gene']
            
            if len(df) < 3:
                return 2.4
            
            # Simple features: GC content and length proxy
            # (In a real scenario we'd use sequence motifs)
            df['length'] = df['locus_tag'].apply(lambda x: 1000) # Placeholder for now
            df['gc'] = 0.3 # Typical low GC for JCVI
            
            X = df[['length', 'gc']]
            y = df['mrna_half_life_min']
            
            model = LinearRegression()
            model.fit(X, y)
            
            # Predict for current gene
            gene_len = len(gene.sequence) if gene.sequence else 1000
            gene_gc = (gene.sequence.count('G') + gene.sequence.count('C')) / gene_len if gene_len > 0 else 0.3
            
            pred = model.predict([[gene_len, gene_gc]])[0]
            # Clip to biological range [0.5, 15.0]
            return round(max(0.5, min(15.0, float(pred))), 2)
            
        except ImportError:
            return 2.4
        except Exception as e:
            print(f"Error in mRNA prediction: {e}")
            return 2.4

    def predict_protein_half_life(self, gene: Gene) -> float:
        """
        AI GAP FILLING: Predict protein half-life using N-End Rule + Features.
        """
        # Default starting point
        half_life = 25.0 
        
        # N-End Rule Lookup (Approximate)
        # Stabilizing: M, G, A, S, T, V
        # Destabilizing: R, K, D, E, F, L, W, Y
        
        if gene.sequence and len(gene.sequence) >= 3:
            # Simple start codon check (usually Met)
            start_codon = gene.sequence[:3].upper()
            if start_codon == 'ATG': # Methionine
                half_life = 40.0
            elif start_codon in ['GTG', 'TTG']: # Valine/Leucine
                half_life = 20.0
        
        # Adjust based on essentiality
        if gene.essentiality_status == "Essential":
            half_life *= 1.5
        elif gene.essentiality_status == "Nonessential":
            half_life *= 0.8
            
        # PEST sequences check (simplified)
        if gene.sequence and 'CCT' in gene.sequence and 'GAA' in gene.sequence:
            half_life *= 0.7 # Likely destabilizing regions
            
        return round(half_life, 2)
