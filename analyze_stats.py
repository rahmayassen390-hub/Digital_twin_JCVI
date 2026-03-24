
import sys
import os
import collections
from Bio.Seq import Seq

# Add project root to path
sys.path.append("/home/rahema/Desktop/bed_files /jcvi_genome_analyzer_updated (6)/")
sys.path.append("/home/rahema/Desktop/bed_files /jcvi_genome_analyzer_updated (6)/jcvi_genome_analyzer_updated/")

from jcvi_genome_analyzer_updated.genome_manager import GenomeManager
from jcvi_genome_analyzer_updated.gui_methods import validate_protein_sequence

def analyze_jcvi_stats():
    # Paths
    fasta_path = "/home/rahema/Desktop/JCVI_3.0.fasta"
    gff_path = "/home/rahema/Desktop/JCVI_3.0.gff3"
    regulation_bed = "/home/rahema/Desktop/bed_files /JCVI_Gene_Regulation_3.0#.bed"
    transcription_bed = "/home/rahema/Desktop/bed_files /JCVI_Transcription_Dynamics_3.0#.bed"
    translation_bed = "/home/rahema/Desktop/bed_files /JCVI_Translation_Dynamics_3.0#.bed"

    # Initialize Manager
    gm = GenomeManager()
    
    print("--- Loading Genome ---")
    # Load FASTA
    seq_name, seq, seq_len = gm.load_fasta(fasta_path)
    gm.sequence = seq
    gm.sequence_name = seq_name
    gm.sequence_length = seq_len
    
    # Load GFF
    genes = gm.load_gff(gff_path)
    gm.genes = genes
    gm.annotate_genes_with_promoters()
    gm.extract_gene_sequences()
    
    print(f"Total Genes Loaded: {len(genes)}")
    
    # --- Analysis 1: Genome Validation ---
    valid_start_count = 0
    missing_stop_count = 0
    invalid_cds_count = 0
    
    valid_starts = {'ATG', 'GTG', 'TTG'}
    valid_stops = {'TAA', 'TAG', 'TGA'}
    
    for gene_id, gene in genes.items():
        if not gene.dna_sequence:
            invalid_cds_count += 1
            continue
            
        seq_str = str(gene.dna_sequence).upper()
        
        # Check start
        if seq_str[:3] in valid_starts:
            valid_start_count += 1
            
        # Check stop
        if seq_str[-3:] not in valid_stops:
            missing_stop_count += 1
            
        # Basic CDS validation (length mainly)
        if len(seq_str) < 90: # < 30 aa
             invalid_cds_count += 1
    
    print(f"\n--- Genome Statistics ---")
    print(f"Genes Identified: {len(genes)} (Expected: ~473)")
    print(f"Genes with Invalid/Short CDS: {invalid_cds_count}")
    print(f"Valid Start Codons: {valid_start_count} ({valid_start_count/len(genes)*100:.1f}%)")
    print(f"Missing Stop Codons: {missing_stop_count}")

    # --- Analysis 2: Regulation Distribution ---
    print(f"\n--- Gene Regulation Distribution ---")
    print(f"Source File: {regulation_bed}")
    
    reg_counts = collections.defaultdict(int)
    
    with open(regulation_bed, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                # Format: gene_id|type|...
                name_col = parts[3]
                sub_parts = name_col.split('|')
                if len(sub_parts) > 1:
                    reg_type = sub_parts[1]
                    reg_counts[reg_type] += 1
    
    for r_type, count in reg_counts.items():
        print(f"{r_type}: {count}")

    # --- Analysis 3: Expression Levels ---
    print(f"\n--- Expression Levels ---")
    
    # mRNA from Transcription BED
    mrna_counts = []
    
    with open(transcription_bed, 'r') as f:
        for line in f:
            if line.startswith('Transcription_Level'): # Track 1
                continue 
            if line.startswith('#') or not line.strip(): continue
                
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                name = parts[3]
                # gene_id|Level:XX%
                if "Level:" in name:
                    try:
                        level_str = name.split("Level:")[1].replace('%', '')
                        mrna_counts.append(float(level_str))
                    except: pass

    if mrna_counts:
        # Convert % to approximate count (assuming max ~50 copies for single cell or based on scaling)
        # The prompt asks for "observed", so I'll report the % value first, 
        # but also check how the code interprets it.
        # In the code: rate = 0.1 + (level / 100.0) * 0.2. 
        # Actually let's just report the raw numbers found in the file first.
        print(f"mRNA Level (from BED %): Min {min(mrna_counts)}%, Max {max(mrna_counts)}%")
    
    # Protein from Translation BED
    protein_counts = []
    ribosomal_examples = []
    
    with open(translation_bed, 'r') as f:
        for line in f:
             if line.startswith('Translation_Level'): continue
             if line.startswith('#') or not line.strip(): continue
             
             parts = line.strip().split('\t')
             if len(parts) >= 4:
                 name = parts[3]
                 # gene_id|Protein:XX%|RBS:strength
                 if "Protein:" in name:
                     try:
                         p_level = float(name.split("Protein:")[1].split('|')[0].replace('%', ''))
                         protein_counts.append(p_level)
                         
                         gene_id = name.split('|')[0]
                         if gene_id.startswith('rps') or gene_id.startswith('rpl'):
                             ribosomal_examples.append((gene_id, p_level))
                     except: pass
                     
    if protein_counts:
        print(f"Protein Level (from BED %): Min {min(protein_counts)}%, Max {max(protein_counts)}%")
        
    print("\nRibosomal Examples (Protein %):")
    for g, p in ribosomal_examples[:5]:
        print(f"{g}: {p}%")

if __name__ == "__main__":
    analyze_jcvi_stats()
