
import sys
import os
import math

# Add project root AND inner package to path to allow direct imports
sys.path.append("/home/rahema/Desktop/bed_files /jcvi_genome_analyzer_updated (6)")
sys.path.append("/home/rahema/Desktop/bed_files /jcvi_genome_analyzer_updated (6)/jcvi_genome_analyzer_updated")

from genome_manager import GenomeManager
from cell_simulation_engine import CellSimulationEngine, SimulationMode

def extract_gene_data():
    print("Loading Genome Data...")
    gm = GenomeManager()
    
    # Load all files
    base_dir = "/home/rahema/Desktop/bed_files /"
    gff_path = "/home/rahema/Desktop/JCVI_3.0.gff3"
    
    # 1. Load Genes
    if os.path.exists(gff_path):
        gm.genes = gm.load_gff(gff_path)
    
    # 2. Load BEDs (Regulations, Dynamics)
    bed_files = [f for f in os.listdir(base_dir) if f.endswith(".bed")]
    for bed in bed_files:
        full_path = os.path.join(base_dir, bed)
        try:
            if "Regulation" in bed:
                gm.load_regulation_bed(full_path)
            elif "Transcription" in bed:
                gm.load_transcription_dynamics_bed(full_path)
            elif "Translation" in bed:
                gm.load_translation_dynamics_bed(full_path)
        except: pass

    print(f"Loaded {len(gm.genes)} genes.")
    
    # Define genes of interest (Manual map or search)
    # Based on grep results or common knowledge for JCVI syn3.0
    # dnaA is usually 0001
    # dnaN is usually 0002
    # rpoB found as 0804
    # gap found as 0607
    
    target_genes = [
        "JCVISYN3_0001", # dnaA
        "JCVISYN3_0002", # dnaN
        "JCVISYN3_0804", # rpoB
        "JCVISYN3_0607", # gap (GAPDH)
        "JCVISYN3_0004", # gyrB
        "JCVISYN3_0216", # rpoA? check
        "JCVISYN3_0418", # ftsZ? check
    ]
    
    print("\n| Gene ID | Name | Category | Regulation | Phase | RBS Strength | Final mRNA | Final Protein |")
    print("|---|---|---|---|---|---|---|---|")
    
    for gid in target_genes:
        # Find gene object
        gene_obj = next((g for g in gm.genes if g.id == gid), None)
        
        # Get BED Data
        reg_data = gm.regulation_data.get(gid)
        tx_data = gm.transcription_dynamics.get(gid)
        tl_data = gm.translation_dynamics.get(gid)
        
        name = gene_obj.name if gene_obj else gid
        # Clean name (remove JCVISYN3_ prefix if redundant)
        if name == gid: 
            # Try to get product from BED data or gene object
            if gene_obj and gene_obj.product:
                name = gene_obj.product.split(':')[0]
        
        category = "Unknown"
        if reg_data and reg_data.category: category = reg_data.category
        elif tx_data and tx_data.category: category = tx_data.category
        
        regulation = reg_data.regulation_type if reg_data else "Unknown"
        
        phase = "Unphased"
        if tx_data and tx_data.time_range:
             # Extract phase from time range or name if stored
             # The loader stores 'phase' in time_range sometimes or we have to parse
             # BED loader puts it in `time_range`? no, `phase` usually.
             pass
        
        # In current GenomeManager, `phase` isn't explicitly stored on TranscriptionDynamics object 
        # but `temporal_order` is. The BED reader logic had phase parsing.
        
        rbs = tl_data.rbs_strength if tl_data else "Unknown"
        
        # Calculate Steady State / Final Values
        # Using CellSimulationEngine logic:
        # mRNA_ss = tx_rate / k_deg_mrna
        # Protein_ss = tl_rate * mRNA_ss / k_deg_prot
        
        # Rates
        tx_rate = 0.0 # derived from level?
        if tx_data:
             # In gui_methods, rate = 0.2 / duration ... complex
             # But simplified SS: level * scalar
             # Let's use the stored "expression_level" or "transcription_level" as proxy for relative count
             # In simulation engine, default tx_rate=1.0.
             # The BED data has 'transcription_level' (0-100)
             # Let's map Level -> Rate: 100% = 2.0 copies/sec?
             rate_val = (tx_data.transcription_level / 100.0) * 2.0
             mrna_hl = 300.0
             k_m = math.log(2)/mrna_hl
             final_mrna = rate_val / k_m
        else:
             final_mrna = 0.0
             
        if tl_data:
             # TL Rate
             rbs_map = {'strong': 0.3, 'medium': 0.2, 'weak': 0.1}
             tl_rate = rbs_map.get(str(rbs).lower(), 0.15)
             prot_hl = 3600.0
             k_p = math.log(2)/prot_hl
             final_prot = (tl_rate * final_mrna) / k_p
        else:
             final_prot = 0.0
             
        # Format Phase
        # We can try to infer phase from order
        order = getattr(tx_data, 'order', getattr(tx_data, 'temporal_order', 0))
        if order > 0:
             if order < 150: phase = "1 (Early)"
             elif order < 300: phase = "2 (Mid)"
             else: phase = "3 (Late)"
             
        print(f"| {gid} | {name} | {category} | {regulation} | {phase} | {rbs} | {final_mrna:.1f} | {final_prot:.0f} |")

def analyze_errors():
    print("\n=== Error Analysis ===")
    gm = GenomeManager()
    gff_path = "/home/rahema/Desktop/JCVI_3.0.gff3"
    
    if os.path.exists(gff_path):
        gm.genes = gm.load_gff(gff_path)
        
    print(f"Genes successfully parsed: {len(gm.genes)}")
    
    # Check for warnings in BED files
    # We can't easily capture print output of the loader here without redirecting stdout
    # But we can check internal consistency
    
    base_dir = "/home/rahema/Desktop/bed_files /"
    bed_files = [f for f in os.listdir(base_dir) if f.endswith(".bed")]
    
    for bed in bed_files:
        full_path = os.path.join(base_dir, bed)
        try:
             # Load and count "bed fallback" cases manually
             # The loader prints them, so user will see them in output
             if "Translation" in bed:
                 print(f"Re-loading {bed} to check fallbacks...")
                 gm.load_translation_dynamics_bed(full_path)
        except: pass

if __name__ == "__main__":
    extract_gene_data()
    analyze_errors()
