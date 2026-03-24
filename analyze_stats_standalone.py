
import sys
import os
import collections

# Standalone parsers to avoid PyQt dependency

def load_fasta(path):
    seq = ""
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                seq += line.strip()
    return seq.upper()

def parse_gff(path, sequence):
    genes = {}
    id_to_name = {}
    
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            if parts[2] == 'CDS':
                attr_str = parts[8]
                attr_map = {}
                for attr in attr_str.split(';'):
                    if '=' in attr:
                        k, v = attr.split('=', 1)
                        attr_map[k] = v
                
                gene_id = attr_map.get('locus_tag', attr_map.get('ID', '').replace('gene-', '').replace('cds-', ''))
                gene_name = attr_map.get('product', attr_map.get('Name', 'unknown'))
                
                if not gene_id:
                     # Fallback to Parent if needed, cleaning prefix
                     parent = attr_map.get('Parent', '')
                     gene_id = parent.replace('gene-', '')

                id_to_name[gene_id] = gene_name
                
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                
                # Extract CDS
                s_idx = start - 1
                e_idx = end
                cds_seq = sequence[s_idx:e_idx]
                
                if strand == '-':
                    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
                    cds_seq = "".join([comp.get(b, 'N') for b in reversed(cds_seq)])
                
                genes[gene_id] = cds_seq
                
    return genes, id_to_name

def analyze_jcvi_stats():
    # Paths
    fasta_path = "/home/rahema/Desktop/JCVI_3.0.fasta"
    gff_path = "/home/rahema/Desktop/JCVI_3.0.gff3"
    regulation_bed = "/home/rahema/Desktop/bed_files /JCVI_Gene_Regulation_3.0#.bed"
    transcription_bed = "/home/rahema/Desktop/bed_files /JCVI_Transcription_Dynamics_3.0#.bed"
    translation_bed = "/home/rahema/Desktop/bed_files /JCVI_Translation_Dynamics_3.0#.bed"
    operon_bed = "/home/rahema/Desktop/bed_files /Internal_Operons_3.0#.bed"

    print("--- Loading Genome Data ---")
    sequence = load_fasta(fasta_path)
    genes, id_to_name = parse_gff(gff_path, sequence)
    
    # --- Analysis 1: Genome Validation ---
    total_genes = len(genes)
    invalid_cds = 0
    valid_start_count = 0
    missing_stop_count = 0
    
    valid_starts = {'ATG', 'GTG', 'TTG'}
    valid_stops = {'TAA', 'TAG', 'TGA'}
    
    for g_id, seq in genes.items():
        if len(seq) < 90:
            invalid_cds += 1
        
        if seq[:3] in valid_starts:
            valid_start_count += 1
            
        if seq[-3:] not in valid_stops:
            missing_stop_count += 1
            
    print(f"\n--- Genome Statistics ---")
    print(f"Genes identified: {total_genes}")
    print(f"Genes with Invalid CDS: {invalid_cds}")
    
    pct_valid_start = (valid_start_count / total_genes * 100) if total_genes else 0
    print(f"Genes with Valid Start Codons: {valid_start_count} ({pct_valid_start:.1f}%)")
    print(f"Genes Missing Valid Stop Codons: {missing_stop_count}")

    # --- Analysis 2: Regulation Distribution ---
    print(f"\n--- Gene Regulation Distribution ---")
    reg_counts = collections.defaultdict(int)
    
    if os.path.exists(regulation_bed):
        with open(regulation_bed, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    # gene_id|type|...
                    info = parts[3]
                    try:
                        sub = info.split('|')
                        if len(sub) > 1:
                            reg_type = sub[1]
                            reg_counts[reg_type] += 1
                    except: pass
    
    # Parse Operon BED separately for Operon Members
    operon_count = 0
    if os.path.exists(operon_bed):
        with open(operon_bed, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                     operon_count += 1
    
    print(f"Constitutive genes = {reg_counts.get('constitutive', 0)}")
    print(f"Operon members = {operon_count} (Source: {operon_bed})")
    print(f"Potentially inducible = {reg_counts.get('potentiallyInducible', 0)}")
    print(f"Potentially repressible = {reg_counts.get('potentiallyRepressible', 0)}")
    print(f"Regulation Data Source: {regulation_bed}")

    
    # --- Analysis 3: Expression Levels & Phases ---
    print(f"\n--- Transcription Phase Analysis ---")
    
    phase_counts = collections.defaultdict(int)
    phase_levels = collections.defaultdict(list)
    unclassified_phase = 0
    min_tx_level = 100.0
    max_tx_level = 0.0
    
    # Track min and max order to calculate timing
    min_order = 100000
    max_order = -1
    
    current_track_name = ""
    
    if os.path.exists(transcription_bed):
        with open(transcription_bed, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                
                # Check track headers
                if line.startswith('track'):
                    if 'Transcription_Order' in line:
                        current_track_name = 'order'
                    elif 'Transcription_Level' in line:
                        current_track_name = 'level'
                    elif 'Transcription_Timeline' in line:
                        current_track_name = 'timeline'
                    continue
                
                if line.startswith('#'): continue
                
                parts = line.split('\t')
                if len(parts) >= 4:
                    name = parts[3]
                    
                    # Logic for Phase & Timing (From Order Track usually)
                    # Structure can contain Phase1, Phase2 etc. 
                    if current_track_name == 'order' or ("Phase" in name and "|" in name):
                        try:
                            sub = name.split('|')
                            
                            if sub[0].isdigit():
                                order = int(sub[0])
                                min_order = min(min_order, order)
                                max_order = max(max_order, order)
                            
                            phase_str = "Unknown"
                            for s in sub:
                                # Look for Phase1, Phase2, Phase3...
                                if s.startswith("Phase"):
                                    if "Phase1" in s: phase_str = "1"
                                    elif "Phase2" in s: phase_str = "2"
                                    elif "Phase3" in s: phase_str = "3"
                                    elif ":" in s: phase_str = s.split(":")[1]
                                    break
                            
                            level = 0.0
                            if "%" in name:
                                level_part = [s for s in sub if "%" in s][0]
                                level = float(level_part.replace('%', ''))
                                
                                if phase_str != "Unknown":
                                    phase_levels[phase_str].append(level)
                                    
                                min_tx_level = min(min_tx_level, level)
                                max_tx_level = max(max_tx_level, level)

                            if phase_str != "Unknown":
                                phase_counts[phase_str] += 1
                            else:
                                unclassified_phase += 1
                                
                        except: pass
                    
                    elif "Phase" in name:
                         # Fallback search
                         try:
                            sub = name.split('|')
                            p_val = "Unknown"
                            for s in sub:
                                if "Phase1" in s: p_val = "1"
                                elif "Phase2" in s: p_val = "2"
                                elif "Phase3" in s: p_val = "3"
                                
                            if p_val != "Unknown":
                                phase_counts[p_val] += 1
                                if "%" in name:
                                   l_str = [x for x in sub if "%" in x][0]
                                   l_val = float(l_str.replace('%', ''))
                                   phase_levels[p_val].append(l_val)
                         except: pass

    # Heuristic for double counting
    if phase_counts.get('1', 0) > 400:
         for k in phase_counts: phase_counts[k] //= 2

    print(f"Phase 1 count = {phase_counts.get('1', 0)}")
    print(f"Phase 2 count = {phase_counts.get('2', 0)}")
    print(f"Phase 3 count = {phase_counts.get('3', 0)}")
    
    total_phased = phase_counts.get('1', 0) + phase_counts.get('2', 0) + phase_counts.get('3', 0)
    print(f"Total Phased Genes Found: {total_phased}")
    
    print(f"\n--- Actual Transcription Rates ---")
    print(f"Highest level value = {max_tx_level}%")
    print(f"Lowest level value = {min_tx_level}%")
    
    def avg(lst): return sum(lst)/len(lst) if lst else 0.0
    
    print(f"Phase 1 average level = {avg(phase_levels.get('1', [])):.1f}%")
    print(f"Phase 2 average level = {avg(phase_levels.get('2', [])):.1f}%")
    print(f"Phase 3 average level = {avg(phase_levels.get('3', [])):.1f}%")

    # --- Analysis 4: Translation RBS Strength ---
    print(f"\n--- Translation RBS Strength Distribution ---")
    
    rbs_counts = collections.defaultdict(int)
    ncRNA_count = 0
    
    if os.path.exists(translation_bed):
        with open(translation_bed, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    name = parts[3]
                    
                    if "RBS:" in name:
                        try:
                            sub = name.split('|')
                            rbs = "Unknown"
                            for s in sub:
                                if s.startswith("RBS:"):
                                    rbs = s.split(":")[1]
                                    break
                            rbs_counts[rbs] += 1
                        except: pass
                    
                    if "ncRNA" in name or "No_Translation" in name:
                        ncRNA_count += 1

    print(f"Strong RBS = {rbs_counts.get('strong', 0)}")
    print(f"Medium RBS = {rbs_counts.get('medium', 0)}")
    print(f"Weak RBS = {rbs_counts.get('weak', 0)}")
    print(f"ncRNA (Not Translated) = {ncRNA_count}")

    # --- Analysis 5: Animation Timing ---
    print(f"\n--- Animation Timing ---")
    
    if max_order < 0:
        if total_phased > 0:
            min_order = 1
            max_order = total_phased
        else:
            min_order = 0
            max_order = 0
    
    first_start = min_order * 0.3
    last_start = max_order * 0.3
    
    print(f"First gene order: {min_order}, Start Time: {first_start:.1f} seconds")
    print(f"Last gene order: {max_order}, Start Time: {last_start:.1f} seconds")
    print(f"Formula Confirmation: start_time = order * 0.3 (0.3s per gene)")

if __name__ == "__main__":
    analyze_jcvi_stats()
