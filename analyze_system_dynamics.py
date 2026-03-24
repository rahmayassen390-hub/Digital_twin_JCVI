
import os
import math
import collections

def analyze_operons():
    operon_bed = "/home/rahema/Desktop/bed_files /Internal_Operons_3.0#.bed"
    translation_bed = "/home/rahema/Desktop/bed_files /JCVI_Translation_Dynamics_3.0#.bed"
    
    # 1. Parse Operons
    # Map Operon_Name -> List of Gene_IDs
    operons = collections.defaultdict(list)
    
    if os.path.exists(operon_bed):
        with open(operon_bed, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    name = parts[3]
                    # Name format: JCVISYN3_xxxx_OperonY or just OperonY?
                    # Previous head showed: JCVISYN3_0003_Operon1
                    if "_Operon" in name:
                        # Extract the Operon Tag
                        # JCVISYN3_0003_Operon1 -> Operon1
                        # JCVISYN3_0106_Operon3 -> Operon3
                        op_tag = name.split('_')[-1]
                        
                        # Extract Gene ID if present (e.g. JCVISYN3_0003)
                        # Assumes format JCVISYN3_XXXX_OperonY
                        if name.startswith("JCVISYN3"):
                            gene_id = "_".join(name.split('_')[:2])
                            operons[op_tag].append(gene_id)
    
    total_operons = len(operons)
    print(f"Total Operons Found: {total_operons}")
    
    # 2. Parse Translation Orders
    # Map Gene_ID -> Order
    gene_orders = {}
    
    if os.path.exists(translation_bed):
        with open(translation_bed, 'r') as f:
            for line in f:
                if line.startswith('track') or line.startswith('#') or not line.strip(): continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    name = parts[3]
                    # Format: 7|JCVISYN3_0452|Phase1...
                    try:
                        sub = name.split('|')
                        if sub[0].isdigit():
                            order = int(sub[0])
                            # Find gene id? usually 2nd part
                            # 7|JCVISYN3_0452|...
                            gid = sub[1]
                            gene_orders[gid] = order
                    except: pass
    
    # 3. Analyze Operon Synchronization
    print("\n--- Operon Synchronization ---")
    simultaneous_count = 0
    sequential_count = 0
    lags = []
    
    for op, members in operons.items():
        if len(members) < 2: continue
        
        # Get orders
        member_orders = []
        for m in members:
            if m in gene_orders:
                member_orders.append(gene_orders[m])
        
        member_orders.sort()
        
        if len(member_orders) < 2: continue
        
        # Check diffs
        is_simul = True
        for i in range(len(member_orders)-1):
            diff = member_orders[i+1] - member_orders[i]
            if diff != 0:
                is_simul = False
            lags.append(diff * 0.3) # 0.3s per order unit
            
        if is_simul:
            simultaneous_count += 1
        else:
            sequential_count += 1
            
    print(f"Operons with Simultaneous Start (Same Order): {simultaneous_count}")
    print(f"Operons with Sequential Start (Diff Order): {sequential_count}")
    
    if lags:
        avg_lag = sum(lags) / len(lags)
        print(f"Average Time Lag between operon genes: {avg_lag:.2f} seconds")
    else:
        print("Could not calculate lags (no matching orders found).")

def analyze_pause_decay():
    print("\n--- Pause Experiment Simulation ---")
    # Constants
    mrna_hl = 300.0 # 5 mins
    prot_hl = 3600.0 # 60 mins
    
    mrna_k = math.log(2) / mrna_hl
    prot_k = math.log(2) / prot_hl
    
    print(f"mRNA Decay Constant (k): {mrna_k:.6f} / sec")
    
    # Time to 50%
    # [A] = [A]0 * e^(-kt)
    # 0.5 = e^(-kt) -> t = ln(2)/k = HL
    print(f"Time to drop to 50% mRNA: {mrna_hl/60:.1f} minutes")
    
    # Protein Persistence
    # Simulating protein decay after mRNA depletion
    # Assuming mRNA depletes instantly for this question ("after mRNA is depleted")
    # Protein decays with half-life 60 mins.
    
    # Time to 10%
    # 0.1 = e^(-k_prot * t) -> ln(0.1) = -k * t -> t = -ln(0.1)/k
    t_10 = -math.log(0.1) / prot_k
    print(f"Time for Protein to drop to 10%: {t_10/60:.1f} minutes")
    
    # Time to 1%
    t_1 = -math.log(0.01) / prot_k
    print(f"Time for Protein to drop to 1%: {t_1/60:.1f} minutes")

def check_ribosome_usage():
    print("\n--- Ribosome Saturation Logic ---")
    # Based on code analysis, ribosomes_in_use is initialized to 0
    # and only modified by allocate_ribosomes().
    # Using 'grep' we found allocate_ribosomes() is defined but NEVER CALLED.
    # Therefore:
    print("Max Ribosome Saturation Reached: 0.0%")
    print("Reason: Resource allocation logic is present but not connected to the simulation loop.")
    print("Consequence: Translation rates are effectively unconstrained (100% capacity always).")

if __name__ == "__main__":
    analyze_operons()
    analyze_pause_decay()
    check_ribosome_usage()
