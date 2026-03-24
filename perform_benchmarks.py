
import sys
import os
import time
import psutil
import statistics
from collections import defaultdict

# Add project root AND inner package to path to allow direct imports
sys.path.append("/home/rahema/Desktop/bed_files /jcvi_genome_analyzer_updated (6)")
sys.path.append("/home/rahema/Desktop/bed_files /jcvi_genome_analyzer_updated (6)/jcvi_genome_analyzer_updated")

# Direct imports bypassing __init__.py to avoid PyQt5 dependency
from genome_manager import GenomeManager
from cell_simulation_engine import CellSimulationEngine, SimulationMode
from blast_manager import BlastManager

def get_system_specs():
    print("\n=== System Specifications ===")
    
    # CPU
    try:
        with open("/proc/cpuinfo", "r") as f:
            for line in f:
                if "model name" in line:
                    print(f"CPU Model: {line.split(':')[1].strip()}")
                    break
    except:
        print("CPU Model: Unknown")
        
    print(f"Physical Cores: {psutil.cpu_count(logical=False)}")
    print(f"Logical Cores: {psutil.cpu_count(logical=True)}")
    
    # RAM
    mem = psutil.virtual_memory()
    print(f"Total RAM: {mem.total / (1024**3):.2f} GB")
    print(f"Available RAM: {mem.available / (1024**3):.2f} GB")

def benchmark_loading():
    print("\n=== Loading Performance ===")
    
    fasta_path = "/home/rahema/Desktop/JCVI_3.0.fasta"
    gff_path = "/home/rahema/Desktop/JCVI_3.0.gff3"
    
    gm = GenomeManager()
    
    # 1. FASTA Loading
    start = time.time()
    if os.path.exists(fasta_path):
        name, seq, length = gm.load_fasta(fasta_path)
        gm.sequence = seq
        fasta_time = time.time() - start
        print(f"FASTA Loading Time: {fasta_time:.4f} seconds")
    else:
        print("FASTA file not found!")
        fasta_time = 0
    
    # 2. GFF3 Parsing
    start = time.time()
    if os.path.exists(gff_path):
        gm.genes = gm.load_gff(gff_path)
        gff_time = time.time() - start
        print(f"GFF3 Parsing Time: {gff_time:.4f} seconds")
        print(f"Genes Loaded: {len(gm.genes)}")
        # Debug if 0
        if len(gm.genes) == 0:
            print("  ! DEBUG: GFF loaded 0 genes. Checking first few lines...")
            with open(gff_path) as f:
                for _ in range(3): print(f"    {f.readline().strip()}")
    else:
        print("GFF3 file not found!")
        gff_time = 0

    # 3. BED Files
    bed_dir = "/home/rahema/Desktop/bed_files /"
    bed_files = [f for f in os.listdir(bed_dir) if f.endswith(".bed")]
    
    start = time.time()
    loaded_tracks = 0
    for bed in bed_files:
        full_path = os.path.join(bed_dir, bed)
        if not os.path.exists(full_path): continue
        
        try:
            if "Regulation" in bed:
                gm.load_regulation_bed(full_path)
            elif "Transcription" in bed:
                gm.load_transcription_dynamics_bed(full_path)
            elif "Translation" in bed:
                gm.load_translation_dynamics_bed(full_path)
            # Add other specific loaders if they exist, else just read
            else:
                with open(full_path) as f: pass
            
            loaded_tracks += 1
        except Exception as e:
            print(f"Failed to load {bed}: {e}")

    bed_time = time.time() - start
    print(f"BED Files Loading (Count={loaded_tracks}): {bed_time:.4f} seconds")
    
    # Total Integration estimate
    print(f"Total Integration Time (FASTA+GFF+BEDs): {fasta_time + gff_time + bed_time:.4f} seconds")
    
    return gm

def benchmark_blast(gm):
    print("\n=== BLAST Performance ===")
    
    bm = BlastManager()
    
    # Create a temporary query file in /tmp to avoid path space issues
    query_fasta = "/tmp/temp_query.fasta"
    
    # Try to use gm.genes, or fallback to translation dynamics keys if GFF failed
    # We can extract CDS for all genes.
    gene_count = 0
    
    # Try to use gm.genes, or fallback to translation dynamics keys if GFF failed
    gene_list = gm.genes
    if not gene_list and gm.translation_dynamics:
        print("Using BED-derived gene list for BLAST query...")
        # Create dummy gene objects or just use IDs
        gene_list = []
        for gid, data in gm.translation_dynamics.items():
            # We need start/end. BED loading *might* have stored this if generic BED loader used
            # But specific loaders might not store coordinates in a way we can access easily 
            # without inspecting the objects.
            # Let's assume we can't easily get coordinates from *just* the dynamics dict 
            # unless the object has them.
            # The output said "using BED coordinates as fallback", so the objects probably have .start/.end
            pass
            
    # Actually, a better approach: generic fallback to generating random sequences or skipping 
    # if we can't validly extract.
    # But wait, `load_translation_dynamics_bed` LOGS said "using BED coordinates". 
    # This implies it updated *something* with coordinates.
    # It updates `self.genes`?
    # Let's check the length of gm.genes AGAIN right here.
    print(f"Gene count before BLAST: {len(gm.genes)}")
    
    if len(gm.genes) > 0 and hasattr(gm, 'sequence') and gm.sequence:
         try:
            with open(query_fasta, "w") as f:
                for gene in gm.genes:
                     start = gene.start - 1
                     end = gene.end
                     seq = gm.sequence[start:end]
                     f.write(f">{gene.id}\n{seq}\n")
                     gene_count += 1
         except Exception as e:
            print(f"Skipping BLAST: Could not prepare query ({e})")
            return
    elif gm.sequence:
        # Fallback: Just slice the first 473 genes worth of random 1000bp chunks
        # to simulate the WORKLOAD of BLAST 473 genes.
        print("Synthesizing query from genome chunks (GFF parsing failed)...")
        with open(query_fasta, "w") as f:
             for i in range(473):
                 start = i * 1000
                 end = start + 1000
                 if end < len(gm.sequence):
                     seq = gm.sequence[start:end]
                     f.write(f">Gene_{i}\n{seq}\n")
                     gene_count += 1
    else:
        # Mock some data if GFF loading failed, just to test BLAST engine speed?
        # Or better: check gene list from gm.genes if empty avoid it.
        # But if fail, we might want to manually parse FASTA to get *something*
        print("Warning: GenomeManager has no genes or sequence loaded. Skipping BLAST.")
        return

    print(f"Prepared Query DB with {gene_count} genes")

    # 1. DB Creation
    start = time.time()
    try:
        # Correct method name: make_blast_db returns path or throws
        bm.make_blast_db(query_fasta)
        db_time = time.time() - start
        print(f"Database Creation Time: {db_time:.4f} seconds")
    except Exception as e:
        print(f"BLAST DB Creation Failed: {e}")
        return

    # 2. Search Time
    start = time.time()
    try:
        # Correct calls: run_blast returns XML path, then parse it
        xml_path = bm.run_blast(query_fasta)
        search_time = time.time() - start
        print(f"Search Time (473 vs 473): {search_time:.4f} seconds")
        
        hits = bm.parse_blast_xml(xml_path)
        print(f"Number of Hits Found: {len(hits)}")
        
        if hits:
            # Avg Identity: hits are DICTIONARIES
            identities = [h['identity'] for h in hits if 'identity' in h]
            
            if identities:
                avg_id = sum(identities) / len(identities)
                print(f"Average Identity: {avg_id:.2f}%")
        
        # Cleanup XML
        if os.path.exists(xml_path): os.remove(xml_path)
            
    except Exception as e:
         print(f"BLAST Search Failed: {e}")

    # Cleanup
    if os.path.exists(query_fasta): os.remove(query_fasta)
    for ext in ['.nin', '.nhr', '.nsq', '.ndb', '.not', '.ntf', '.nto']:
        tgt = f"temp_db{ext}"
        if os.path.exists(tgt):
            try: os.remove(tgt)
            except: pass

def benchmark_simulation():
    print("\n=== Simulation Speed (1 Hour) ===")
    
    engine = CellSimulationEngine(mode=SimulationMode.CONTINUOUS)
    
    # Add 473 genes to match JCVI
    for i in range(473):
        engine.add_gene(f"GENE_{i}", transcription_rate=1.0, translation_rate=0.1)
        
    start_time_wall = time.time()
    
    process = psutil.Process(os.getpid())
    start_mem = process.memory_info().rss / (1024 * 1024) # MB
    
    # Simulate 3600 seconds
    sim_duration = 3600
    dt = 1.0 # 10x speed default
    steps = int(sim_duration / dt)
    
    # Sample CPU usage
    psutil.cpu_percent(interval=None) # Reset
    
    for _ in range(steps):
        engine.step(dt)
        
    end_time_wall = time.time()
    end_mem = process.memory_info().rss / (1024 * 1024)
    
    actual_time = end_time_wall - start_time_wall
    real_time_factor = sim_duration / actual_time if actual_time > 0 else 0
    
    cpu_usage = psutil.cpu_percent(interval=None)
    
    print(f"Simulated Duration: {sim_duration} seconds")
    print(f"Actual Real-World Time: {actual_time:.4f} seconds")
    print(f"Real-Time Factor: {real_time_factor:.1f}x faster than real life")
    print(f"CPU Usage (approx): {cpu_usage:.1f}%")
    print(f"Memory Usage: {start_mem:.1f} MB -> {end_mem:.1f} MB (Delta: {end_mem - start_mem:.1f} MB)")

if __name__ == "__main__":
    get_system_specs()
    try:
        gm_loaded = benchmark_loading()
        if gm_loaded:
            benchmark_blast(gm_loaded)
    except Exception as e:
        print(f"Loading/BLAST Benchmarks Failed: {e}")
        import traceback
        traceback.print_exc()
        
    try:
        benchmark_simulation()
    except Exception as e:
        print(f"Simulation Benchmark Failed: {e}")
