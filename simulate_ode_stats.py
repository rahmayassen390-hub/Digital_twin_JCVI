
import math
import random
import statistics
import collections

# === SIMULATION ENGINE REPLICA ===
# Replicating logic from cell_simulation_engine.py / gui_methods.py

class GeneState:
    def __init__(self, gene_id, tx_rate=1.0, tl_rate=0.1, mrna_hl=300.0, prot_hl=3600.0, noise_level=0.05):
        self.gene_id = gene_id
        
        # Parameters
        self.tx_rate = tx_rate
        self.tl_rate = tl_rate
        self.mrna_hl = mrna_hl
        self.prot_hl = prot_hl
        self.noise_level = noise_level
        
        # State
        self.mrna = 0.0
        self.protein = 0.0
        
        # Decay constants
        self.mrna_k = math.log(2) / mrna_hl if mrna_hl > 0 else 0
        self.prot_k = math.log(2) / prot_hl if prot_hl > 0 else 0
        
        # Steady states
        self.mrna_ss = tx_rate / self.mrna_k if self.mrna_k else 0
        self.prot_ss = (tl_rate * self.mrna_ss) / self.prot_k if self.prot_k else 0
        
        # History for convergence check
        self.history = collections.deque(maxlen=60) # Store last 60 steps (seconds) to check stability
        
        self.steady_state_reached_time = None

    def apply_noise(self, val, scale, enable_noise):
        if not enable_noise or self.noise_level <= 0:
            return val
        
        # Logic from cell_simulation_engine.py:
        # std_dev = effective_noise * value * noise_scale
        std_dev = self.noise_level * val * scale
        noise = random.gauss(0, std_dev) if std_dev > 0 else 0
        return max(0.0, val + noise)

    def step(self, dt, time, enable_noise=True):
        # mRNA Logic
        # Synthesis
        synthesis = self.tx_rate * dt
        synthesis = self.apply_noise(synthesis, 0.5, enable_noise)
        
        # Degradation
        decay = self.mrna_k * self.mrna * dt
        decay = self.apply_noise(decay, 0.3, enable_noise)
        
        self.mrna = max(0.0, self.mrna + synthesis - decay)
        
        # Protein Logic
        # Synthesis (TL = 0 if mRNA = 0)
        if self.mrna > 0:
            prot_syn_base = self.tl_rate * self.mrna * dt
            prot_syn = self.apply_noise(prot_syn_base, 0.5, enable_noise)
        else:
            prot_syn = 0.0
            
        # Degradation
        prot_decay_base = self.prot_k * self.protein * dt
        prot_decay = self.apply_noise(prot_decay_base, 0.3, enable_noise)
        
        self.protein = max(0.0, self.protein + prot_syn - prot_decay)
        
        # Check Steady State Convergence
        # "Steady: production ≈ degradation (within 1-5%)" - from code comments
        # Here we define metric: consistent values close to theoretical SS?
        # Or net balance approach from cell_simulation_engine.py
        
        net_prot = prot_syn - prot_decay
        avg_prot = self.protein
        
        threshold = 0.01 * avg_prot + 0.001
        
        if self.steady_state_reached_time is None:
            if avg_prot > 0.95 * self.prot_ss and abs(net_prot) < threshold:
                # Simple check: reached 95% of SS and net change is low
                self.steady_state_reached_time = time

def run_simulation(dt=1.0, duration=20000, noise=0.05, enable_noise=True, seed=None):
    if seed is not None:
        random.seed(seed)
    
    # Create representative genes
    # 1. Average Gene (TX=1.0, TL=0.1)
    # 2. Fast Gene (TX=2.0)
    # 3. Slow Gene (TX=0.5)
    
    genes = [
        GeneState("Avg", tx_rate=1.0, tl_rate=0.1, noise_level=noise),
        GeneState("Fast", tx_rate=2.0, tl_rate=0.1, noise_level=noise),
        GeneState("Slow", tx_rate=0.5, tl_rate=0.1, noise_level=noise)
    ]
    
    time = 0.0
    steps = int(duration / dt)
    
    steady_times = {}
    
    # Trace for Avg gene for variation analysis
    avg_gene_trace = []
    
    for _ in range(steps):
        time += dt
        for g in genes:
            g.step(dt, time, enable_noise)
            if g.steady_state_reached_time and g.gene_id not in steady_times:
                steady_times[g.gene_id] = g.steady_state_reached_time
        
        avg_gene_trace.append(genes[0].protein)
        
    return {
        'steady_times': steady_times,
        'final_states': {g.gene_id: g.protein for g in genes},
        'trace': avg_gene_trace
    }

def analyze_stats():
    print("--- ODE Simulation Analysis ---")
    
    # 1. Simulation Settings & Determinism
    print("\n[1. Simulation Settings & Determinism]")
    print("Checking variation between two identical runs (Noise=0)...")
    
    # Run 1
    run1 = run_simulation(dt=1.0, noise=0.0, enable_noise=False, seed=42)
    # Run 2
    run2 = run_simulation(dt=1.0, noise=0.0, enable_noise=False, seed=42)
    
    diffs = []
    for i in range(len(run1['trace'])):
        diffs.append(abs(run1['trace'][i] - run2['trace'][i]))
    
    max_diff = max(diffs)
    print(f"Max difference between runs (Noise=0): {max_diff:.10f}")
    if max_diff < 1e-9:
        print(">> Results are EXACTLY the same.")
    else:
        print(f">> Results differ by {max_diff}")
        
    print(f"\nExact Time Step Used: The code uses `dt = 0.1 * speed`. Default speed is usually 1.0 or 10.0.")
    print(f"If you see 1.0s, the time step is 1.0s.")

    # 2. Steady State Convergence
    print("\n[2. Steady-State Convergence]")
    # Run for long duration to ensure SS
    long_run = run_simulation(dt=1.0, duration=50000, noise=0.0, enable_noise=False)
    
    times = long_run['steady_times']
    print("Time to reach steady state (95% + balance check):")
    for gid, t in times.items():
        print(f"  {gid} Gene: {t:.1f} seconds ({t/60:.1f} minutes)")
    
    if not times:
        print("  (Did not reach steady state in timeframe)")
    
    avg_mins = sum(times.values()) / len(times) / 60 if times else 0
    print(f"Average time: {avg_mins:.1f} minutes")

    # 3. Noise Level Testing
    print("\n[3. Noise Level Testing]")
    noise_levels = [0.0, 0.05, 0.10]
    
    for n in noise_levels:
        print(f"\nTesting Noise Level = {n}")
        # Run 3 trials
        trials = []
        for i in range(3):
            res = run_simulation(dt=1.0, duration=10000, noise=n, enable_noise=True, seed=None) # Random seed
            trials.append(res['final_states']['Avg'])
        
        mean_val = statistics.mean(trials)
        stdev_val = statistics.stdev(trials) if len(trials) > 1 else 0
        cv = (stdev_val / mean_val * 100) if mean_val > 0 else 0
        
        print(f"  Trials (Final Protein Count): {[f'{x:.1f}' for x in trials]}")
        print(f"  Mean: {mean_val:.1f}, StdDev: {stdev_val:.1f}")
        print(f"  Variation between runs (CV): {cv:.2f}%")

if __name__ == "__main__":
    analyze_stats()
