import metabolic_solver
import numpy as np

def probe_network():
    solver = metabolic_solver.MetabolicSolver('../Data/metabolic_reconstruction.xlsx')
    
    dead_ends = []
    no_sinks = []
    
    for i, m_id in enumerate(solver.metabolites):
        row = solver.s_matrix[i, :]
        producers = np.where(row > 0)[0]
        consumers = np.where(row < 0)[0]
        
        # Dead end: no way to produce it
        if len(producers) == 0:
            dead_ends.append(m_id)
            
        # No sink: no way to consume it
        if len(consumers) == 0:
            no_sinks.append(m_id)
            
    print(f"Total dead-ends (no producers): {len(dead_ends)}")
    print(f"Total no-sinks (no consumers): {len(no_sinks)}")
    
    # Check if any dead-end is in BIOMASS
    if 'BIOMASS' in solver.reactions:
        idx = solver.reactions.index('BIOMASS')
        col = solver.s_matrix[:, idx]
        biomass_metabs = [solver.metabolites[m_i] for m_i in np.where(col != 0)[0]]
        
        biomass_dead_ends = [m for m in biomass_metabs if m in dead_ends]
        print(f"\nBiomass metabolites that are dead-ends: {len(biomass_dead_ends)}")
        if biomass_dead_ends:
            print(f"  {biomass_dead_ends[:10]}")
            
    # Check for metabolites that are produced but NEVER consumed (must be exchanges or products)
    non_exchange_no_sinks = [m for m in no_sinks if not any(r.startswith('EX_') for r in [solver.reactions[p_i] for p_i in np.where(solver.s_matrix[solver.metabolites.index(m), :] > 0)[0]])]
    # Wait, simpler logic:
    real_no_sinks = []
    for m in no_sinks:
        m_idx = solver.metabolites.index(m)
        is_exchange = False
        for r_idx in np.where(solver.s_matrix[m_idx, :] != 0)[0]:
            if solver.reactions[r_idx].startswith('EX_'):
                is_exchange = True
                break
        if not is_exchange:
            real_no_sinks.append(m)
            
    print(f"\nIntracellular metabolites with NO sinks: {len(real_no_sinks)}")
    if real_no_sinks:
        print(f"  {real_no_sinks[:20]}")

if __name__ == "__main__":
    probe_network()
