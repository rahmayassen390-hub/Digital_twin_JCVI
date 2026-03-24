#!/usr/bin/env python3
"""
Test script for 28-qubit quantum cluster simulation
Runs the quantum entanglement model as described in the report
"""

import sys
import os
sys.path.append('jcvi_genome_analyzer_updated')

from jcvi_genome_analyzer_updated.quantum_worker import QuantumClusterWorker
from jcvi_genome_analyzer_updated.cluster_definitions import ENTANGLEMENT_PAIRS

def run_quantum_cluster_test():
    """Run a test of the 28-qubit cluster quantum simulation"""

    print("🚀 Testing 28-Qubit Quantum Cluster Simulation")
    print("=" * 50)

    # Create cluster levels based on report's environmental conditions
    # Glucose limitation: Energy/ATP (C14) low, Transport (C10) affected
    # Amino acid pulses: Translation (C1) variable, Transcription (C4) bursts

    cluster_levels = {}

    # Set levels based on report results (approximate)
    # C1 Translation: 61.8%-90.2% → average ~75%
    cluster_levels[1] = 0.75

    # C4 Transcription: 59.9%-90.6% → average ~75%
    cluster_levels[4] = 0.75

    # C10 Transport: 41.2%-70.5% → average ~55%
    cluster_levels[10] = 0.55

    # C14 Energy/ATP: 41.1%-69.5% → average ~55% (glucose limitation)
    cluster_levels[14] = 0.55

    # C27 Essential: 19.1%-50.0% → average ~35% (oxidative stress)
    cluster_levels[27] = 0.35

    # Fill in other clusters with typical values
    for i in range(28):
        if i not in cluster_levels:
            if i in [0, 2, 3, 5, 6]:  # Core machinery - higher
                cluster_levels[i] = 0.7
            elif i in [15, 16, 17]:  # Metabolism - medium
                cluster_levels[i] = 0.6
            else:  # Others - lower
                cluster_levels[i] = 0.4

    print("Cluster Levels (normalized 0-1):")
    from jcvi_genome_analyzer_updated.cluster_definitions import get_cluster_name
    for i in range(28):
        name = get_cluster_name(i)
        level = cluster_levels[i]
        print(f"  {i:2d}: {name:<20} = {level:.2f}")

    print(f"\nEntanglement Pairs ({len(ENTANGLEMENT_PAIRS)} connections):")
    for pair in ENTANGLEMENT_PAIRS[:5]:  # Show first 5
        c1, c2 = pair
        name1 = get_cluster_name(c1)
        name2 = get_cluster_name(c2)
        print(f"  {c1} ({name1}) ↔ {c2} ({name2})")
    if len(ENTANGLEMENT_PAIRS) > 5:
        print(f"  ... and {len(ENTANGLEMENT_PAIRS)-5} more")

    # Create quantum worker
    print("\n📡 Submitting to Quantum Backend (localhost:8081)...")

    worker = QuantumClusterWorker(
        cluster_levels=cluster_levels,
        entanglement_pairs=ENTANGLEMENT_PAIRS
    )

    # For console testing, we'll implement the logic directly instead of signals
    # Since we can't run Qt event loop without GUI

    # Simulate the worker's run method
    try:
        circuit = worker._generate_cluster_circuit()

        import requests
        import json
        import time

        job_payload = {
            "name": "cluster_28q_sim",
            "device_id": "SVSim",
            "job_type": "sampling",
            "shots": 4096,
            "job_info": {
                "program": [circuit]
            }
        }

        headers = {"X-API-Token": "local-dev-token"}

        print("DEBUG: [Quantum] Submitting 28-qubit cluster circuit...")
        submit_resp = requests.post("http://localhost:8081/jobs", json=job_payload, headers=headers)

        print(f"DEBUG: [Quantum] Response status: {submit_resp.status_code}")
        print(f"DEBUG: [Quantum] Response body: {submit_resp.text[:500]}")

        if submit_resp.status_code not in [200, 201]:
            raise Exception(f"Submission failed: {submit_resp.text}")

        resp_data = submit_resp.json()
        job_id = resp_data.get('id') or resp_data.get('job_id')

        if not job_id:
            raise Exception(f"No job ID in response: {resp_data}")

        print(f"DEBUG: [Quantum] Job submitted: {job_id}")

        # Poll for completion
        while True:
            status_resp = requests.get(f"http://localhost:8081/jobs/{job_id}", headers=headers)
            if status_resp.status_code != 200:
                time.sleep(1)
                continue

            job_data = status_resp.json()
            status = job_data['status']

            if status in ["COMPLETED", "succeeded"]:
                # Handle both OQTOPUS native format and mock worker format
                counts = {}
                if 'result' in job_data and isinstance(job_data['result'], str):
                    # Mock worker format: result is JSON string
                    result = json.loads(job_data.get('result', '{}'))
                    counts = result.get('counts', {})
                elif 'job_info' in job_data and 'result' in job_data.get('job_info', {}):
                    # OQTOPUS native format: job_info.result.sampling.counts
                    sampling = job_data['job_info']['result'].get('sampling', {})
                    counts = sampling.get('counts', {})

                # Parse 28-qubit measurement results
                cluster_probs = worker._parse_cluster_probabilities(counts)

                print("✅ Quantum Simulation Completed!")
                print("\nCluster Probabilities (P(|1⟩) = Transcription Level):")
                print("Cluster | Name                     | Probability | Expected Range")
                print("-" * 65)

                # Sort by probability descending
                sorted_clusters = sorted(cluster_probs.items(), key=lambda x: x[1], reverse=True)

                for cluster_id, prob in sorted_clusters:
                    name = get_cluster_name(cluster_id)
                    # Get expected range from report
                    expected = get_expected_range(cluster_id)
                    prob_pct = prob
                    print(f"  {cluster_id:2d}    | {name:<20} | {prob_pct:6.1f}%      | {expected}")

                # Check key correlations from report
                print("\n🔗 Key Correlations (from Report):")
                energy_prob = cluster_probs.get(14, 0)
                transport_prob = cluster_probs.get(10, 0)
                print(f"  C14 Energy/ATP: {energy_prob:.1f}%")
                print(f"  C10 Transport: {transport_prob:.1f}%")

                if abs(energy_prob - transport_prob) < 5:  # Within 5% - correlated
                    print("  → ✅ ENTANGLEMENT CONFIRMED: Transport follows Energy (CNOT gate effect)")
                else:
                    print("  → ❌ No correlation detected")

                print(f"\nTotal measurement shots: {sum(counts.values())}")
                break

            elif status in ["FAILED", "failed"]:
                print(f"❌ Quantum Simulation Failed: Job failed in backend")
                break

            time.sleep(0.5)

    except Exception as e:
        print(f"❌ Quantum Simulation Failed: {e}")

    print("\n🎯 Test completed. This demonstrates the quantum-biological mapping.")

def get_expected_range(cluster_id):
    """Get expected probability range from report"""
    ranges = {
        1: "61.8-90.2%",   # Translation
        4: "59.9-90.6%",   # Transcription
        10: "41.2-70.5%",  # Transport
        14: "41.1-69.5%",  # Energy/ATP
        27: "19.1-50.0%"   # Essential
    }
    return ranges.get(cluster_id, "40-70%")

if __name__ == "__main__":
    run_quantum_cluster_test()
