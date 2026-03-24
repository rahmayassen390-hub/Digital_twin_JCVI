import requests
import time
import json
import random
import sys

# Configuration
USER_API_URL = "http://localhost:8081/jobs"      # For reading jobs
PROVIDER_API_URL = "http://localhost:8888/jobs"  # For updating jobs
POLL_INTERVAL = 1.0  # seconds

def get_queued_jobs():
    try:
        response = requests.get(f"{USER_API_URL}?status=submitted")
        if response.status_code == 200:
            return response.json()
        return []
    except Exception as e:
        print(f"Error polling jobs: {e}")
        return []

def update_job_status(job_id, status, result=None):
    """Update job status via Provider API"""
    try:
        # Update status
        status_resp = requests.patch(
            f"{PROVIDER_API_URL}/{job_id}/status",
            json={"status": status}
        )
        if status_resp.status_code not in [200, 204]:
            print(f"Warning: Status update returned {status_resp.status_code}")
        
        # If we have results, update job_info
        if result and status == "succeeded":
            # Put result in job_info.result.sampling format
            result_data = json.loads(result) if isinstance(result, str) else result
            job_info = {
                "result": {
                    "sampling": {
                        "counts": result_data.get("counts", {})
                    }
                }
            }
            info_resp = requests.patch(
                f"{PROVIDER_API_URL}/{job_id}/job_info",
                json=job_info
            )
            if info_resp.status_code not in [200, 204]:
                print(f"Warning: job_info update returned {info_resp.status_code}")
                
    except Exception as e:
        print(f"Error updating job {job_id}: {e}")

def run_mock_worker():
    print("🤖 Mock Simulator Worker Started (28-qubit capable)")
    print(f"Polling {USER_API_URL} every {POLL_INTERVAL}s...")
    print(f"Updating via Provider API: {PROVIDER_API_URL}")
    
    while True:
        jobs = get_queued_jobs()
        
        if jobs:
            print(f"Found {len(jobs)} submitted jobs.")
            for job in jobs:
                job_id = job.get('job_id') or job.get('id')
                job_name = job.get('name', 'unknown')
                print(f"Processing job {job_id} ({job_name})...")
                
                # 1. Mark RUNNING
                update_job_status(job_id, "running")
                
                # 2. Simulate processing delay
                time.sleep(0.5)
                
                # 3. Generate Result based on job type
                if "cluster" in job_name.lower() or "28q" in job_name.lower():
                    # 28-qubit cluster simulation
                    result_data = generate_28qubit_result()
                    print(f"Job {job_id}: 28-qubit cluster result generated")
                else:
                    # Single qubit job (legacy)
                    prob = random.uniform(70.0, 99.0)
                    counts = {"0": int(1000 * (100-prob)/100), "1": int(1000 * prob/100)}
                    result_data = {"probability": prob, "counts": counts}
                    print(f"Job {job_id}: single qubit p={prob:.1f}%")
                
                # 4. Mark COMPLETED (use 'succeeded' for OQTOPUS)
                update_job_status(job_id, "succeeded", result=json.dumps(result_data))
        
        time.sleep(POLL_INTERVAL)


def generate_28qubit_result(num_qubits=28, num_shots=4096):
    """Generate realistic 28-qubit measurement results.
    
    Simulates a 28-qubit circuit where each qubit has varying probability
    of being measured as |1⟩ (representing cluster transcription level).
    """
    counts = {}
    
    # Generate per-qubit probabilities (simulating cluster transcription levels)
    # Mix of high (ribosome, translation) and low (stress, hypothetical) activity
    qubit_probs = []
    for i in range(num_qubits):
        if i < 5:  # Core machinery clusters - high activity
            p = random.uniform(0.6, 0.9)
        elif i < 15:  # Biosynthesis clusters - medium activity
            p = random.uniform(0.4, 0.7)
        else:  # Support clusters - lower activity
            p = random.uniform(0.2, 0.5)
        qubit_probs.append(p)
    
    # Generate shots
    for _ in range(num_shots):
        bitstring = ""
        for p in qubit_probs:
            bit = "1" if random.random() < p else "0"
            bitstring += bit
        
        # Reverse bitstring (q[0] is rightmost)
        bitstring = bitstring[::-1]
        counts[bitstring] = counts.get(bitstring, 0) + 1
    
    # Calculate per-qubit probabilities for backward compatibility
    prob = sum(qubit_probs) / num_qubits * 100  # Average probability
    
    return {"probability": prob, "counts": counts}


if __name__ == "__main__":
    try:
        run_mock_worker()
    except KeyboardInterrupt:
        print("Worker stopped.")
