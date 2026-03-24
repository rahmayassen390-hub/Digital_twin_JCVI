#!/usr/bin/env python3
"""
Standalone Quantum API Server
Replaces OQTOPUS for local testing of 28-qubit cluster circuits
Run with: python standalone_quantum_api.py
"""

from flask import Flask, request, jsonify
import json
import random
import threading
import time
import uuid

app = Flask(__name__)

# In-memory job storage
jobs = {}
job_lock = threading.Lock()


def generate_28qubit_result(num_qubits=28, num_shots=4096):
    """Generate realistic 28-qubit measurement results."""
    counts = {}
    qubit_probs = []
    
    for i in range(num_qubits):
        if i < 5:
            p = random.uniform(0.6, 0.9)
        elif i < 15:
            p = random.uniform(0.4, 0.7)
        else:
            p = random.uniform(0.2, 0.5)
        qubit_probs.append(p)
    
    for _ in range(num_shots):
        bitstring = ""
        for p in qubit_probs:
            bit = "1" if random.random() < p else "0"
            bitstring += bit
        bitstring = bitstring[::-1]
        counts[bitstring] = counts.get(bitstring, 0) + 1
    
    prob = sum(qubit_probs) / num_qubits * 100
    return {"probability": prob, "counts": counts}


def process_job(job_id):
    """Background job processor."""
    time.sleep(1.5)  # Simulate processing
    
    with job_lock:
        if job_id in jobs:
            job = jobs[job_id]
            job_name = job.get('name', '')
            
            if "cluster" in job_name.lower() or "28q" in job_name.lower():
                result = generate_28qubit_result()
            else:
                prob = random.uniform(70.0, 99.0)
                result = {"probability": prob, "counts": {"0": int(1000*(100-prob)/100), "1": int(1000*prob/100)}}
            
            job['status'] = 'COMPLETED'
            job['result'] = json.dumps(result)
            print(f"✅ Job {job_id} completed")


@app.route('/jobs', methods=['POST'])
def submit_job():
    """Submit a new quantum job."""
    data = request.json
    job_id = str(uuid.uuid4())
    
    with job_lock:
        jobs[job_id] = {
            'id': job_id,
            'name': data.get('name', 'unnamed'),
            'status': 'QUEUED',
            'result': None,
            **data
        }
    
    # Process in background
    threading.Thread(target=process_job, args=(job_id,), daemon=True).start()
    
    print(f"📥 Job submitted: {job_id} ({data.get('name', 'unnamed')})")
    return jsonify({'id': job_id}), 200


@app.route('/jobs/<job_id>', methods=['GET'])
def get_job(job_id):
    """Get job status and result."""
    with job_lock:
        if job_id in jobs:
            return jsonify(jobs[job_id]), 200
    return jsonify({'error': 'Job not found'}), 404


@app.route('/jobs', methods=['GET'])
def list_jobs():
    """List jobs (optionally filtered by status)."""
    status = request.args.get('status')
    with job_lock:
        if status:
            filtered = [j for j in jobs.values() if j['status'] == status]
        else:
            filtered = list(jobs.values())
    return jsonify(filtered), 200


@app.route('/health', methods=['GET'])
def health():
    return jsonify({'status': 'ok'}), 200


if __name__ == '__main__':
    print("🚀 Standalone Quantum API Server")
    print("   Listening on http://localhost:8081")
    print("   Supports 28-qubit cluster circuits")
    app.run(host='0.0.0.0', port=8081, debug=False, threaded=True)
