import json
import time
import requests
from PyQt5.QtCore import QThread, pyqtSignal

class QuantumJobWorker(QThread):
    """
    Worker to submit quantum jobs to OQTOPUS and poll for results.
    Running in a separate thread to keep GUI responsive.
    """
    job_finished = pyqtSignal(str, float, dict)  # gene_id, probability, counts
    job_failed = pyqtSignal(str, str)            # gene_id, error_message
    
    # OQTOPUS API Config
    API_URL = "http://localhost:8081/jobs"
    
    def __init__(self, gene_id, transcription_rate):
        super().__init__()
        self.gene_id = gene_id
        self.transcription_rate = transcription_rate
        self.running = True
        
    def run(self):
        try:
            # 1. Submit Job
            job_payload = {
                "name": f"gene_{self.gene_id}_sim",
                "device_id": "SVSim",
                "job_type": "sampling",
                "job_info": {
                    "program": self._generate_circuit(self.transcription_rate),
                    "shots": 1000
                }
            }
            
            # Auto-auth header (simulated for local mode)
            headers = {"X-API-Token": "local-dev-token"}
            
            submit_resp = requests.post(self.API_URL, json=job_payload, headers=headers, timeout=5)
            if submit_resp.status_code != 200:
                raise Exception(f"Submission failed: {submit_resp.text}")
                
            job_id = submit_resp.json()['id']
            
            while self.running:
                status_resp = requests.get(f"{self.API_URL}/{job_id}", headers=headers, timeout=5)
                if status_resp.status_code != 200:
                    time.sleep(1)
                    continue
                    
                job_data = status_resp.json()
                status = job_data['status']
                
                if status in ["COMPLETED", "succeeded"]:
                    # Handle both OQTOPUS native format and mock worker format
                    result = {}
                    if 'result' in job_data and isinstance(job_data['result'], str):
                        # Mock worker format: result is JSON string
                        result = json.loads(job_data.get('result', '{}'))
                    elif 'job_info' in job_data and 'result' in job_data.get('job_info', {}):
                        # OQTOPUS native format: job_info.result.sampling.counts
                        sampling = job_data['job_info']['result'].get('sampling', {})
                        counts = sampling.get('counts', {})
                        total = sum(counts.values()) if counts else 1
                        prob = (counts.get('1', 0) / total) * 100 if total > 0 else 0
                        result = {'probability': prob, 'counts': counts}
                    
                    prob = result.get('probability', 0.0)
                    counts = result.get('counts', {})
                    
                    self.job_finished.emit(self.gene_id, prob, counts)
                    return
                    
                elif status in ["FAILED", "failed"]:
                    self.job_failed.emit(self.gene_id, "Job failed in backend")
                    return
                
                time.sleep(0.5)
                
        except Exception as e:
            if self.running:
                self.job_failed.emit(self.gene_id, str(e))
    
    def stop(self):
        self.running = False
        
    def _generate_circuit(self, rate):
        """Generate a simple OpenQASM circuit based on transcription rate"""
        # Mapping: rate -> rotation angle theta
        # Max rate ~ 2.0 -> pi
        import math
        theta = min(math.pi, rate * math.pi)
        
        qasm = f"""
        OPENQASM 3.0;
        include "stdgates.inc";
        bit[1] c;
        qubit[1] q;
        
        ry({theta}) q[0];
        c[0] = measure q[0];
        """
        # Return base64 encoded string if needed, but our mock accepts raw string or base64
        # For this simple integration, we send the string.
        # Ensure it's compatible with what the backend expects (often base64)
        import base64
        return base64.b64encode(qasm.encode('utf-8')).decode('utf-8')


class QuantumClusterWorker(QThread):
    """
    Worker for 28-qubit cluster-based quantum simulation.
    
    Each qubit represents a functional gene cluster:
    - RY rotation encodes cluster's transcription level
    - CNOT gates create entanglement between related clusters
    - Measurement gives probability distribution for each cluster
    """
    # Emits: dict of cluster_id -> probability, full counts dict
    cluster_finished = pyqtSignal(dict, dict)
    cluster_failed = pyqtSignal(str)
    
    API_URL = "http://localhost:8081/jobs"
    NUM_CLUSTERS = 28
    
    def __init__(self, cluster_levels: dict, entanglement_pairs: list = None):
        """
        Args:
            cluster_levels: Dict of cluster_id (0-27) -> normalized level [0, 1]
            entanglement_pairs: List of (cluster_a, cluster_b) tuples for CNOT gates
        """
        super().__init__()
        self.cluster_levels = cluster_levels
        self.entanglement_pairs = entanglement_pairs or []
        self.running = True
    
    def run(self):
        try:
            circuit = self._generate_cluster_circuit()
            
            job_payload = {
                "name": "cluster_28q_sim",
                "device_id": "SVSim",
                "job_type": "sampling",
                "shots": 4096,  # Shots at top level per API spec
                "job_info": {
                    "program": [circuit]  # Program as list per API spec
                }
            }
            
            headers = {"X-API-Token": "local-dev-token"}
            
            print(f"DEBUG: [Quantum] Submitting 28-qubit cluster circuit...")
            submit_resp = requests.post(self.API_URL, json=job_payload, headers=headers, timeout=5)
            
            # Debug: print full response
            print(f"DEBUG: [Quantum] Response status: {submit_resp.status_code}")
            print(f"DEBUG: [Quantum] Response body: {submit_resp.text[:500]}")
            
            if submit_resp.status_code not in [200, 201]:
                raise Exception(f"Submission failed: {submit_resp.text}")
            
            resp_data = submit_resp.json()
            # Handle different response formats
            if isinstance(resp_data, dict):
                job_id = resp_data.get('id') or resp_data.get('job_id') or resp_data.get('data', {}).get('id')
            else:
                job_id = resp_data
            
            if not job_id:
                raise Exception(f"No job ID in response: {resp_data}")
            
            print(f"DEBUG: [Quantum] Job submitted: {job_id}")
            
            # Poll for completion
            while self.running:
                status_resp = requests.get(f"{self.API_URL}/{job_id}", headers=headers, timeout=5)
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
                    cluster_probs = self._parse_cluster_probabilities(counts)
                    
                    print(f"DEBUG: [Quantum] 28-qubit results parsed successfully")
                    self.cluster_finished.emit(cluster_probs, counts)
                    return
                
                elif status in ["FAILED", "failed"]:
                    self.cluster_failed.emit("Cluster job failed in backend")
                    return
                
                time.sleep(0.5)
                
        except Exception as e:
            if self.running:
                self.cluster_failed.emit(str(e))
    
    def stop(self):
        self.running = False
    
    def _generate_cluster_circuit(self) -> str:
        """
        Generate 28-qubit OpenQASM circuit.
        
        Structure:
        1. Declare 28 qubits and 28 classical bits
        2. Apply RY(θ) to each qubit based on cluster level
        3. Apply CNOT gates for entanglement topology
        4. Measure all qubits
        """
        import math
        import base64
        
        lines = [
            "OPENQASM 3.0;",
            'include "stdgates.inc";',
            f"qubit[{self.NUM_CLUSTERS}] q;",
            f"bit[{self.NUM_CLUSTERS}] c;",
            "",
            "// Initialize each cluster qubit based on transcription level"
        ]
        
        # Layer 1: RY gates for each cluster
        for cluster_id in range(self.NUM_CLUSTERS):
            level = self.cluster_levels.get(cluster_id, 0.0)
            # Map level [0,1] to theta [0, π]
            # level=0 → |0⟩, level=1 → |1⟩, level=0.5 → superposition
            theta = level * math.pi
            lines.append(f"ry({theta:.6f}) q[{cluster_id}];")
        
        lines.append("")
        lines.append("// Entanglement layer: CNOT between related clusters")
        
        # Layer 2: CNOT gates for entanglement
        for cluster_a, cluster_b in self.entanglement_pairs:
            if 0 <= cluster_a < self.NUM_CLUSTERS and 0 <= cluster_b < self.NUM_CLUSTERS:
                lines.append(f"cx q[{cluster_a}], q[{cluster_b}];")
        
        lines.append("")
        lines.append("// Measure all clusters")
        
        # Layer 3: Measurement
        for i in range(self.NUM_CLUSTERS):
            lines.append(f"c[{i}] = measure q[{i}];")
        
        qasm = "\n".join(lines)
        return base64.b64encode(qasm.encode('utf-8')).decode('utf-8')
    
    def _parse_cluster_probabilities(self, counts: dict) -> dict:
        """
        Parse measurement counts to get per-cluster probabilities.
        
        For each qubit position, calculate P(|1⟩) across all measurement outcomes.
        
        Args:
            counts: Dict of bitstring -> count (e.g., {"0101...": 123, ...})
        
        Returns:
            Dict of cluster_id -> probability [0, 1]
        """
        total_shots = sum(counts.values())
        if total_shots == 0:
            return {i: 0.0 for i in range(self.NUM_CLUSTERS)}
        
        # Count how many times each qubit measured |1⟩
        one_counts = {i: 0 for i in range(self.NUM_CLUSTERS)}
        
        for bitstring, count in counts.items():
            # Bitstring may be in different formats
            # Standard: rightmost bit is q[0]
            bits = bitstring.replace(" ", "")
            
            for i, bit in enumerate(reversed(bits)):
                if i >= self.NUM_CLUSTERS:
                    break
                if bit == '1':
                    one_counts[i] += count
        
        # Calculate probabilities
        probabilities = {
            cluster_id: (one_counts[cluster_id] / total_shots) * 100.0
            for cluster_id in range(self.NUM_CLUSTERS)
        }
        
        return probabilities
