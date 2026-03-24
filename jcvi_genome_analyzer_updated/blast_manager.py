"""
BLAST Manager for JCVI Genome Analyzer
Handles all BLAST operations including database creation and search
"""

import os
import subprocess
import tempfile
from typing import List, Dict, Tuple

try:
    from Bio.Blast import NCBIXML
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


class BlastManager:
    """Manages BLAST operations"""
    
    def __init__(self):
        self.blast_path = self.find_blast()
        self.temp_dir = tempfile.mkdtemp()
        self.db_path = None
    
    def find_blast(self) -> str:
        """Find BLAST executable"""
        possible_paths = ['blastn', '/usr/bin/blastn', '/usr/local/bin/blastn']
        
        for path in possible_paths:
            try:
                result = subprocess.run([path, '-version'], 
                                      capture_output=True, text=True, timeout=2)
                if result.returncode == 0:
                    return path
            except:
                continue
        return None
    
    def check_blast(self) -> Tuple[bool, str]:
        """Check if BLAST is available"""
        if self.blast_path:
            try:
                result = subprocess.run([self.blast_path, '-version'],
                                      capture_output=True, text=True)
                version = result.stdout.strip().split('\n')[0]
                return True, f"BLAST found: {version}"
            except:
                return False, "BLAST found but not working"
        return False, "BLAST NOT found!"
    
    def make_blast_db(self, fasta_file: str) -> str:
        """Create BLAST database"""
        if not os.path.exists(fasta_file):
            raise Exception(f"FASTA file not found: {fasta_file}")
            
        if os.path.isdir(fasta_file):
            raise Exception(f"FASTA path is a directory, not a file: {fasta_file}\n"
                            f"Please select the actual .fasta file instead of the folder.")
        
        db_path = os.path.join(self.temp_dir, "blast_db")
        makeblastdb_path = self.blast_path.replace('blastn', 'makeblastdb')
        
        # Wrap paths in quotes to help BLAST handle spaces
        fasta_abs = os.path.abspath(fasta_file)
        
        cmd = [
            makeblastdb_path,
            '-in', fasta_abs,
            '-dbtype', 'nucl',
            '-out', db_path
        ]
        
        # If there are spaces, some BLAST versions require the path to be quoted EVEN when passed as a list
        if " " in fasta_abs:
            cmd[2] = f'"{fasta_abs}"'
        if " " in db_path:
            cmd[6] = f'"{db_path}"'
            
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0:
            self.db_path = db_path
            return db_path
        else:
            stderr = result.stderr
            # Catch common error where path with spaces is not quoted internally by BLAST
            if "does not exist" in stderr and " " in fasta_file:
                stderr += f"\n\nTIP: The path contains spaces. Ensure the FASTA file exists at: {fasta_file}"
                
            raise Exception(f"makeblastdb failed:\n{stderr}")
    
    def run_blast(self, query_file: str) -> str:
        """Run BLASTN"""
        output_file = os.path.join(self.temp_dir, "blast_results.xml")
        
        query_abs = os.path.abspath(query_file)
        db_abs = os.path.abspath(self.db_path)
        output_abs = os.path.abspath(output_file)
        
        cmd = [
            self.blast_path,
            '-query', query_abs,
            '-db', db_abs,
            '-out', output_abs,
            '-outfmt', '5',
            '-num_alignments', '1000',
            '-evalue', '10',
            '-word_size', '11',
            '-num_threads', '4',
        ]
        
        # Wrap in quotes if spaces are present
        if " " in query_abs: cmd[2] = f'"{query_abs}"'
        if " " in db_abs: cmd[4] = f'"{db_abs}"'
        if " " in output_abs: cmd[6] = f'"{output_abs}"'
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        
        if result.returncode == 0:
            return output_file
        else:
            raise Exception(f"BLAST error: {result.stderr}")
    
    def parse_blast_xml(self, xml_file: str) -> List[Dict]:
        """Parse BLAST XML output"""
        hits = []
        
        with open(xml_file, 'r') as f:
            blast_records = NCBIXML.parse(f)
            
            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        identity = (hsp.identities / hsp.align_length * 100)
                        
                        hit = {
                            'query_id': record.query,
                            'target_id': alignment.title,
                            'identity': identity,
                            'alignment_length': hsp.align_length,
                            'mismatches': hsp.align_length - hsp.identities,
                            'gap_opens': hsp.gaps,
                            'query_start': hsp.query_start,
                            'query_end': hsp.query_end,
                            'target_start': hsp.sbjct_start,
                            'target_end': hsp.sbjct_end,
                            'evalue': hsp.expect,
                            'bit_score': hsp.score,
                            'query_seq': hsp.query,
                            'target_seq': hsp.sbjct,
                        }
                        hits.append(hit)
        
        return hits
