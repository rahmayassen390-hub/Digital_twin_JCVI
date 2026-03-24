"""
AI Transcription Agent for JCVI Genome Analyzer
Displays transcription data DIRECTLY from Transcription Dynamics BED file
NO COMPUTATION - Direct data display from BED file
"""

import re
from typing import Dict, List, Tuple
from dataclasses import dataclass


@dataclass
class TranscriptionData:
    """Data from Transcription Dynamics BED file"""
    gene_id: str
    level: float  # 0-100%
    start: int
    end: int
    strand: str
    score: float
    color: Tuple[int, int, int]


class TranscriptionDynamicsAgent:
    """
    Agent #2: Transcription Dynamics Analysis
    
    Displays transcription data DIRECTLY from Transcription Dynamics BED file.
    This agent does NOT compute anything - it reads and displays pre-computed data.
    
    Input: JCVI_Transcription_Dynamics_3_0_.bed
    Format: gene_id|Level:XX%
    """
    
    def __init__(self):
        self.transcription_data: Dict[str, TranscriptionData] = {}
    
    def load_transcription_bed(self, bed_path: str) -> int:
        """
        Load Transcription Dynamics BED file directly
        
        Format: JCVISYN3_XXXX|Level:XX%
        Example: JCVISYN3_0253|Level:90%
        
        Returns: number of entries loaded
        """
        self.transcription_data.clear()
        count = 0
        
        with open(bed_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('track') or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 6:
                    continue
                
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3]
                score = float(parts[4])
                strand = parts[5]
                
                # Parse RGB color
                color = (128, 128, 128)
                if len(parts) >= 9:
                    rgb_str = parts[8]
                    rgb_parts = rgb_str.split(',')
                    if len(rgb_parts) == 3:
                        try:
                            color = (int(rgb_parts[0]), int(rgb_parts[1]), int(rgb_parts[2]))
                        except:
                            pass
                
                # Parse name: gene_id|Level:XX%
                name_parts = name.split('|')
                gene_id = name_parts[0]
                
                level = 0.0
                if len(name_parts) >= 2:
                    level_match = re.search(r'Level:(\d+)%', name_parts[1])
                    if level_match:
                        level = float(level_match.group(1))
                
                # Keep highest level for each gene
                if gene_id in self.transcription_data:
                    if level <= self.transcription_data[gene_id].level:
                        continue
                
                self.transcription_data[gene_id] = TranscriptionData(
                    gene_id=gene_id,
                    level=level,
                    start=start,
                    end=end,
                    strand=strand,
                    score=score,
                    color=color
                )
                count += 1
        
        print(f"DEBUG: Loaded {len(self.transcription_data)} transcription entries from BED")
        return len(self.transcription_data)
    
    def get_transcription_data(self) -> Dict[str, TranscriptionData]:
        """Get all transcription data"""
        return self.transcription_data
    
    def get_statistics(self) -> Dict:
        """Get transcription statistics"""
        if not self.transcription_data:
            return {}
        
        # Count by level ranges
        high_level = 0  # >= 70%
        medium_level = 0  # 40-70%
        low_level = 0  # < 40%
        
        levels = []
        for gene_id, data in self.transcription_data.items():
            levels.append(data.level)
            if data.level >= 70:
                high_level += 1
            elif data.level >= 40:
                medium_level += 1
            else:
                low_level += 1
        
        return {
            'total_genes': len(self.transcription_data),
            'high_level_count': high_level,
            'medium_level_count': medium_level,
            'low_level_count': low_level,
            'max_level': max(levels) if levels else 0,
            'min_level': min(levels) if levels else 0,
            'avg_level': sum(levels) / len(levels) if levels else 0
        }
    
    def get_sorted_by_level(self, descending: bool = True) -> List[TranscriptionData]:
        """Get all genes sorted by transcription level"""
        sorted_data = sorted(
            self.transcription_data.values(),
            key=lambda x: x.level,
            reverse=descending
        )
        return sorted_data
    
    def get_animation_data(self) -> List[Dict]:
        """
        Prepare data for animation display
        Returns list of dicts with animation parameters
        """
        animation_data = []
        
        # Sort by level (highest first = earliest transcription)
        sorted_data = self.get_sorted_by_level(descending=True)
        
        for order, data in enumerate(sorted_data):
            # Start time based on order
            start_time = order * 0.3  # 0.3 seconds between each gene
            
            # Rate based on level
            rate = 0.1 + (data.level / 100.0) * 0.2
            
            animation_data.append({
                'gene_id': data.gene_id,
                'level': data.level,
                'strand': data.strand,
                'start': data.start,
                'end': data.end,
                'temporal_order': order + 1,
                'start_time': start_time,
                'rate': rate,
                'current_level': 0.0,
                'status': 'waiting'
            })
        
        return animation_data
    
    def predict_transcription(self, genes, promoter_results=None, sequence=""):
        """
        For compatibility - returns transcription data for display
        This method is kept for backward compatibility but uses BED data
        """
        return self.transcription_data
