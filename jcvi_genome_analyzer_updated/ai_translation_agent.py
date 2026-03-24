"""
AI Translation Agent for JCVI Genome Analyzer
Displays translation data DIRECTLY from Translation Dynamics BED file
NO COMPUTATION - Direct data display from BED file
"""

import re
from typing import Dict, List, Tuple
from dataclasses import dataclass


@dataclass
class TranslationData:
    """Data from Translation Dynamics BED file"""
    gene_id: str
    protein_level: float  # 0-100%
    rbs_strength: str  # strong, medium
    start: int
    end: int
    strand: str
    score: float
    color: Tuple[int, int, int]


class AITranslationAgent:
    """
    Agent #3: Translation Dynamics Analysis
    
    Displays translation data DIRECTLY from Translation Dynamics BED file.
    This agent does NOT compute anything - it reads and displays pre-computed data.
    
    Input: JCVI_Translation_Dynamics_3_0_.bed
    Format: gene_id|Protein:XX%|RBS:strength
    """
    
    def __init__(self):
        self.translation_data: Dict[str, TranslationData] = {}
    
    def load_translation_bed(self, bed_path: str) -> int:
        """
        Load Translation Dynamics BED file directly
        
        Format: JCVISYN3_XXXX|Protein:XX%|RBS:strength
        Example: JCVISYN3_0294|Protein:95%|RBS:strong
        
        Returns: number of entries loaded
        """
        self.translation_data.clear()
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
                
                # Parse name: gene_id|Protein:XX%|RBS:strength
                name_parts = name.split('|')
                gene_id = name_parts[0]
                
                protein_level = 0.0
                rbs_strength = ""
                
                for part in name_parts[1:]:
                    # Protein:XX%
                    protein_match = re.search(r'Protein:(\d+)%', part)
                    if protein_match:
                        protein_level = float(protein_match.group(1))
                    
                    # RBS:strength
                    rbs_match = re.search(r'RBS:(\w+)', part)
                    if rbs_match:
                        rbs_strength = rbs_match.group(1)
                
                # Keep highest protein level for each gene
                if gene_id in self.translation_data:
                    if protein_level <= self.translation_data[gene_id].protein_level:
                        continue
                
                self.translation_data[gene_id] = TranslationData(
                    gene_id=gene_id,
                    protein_level=protein_level,
                    rbs_strength=rbs_strength,
                    start=start,
                    end=end,
                    strand=strand,
                    score=score,
                    color=color
                )
                count += 1
        
        print(f"DEBUG: Loaded {len(self.translation_data)} translation entries from BED")
        return len(self.translation_data)
    
    def get_translation_data(self) -> Dict[str, TranslationData]:
        """Get all translation data"""
        return self.translation_data
    
    def get_statistics(self) -> Dict:
        """Get translation statistics"""
        if not self.translation_data:
            return {}
        
        # Count by RBS strength
        rbs_counts = {'strong': 0, 'medium': 0, 'unknown': 0}
        
        # Count by protein level ranges
        high_level = 0  # >= 80%
        medium_level = 0  # 50-80%
        low_level = 0  # < 50%
        
        levels = []
        for gene_id, data in self.translation_data.items():
            levels.append(data.protein_level)
            
            # RBS
            rbs = data.rbs_strength.lower() if data.rbs_strength else 'unknown'
            if rbs in rbs_counts:
                rbs_counts[rbs] += 1
            else:
                rbs_counts['unknown'] += 1
            
            # Level
            if data.protein_level >= 80:
                high_level += 1
            elif data.protein_level >= 50:
                medium_level += 1
            else:
                low_level += 1
        
        return {
            'total_genes': len(self.translation_data),
            'rbs_distribution': rbs_counts,
            'high_level_count': high_level,
            'medium_level_count': medium_level,
            'low_level_count': low_level,
            'max_level': max(levels) if levels else 0,
            'min_level': min(levels) if levels else 0,
            'avg_level': sum(levels) / len(levels) if levels else 0
        }
    
    def get_sorted_by_level(self, descending: bool = True) -> List[TranslationData]:
        """Get all genes sorted by protein level"""
        sorted_data = sorted(
            self.translation_data.values(),
            key=lambda x: x.protein_level,
            reverse=descending
        )
        return sorted_data
    
    def get_animation_data(self) -> List[Dict]:
        """
        Prepare data for animation display
        Returns list of dicts with animation parameters
        """
        animation_data = []
        
        # Sort by protein level (highest first)
        sorted_data = self.get_sorted_by_level(descending=True)
        
        for order, data in enumerate(sorted_data):
            # Start time based on order
            start_time = order * 0.3  # 0.3 seconds between each gene
            
            # Rate based on RBS strength
            rate_map = {'strong': 0.3, 'medium': 0.2}
            rate = rate_map.get(data.rbs_strength.lower(), 0.15) if data.rbs_strength else 0.15
            
            animation_data.append({
                'gene_id': data.gene_id,
                'protein_level': data.protein_level,
                'rbs_strength': data.rbs_strength or 'unknown',
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
