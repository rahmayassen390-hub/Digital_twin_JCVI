"""
AI Promoter Agent for JCVI Genome Analyzer
Displays regulation data DIRECTLY from Regulation Analysis BED file
NO COMPUTATION - Direct data display from BED file
"""

import re
from typing import Dict, List, Tuple
from dataclasses import dataclass


@dataclass
class RegulationData:
    """Data from Regulation BED file"""
    gene_id: str
    regulation_type: str  # constitutive, operonMember, potentiallyInducible, potentiallyRepressible, unknown
    category: str  # Housekeeping, Inducible, Repressible, Unknown
    strength: str  # high, medium, low
    start: int
    end: int
    strand: str
    score: float
    color: Tuple[int, int, int]


class PromoterAnalyzerAgent:
    """
    Agent #1: Promoter/Regulation Analysis
    
    Displays regulation data DIRECTLY from Regulation Analysis BED file.
    This agent does NOT compute anything - it reads and displays pre-computed data.
    
    Input: JCVI_syn3_Gene_Regulation_3_0_.bed
    Format: gene_id|regulation_type|category|strength
    """
    
    def __init__(self):
        self.regulation_data: Dict[str, RegulationData] = {}
    
    def load_regulation_bed(self, bed_path: str) -> int:
        """
        Load Regulation Analysis BED file directly
        
        Format: JCVISYN3_XXXX|regulation_type|category|strength
        Example: JCVISYN3_0002|constitutive|Housekeeping|high
        
        Returns: number of entries loaded
        """
        self.regulation_data.clear()
        count = 0
        
        with open(bed_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('track') or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 9:
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
                
                # Parse name: gene_id|regulation_type|category|strength
                name_parts = name.split('|')
                if len(name_parts) >= 4:
                    gene_id = name_parts[0]
                    regulation_type = name_parts[1]
                    category = name_parts[2]
                    strength = name_parts[3]
                    
                    self.regulation_data[gene_id] = RegulationData(
                        gene_id=gene_id,
                        regulation_type=regulation_type,
                        category=category,
                        strength=strength,
                        start=start,
                        end=end,
                        strand=strand,
                        score=score,
                        color=color
                    )
                    count += 1
        
        print(f"DEBUG: Loaded {count} regulation entries from BED")
        return count
    
    def get_regulation_data(self) -> Dict[str, RegulationData]:
        """Get all regulation data"""
        return self.regulation_data
    
    def get_statistics(self) -> Dict:
        """Get regulation statistics"""
        if not self.regulation_data:
            return {}
        
        # Count by regulation type
        type_counts = {}
        for gene_id, data in self.regulation_data.items():
            reg_type = data.regulation_type
            type_counts[reg_type] = type_counts.get(reg_type, 0) + 1
        
        # Count by category
        category_counts = {}
        for gene_id, data in self.regulation_data.items():
            cat = data.category
            category_counts[cat] = category_counts.get(cat, 0) + 1
        
        # Count by strength
        strength_counts = {}
        for gene_id, data in self.regulation_data.items():
            strength = data.strength
            strength_counts[strength] = strength_counts.get(strength, 0) + 1
        
        return {
            'total_genes': len(self.regulation_data),
            'regulation_types': type_counts,
            'categories': category_counts,
            'strengths': strength_counts
        }
    
    def get_genes_by_type(self, regulation_type: str) -> List[RegulationData]:
        """Get all genes with a specific regulation type"""
        return [
            data for data in self.regulation_data.values()
            if data.regulation_type == regulation_type
        ]
    
    def analyze_genes(self, genes, sequence="", likely_promoters=None, internal_operons=None):
        """
        For compatibility - returns regulation data for display
        This method is kept for backward compatibility but uses BED data
        """
        # Just return the loaded regulation data
        return self.regulation_data
