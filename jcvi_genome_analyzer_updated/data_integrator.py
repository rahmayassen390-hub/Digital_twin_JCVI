"""
Data Integrator for JCVI Genome Analyzer
Integrates and interprets data from all BED files and GFF3

This module understands:
- Color meanings in BED files
- Relationships between files
- How to extract complete gene information
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class IntegratedGeneData:
    """Complete integrated data for a gene from all sources"""
    gene_id: str
    
    # From GFF3
    function: str = ""
    product: str = ""
    start: int = 0
    end: int = 0
    strand: str = "+"
    
    # From Regulation BED
    regulation_type: str = ""  # constitutive, operonMember, potentiallyInducible, etc.
    regulation_category: str = ""  # Housekeeping, Inducible, etc. (or function from GFF3)
    regulation_strength: str = ""  # high, medium, low
    
    # From Transcription BED
    transcription_level: float = 0.0
    transcription_phase: str = ""  # Phase1, Phase2, etc.
    transcription_category: str = ""  # DNA/Transcription, etc.
    transcription_order: int = 0
    transcription_time: str = ""  # 0-5min, etc.
    transcription_status: str = ""  # active, delayed, etc.
    
    # From Translation BED
    protein_level: float = 0.0
    translation_phase: str = ""  # Phase2, etc.
    translation_category: str = ""  # Ribosome Assembly, etc.
    translation_order: int = 0
    rbs_strength: str = ""  # strong, medium, weak
    protein_stability: str = ""  # very stable, stable, unstable
    is_ncrna: bool = False  # Non-coding RNA - no translation
    tx_vs_tl_efficiency: float = 0.0
    
    # Interpreted status
    expression_status: str = ""  # "Active", "Stopped at Transcription", "No Translation (ncRNA)"
    status_reason: str = ""  # Explanation of the status


class DataIntegrator:
    """
    Integrates data from all sources and interprets meanings
    
    Files used:
    - GFF3: Gene functions and products
    - Regulation BED: Gene regulation types
    - Transcription BED: Transcription dynamics (3 tracks)
    - Translation BED: Translation dynamics (6 tracks)
    """
    
    # Color interpretation for Transcription BED
    TRANSCRIPTION_COLORS = {
        (139, 0, 0): ("Phase1", "Early transcription - Essential genes"),
        (0, 0, 139): ("Phase2", "Secondary transcription"),
        (34, 139, 34): ("Phase3", "Metabolic genes"),
        (255, 165, 0): ("Phase4", "Growth-dependent"),
        (128, 0, 128): ("Phase5", "Late phase"),
        (0, 128, 128): ("Phase6", "Final phase"),
    }
    
    # Color interpretation for Translation BED
    TRANSLATION_COLORS = {
        (0, 100, 0): ("Strong RBS", "High protein production"),
        (0, 247, 0): ("Very high protein", "Efficient translation"),
        (255, 0, 255): ("Phase2 translation", "Ribosome assembly"),
        (255, 0, 0): ("Low efficiency", "Transcription active but translation reduced"),
        (255, 69, 0): ("ncRNA", "Non-coding - Transcribed but NOT translated"),
        (0, 0, 255): ("Stable protein", "Long half-life"),
        (255, 165, 0): ("Unstable protein", "Short half-life"),
    }
    
    def __init__(self):
        self.genes: Dict[str, IntegratedGeneData] = {}
        self.gene_functions: Dict[str, str] = {}  # From GFF3
    
    def load_gff3_functions(self, gff_path: str):
        """Load gene functions from GFF3 file"""
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature_type = parts[2]
                if feature_type != 'CDS':
                    continue
                
                attributes = parts[8]
                gene_id = ""
                product = ""
                
                for attr in attributes.split(';'):
                    if attr.startswith('locus_tag='):
                        gene_id = attr.split('=')[1]
                    elif attr.startswith('product='):
                        product = attr.split('=')[1]
                
                if gene_id and product:
                    self.gene_functions[gene_id] = product
        
        print(f"DEBUG: Loaded {len(self.gene_functions)} gene functions from GFF3")
    
    def integrate_regulation_data(self, regulation_data: dict):
        """Integrate regulation BED data and replace Unknown with GFF3 functions"""
        for gene_id, data in regulation_data.items():
            if gene_id not in self.genes:
                self.genes[gene_id] = IntegratedGeneData(gene_id=gene_id)
            
            gene = self.genes[gene_id]
            gene.regulation_type = data.regulation_type
            gene.regulation_strength = data.strength
            gene.start = data.start
            gene.end = data.end
            gene.strand = data.strand
            
            # Replace "Unknown" category with function from GFF3
            if data.category == "Unknown" and gene_id in self.gene_functions:
                # Extract short category from product
                product = self.gene_functions[gene_id]
                gene.regulation_category = self._extract_category(product)
                gene.function = product
            else:
                gene.regulation_category = data.category
    
    def integrate_transcription_data(self, transcription_data: dict):
        """Integrate transcription BED data with color interpretation"""
        for gene_id, data in transcription_data.items():
            if gene_id not in self.genes:
                self.genes[gene_id] = IntegratedGeneData(gene_id=gene_id)
            
            gene = self.genes[gene_id]
            gene.transcription_level = data.transcription_level
            gene.transcription_phase = data.phase
            gene.transcription_category = data.category
            gene.transcription_order = data.order
            gene.transcription_time = data.time_range
            
            # Interpret color
            if hasattr(data, 'color') and data.color:
                status = self._interpret_transcription_color(data.color)
                gene.transcription_status = status
    
    def integrate_translation_data(self, translation_data: dict):
        """Integrate translation BED data with color interpretation"""
        for gene_id, data in translation_data.items():
            if gene_id not in self.genes:
                self.genes[gene_id] = IntegratedGeneData(gene_id=gene_id)
            
            gene = self.genes[gene_id]
            gene.protein_level = data.protein_level
            gene.translation_phase = data.phase
            gene.translation_category = data.category
            gene.translation_order = data.order
            gene.rbs_strength = data.rbs_strength
            gene.protein_stability = data.stability
            gene.is_ncrna = data.is_ncrna
            gene.tx_vs_tl_efficiency = data.efficiency
            
            # Determine expression status
            self._determine_expression_status(gene, data)
    
    def _determine_expression_status(self, gene: IntegratedGeneData, data):
        """Determine the expression status based on all data"""
        
        # Check for ncRNA (transcribed but not translated)
        if data.is_ncrna:
            gene.expression_status = "🔴 No Translation (ncRNA)"
            gene.status_reason = "Non-coding RNA - Transcribed but NOT translated"
            return
        
        # Check for low translation efficiency (red color = stopped at transcription)
        if hasattr(data, 'color') and data.color == (255, 0, 0):
            gene.expression_status = "🟡 Low Translation"
            gene.status_reason = f"TX:{data.tx_level}% but TL:{data.tl_level}% - Translation bottleneck"
            return
        
        # Check for high efficiency
        if data.efficiency > 100:
            gene.expression_status = "🟢 High Expression"
            gene.status_reason = f"Efficient translation (Eff:{data.efficiency}%)"
            return
        
        # Normal active translation
        if data.protein_level > 80:
            gene.expression_status = "🟢 Active"
            gene.status_reason = f"High protein production ({data.protein_level}%)"
        elif data.protein_level > 50:
            gene.expression_status = "🟡 Moderate"
            gene.status_reason = f"Moderate protein production ({data.protein_level}%)"
        else:
            gene.expression_status = "🟠 Low"
            gene.status_reason = f"Low protein production ({data.protein_level}%)"
    
    def _interpret_transcription_color(self, color: tuple) -> str:
        """Interpret transcription BED color"""
        # Check exact match first
        if color in self.TRANSCRIPTION_COLORS:
            return self.TRANSCRIPTION_COLORS[color][1]
        
        # Check by darkness (grayscale indicates level)
        r, g, b = color
        if r == g == b:  # Grayscale
            if r < 50:
                return "Very high expression"
            elif r < 100:
                return "High expression"
            elif r < 150:
                return "Medium expression"
            else:
                return "Low expression"
        
        return "Active"
    
    def _extract_category(self, product: str) -> str:
        """Extract a short category from product description"""
        product_lower = product.lower()
        
        # Map common keywords to categories
        if any(kw in product_lower for kw in ['ribosom', 'rrna', 'trna']):
            return "Ribosome/Translation"
        elif any(kw in product_lower for kw in ['dna', 'replication', 'gyrase', 'helicase']):
            return "DNA/Replication"
        elif any(kw in product_lower for kw in ['transcription', 'rna polymerase', 'sigma']):
            return "Transcription"
        elif any(kw in product_lower for kw in ['membrane', 'transport', 'permease']):
            return "Membrane/Transport"
        elif any(kw in product_lower for kw in ['metabol', 'synthase', 'kinase', 'dehydrogenase']):
            return "Metabolism"
        elif any(kw in product_lower for kw in ['chaperone', 'protease', 'heat shock']):
            return "Protein Processing"
        elif any(kw in product_lower for kw in ['cell division', 'fts']):
            return "Cell Division"
        elif 'hypothetical' in product_lower:
            return "Hypothetical"
        else:
            # Return first part of product (usually the gene name)
            return product.split(':')[0] if ':' in product else product[:30]
    
    def get_integrated_data(self) -> Dict[str, IntegratedGeneData]:
        """Get all integrated gene data"""
        return self.genes
    
    def get_genes_by_status(self, status: str) -> List[IntegratedGeneData]:
        """Get genes filtered by expression status"""
        return [g for g in self.genes.values() if status in g.expression_status]
    
    def get_statistics(self) -> Dict:
        """Get statistics about integrated data"""
        stats = {
            'total_genes': len(self.genes),
            'with_function': sum(1 for g in self.genes.values() if g.function),
            'regulation_types': {},
            'expression_statuses': {},
            'ncrna_count': sum(1 for g in self.genes.values() if g.is_ncrna),
        }
        
        for gene in self.genes.values():
            # Count regulation types
            if gene.regulation_type:
                stats['regulation_types'][gene.regulation_type] = \
                    stats['regulation_types'].get(gene.regulation_type, 0) + 1
            
            # Count expression statuses
            if gene.expression_status:
                stats['expression_statuses'][gene.expression_status] = \
                    stats['expression_statuses'].get(gene.expression_status, 0) + 1
        
        return stats
    
    def generate_files_used_report(self) -> str:
        """Generate report of which files are used for each data type"""
        report = """
=== FILES USED FOR EACH DATA TYPE ===

📁 Promoter/Regulation Analysis:
   • Source: Regulation Analysis BED (JCVI_syn3_Gene_Regulation_3_0_.bed)
   • Fields: regulation_type, category, strength
   • Category "Unknown" → Replaced with function from GFF3

📁 Transcription Dynamics:
   • Source: Transcription Dynamics BED (JCVI_Transcription_Dynamics_3_0_.bed)
   • Track 1 - Transcription_Level: Level percentage
   • Track 2 - Transcription_Order: Phase, Category, Order
   • Track 3 - Transcription_Timeline: Time range (0-5min, etc.)
   • Colors: Grayscale = expression level (darker = higher)

📁 Translation Dynamics:
   • Source: Translation Dynamics BED (JCVI_Translation_Dynamics_3_0_.bed)
   • Track 1 - Translation_Level: Protein percentage
   • Track 2 - Translation_Order: Phase, Category, Order
   • Track 3 - RBS_Strength: RBS strength, Stability
   • Track 4 - Protein_Stability: Half-life
   • Track 5 - ncRNA_No_Translation: Non-coding RNAs (NO TRANSLATION)
   • Track 6 - TX_vs_TL: Efficiency comparison
   
   COLOR MEANINGS:
   • Green (0,100,0) / (0,247,0): High protein, efficient translation
   • Red (255,0,0): LOW EFFICIENCY - Transcription OK but Translation reduced
   • Orange (255,69,0): ncRNA - Transcribed but NOT Translated
   • Purple (255,0,255): Ribosome assembly phase
   • Blue: Stable protein

📁 Gene Functions:
   • Source: GFF3 file (JCVI_3_0.gff3)
   • Used to: Replace "Unknown" categories with actual functions
"""
        return report
