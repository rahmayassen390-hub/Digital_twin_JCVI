"""
Smart Data Engine for JCVI Genome Analyzer
===========================================

A THINKING data extraction system that:
1. Knows ALL available files and what each contains
2. Extracts data intelligently - if not found in one file, searches others
3. Links data across files to build complete gene profiles
4. Never leaves "Unknown" if information exists somewhere

Author: JCVI Genome Analysis System
Version: 2.0 - Smart Extraction

FILE MAP:
=========
┌─────────────────────────────────────────────────────────────────────────────┐
│ FILE TYPE          │ WHAT IT CONTAINS                  │ KEY FIELDS         │
├─────────────────────────────────────────────────────────────────────────────┤
│ GFF3               │ Gene definitions & functions      │ product, locus_tag │
│ Regulation BED     │ Gene regulation types             │ type, category     │
│ Transcription BED  │ Transcription dynamics (3 tracks) │ level, phase, time │
│ Translation BED    │ Translation dynamics (6 tracks)   │ protein%, RBS, etc │
│ Promoters BED      │ Promoter locations                │ -35, -10 boxes     │
│ Operons BED        │ Operon structure                  │ operon members     │
│ Complex Cases BED  │ Special regulation cases          │ case details       │
└─────────────────────────────────────────────────────────────────────────────┘

THINKING PROCESS:
=================
1. Load ALL files first → Know what's available
2. For each gene, collect ALL information from ALL files
3. If category is "Unknown" → Search in:
   - GFF3 (product field) → Extract function category
   - Transcription BED (category field)
   - Translation BED (category field)
   - Operons BED (maybe part of operon)
4. Build complete IntegratedGene with NO unknowns if possible
"""

import os
import re
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
from collections import defaultdict


@dataclass
class IntegratedGene:
    """Complete gene profile from ALL sources (Symmetric with Gene class)"""
    id: str
    name: str = ""              # Short name (e.g. dnaA)
    product: str = ""           # Full product description
    function: str = "Unknown"   # Broad function
    function_category: str = "" # Detailed category (e.g., "Ribosome/Translation")
    type: str = "gene"          # gene, CDS, rRNA, tRNA
    start: int = 0
    end: int = 0
    strand: str = "+"
    sequence: str = ""
    
    # === REGULATION & OPERONS ===
    has_promoter: bool = False
    is_operon_internal: bool = False
    operon_id: Optional[str] = None
    regulation_type: str = ""   # constitutive, operonMember, etc.
    regulation_category: str = ""
    regulation_strength: str = ""
    
    # === TRANSCRIPTION & TRANSLATION DYNAMICS ===
    transcription_level: float = 0.0
    transcription_phase: str = ""
    transcription_order: int = 0
    protein_level: float = 0.0
    translation_phase: str = ""
    translation_order: int = 0
    rbs_category: str = ""          # categorical: strong, medium, weak
    rbs_strength: float = 1.0       # numeric index (for simulation)
    protein_stability: str = ""
    is_ncrna: bool = False
    
    # === KINETICS & ENRICHMENT (Sync with Gene class) ===
    mrna_expression: float = 0.0
    protein_abundance: float = 0.0
    mrna_half_life: float = 5.0     # minutes
    protein_half_life: float = 25.0 # hours
    essentiality_status: str = "Unknown"
    function_class: str = "Unknown"
    ortholog_id: str = ""
    kegg_id: str = ""
    kinetics: Dict = field(default_factory=dict)
    
    # === UNIVERSAL REPRESENTATIONS ===
    molecular_ur: Dict = field(default_factory=dict)
    cellular_ur: Optional[np.ndarray] = None
    
    # === ADDITIONAL METADATA ===
    minus35_box: str = ""
    minus10_box: str = ""
    gc_content: float = 0.0
    experimental_data: Dict = field(default_factory=dict)
    notes: str = ""
    
    # === OPERON (from Operons BED) ===
    operon_id: str = ""
    operon_position: int = 0        # Position in operon (1, 2, 3, etc.)
    operon_members: List[str] = field(default_factory=list)
    
    # === METADATA ===
    data_sources: List[str] = field(default_factory=list)  # Which files provided data
    confidence: float = 1.0         # Data confidence score
    notes: str = ""                 # Any special notes


class SmartDataEngine:
    """
    SMART Data Engine - Thinks and searches for complete information
    
    PHILOSOPHY:
    - Never accept "Unknown" if information exists somewhere
    - Always try multiple sources
    - Build the most complete picture possible
    """
    
    # === CATEGORY MAPPING ===
    # Maps keywords in product descriptions to categories
    CATEGORY_KEYWORDS = {
        "Ribosome/Translation": [
            'ribosom', 'rrna', 'trna', '30s', '50s', 'translation', 
            'elongation factor', 'initiation factor', 'release factor',
            'aminoacyl-trna', 'peptidyl', 'ribosomal protein'
        ],
        "DNA/Replication": [
            'dna polymerase', 'replication', 'gyrase', 'helicase', 
            'topoisomerase', 'primase', 'ligase', 'dnaa', 'dnab',
            'dna repair', 'recombination'
        ],
        "Transcription": [
            'rna polymerase', 'sigma', 'transcription', 'rnase',
            'anti-sigma', 'nusa', 'nus', 'rho'
        ],
        "Membrane/Transport": [
            'membrane', 'transport', 'permease', 'atp-binding cassette',
            'abc transporter', 'channel', 'porin', 'secretion',
            'efflux', 'import', 'export'
        ],
        "Metabolism": [
            'synthase', 'kinase', 'dehydrogenase', 'reductase', 'oxidase',
            'transferase', 'isomerase', 'hydrolase', 'lyase', 'ligase',
            'metabol', 'glycolysis', 'tca', 'pentose'
        ],
        "Protein Processing": [
            'chaperone', 'protease', 'heat shock', 'clp', 'groel', 'groes',
            'dnak', 'dnaj', 'hsp', 'foldase', 'peptidase'
        ],
        "Cell Division": [
            'cell division', 'fts', 'division', 'septum', 'min',
            'chromosome partition', 'segregation'
        ],
        "Lipid Metabolism": [
            'lipid', 'fatty acid', 'phospholipid', 'acyl', 'lipase',
            'phosphatase', 'cardiolipin'
        ],
        "Cell Wall": [
            'peptidoglycan', 'murein', 'mur', 'lipopolysaccharide',
            'lps', 'cell wall'
        ],
        "Stress Response": [
            'stress', 'oxidative', 'heat shock', 'cold shock',
            'osmotic', 'acid', 'alkali'
        ],
        "Signaling": [
            'kinase', 'phosphatase', 'two-component', 'sensor',
            'regulator', 'signal'
        ],
        "Energy": [
            'atp synthase', 'atpase', 'electron transport', 'cytochrome',
            'nadh', 'fadh', 'quinone'
        ],
        "Nucleotide Metabolism": [
            'nucleotide', 'purine', 'pyrimidine', 'nucleoside',
            'thymidylate', 'guanylate', 'adenylate'
        ],
        "Amino Acid Metabolism": [
            'amino acid', 'aminotransferase', 'biosynthesis',
            'glutamate', 'aspartate', 'serine', 'glycine'
        ],
        "Cofactor Biosynthesis": [
            'cofactor', 'vitamin', 'coenzyme', 'nad', 'fad',
            'folate', 'biotin', 'thiamine'
        ],
    }
    
    def __init__(self):
        """Initialize the Smart Data Engine"""
        self.genes: Dict[str, IntegratedGene] = {}
        
        # === RAW DATA STORAGE ===
        self.gff3_data: Dict[str, dict] = {}          # gene_id -> {product, start, end, strand}
        self.regulation_data: Dict[str, dict] = {}     # gene_id -> regulation info
        self.transcription_data: Dict[str, dict] = {}  # gene_id -> transcription info
        self.translation_data: Dict[str, dict] = {}    # gene_id -> translation info
        self.promoter_data: Dict[str, dict] = {}       # gene_id -> promoter info
        self.operon_data: Dict[str, dict] = {}         # gene_id -> operon info
        self.complex_cases: Dict[str, dict] = {}       # gene_id -> complex case info
        
        # === FILE TRACKING ===
        self.loaded_files: List[str] = []
        self.file_gene_counts: Dict[str, int] = {}
        
        # === THINKING LOG ===
        self.thinking_log: List[str] = []
    
    def log_thinking(self, message: str):
        """Log the engine's thinking process"""
        self.thinking_log.append(message)
        print(f"🧠 THINKING: {message}")
    
    # =========================================================================
    # FILE LOADING METHODS
    # =========================================================================
    
    def load_gff3(self, gff_path: str) -> int:
        """
        Load GFF3 file - PRIMARY source for gene functions
        
        Extracts:
        - Gene ID (locus_tag)
        - Product description (function)
        - Start, End, Strand
        """
        if not gff_path or not os.path.exists(gff_path):
            self.log_thinking(f"GFF3 file not found: {gff_path}")
            return 0
        
        self.log_thinking(f"Loading GFF3: {gff_path}")
        count = 0
        
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature_type = parts[2]
                if feature_type not in ['CDS', 'gene']:
                    continue
                
                # Parse attributes
                attributes = {}
                for attr in parts[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key] = value
                
                gene_id = attributes.get('locus_tag', '')
                if not gene_id:
                    continue
                
                # Store data
                if gene_id not in self.gff3_data:
                    self.gff3_data[gene_id] = {}
                
                self.gff3_data[gene_id].update({
                    'product': attributes.get('product', ''),
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6],
                    'feature_type': feature_type
                })
                count += 1
        
        self.loaded_files.append(gff_path)
        self.file_gene_counts['GFF3'] = len(self.gff3_data)
        self.log_thinking(f"Loaded {len(self.gff3_data)} genes from GFF3")
        
        return len(self.gff3_data)
    
    def load_regulation_bed(self, bed_path: str) -> int:
        """
        Load Regulation BED file
        
        Format: gene_id|regulation_type|category|strength
        Example: JCVISYN3_0002|constitutive|Housekeeping|high
        """
        if not bed_path or not os.path.exists(bed_path):
            self.log_thinking(f"Regulation BED not found: {bed_path}")
            return 0
        
        self.log_thinking(f"Loading Regulation BED: {bed_path}")
        count = 0
        
        with open(bed_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('track') or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                
                # Parse name field: gene_id|type|category|strength
                name_parts = parts[3].split('|')
                if len(name_parts) < 4:
                    continue
                
                gene_id = name_parts[0]
                
                self.regulation_data[gene_id] = {
                    'regulation_type': name_parts[1],
                    'category': name_parts[2],
                    'strength': name_parts[3],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'strand': parts[5] if len(parts) > 5 else '+',
                    'score': float(parts[4]) if len(parts) > 4 else 0,
                    'color': self._parse_color(parts[8] if len(parts) > 8 else '')
                }
                count += 1
        
        self.loaded_files.append(bed_path)
        self.file_gene_counts['Regulation BED'] = count
        self.log_thinking(f"Loaded {count} regulation entries")
        
        return count
    
    def load_transcription_bed(self, bed_path: str) -> int:
        """
        Load Transcription Dynamics BED file (3 tracks)
        
        Track 1: Transcription_Level - gene_id|Level:XX%
        Track 2: Transcription_Order - order|gene_id|Phase:category|level%
        Track 3: Transcription_Timeline - gene_id|time_range|level%
        """
        if not bed_path or not os.path.exists(bed_path):
            self.log_thinking(f"Transcription BED not found: {bed_path}")
            return 0
        
        self.log_thinking(f"Loading Transcription BED: {bed_path}")
        current_track = None
        
        with open(bed_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Detect track
                if line.startswith('track'):
                    if 'Transcription_Level' in line:
                        current_track = 'level'
                    elif 'Transcription_Order' in line:
                        current_track = 'order'
                    elif 'Transcription_Timeline' in line:
                        current_track = 'timeline'
                    continue
                
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                
                name = parts[3]
                color = self._parse_color(parts[8] if len(parts) > 8 else '')
                
                if current_track == 'level':
                    # Format: gene_id|Level:XX%
                    match = re.match(r'([^|]+)\|Level:(\d+)%', name)
                    if match:
                        gene_id = match.group(1)
                        level = float(match.group(2))
                        
                        if gene_id not in self.transcription_data:
                            self.transcription_data[gene_id] = {}
                        
                        self.transcription_data[gene_id].update({
                            'level': level,
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'strand': parts[5] if len(parts) > 5 else '+',
                            'color': color
                        })
                
                elif current_track == 'order':
                    # Format: order|gene_id|Phase:category|level%
                    match = re.match(r'(\d+)\|([^|]+)\|([^:]+):([^|]+)\|(\d+)%', name)
                    if match:
                        order = int(match.group(1))
                        gene_id = match.group(2)
                        phase = match.group(3)
                        category = match.group(4)
                        level = float(match.group(5))
                        
                        if gene_id not in self.transcription_data:
                            self.transcription_data[gene_id] = {}
                        
                        self.transcription_data[gene_id].update({
                            'order': order,
                            'phase': phase,
                            'category': category,
                            'color': color
                        })
                
                elif current_track == 'timeline':
                    # Format: gene_id|time_range|level%
                    match = re.match(r'([^|]+)\|([^|]+)\|(\d+)%', name)
                    if match:
                        gene_id = match.group(1)
                        time_range = match.group(2)
                        
                        if gene_id not in self.transcription_data:
                            self.transcription_data[gene_id] = {}
                        
                        self.transcription_data[gene_id]['time_range'] = time_range
        
        self.loaded_files.append(bed_path)
        self.file_gene_counts['Transcription BED'] = len(self.transcription_data)
        self.log_thinking(f"Loaded {len(self.transcription_data)} transcription entries")
        
        return len(self.transcription_data)
    
    def load_translation_bed(self, bed_path: str) -> int:
        """
        Load Translation Dynamics BED file (6 tracks)
        
        Track 1: Translation_Level - gene_id|Protein:XX%|RBS:strength
        Track 2: Translation_Order - order|gene_id|Phase:category|level%
        Track 3: RBS_Strength - gene_id|RBS:strength|Stability:stability
        Track 4: Protein_Stability - gene_id|Stability:stability|HalfLife:time
        Track 5: ncRNA_No_Translation - gene_id|ncRNA|No_Translation
        Track 6: TX_vs_TL - gene_id|TX:XX%|TL:XX%|Eff:XX%
        """
        if not bed_path or not os.path.exists(bed_path):
            self.log_thinking(f"Translation BED not found: {bed_path}")
            return 0
        
        self.log_thinking(f"Loading Translation BED: {bed_path}")
        current_track = None
        
        with open(bed_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Detect track
                if line.startswith('track'):
                    if 'Translation_Level' in line:
                        current_track = 'level'
                    elif 'Translation_Order' in line:
                        current_track = 'order'
                    elif 'RBS_Strength' in line:
                        current_track = 'rbs'
                    elif 'Protein_Stability' in line:
                        current_track = 'stability'
                    elif 'ncRNA_No_Translation' in line:
                        current_track = 'ncrna'
                    elif 'Transcription_vs_Translation' in line or 'TX_vs_TL' in line:
                        current_track = 'efficiency'
                    continue
                
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                
                name = parts[3]
                color = self._parse_color(parts[8] if len(parts) > 8 else '')
                
                if current_track == 'level':
                    # Format: gene_id|Protein:XX%|RBS:strength
                    match = re.match(r'([^|]+)\|Protein:(\d+)%\|RBS:(\w+)', name)
                    if match:
                        gene_id = match.group(1)
                        protein = float(match.group(2))
                        rbs = match.group(3)
                        
                        if gene_id not in self.translation_data:
                            self.translation_data[gene_id] = {}
                        
                        self.translation_data[gene_id].update({
                            'protein_level': protein,
                            'rbs_strength': rbs,
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'strand': parts[5] if len(parts) > 5 else '+',
                            'color': color
                        })
                
                elif current_track == 'order':
                    # Format: order|gene_id|Phase:category|level%
                    match = re.match(r'(\d+)\|([^|]+)\|([^:]+):([^|]+)\|(\d+)%', name)
                    if match:
                        order = int(match.group(1))
                        gene_id = match.group(2)
                        phase = match.group(3)
                        category = match.group(4)
                        
                        if gene_id not in self.translation_data:
                            self.translation_data[gene_id] = {}
                        
                        self.translation_data[gene_id].update({
                            'order': order,
                            'phase': phase,
                            'category': category
                        })
                
                elif current_track == 'rbs':
                    # Format: gene_id|RBS:strength|Stability:stability
                    match = re.match(r'([^|]+)\|RBS:(\w+)\|Stability:(.+)', name)
                    if match:
                        gene_id = match.group(1)
                        rbs = match.group(2)
                        stability = match.group(3)
                        
                        if gene_id not in self.translation_data:
                            self.translation_data[gene_id] = {}
                        
                        self.translation_data[gene_id].update({
                            'rbs_strength': rbs,
                            'stability': stability
                        })
                
                elif current_track == 'stability':
                    # Format: gene_id|Stability:stability|HalfLife:time
                    match = re.match(r'([^|]+)\|Stability:([^|]+)', name)
                    if match:
                        gene_id = match.group(1)
                        stability = match.group(2)
                        
                        if gene_id not in self.translation_data:
                            self.translation_data[gene_id] = {}
                        
                        self.translation_data[gene_id]['stability'] = stability
                
                elif current_track == 'ncrna':
                    # Format: gene_id|ncRNA|...
                    if '|ncRNA' in name or 'ncRNA' in name:
                        gene_id = name.split('|')[0]
                        
                        if gene_id not in self.translation_data:
                            self.translation_data[gene_id] = {}
                        
                        self.translation_data[gene_id]['is_ncrna'] = True
                
                elif current_track == 'efficiency':
                    # Format: gene_id|TX:XX%|TL:XX%|Eff:XX%
                    match = re.match(r'([^|]+)\|TX:(\d+)%\|TL:(\d+)%\|Eff:(\d+)%', name)
                    if match:
                        gene_id = match.group(1)
                        tx = float(match.group(2))
                        tl = float(match.group(3))
                        eff = float(match.group(4))
                        
                        if gene_id not in self.translation_data:
                            self.translation_data[gene_id] = {}
                        
                        self.translation_data[gene_id].update({
                            'tx_level': tx,
                            'tl_level': tl,
                            'efficiency': eff,
                            'color': color
                        })
        
        self.loaded_files.append(bed_path)
        self.file_gene_counts['Translation BED'] = len(self.translation_data)
        self.log_thinking(f"Loaded {len(self.translation_data)} translation entries")
        
        return len(self.translation_data)
    
    def load_promoters_bed(self, bed_path: str) -> int:
        """Load Promoters BED file"""
        if not bed_path or not os.path.exists(bed_path):
            return 0
        
        self.log_thinking(f"Loading Promoters BED: {bed_path}")
        count = 0
        
        with open(bed_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('track') or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                
                name = parts[3]
                # Extract gene_id from various formats
                gene_id = name.split('|')[0] if '|' in name else name.split('_promoter')[0]
                
                self.promoter_data[gene_id] = {
                    'has_promoter': True,
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'name': name
                }
                count += 1
        
        self.loaded_files.append(bed_path)
        self.file_gene_counts['Promoters BED'] = count
        
        return count
    
    def load_operons_bed(self, bed_path: str) -> int:
        """Load Operons BED file"""
        if not bed_path or not os.path.exists(bed_path):
            return 0
        
        self.log_thinking(f"Loading Operons BED: {bed_path}")
        count = 0
        
        with open(bed_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('track') or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                
                name = parts[3]
                # Parse operon format: gene1,gene2,gene3 or similar
                gene_ids = re.findall(r'JCVISYN3_\d+', name)
                
                for i, gene_id in enumerate(gene_ids):
                    self.operon_data[gene_id] = {
                        'operon_id': name,
                        'position': i + 1,
                        'members': gene_ids
                    }
                    count += 1
        
        self.loaded_files.append(bed_path)
        self.file_gene_counts['Operons BED'] = count
        
        return count
    
    # =========================================================================
    # SMART CATEGORY DERIVATION
    # =========================================================================
    
    def _derive_category(self, product: str) -> str:
        """
        SMART: Derive category from product description
        Uses keyword matching to find the best category
        """
        if not product:
            return ""
        
        product_lower = product.lower()
        
        # Check each category's keywords
        for category, keywords in self.CATEGORY_KEYWORDS.items():
            for keyword in keywords:
                if keyword in product_lower:
                    return category
        
        # If no keyword match, try to extract from product format
        # Many products are formatted as "GeneName: description"
        if ':' in product:
            gene_name = product.split(':')[0].strip()
            if len(gene_name) <= 10:  # Short gene names
                return gene_name
        
        # Check for hypothetical
        if 'hypothetical' in product_lower or 'unknown function' in product_lower:
            return "Hypothetical"
        
        # Return shortened product as last resort
        return product[:30] + "..." if len(product) > 30 else product
    
    def _parse_color(self, color_str: str) -> Tuple[int, int, int]:
        """Parse RGB color from string"""
        if not color_str:
            return (128, 128, 128)
        
        try:
            parts = color_str.split(',')
            if len(parts) >= 3:
                return (int(parts[0]), int(parts[1]), int(parts[2]))
        except:
            pass
        
        return (128, 128, 128)
    
    # =========================================================================
    # INTEGRATION - THE SMART PART
    # =========================================================================
    
    def integrate_all_data(self) -> Dict[str, IntegratedGene]:
        """
        MAIN INTEGRATION METHOD
        
        Combines all loaded data into complete IntegratedGene objects.
        Uses THINKING to fill gaps and eliminate "Unknown" values.
        """
        self.log_thinking("=== STARTING SMART INTEGRATION ===")
        
        # Collect ALL gene IDs from ALL sources
        all_gene_ids: Set[str] = set()
        all_gene_ids.update(self.gff3_data.keys())
        all_gene_ids.update(self.regulation_data.keys())
        all_gene_ids.update(self.transcription_data.keys())
        all_gene_ids.update(self.translation_data.keys())
        all_gene_ids.update(self.promoter_data.keys())
        all_gene_ids.update(self.operon_data.keys())
        
        self.log_thinking(f"Found {len(all_gene_ids)} unique genes across all files")
        
        # Process each gene
        for gene_id in sorted(all_gene_ids):
            gene = self._integrate_gene(gene_id)
            self.genes[gene_id] = gene
        
        # Statistics
        unknown_count = sum(1 for g in self.genes.values() 
                          if g.function_category in ['', 'Unknown', 'Hypothetical'])
        
        self.log_thinking(f"Integration complete: {len(self.genes)} genes")
        self.log_thinking(f"Genes with known categories: {len(self.genes) - unknown_count}")
        self.log_thinking(f"Genes still unknown: {unknown_count}")
        
        return self.genes
    
    def _integrate_gene(self, gene_id: str) -> IntegratedGene:
        """
        Integrate all data for a single gene
        THINKS about where to find missing information
        """
        gene = IntegratedGene(id=gene_id)
        
        # === STEP 1: Load GFF3 data (primary source for function) ===
        if gene_id in self.gff3_data:
            gff = self.gff3_data[gene_id]
            gene.product = gff.get('product', '')
            gene.start = gff.get('start', 0)
            gene.end = gff.get('end', 0)
            gene.strand = gff.get('strand', '+')
            gene.data_sources.append('GFF3')
            
            # Derive category from product
            gene.function_category = self._derive_category(gene.product)
        
        # === STEP 2: Load Regulation data ===
        if gene_id in self.regulation_data:
            reg = self.regulation_data[gene_id]
            gene.regulation_type = reg.get('regulation_type', '')
            gene.regulation_category = reg.get('category', '')
            gene.regulation_strength = reg.get('strength', '')
            gene.data_sources.append('Regulation BED')
            
            # THINK: If regulation category is "Unknown", use function_category
            if gene.regulation_category in ['', 'Unknown']:
                if gene.function_category and gene.function_category != 'Hypothetical':
                    gene.regulation_category = gene.function_category
        
        # === STEP 3: Load Transcription data ===
        if gene_id in self.transcription_data:
            trans = self.transcription_data[gene_id]
            gene.transcription_level = trans.get('level', 0)
            gene.transcription_phase = trans.get('phase', '')
            gene.transcription_category = trans.get('category', '')
            gene.transcription_order = trans.get('order', 0)
            gene.transcription_time = trans.get('time_range', '')
            gene.data_sources.append('Transcription BED')
            
            # THINK: Use transcription category if regulation category is still unknown
            if gene.regulation_category in ['', 'Unknown'] and gene.transcription_category:
                gene.regulation_category = gene.transcription_category
        
        # === STEP 4: Load Translation data ===
        if gene_id in self.translation_data:
            transl = self.translation_data[gene_id]
            gene.protein_level = transl.get('protein_level', 0)
            gene.translation_phase = transl.get('phase', '')
            gene.translation_category = transl.get('category', '')
            gene.translation_order = transl.get('order', 0)
            gene.rbs_category = transl.get('rbs_strength', '')
            gene.protein_stability = transl.get('stability', '')
            gene.is_ncrna = transl.get('is_ncrna', False)
            gene.tx_level = transl.get('tx_level', 0)
            gene.tl_level = transl.get('tl_level', 0)
            gene.efficiency = transl.get('efficiency', 0)
            gene.data_sources.append('Translation BED')
            
            # THINK: Use translation category if still unknown
            if gene.regulation_category in ['', 'Unknown'] and gene.translation_category:
                gene.regulation_category = gene.translation_category
        
        # === STEP 5: Load Promoter data ===
        if gene_id in self.promoter_data:
            prom = self.promoter_data[gene_id]
            gene.has_promoter = True
            gene.promoter_type = prom.get('name', '')
            gene.data_sources.append('Promoters BED')
        
        # === STEP 6: Load Operon data ===
        if gene_id in self.operon_data:
            op = self.operon_data[gene_id]
            gene.operon_id = op.get('operon_id', '')
            gene.operon_position = op.get('position', 0)
            gene.operon_members = op.get('members', [])
            gene.data_sources.append('Operons BED')
            
            # THINK: If still unknown and part of operon, check operon members
            if gene.regulation_category in ['', 'Unknown']:
                for member_id in gene.operon_members:
                    if member_id != gene_id and member_id in self.gff3_data:
                        member_product = self.gff3_data[member_id].get('product', '')
                        member_category = self._derive_category(member_product)
                        if member_category and member_category not in ['', 'Unknown', 'Hypothetical']:
                            gene.regulation_category = f"Operon: {member_category}"
                            gene.notes = f"Category derived from operon member {member_id}"
                            break
        
        # === FINAL CHECK: Set confidence based on data sources ===
        gene.confidence = len(gene.data_sources) / 6.0  # 6 possible sources
        
        return gene
    
    # =========================================================================
    # REPORTING
    # =========================================================================
    
    def get_statistics(self) -> dict:
        """Get comprehensive statistics"""
        stats = {
            'total_genes': len(self.genes),
            'files_loaded': len(self.loaded_files),
            'file_counts': self.file_gene_counts.copy(),
            'categories': defaultdict(int),
            'regulation_types': defaultdict(int),
            'transcription_phases': defaultdict(int),
            'rbs_strengths': defaultdict(int),
            'ncrna_count': 0,
            'unknown_count': 0,
            'hypothetical_count': 0
        }
        
        for gene in self.genes.values():
            # Categories
            cat = gene.function_category or gene.regulation_category or 'Unknown'
            stats['categories'][cat] += 1
            
            if cat == 'Unknown':
                stats['unknown_count'] += 1
            elif cat == 'Hypothetical':
                stats['hypothetical_count'] += 1
            
            # Regulation types
            if gene.regulation_type:
                stats['regulation_types'][gene.regulation_type] += 1
            
            # Transcription phases
            if gene.transcription_phase:
                stats['transcription_phases'][gene.transcription_phase] += 1
            
            # RBS strengths
            if gene.rbs_strength:
                stats['rbs_strengths'][gene.rbs_strength] += 1
            
            # ncRNA
            if gene.is_ncrna:
                stats['ncrna_count'] += 1
        
        return stats
    
    def generate_report(self) -> str:
        """Generate comprehensive integration report"""
        stats = self.get_statistics()
        
        report = """
╔══════════════════════════════════════════════════════════════════════════════╗
║                    SMART DATA ENGINE - INTEGRATION REPORT                     ║
╚══════════════════════════════════════════════════════════════════════════════╝

📁 FILES LOADED:
"""
        for file in self.loaded_files:
            fname = os.path.basename(file)
            report += f"   ✅ {fname}\n"
        
        report += f"""
📊 GENE COUNTS BY SOURCE:
"""
        for source, count in stats['file_counts'].items():
            report += f"   • {source}: {count} genes\n"
        
        report += f"""
🧬 TOTAL INTEGRATED GENES: {stats['total_genes']}

📈 CATEGORY DISTRIBUTION:
"""
        # Sort categories by count
        sorted_cats = sorted(stats['categories'].items(), key=lambda x: -x[1])
        for cat, count in sorted_cats[:15]:  # Top 15
            pct = count / stats['total_genes'] * 100
            bar = '█' * int(pct / 5)
            report += f"   {cat[:25]:25} {count:4} ({pct:5.1f}%) {bar}\n"
        
        report += f"""
🔬 REGULATION TYPES:
"""
        for reg_type, count in sorted(stats['regulation_types'].items(), key=lambda x: -x[1]):
            pct = count / stats['total_genes'] * 100
            report += f"   • {reg_type}: {count} ({pct:.1f}%)\n"
        
        report += f"""
⏱️ TRANSCRIPTION PHASES:
"""
        for phase, count in sorted(stats['transcription_phases'].items()):
            report += f"   • {phase}: {count} genes\n"
        
        report += f"""
🧪 RBS STRENGTH DISTRIBUTION:
"""
        for rbs, count in sorted(stats['rbs_strengths'].items(), key=lambda x: -x[1]):
            report += f"   • {rbs}: {count} genes\n"
        
        report += f"""
⚠️ DATA QUALITY:
   • ncRNA (No Translation): {stats['ncrna_count']}
   • Unknown Category: {stats['unknown_count']}
   • Hypothetical: {stats['hypothetical_count']}
   • Known Function: {stats['total_genes'] - stats['unknown_count'] - stats['hypothetical_count']}

🧠 THINKING LOG (Last 10):
"""
        for log in self.thinking_log[-10:]:
            report += f"   {log}\n"
        
        return report
