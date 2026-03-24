"""
Genome Manager for JCVI Genome Analyzer
Handles genome data loading and annotation management
"""

from typing import List, Tuple

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

from data_structures import Gene, PromoterAnnotation, RegulationAnnotation, TranscriptionDynamics, TranslationDynamics


class GenomeManager:
    """Manages genome data and annotations"""
    
    def __init__(self):
        self.genes = []
        self.sequence = ""
        self.sequence_name = ""
        self.sequence_length = 0
        self.promoter_annotations = {}
        self.query_promoter_annotations = {}
        self.regulation_data = {}
        self.transcription_dynamics = {}
        self.translation_dynamics = {}
        
        # === GFF AS SINGLE SOURCE OF TRUTH FOR CDS ===
        # CDS features extracted from GFF - defines translatable regions
        self.cds_features = {}  # gene_id -> {start, end, strand, phase, sequence, is_valid}
    
    def load_fasta(self, fasta_path: str) -> Tuple[str, str, int]:
        """Load FASTA file"""
        with open(fasta_path, 'r') as f:
            record = SeqIO.read(f, "fasta")
            seq = str(record.seq).upper()
            return record.id, seq, len(seq)
    
    def load_gff(self, gff_path: str) -> List[Gene]:
        """Load GFF3 file - properly merge gene and CDS information"""
        
        # First pass: collect all features
        gene_features = {}  # gene features
        cds_features = {}   # CDS features
        
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                feature_type = parts[2]
                if feature_type not in ['gene', 'CDS', 'mRNA']:
                    continue
                
                # Parse attributes - handle different separators and formats
                attrs = {}
                attr_string = parts[8]
                
                # Split by semicolon
                for attr in attr_string.split(';'):
                    attr = attr.strip()
                    if not attr:
                        continue
                    
                    # Handle both '=' and ' ' as separators
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attrs[key.strip()] = value.strip()
                    elif ' ' in attr and '"' in attr:
                        # Handle format like: product "some protein"
                        parts_attr = attr.split(' ', 1)
                        if len(parts_attr) == 2:
                            key = parts_attr[0].strip()
                            value = parts_attr[1].strip().strip('"')
                            attrs[key] = value
                
                # Store feature info
                feature_info = {
                    'type': feature_type,
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6],
                    'attrs': attrs
                }
                
                # Get feature ID
                feature_id = attrs.get('ID', attrs.get('locus_tag', attrs.get('Name', '')))
                parent_id = attrs.get('Parent', '')
                
                if feature_type == 'gene':
                    gene_features[feature_id] = feature_info
                elif feature_type == 'CDS':
                    # Link to parent gene if exists
                    key = parent_id if parent_id else feature_id
                    if key not in cds_features:
                        cds_features[key] = []
                    cds_features[key].append(feature_info)
        
        # Second pass: merge gene and CDS information
        genes = []
        processed = set()
        
        # Process genes with their CDS
        for gene_id, gene_info in gene_features.items():
            if gene_id in processed:
                continue
            
            processed.add(gene_id)
            
            # Get attributes from gene
            attrs = gene_info['attrs'].copy()
            
            # Merge with CDS attributes if available
            if gene_id in cds_features:
                for cds_info in cds_features[gene_id]:
                    cds_attrs = cds_info['attrs']
                    # CDS attributes override gene attributes for specific fields
                    if 'product' in cds_attrs:
                        attrs['product'] = cds_attrs['product']
                    if 'function' in cds_attrs:
                        attrs['function'] = cds_attrs['function']
                    if 'note' in cds_attrs:
                        if 'note' not in attrs:
                            attrs['note'] = cds_attrs['note']
                    if 'gene' in cds_attrs:
                        attrs['gene'] = cds_attrs['gene']
                    if 'protein_id' in cds_attrs:
                        attrs['protein_id'] = cds_attrs['protein_id']
                    if 'translation' in cds_attrs:
                        attrs['translation'] = cds_attrs['translation']
            
            # Extract gene name (try multiple fields)
            gene_name = (attrs.get('Name') or 
                        attrs.get('gene') or 
                        attrs.get('locus_tag') or 
                        attrs.get('gene_name') or
                        gene_id)
            
            # Extract function (try multiple fields)
            function = (attrs.get('product') or 
                       attrs.get('function') or 
                       attrs.get('note') or 
                       attrs.get('description') or
                       'Unknown function')
            
            # Extract product
            product = (attrs.get('product') or 
                      attrs.get('protein') or
                      'Unknown product')
            
            gene = Gene(
                id=gene_id,
                name=gene_name,
                start=gene_info['start'],
                end=gene_info['end'],
                strand=gene_info['strand'],
                type=gene_info['type'],
                function=function,
                product=product
            )
            
            genes.append(gene)
        
        # Process orphan CDS (CDS without gene features)
        for cds_id, cds_list in cds_features.items():
            if cds_id in processed or cds_id in gene_features:
                continue
            
            processed.add(cds_id)
            
            # Use first CDS
            cds_info = cds_list[0]
            attrs = cds_info['attrs']
            
            gene_id = attrs.get('ID', attrs.get('locus_tag', cds_id))
            gene_name = (attrs.get('Name') or 
                        attrs.get('gene') or 
                        attrs.get('locus_tag') or 
                        gene_id)
            
            function = (attrs.get('product') or 
                       attrs.get('function') or 
                       attrs.get('note') or 
                       'Unknown function')
            
            product = attrs.get('product', 'Unknown product')
            
            gene = Gene(
                id=gene_id,
                name=gene_name,
                start=cds_info['start'],
                end=cds_info['end'],
                strand=cds_info['strand'],
                type='CDS',
                function=function,
                product=product
            )
            
            genes.append(gene)
        
        # Sort genes by position
        genes.sort(key=lambda g: g.start)
        
        # === STORE CDS FEATURES (GFF as single source of truth) ===
        # Store all CDS entries for later sequence extraction
        self.cds_features.clear()
        for cds_key, cds_list in cds_features.items():
            for cds_info in cds_list:
                # Determine gene_id - prefer parent gene ID, fallback to CDS ID
                parent_id = cds_info['attrs'].get('Parent', '')
                cds_id = cds_info['attrs'].get('ID', cds_info['attrs'].get('locus_tag', cds_key))
                gene_id = parent_id if parent_id and parent_id in gene_features else cds_key
                
                # Parse phase/frame from GFF column 8 (0, 1, or 2)
                phase = 0
                raw_phase = cds_info.get('phase', '.')
                if raw_phase in ['0', '1', '2']:
                    phase = int(raw_phase)
                
                self.cds_features[gene_id] = {
                    'start': cds_info['start'],
                    'end': cds_info['end'],
                    'strand': cds_info['strand'],
                    'phase': phase,
                    'sequence': '',  # Will be populated by extract_cds_from_gff
                    'protein_length': 0,
                    'is_valid': False,
                    'validation_issues': [],
                    'cds_id': cds_id,
                    'product': cds_info['attrs'].get('product', '')
                }
        
        print(f"\nDEBUG: Loaded {len(genes)} genes from GFF")
        print(f"DEBUG: Found {len(self.cds_features)} CDS features (single source of truth)")
        print(f"DEBUG: Sample genes:")
        for gene in genes[:3]:
            print(f"  - {gene.id}: {gene.name} | {gene.function}")
        
        return genes
    
    def extract_cds_from_gff(self):
        """
        Extract CDS sequences using ONLY GFF coordinates (single source of truth).
        
        This method:
        1. Uses start/end/strand from GFF CDS features
        2. Applies reverse complement for minus strand
        3. Applies phase offset for partial codons
        4. Validates sequences biologically (length, start/stop codons)
        
        Must be called AFTER load_fasta() and load_gff().
        """
        if not self.sequence:
            print("WARNING: No genome sequence loaded, cannot extract CDS sequences")
            return
        
        if not self.cds_features:
            print("WARNING: No CDS features found in GFF")
            return
        
        # Validation constants
        VALID_START_CODONS = {'ATG', 'GTG', 'TTG'}
        VALID_STOP_CODONS = {'TAA', 'TAG', 'TGA'}
        MIN_PROTEIN_LENGTH = 30
        
        valid_count = 0
        invalid_count = 0
        
        for gene_id, cds_info in self.cds_features.items():
            try:
                # Extract using GFF coordinates (1-based, convert to 0-based)
                start = max(0, cds_info['start'] - 1)
                end = min(len(self.sequence), cds_info['end'])
                strand = cds_info['strand']
                phase = cds_info['phase']
                
                if start >= end:
                    cds_info['sequence'] = ''
                    cds_info['is_valid'] = False
                    cds_info['validation_issues'] = ['Invalid coordinates (start >= end)']
                    invalid_count += 1
                    continue
                
                # Extract raw sequence from FASTA
                cds_seq = self.sequence[start:end]
                
                # Apply reverse complement for minus strand
                if strand == '-':
                    if BIOPYTHON_AVAILABLE:
                        cds_seq = str(Seq(cds_seq).reverse_complement())
                    else:
                        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                                     'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
                        cds_seq = ''.join(complement.get(b, b) for b in reversed(cds_seq))
                
                # Apply phase/frame offset (0, 1, or 2 nucleotides)
                if phase > 0:
                    cds_seq = cds_seq[phase:]
                
                cds_info['sequence'] = cds_seq
                cds_info['protein_length'] = len(cds_seq) // 3
                
                # === BIOLOGICAL VALIDATION ===
                issues = []
                
                # Check length
                if len(cds_seq) < 9:
                    issues.append('Ultra-short sequence (<3 codons)')
                elif cds_info['protein_length'] < MIN_PROTEIN_LENGTH:
                    issues.append(f'Short protein ({cds_info["protein_length"]}aa < {MIN_PROTEIN_LENGTH}aa)')
                
                # Check divisible by 3
                if len(cds_seq) % 3 != 0:
                    issues.append('Length not divisible by 3')
                
                # Check start codon
                start_codon = cds_seq[:3].upper() if len(cds_seq) >= 3 else ''
                if start_codon not in VALID_START_CODONS:
                    issues.append(f'Invalid start codon ({start_codon})')
                
                # Check stop codon
                stop_codon = cds_seq[-3:].upper() if len(cds_seq) >= 3 else ''
                if stop_codon not in VALID_STOP_CODONS:
                    issues.append(f'Missing stop codon ({stop_codon})')
                
                cds_info['validation_issues'] = issues
                cds_info['is_valid'] = len(issues) == 0
                
                if cds_info['is_valid']:
                    valid_count += 1
                else:
                    invalid_count += 1
                    
            except Exception as e:
                cds_info['sequence'] = ''
                cds_info['is_valid'] = False
                cds_info['validation_issues'] = [f'Extraction error: {e}']
                invalid_count += 1
        
        print(f"DEBUG: Extracted CDS sequences from GFF coordinates:")
        print(f"  - Valid CDS: {valid_count}")
        print(f"  - Invalid CDS: {invalid_count}")
        print(f"  - Total: {len(self.cds_features)}")
    
    def _parse_promoter_name(self, name: str) -> dict:
        """
        Parse promoter name field from BED file column 4.
        New format: gene_id|strand_info|M10:motif@score|M35:motif@score|GC:XX.XX%
        Old format: simple string
        
        Returns dict with parsed fields.
        """
        import re
        
        result = {
            'gene_id': '',
            'strand_info': '',
            'motif_10_seq': '',
            'motif_10_score': 0,
            'motif_35_seq': '',
            'motif_35_score': 0,
            'gc_content': 0.0
        }
        
        # Check if this is the new structured format (contains pipe delimiters)
        if '|' in name:
            parts = name.split('|')
            
            # Parse each component
            for i, part in enumerate(parts):
                part = part.strip()
                
                if i == 0:
                    # First part is gene_id
                    result['gene_id'] = part
                elif i == 1:
                    # Second part is strand_info (BOTH, FWD, REV)
                    result['strand_info'] = part
                elif part.startswith('M10:'):
                    # -10 box: M10:TATAAT@128
                    match = re.match(r'M10:([A-Z]+)@(\d+)', part)
                    if match:
                        result['motif_10_seq'] = match.group(1)
                        result['motif_10_score'] = int(match.group(2))
                elif part.startswith('M35:'):
                    # -35 box: M35:TTAACA@112
                    match = re.match(r'M35:([A-Z]+)@(\d+)', part)
                    if match:
                        result['motif_35_seq'] = match.group(1)
                        result['motif_35_score'] = int(match.group(2))
                elif part.startswith('GC:'):
                    # GC content: GC:20.67%
                    match = re.match(r'GC:([\d.]+)%?', part)
                    if match:
                        result['gc_content'] = float(match.group(1))
        else:
            # Old format - just use the name as gene_id
            result['gene_id'] = name
        
        return result
    
    def load_bed_file(self, bed_path: str, annotation_type: str, 
                     is_query: bool = False) -> int:
        """Load BED file with promoter annotations"""
        count = 0
        annotations = {}
        
        with open(bed_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3] if len(parts) > 3 else f"{annotation_type}_{count}"
                score = float(parts[4]) if len(parts) > 4 else 1000.0
                strand = parts[5] if len(parts) > 5 else '+'
                
                # Parse the name field to extract new structured data
                parsed = self._parse_promoter_name(name)
                
                annotation = PromoterAnnotation(
                    chrom=chrom,
                    start=start,
                    end=end,
                    name=name,
                    score=score,
                    strand=strand,
                    annotation_type=annotation_type,
                    gene_id=parsed['gene_id'],
                    strand_info=parsed['strand_info'],
                    motif_10_seq=parsed['motif_10_seq'],
                    motif_10_score=parsed['motif_10_score'],
                    motif_35_seq=parsed['motif_35_seq'],
                    motif_35_score=parsed['motif_35_score'],
                    gc_content=parsed['gc_content']
                )
                
                key = (start, end)
                annotations[key] = annotation
                count += 1
        
        if is_query:
            self.query_promoter_annotations[annotation_type] = annotations
        else:
            self.promoter_annotations[annotation_type] = annotations
        
        return count
    
    def extract_gene_sequences(self):
        """Extract gene sequences from loaded FASTA"""
        if not self.sequence or not self.genes:
            return
        
        for gene in self.genes:
            try:
                # Extract sequence
                start = max(0, gene.start - 1)  # 0-indexed
                end = min(len(self.sequence), gene.end)
                
                gene_seq = self.sequence[start:end]
                
                # Handle reverse strand
                if gene.strand == '-':
                    gene_seq = str(Seq(gene_seq).reverse_complement())
                
                gene.sequence = gene_seq
                
            except Exception as e:
                print(f"Warning: Could not extract sequence for {gene.id}: {e}")
                gene.sequence = ""
    
    def annotate_genes_with_promoters(self):
        """Annotate genes with promoter information"""
        likely_promoters = self.promoter_annotations.get('likely_promoters', {})
        internal_operons = self.promoter_annotations.get('internal_operons', {})
        
        for gene in self.genes:
            # Check for likely promoter
            for (start, end), annotation in likely_promoters.items():
                if abs(gene.start - end) < 500 or abs(gene.start - start) < 500:
                    gene.has_promoter = True
                    break
            
            # Check for internal operon
            for (start, end), annotation in internal_operons.items():
                if start <= gene.start <= end or start <= gene.end <= end:
                    gene.is_operon_internal = True
                    break
    
    def populate_translation_cds_sequences(self):
        """
        Populate CDS sequences for all translation dynamics entries.
        
        NEW ARCHITECTURE (GFF as single source of truth):
        1. First call extract_cds_from_gff() to populate cds_features
        2. Copy GFF-derived sequences to translation_dynamics
        3. Only fallback to BED coordinates if no GFF CDS
        """
        if not self.sequence:
            print("Warning: No genome sequence loaded, cannot populate CDS sequences")
            return
        
        # === STEP 1: Extract CDS from GFF (single source of truth) ===
        if not self.cds_features:
            print("DEBUG: No CDS features from GFF - this may indicate GFF loading issue")
        else:
            self.extract_cds_from_gff()
        
        # === STEP 2: Copy GFF CDS sequences to translation dynamics ===
        gff_count = 0
        fallback_count = 0
        
        for gene_id, trans_data in self.translation_dynamics.items():
            # Check if we have GFF-derived CDS sequence
            if gene_id in self.cds_features and self.cds_features[gene_id].get('sequence'):
                gff_cds = self.cds_features[gene_id]
                trans_data.cds_sequence = gff_cds['sequence']
                trans_data.start = gff_cds['start']
                trans_data.end = gff_cds['end']
                trans_data.strand = gff_cds['strand']
                gff_count += 1
            elif not trans_data.cds_sequence:
                # Fallback: extract from BED coordinates (not recommended)
                try:
                    start = max(0, trans_data.start - 1)
                    end = min(len(self.sequence), trans_data.end)
                    
                    if start < end:
                        cds_seq = self.sequence[start:end]
                        
                        if trans_data.strand == '-':
                            if BIOPYTHON_AVAILABLE:
                                cds_seq = str(Seq(cds_seq).reverse_complement())
                            else:
                                complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                                             'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
                                cds_seq = ''.join(complement.get(b, b) for b in reversed(cds_seq))
                        
                        trans_data.cds_sequence = cds_seq
                        fallback_count += 1
                except Exception as e:
                    print(f"Warning: Could not extract CDS for {gene_id}: {e}")
                    trans_data.cds_sequence = ""
        
        print(f"DEBUG: Populated CDS sequences for {len(self.translation_dynamics)} entries")
        print(f"  - From GFF (source of truth): {gff_count}")
        print(f"  - From BED fallback: {fallback_count}")
    
    def load_regulation_bed(self, bed_path: str, is_query: bool = False) -> int:
        """
        Load Regulation Analysis BED file
        Format: gene|reg_type|category|strength
        """
        import re
        count = 0
        data = {}
        
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
                
                # Parse name: gene|reg_type|category|strength
                name_parts = name.split('|')
                if len(name_parts) >= 4:
                    gene_id = name_parts[0]
                    reg_type = name_parts[1]
                    category = name_parts[2]
                    strength = name_parts[3]
                    
                    data[gene_id] = RegulationAnnotation(
                        gene_id=gene_id,
                        regulation_type=reg_type,
                        category=category,
                        strength=strength,
                        start=start,
                        end=end,
                        strand=strand,
                        score=score
                    )
                    count += 1
        
        self.regulation_data = data
        print(f"DEBUG: Loaded {count} regulation annotations")
        return count
    
    def load_transcription_dynamics_bed(self, bed_path: str, is_query: bool = False) -> int:
        """
        Load Transcription Dynamics BED file with multiple tracks
        
        Track 1: Transcription_Level - gene_id|Level:XX%
        Track 2: Transcription_Order - order|gene_id|Phase:category|level%
        Track 3: Transcription_Timeline - gene_id|time_range|level%
        """
        import re
        count = 0
        data = {}
        current_track = ""
        
        with open(bed_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Track header
                if line.startswith('track'):
                    if 'Transcription_Level' in line:
                        current_track = 'level'
                    elif 'Transcription_Order' in line:
                        current_track = 'order'
                    elif 'Transcription_Timeline' in line:
                        current_track = 'timeline'
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
                
                # Parse based on track type
                gene_id = ""
                transcription_level = 0.0
                time_range = ""
                phase = ""
                category = ""
                order = 0
                
                name_parts = name.split('|')
                
                if current_track == 'level':
                    # Format: gene_id|Level:XX%
                    gene_id = name_parts[0]
                    if len(name_parts) >= 2:
                        level_match = re.search(r'Level:(\d+)%', name_parts[1])
                        if level_match:
                            transcription_level = float(level_match.group(1))
                
                elif current_track == 'order':
                    # Format: order|gene_id|Phase:category|level%
                    if len(name_parts) >= 4:
                        try:
                            order = int(name_parts[0])
                        except:
                            order = 0
                        gene_id = name_parts[1]
                        
                        # Parse Phase:category
                        phase_match = re.search(r'Phase(\d+):(.+)', name_parts[2])
                        if phase_match:
                            phase = f"Phase{phase_match.group(1)}"
                            category = phase_match.group(2)
                        
                        # Parse level%
                        level_match = re.search(r'(\d+)%', name_parts[3])
                        if level_match:
                            transcription_level = float(level_match.group(1))
                
                elif current_track == 'timeline':
                    # Format: gene_id|time_range|level%
                    gene_id = name_parts[0]
                    if len(name_parts) >= 2:
                        time_range = name_parts[1]
                    if len(name_parts) >= 3:
                        level_match = re.search(r'(\d+)%', name_parts[2])
                        if level_match:
                            transcription_level = float(level_match.group(1))
                
                if not gene_id:
                    continue
                
                # Create or update entry
                if gene_id not in data:
                    data[gene_id] = TranscriptionDynamics(
                        gene_id=gene_id,
                        transcription_level=transcription_level,
                        time_range=time_range,
                        phase=phase,
                        category=category,
                        order=order,
                        start=start,
                        end=end,
                        strand=strand,
                        score=score,
                        color=color
                    )
                else:
                    # Update existing entry with new info
                    existing = data[gene_id]
                    if transcription_level > existing.transcription_level:
                        existing.transcription_level = transcription_level
                    if time_range and not existing.time_range:
                        existing.time_range = time_range
                    if phase and not existing.phase:
                        existing.phase = phase
                    if category and not existing.category:
                        existing.category = category
                    if order > 0 and existing.order == 0:
                        existing.order = order
                    if color != (128, 128, 128):
                        existing.color = color
                
                count += 1
        
        self.transcription_dynamics = data
        print(f"DEBUG: Loaded {len(data)} transcription dynamics entries from {count} lines")
        return len(data)
    
    def load_translation_dynamics_bed(self, bed_path: str, is_query: bool = False) -> int:
        """
        Load Translation Dynamics BED file with multiple tracks
        
        Track 1: Translation_Level - gene_id|Protein:XX%|RBS:strength
        Track 2: Translation_Order - order|gene_id|Phase:category|level%
        Track 3: RBS_Strength - gene_id|RBS:strength|Stability:stability
        Track 4: Protein_Stability - gene_id|Stability:stability|HalfLife:time
        Track 5: ncRNA_No_Translation - gene_id|ncRNA|No_Translation
        Track 6: Transcription_vs_Translation - gene_id|TX:XX%|TL:XX%|Efficiency:XX%
        """
        import re
        count = 0
        data = {}
        current_track = ""
        
        with open(bed_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Track header
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
                    elif 'Transcription_vs_Translation' in line:
                        current_track = 'efficiency'
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
                
                # Parse based on track type
                gene_id = ""
                protein_level = 0.0
                rbs_strength = ""
                stability = ""
                phase = ""
                category = ""
                order = 0
                tx_level = 0.0
                tl_level = 0.0
                efficiency = 0.0
                is_ncrna = False
                
                name_parts = name.split('|')
                
                if current_track == 'level':
                    # Format: gene_id|Protein:XX%|RBS:strength
                    gene_id = name_parts[0]
                    for part in name_parts[1:]:
                        protein_match = re.search(r'Protein:(\d+)%', part)
                        if protein_match:
                            protein_level = float(protein_match.group(1))
                        rbs_match = re.search(r'RBS:(\w+)', part)
                        if rbs_match:
                            rbs_strength = rbs_match.group(1)
                
                elif current_track == 'order':
                    # Format: order|gene_id|Phase:category|level%
                    if len(name_parts) >= 4:
                        try:
                            order = int(name_parts[0])
                        except:
                            order = 0
                        gene_id = name_parts[1]
                        
                        # Parse Phase:category
                        phase_match = re.search(r'Phase(\d+):(.+)', name_parts[2])
                        if phase_match:
                            phase = f"Phase{phase_match.group(1)}"
                            category = phase_match.group(2)
                        
                        # Parse level%
                        level_match = re.search(r'(\d+)%', name_parts[3])
                        if level_match:
                            protein_level = float(level_match.group(1))
                
                elif current_track == 'rbs':
                    # Format: gene_id|RBS:strength|Stability:stability
                    gene_id = name_parts[0]
                    for part in name_parts[1:]:
                        rbs_match = re.search(r'RBS:(\w+)', part)
                        if rbs_match:
                            rbs_strength = rbs_match.group(1)
                        stab_match = re.search(r'Stability:(.+)', part)
                        if stab_match:
                            stability = stab_match.group(1)
                
                elif current_track == 'stability':
                    # Format: gene_id|Stability:stability|HalfLife:XXmin
                    gene_id = name_parts[0]
                    protein_half_life = 0.0  # Will be set if found
                    for part in name_parts[1:]:
                        stab_match = re.search(r'Stability:(.+)', part)
                        if stab_match:
                            stability = stab_match.group(1).strip()
                        # === DATA-DRIVEN: Extract HalfLife from BED ===
                        halflife_match = re.search(r'HalfLife:(\d+(?:\.\d+)?)\s*(min|hr|h|sec|s)?', part, re.IGNORECASE)
                        if halflife_match:
                            hl_value = float(halflife_match.group(1))
                            hl_unit = halflife_match.group(2).lower() if halflife_match.group(2) else 'min'
                            # Convert to seconds
                            if hl_unit in ['hr', 'h']:
                                protein_half_life = hl_value * 3600
                            elif hl_unit in ['sec', 's']:
                                protein_half_life = hl_value
                            else:  # minutes (default)
                                protein_half_life = hl_value * 60
                
                elif current_track == 'ncrna':
                    # Format: gene_id|ncRNA|No_Translation
                    gene_id = name_parts[0]
                    is_ncrna = True
                
                elif current_track == 'efficiency':
                    # Format: gene_id|TX:XX%|TL:XX%|Efficiency:XX%
                    gene_id = name_parts[0]
                    for part in name_parts[1:]:
                        tx_match = re.search(r'TX:(\d+)%', part)
                        if tx_match:
                            tx_level = float(tx_match.group(1))
                        tl_match = re.search(r'TL:(\d+)%', part)
                        if tl_match:
                            tl_level = float(tl_match.group(1))
                        eff_match = re.search(r'Efficiency:(\d+)%', part)
                        if eff_match:
                            efficiency = float(eff_match.group(1))
                
                if not gene_id:
                    continue
                
                # Create or update entry
                if gene_id not in data:
                    # === DATA-DRIVEN: Derive half-life from stability if not explicit ===
                    derived_half_life = protein_half_life if 'protein_half_life' in dir() and protein_half_life > 0 else 0.0
                    if derived_half_life == 0.0 and stability:
                        # Derive from stability category (JCVI syn3.0 typical values)
                        stability_to_halflife = {
                            'stable': 7200.0,      # 2 hours
                            'very_stable': 14400.0, # 4 hours
                            'unstable': 1800.0,     # 30 minutes
                            'labile': 600.0,        # 10 minutes
                            'medium': 3600.0,       # 1 hour
                        }
                        derived_half_life = stability_to_halflife.get(stability.lower(), 0.0)
                    
                    # Determine parameter source
                    if 'protein_half_life' in dir() and protein_half_life > 0:
                        params_source = "bed_file"
                        uses_defaults = False
                    elif derived_half_life > 0:
                        params_source = "derived_from_stability"
                        uses_defaults = False
                    else:
                        params_source = "default"
                        uses_defaults = True
                    
                    # === GFF AS SINGLE SOURCE OF TRUTH FOR CDS ===
                    # Use GFF coordinates if available, otherwise log warning
                    cds_start = start
                    cds_end = end
                    cds_strand = strand
                    cds_sequence = ''
                    cds_source = 'bed_fallback'
                    has_gff_cds = False
                    
                    if gene_id in self.cds_features:
                        gff_cds = self.cds_features[gene_id]
                        cds_start = gff_cds['start']
                        cds_end = gff_cds['end']
                        cds_strand = gff_cds['strand']
                        cds_sequence = gff_cds.get('sequence', '')
                        cds_source = 'gff'
                        has_gff_cds = True
                    else:
                        # Log warning but continue - BED entry doesn't match GFF CDS
                        print(f"WARNING: BED entry {gene_id} has no CDS in GFF - using BED coordinates as fallback")
                    
                    data[gene_id] = TranslationDynamics(
                        gene_id=gene_id,
                        protein_level=protein_level,
                        rbs_strength=rbs_strength,
                        stability=stability,
                        phase=phase,
                        category=category,
                        order=order,
                        tx_level=tx_level,
                        tl_level=tl_level,
                        efficiency=efficiency,
                        is_ncrna=is_ncrna,
                        start=cds_start,  # Use GFF coordinates
                        end=cds_end,      # Use GFF coordinates
                        strand=cds_strand, # Use GFF strand
                        score=score,
                        color=color,
                        protein_half_life=derived_half_life,
                        kinetic_params_source=params_source,
                        uses_default_params=uses_defaults,
                        cds_sequence=cds_sequence  # Pre-populate from GFF
                    )
                else:
                    # Update existing entry with new info
                    existing = data[gene_id]
                    if protein_level > existing.protein_level:
                        existing.protein_level = protein_level
                    if rbs_strength and not existing.rbs_strength:
                        existing.rbs_strength = rbs_strength
                    if stability and not existing.stability:
                        existing.stability = stability
                    if phase and not existing.phase:
                        existing.phase = phase
                    if category and not existing.category:
                        existing.category = category
                    if order > 0 and existing.order == 0:
                        existing.order = order
                    if tx_level > 0:
                        existing.tx_level = tx_level
                    if tl_level > 0:
                        existing.tl_level = tl_level
                    if efficiency > 0:
                        existing.efficiency = efficiency
                    if is_ncrna:
                        existing.is_ncrna = True
                    if color != (128, 128, 128):
                        existing.color = color
                    # === DATA-DRIVEN: Update half-life if found in this track ===
                    if 'protein_half_life' in dir() and protein_half_life > 0 and existing.protein_half_life == 0:
                        existing.protein_half_life = protein_half_life
                        existing.kinetic_params_source = "bed_file"
                        existing.uses_default_params = False
                
                count += 1
        
        # Count how many have GFF CDS vs BED fallback
        gff_count = sum(1 for gid in data if gid in self.cds_features)
        bed_only_count = len(data) - gff_count
        
        self.translation_dynamics = data
        print(f"DEBUG: Loaded {len(data)} translation dynamics entries from {count} lines")
        print(f"  - With GFF CDS (source of truth): {gff_count}")
        print(f"  - BED-only (fallback): {bed_only_count}")
        return len(data)
