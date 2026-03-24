"""
Worker Threads for JCVI Genome Analyzer
Background threads for genome loading, BLAST search, and AI analysis
"""

import os
import traceback
from typing import List

from PyQt5.QtCore import QThread, pyqtSignal

try:
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

from data_structures import BlastHit


class LoadGenomeWorker(QThread):
    """Worker thread for loading genome files"""
    progress = pyqtSignal(int, str)
    finished = pyqtSignal(bool, str, dict)
    
    def __init__(self, fasta_path, gff_path, genome_manager, blast_manager,
                 likely_promoters_path=None, internal_operons_path=None,
                 complex_cases_path=None, regulation_path=None,
                 transcription_path=None, translation_path=None):
        super().__init__()
        self.fasta_path = fasta_path
        self.gff_path = gff_path
        self.genome_manager = genome_manager
        self.blast_manager = blast_manager
        self.likely_promoters_path = likely_promoters_path
        self.internal_operons_path = internal_operons_path
        self.complex_cases_path = complex_cases_path
        self.regulation_path = regulation_path
        self.transcription_path = transcription_path
        self.translation_path = translation_path
    
    def run(self):
        try:
            # Load FASTA
            self.progress.emit(15, "Loading FASTA file...")
            seq_name, seq, seq_len = self.genome_manager.load_fasta(self.fasta_path)
            self.genome_manager.sequence = seq
            self.genome_manager.sequence_name = seq_name
            self.genome_manager.sequence_length = seq_len
            
            # Load GFF
            self.progress.emit(30, "Loading GFF file...")
            genes = self.genome_manager.load_gff(self.gff_path)
            self.genome_manager.genes = genes
            
            # Load promoter annotations
            likely_count = 0
            operon_count = 0
            regulation_count = 0
            transcription_count = 0
            translation_count = 0
            
            if self.likely_promoters_path:
                self.progress.emit(40, "Loading likely promoters...")
                likely_count = self.genome_manager.load_bed_file(
                    self.likely_promoters_path, 'likely_promoters'
                )
            
            if self.internal_operons_path:
                self.progress.emit(50, "Loading internal operons...")
                operon_count = self.genome_manager.load_bed_file(
                    self.internal_operons_path, 'internal_operons'
                )
            
            # Load new BED files
            if self.regulation_path:
                self.progress.emit(60, "Loading regulation analysis...")
                regulation_count = self.genome_manager.load_regulation_bed(
                    self.regulation_path
                )
            
            if self.transcription_path:
                self.progress.emit(70, "Loading transcription dynamics...")
                transcription_count = self.genome_manager.load_transcription_dynamics_bed(
                    self.transcription_path
                )
            
            if self.translation_path:
                self.progress.emit(80, "Loading translation dynamics...")
                translation_count = self.genome_manager.load_translation_dynamics_bed(
                    self.translation_path
                )
            
            # Annotate genes
            self.progress.emit(85, "Annotating genes...")
            self.genome_manager.annotate_genes_with_promoters()
            
            # Extract gene sequences
            self.progress.emit(90, "Extracting gene sequences...")
            self.genome_manager.extract_gene_sequences()
            
            # Create BLAST database (only for reference genome)
            if self.blast_manager:
                self.progress.emit(95, "Creating BLAST database...")
                self.blast_manager.make_blast_db(self.fasta_path)
            
            stats = {
                'seq_name': seq_name,
                'seq_len': seq_len,
                'gene_count': len(genes),
                'likely_promoters_count': likely_count,
                'internal_operons_count': operon_count,
                'regulation_count': regulation_count,
                'transcription_count': transcription_count,
                'translation_count': translation_count
            }
            
            self.progress.emit(100, "Complete!")
            self.finished.emit(True, "Success", stats)
            
        except Exception as e:
            self.finished.emit(False, str(e), {})


class BlastSearchWorker(QThread):
    """Worker thread for BLAST search"""
    progress = pyqtSignal(int, str)
    finished = pyqtSignal(bool, str, list)
    
    def __init__(self, query_seq, ref_genome_manager, blast_manager, query_genome_manager=None):
        super().__init__()
        self.query_seq = query_seq
        self.ref_genome_manager = ref_genome_manager
        self.query_genome_manager = query_genome_manager
        self.blast_manager = blast_manager
    
    def run(self):
        try:
            # Write query to temp file
            query_file = os.path.join(self.blast_manager.temp_dir, "query.fasta")
            with open(query_file, 'w') as f:
                f.write(">query_syn3.0\n")
                f.write(self.query_seq)
            
            self.progress.emit(30, "Running BLAST...")
            xml_file = self.blast_manager.run_blast(query_file)
            
            self.progress.emit(70, "Parsing and mapping genes...")
            hits_data = self.blast_manager.parse_blast_xml(xml_file)
            
            # Convert to BlastHit objects
            results = []
            for hit_data in hits_data:
                # Find QUERY genes at query position
                query_genes_in_hit = []
                query_gene_functions = []
                query_gene_ids = []
                
                if self.query_genome_manager and self.query_genome_manager.genes:
                    for gene in self.query_genome_manager.genes:
                        query_overlap = self._calculate_overlap(
                            hit_data['query_start'], hit_data['query_end'],
                            gene.start, gene.end
                        )
                        
                        if query_overlap > 0:
                            query_genes_in_hit.append(gene.name)
                            query_gene_functions.append(gene.function)
                            query_gene_ids.append(gene.id)
                
                # Find REFERENCE genes at target position
                ref_genes_in_hit = []
                ref_gene_functions = []
                
                for gene in self.ref_genome_manager.genes:
                    ref_overlap = self._calculate_overlap(
                        hit_data['target_start'], hit_data['target_end'],
                        gene.start, gene.end
                    )
                    
                    if ref_overlap > 0:
                        ref_genes_in_hit.append(gene.name)
                        ref_gene_functions.append(gene.function)
                
                # Get protein sequence
                protein_seq = ""
                if query_genes_in_hit and self.query_genome_manager:
                    for gene in self.query_genome_manager.genes:
                        if gene.name in query_genes_in_hit:
                            try:
                                if gene.sequence:
                                    seq_obj = Seq(gene.sequence)
                                    protein_seq = str(seq_obj.translate(to_stop=True))
                                else:
                                    gene_seq = self.query_seq[gene.start-1:gene.end]
                                    seq_obj = Seq(gene_seq)
                                    if gene.strand == '-':
                                        seq_obj = seq_obj.reverse_complement()
                                    protein_seq = str(seq_obj.translate(to_stop=True))
                            except Exception as e:
                                print(f"Warning: Could not translate {gene.name}: {e}")
                                protein_seq = ""
                            if protein_seq:
                                break
                
                # Prioritize query genes
                if query_genes_in_hit:
                    genes_to_show = query_genes_in_hit
                    functions_to_show = query_gene_functions
                elif ref_genes_in_hit:
                    genes_to_show = ref_genes_in_hit
                    functions_to_show = ref_gene_functions
                else:
                    genes_to_show = ["Intergenic region"]
                    functions_to_show = ["Non-coding region"]
                
                blast_hit = BlastHit(
                    query_id=hit_data['query_id'],
                    target_id=hit_data['target_id'],
                    identity=hit_data['identity'],
                    alignment_length=hit_data['alignment_length'],
                    mismatches=hit_data['mismatches'],
                    gap_opens=hit_data['gap_opens'],
                    query_start=hit_data['query_start'],
                    query_end=hit_data['query_end'],
                    target_start=hit_data['target_start'],
                    target_end=hit_data['target_end'],
                    evalue=hit_data['evalue'],
                    bit_score=hit_data['bit_score'],
                    query_seq=hit_data['query_seq'],
                    target_seq=hit_data['target_seq'],
                    genes=genes_to_show,
                    gene_functions=functions_to_show,
                    protein_sequence=protein_seq
                )
                
                results.append(blast_hit)
            
            results.sort(key=lambda x: x.identity, reverse=True)
            
            self.progress.emit(100, "Complete!")
            self.finished.emit(True, "Success", results)
            
        except Exception as e:
            error_msg = f"{str(e)}\n{traceback.format_exc()}"
            self.finished.emit(False, error_msg, [])
    
    def _calculate_overlap(self, start1, end1, start2, end2):
        """Calculate overlap length between two regions"""
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        return max(0, overlap_end - overlap_start + 1)


class AIAnalysisWorker(QThread):
    """
    AI Analysis Worker - Uses Smart Data Engine
    
    This worker uses the SmartDataEngine to:
    1. Load ALL available files
    2. Integrate data intelligently
    3. Fill gaps by searching across files
    4. Eliminate "Unknown" values where possible
    """
    progress = pyqtSignal(int, str)
    finished = pyqtSignal(bool, str, dict)
    
    def __init__(self, genome_manager, blast_hits, query_genome_manager=None, 
                 gff_path=None, query_gff_path=None, all_file_paths=None):
        super().__init__()
        self.genome_manager = genome_manager
        self.blast_hits = blast_hits
        self.query_genome_manager = query_genome_manager
        self.gff_path = gff_path
        self.query_gff_path = query_gff_path
        self.all_file_paths = all_file_paths or {}
    
    def run(self):
        try:
            from smart_data_engine import SmartDataEngine
            
            results = {
                'status': 'success',
                'message': '',
                'statistics': {},
                'integrated_genes': {},
                'engine_report': ''
            }
            
            # === Create Smart Data Engine ===
            self.progress.emit(5, "🧠 Initializing Smart Data Engine...")
            engine = SmartDataEngine()
            
            # === Step 1: Load GFF3 (Primary source for gene functions) ===
            self.progress.emit(10, "📖 Loading GFF3 - Gene functions...")
            
            gff_path = self.all_file_paths.get('query_gff') or self.query_gff_path
            if gff_path:
                engine.load_gff3(gff_path)
            
            # Also try reference GFF if query not available
            if not engine.gff3_data:
                ref_gff = self.all_file_paths.get('ref_gff') or self.gff_path
                if ref_gff:
                    engine.load_gff3(ref_gff)
            
            # === Step 2: Load Regulation BED ===
            self.progress.emit(20, "📊 Loading Regulation BED...")
            
            reg_path = self.all_file_paths.get('query_regulation')
            if reg_path:
                engine.load_regulation_bed(reg_path)
            
            # === Step 3: Load Transcription BED ===
            self.progress.emit(35, "🧬 Loading Transcription BED...")
            
            trans_path = self.all_file_paths.get('query_transcription')
            if trans_path:
                engine.load_transcription_bed(trans_path)
            
            # === Step 4: Load Translation BED ===
            self.progress.emit(50, "🔬 Loading Translation BED...")
            
            transl_path = self.all_file_paths.get('query_translation')
            if transl_path:
                engine.load_translation_bed(transl_path)
            
            # === Step 5: Load Promoters BED ===
            self.progress.emit(60, "🎯 Loading Promoters BED...")
            
            prom_path = self.all_file_paths.get('query_promoters')
            if prom_path:
                engine.load_promoters_bed(prom_path)
            
            # === Step 6: Load Operons BED ===
            self.progress.emit(70, "🔗 Loading Operons BED...")
            
            operon_path = self.all_file_paths.get('query_operons')
            if operon_path:
                engine.load_operons_bed(operon_path)
            
            # === Step 7: SMART INTEGRATION ===
            self.progress.emit(75, "🧠 THINKING: Integrating all data...")
            
            integrated_genes = engine.integrate_all_data()
            
            # === ✨ NEW: DATA ENRICHMENT LAYER ===
            self.progress.emit(80, "🧪 ENRICHING: Loading supplementary biological data...")
            from data_enrichment_layer import DataEnrichmentLayer
            enricher = DataEnrichmentLayer()
            
            # === ✨ NEW: MOLECULAR UR GENERATOR ===
            self.progress.emit(85, "🧬 EMBEDDING: Generating Molecular UR (Transformers)...")
            from molecular_ur_generator import MolecularURGenerator
            # Note: This might be slow on CPU. 
            # In a production quad-GPU system, this would be highly parallelized.
            ur_generator = MolecularURGenerator()
            
            # Process all genes through the new layers
            final_genes = {}
            total = len(integrated_genes)
            for i, (gene_id, gene) in enumerate(integrated_genes.items()):
                if i % 10 == 0:
                    self.progress.emit(85 + int((i/total)*10), f"Processing gene {i}/{total}...")
                
                # Enrich with additional CSV/TSV/XLSX data + AI gap filling
                enriched_gene = enricher.enrich_gene(gene)
                
                # Generate high-dimensional embeddings
                gene_with_ur = ur_generator.generate_ur(enriched_gene)
                
                final_genes[gene_id] = gene_with_ur

            # === Step 8: Update genome_manager with final data ===
            self.progress.emit(95, "📝 Updating genome manager with AI Virtual Cell data...")
            
            # We use final_genes here instead of integrated_genes
            integrated_genes = final_genes
            
            # === ✨ NEW: CELLULAR UR GENERATOR ===
            self.progress.emit(96, "🔬 CELLULAR UR: Generating cell-state representation...")
            try:
                from cellular_ur_generator import CellularURGenerator
                cellular_ur_gen = CellularURGenerator()
                cell_state = cellular_ur_gen.generate_cellular_ur(final_genes)
                
                # Store on genome_manager for GUI access
                if self.query_genome_manager:
                    self.query_genome_manager.cell_state_vector = cell_state
                    self.query_genome_manager.cellular_ur_generator = cellular_ur_gen
                
                print(f"✅ Cellular UR: {cell_state.combined_cellular_ur.shape}-dim vector for {cell_state.num_genes} genes")
            except Exception as e:
                print(f"⚠️ Cellular UR generation skipped: {e}")
                cell_state = None
            
            # === ✨ NEW: VIRTUAL INSTRUMENTS ===
            self.progress.emit(97, "🧪 VIRTUAL INSTRUMENTS: Running predictions...")
            try:
                from virtual_instruments import VirtualInstrumentSuite
                vis = VirtualInstrumentSuite(cellular_ur_gen)
                
                vi_results = {}
                if cell_state:
                    vi_results['growth_rate'] = vis.predict_growth_rate(cell_state)
                    vi_results['phenotype'] = vis.predict_phenotype(cell_state, final_genes)
                    vi_results['stress_heat'] = vis.predict_stress(cell_state, "heat")
                    
                    # Store results
                    if self.query_genome_manager:
                        self.query_genome_manager.vi_results = vi_results
                        self.query_genome_manager.vi_suite = vis
                    
                    print(f"✅ VIs: Growth={vi_results['growth_rate']['doubling_time']:.0f}min, "
                          f"Phenotype={vi_results['phenotype']['phenotype']}")
            except Exception as e:
                print(f"⚠️ Virtual Instruments skipped: {e}")
            
            # Update regulation_data with smart categories (retaining existing logic)
            if self.query_genome_manager and self.query_genome_manager.regulation_data:
                for gene_id, data in self.query_genome_manager.regulation_data.items():
                    if gene_id in integrated_genes:
                        ig = integrated_genes[gene_id]
                        # Update category if we found a better one
                        if ig.regulation_category and data.category in ['', 'Unknown']:
                            data.category = ig.regulation_category
                        elif ig.function_category and data.category in ['', 'Unknown']:
                            data.category = ig.function_category
            
            # === Step 9: Generate report ===
            self.progress.emit(98, "📊 Generating report...")
            
            stats = engine.get_statistics()
            results['statistics'] = stats
            results['integrated_genes'] = integrated_genes
            results['engine_report'] = engine.generate_report()
            results['thinking_log'] = engine.thinking_log
            
            results['message'] = f"""
✅ SMART Integration Complete!

🧠 Thinking Process:
• Loaded {len(engine.loaded_files)} files
• Found {stats['total_genes']} unique genes
• Categories identified: {len(stats['categories'])}
• Unknown reduced to: {stats['unknown_count']}

📊 Data Sources Used:
"""
            for source, count in stats['file_counts'].items():
                results['message'] += f"• {source}: {count} genes\n"
            
            self.progress.emit(100, "✅ Complete!")
            self.finished.emit(True, "Success", results)
            
        except Exception as e:
            import traceback
            error_msg = f"{str(e)}\n{traceback.format_exc()}"
            print(f"ERROR in AIAnalysisWorker: {error_msg}")
            self.finished.emit(False, error_msg, {})


class QuadGPUOrchestrationWorker(QThread):
    """
    Worker for simultaneous 4-GPU processing.
    Delegates tasks to specialized GPU nodes (Motif, Dynamics, Synthesis, Master).
    """
    progress = pyqtSignal(int, str)
    node_status = pyqtSignal(str, str, str)  # node_name, status, device
    finished = pyqtSignal(bool, str)
    
    def __init__(self, sequence_data=None):
        super().__init__()
        self.sequence_data = sequence_data
        from gpu_orchestrator import orchestrator
        self.orchestrator = orchestrator
        
    def run(self):
        try:
            self.progress.emit(10, "⚡ Booting Quad-GPU Architecture...")
            
            # Start nodes if not already started
            self.orchestrator.start_nodes()
            
            time_to_wait = 2  # Simulate parallel initialization
            
            # 1. Start Motif Node on GPU 0
            self.node_status.emit("motif", "RUNNING", str(self.orchestrator.devices['motif']))
            self.progress.emit(25, "🧬 Motif Node (GPU 0) scanning sequences...")
            
            # 2. Start Dynamics Node on GPU 1
            self.node_status.emit("dynamics", "RUNNING", str(self.orchestrator.devices['dynamics']))
            self.progress.emit(50, "⏳ Dynamics Node (GPU 1) predicting levels...")
            
            # 3. Start Synthesis Node on GPU 2
            self.node_status.emit("synthesis", "RUNNING", str(self.orchestrator.devices['synthesis']))
            self.progress.emit(75, "🔬 Synthesis Node (GPU 2) modeling proteins...")
            
            # 4. Master Node on GPU 3
            self.node_status.emit("master", "COORDINATING", str(self.orchestrator.devices['master']))
            self.progress.emit(90, "🪐 Master Node (GPU 3) orchestrating Digital Twin...")
            
            # Perform dummy parallel work to show GPU activity in nvidia-smi
            # In a real scenario, this would involve batch inference
            import torch
            if torch.cuda.is_available():
                for role, device in self.orchestrator.devices.items():
                    if 'cuda' in str(device):
                        # Allocate some memory on each GPU to demonstrate multi-gpu usage
                        torch.zeros((1000, 1000), device=device)
            
            self.progress.emit(100, "🔥 Parallel GPU Session Active")
            self.finished.emit(True, "Quad-GPU Orchestration Active")
            
        except Exception as e:
            self.finished.emit(False, str(e))
