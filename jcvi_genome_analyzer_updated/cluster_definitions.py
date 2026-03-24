"""
Cluster Definitions for 28-Qubit Quantum Network Model
=======================================================

Defines 28 functional gene clusters for JCVI-syn3.0 minimal genome.
Each cluster maps to one qubit in the quantum circuit.

Clusters are based on biological function categories:
- Superposition encodes intermediate transcription states
- Entanglement models co-regulation between related clusters
- Measurement probabilities reflect transcription distribution
"""

from typing import Dict, List, Tuple, Optional
import re


# =============================================================================
# 28 FUNCTIONAL CLUSTER DEFINITIONS
# =============================================================================
# Each cluster has: name, gene patterns (regex), category keywords

CLUSTER_DEFINITIONS: Dict[int, Dict] = {
    0: {
        "name": "Ribosome Assembly",
        "patterns": [r"^rps[A-Z]?", r"^rpl[A-Z]?", r"^rpm[A-Z]?"],
        "categories": ["Ribosome", "ribosomal"],
        "description": "Ribosomal proteins (30S + 50S subunits)"
    },
    1: {
        "name": "Translation",
        "patterns": [r"^tuf", r"^fus[A-Z]?", r"^prf[A-Z]?", r"^inf[A-Z]?", r"^efp"],
        "categories": ["Translation"],
        "description": "Translation factors and elongation"
    },
    2: {
        "name": "tRNA Synthetases",
        "patterns": [r"[a-z]{3}S$", r"^gat[A-Z]", r"^fmt"],
        "categories": ["tRNA", "aminoacyl"],
        "description": "Aminoacyl-tRNA synthetases"
    },
    3: {
        "name": "DNA Replication",
        "patterns": [r"^dna[A-Z]", r"^pol[A-Z]?", r"^ssb", r"^gyr[A-Z]"],
        "categories": ["Replication", "DNA replication"],
        "description": "DNA replication machinery"
    },
    4: {
        "name": "Transcription",
        "patterns": [r"^rpo[A-Z]", r"^nus[A-Z]", r"^greA"],
        "categories": ["Transcription", "RNA polymerase"],
        "description": "RNA polymerase and transcription factors"
    },
    5: {
        "name": "RNA Processing",
        "patterns": [r"^rnc", r"^rnh[A-Z]?", r"^pnp", r"^rnr"],
        "categories": ["RNA processing", "RNase"],
        "description": "RNA processing and degradation"
    },
    6: {
        "name": "DNA Repair",
        "patterns": [r"^rec[A-Z]", r"^ung", r"^mut[A-Z]", r"^lig[A-Z]?"],
        "categories": ["Repair", "recombination"],
        "description": "DNA repair and recombination"
    },
    7: {
        "name": "Cell Division",
        "patterns": [r"^fts[A-Z]", r"^min[A-Z]", r"^sep[A-Z]"],
        "categories": ["Division", "septum"],
        "description": "Cell division and septum formation"
    },
    8: {
        "name": "Membrane Biogenesis",
        "patterns": [r"^mra[A-Z]", r"^mur[A-Z]", r"^mrc[A-Z]"],
        "categories": ["Membrane", "peptidoglycan"],
        "description": "Membrane and peptidoglycan synthesis"
    },
    9: {
        "name": "Lipid Synthesis",
        "patterns": [r"^pls[A-Z]", r"^fab[A-Z]", r"^acc[A-Z]"],
        "categories": ["Lipid", "fatty acid"],
        "description": "Lipid and fatty acid biosynthesis"
    },
    10: {
        "name": "Membrane Transport",
        "patterns": [r"^sec[A-Z]", r"^yid[A-Z]", r"^opp[A-Z]"],
        "categories": ["Transport", "permease"],
        "description": "Membrane transport systems"
    },
    11: {
        "name": "Protein Folding",
        "patterns": [r"^gro[ELS]", r"^dna[KJ]", r"^hsp", r"^tig"],
        "categories": ["Chaperone", "folding", "heat shock"],
        "description": "Chaperones and protein folding"
    },
    12: {
        "name": "Protein Degradation",
        "patterns": [r"^lon", r"^clp[A-Z]", r"^ftsH", r"^hsl"],
        "categories": ["Protease", "degradation"],
        "description": "Proteases and protein degradation"
    },
    13: {
        "name": "Nucleotide Metabolism",
        "patterns": [r"^pyr[A-Z]", r"^pur[A-Z]", r"^nrd[A-Z]", r"^adk"],
        "categories": ["Nucleotide", "purine", "pyrimidine"],
        "description": "Nucleotide biosynthesis"
    },
    14: {
        "name": "Energy/ATP",
        "patterns": [r"^atp[A-Z]", r"^pgi", r"^pfk[A-Z]?"],
        "categories": ["ATP", "energy", "synthase"],
        "description": "ATP synthesis and energy metabolism"
    },
    15: {
        "name": "Glycolysis",
        "patterns": [r"^gap[A-Z]?", r"^pgk", r"^eno", r"^pyk", r"^fba"],
        "categories": ["Glycolysis", "glycolytic"],
        "description": "Glycolysis enzymes"
    },
    16: {
        "name": "Amino Acid Metabolism",
        "patterns": [r"^ilv[A-Z]", r"^leu[A-Z]", r"^ser[A-Z]", r"^met[A-Z]"],
        "categories": ["Amino acid", "biosynthesis"],
        "description": "Amino acid biosynthesis"
    },
    17: {
        "name": "Cofactor Synthesis",
        "patterns": [r"^coa[A-Z]", r"^nad[A-Z]", r"^fol[A-Z]", r"^thi[A-Z]"],
        "categories": ["Cofactor", "vitamin", "coenzyme"],
        "description": "Cofactor and vitamin synthesis"
    },
    18: {
        "name": "Stress Response",
        "patterns": [r"^csp[A-Z]", r"^grpE", r"^hfq", r"^rpo[SH]"],
        "categories": ["Stress", "cold shock", "heat shock"],
        "description": "Stress response proteins"
    },
    19: {
        "name": "Signal Transduction",
        "patterns": [r"^pho[A-Z]", r"kinase$", r"^che[A-Z]"],
        "categories": ["Signal", "kinase", "regulation"],
        "description": "Signal transduction and kinases"
    },
    20: {
        "name": "Metal Homeostasis",
        "patterns": [r"^zur", r"^znu[A-Z]", r"^mnt[A-Z]", r"^fer"],
        "categories": ["Metal", "iron", "zinc"],
        "description": "Metal ion homeostasis"
    },
    21: {
        "name": "Cell Wall",
        "patterns": [r"^ddl", r"^dap[A-Z]", r"^alr"],
        "categories": ["Cell wall", "peptidoglycan precursor"],
        "description": "Cell wall precursor synthesis"
    },
    22: {
        "name": "Secretion",
        "patterns": [r"^ffh", r"^yaj[A-Z]", r"^lsp[A-Z]", r"^lgt"],
        "categories": ["Secretion", "signal peptide"],
        "description": "Protein secretion and export"
    },
    23: {
        "name": "Quality Control",
        "patterns": [r"^ssrA", r"^smpB", r"^arfA"],
        "categories": ["tmRNA", "ribosome rescue"],
        "description": "Translation quality control"
    },
    24: {
        "name": "Gene Regulation",
        "patterns": [r"^sig[A-Z]", r"^hrcA", r"^lexA"],
        "categories": ["Transcription factor", "regulator", "sigma"],
        "description": "Transcriptional regulators"
    },
    25: {
        "name": "ncRNA",
        "patterns": [r"^trn[A-Z]", r"^rrf", r"^rrl", r"^rrs"],
        "categories": ["tRNA", "rRNA", "ncRNA"],
        "description": "Non-coding RNAs (tRNA, rRNA)"
    },
    26: {
        "name": "Hypothetical Essential",
        "patterns": [r"^hyp", r"^y[a-z]{2}[A-Z]"],
        "categories": ["Hypothetical", "unknown function"],
        "description": "Essential genes of unknown function"
    },
    27: {
        "name": "Other Essential",
        "patterns": [],
        "categories": [],
        "description": "Other essential genes not in above categories"
    }
}


# =============================================================================
# ENTANGLEMENT TOPOLOGY (Chain Structure)
# =============================================================================
# Biologically motivated connections between related clusters
# Format: (cluster_a, cluster_b) - bidirectional entanglement

ENTANGLEMENT_PAIRS: List[Tuple[int, int]] = [
    # Central Dogma chain
    (0, 1),   # Ribosome ↔ Translation
    (1, 2),   # Translation ↔ tRNA Synthetases
    (4, 5),   # Transcription ↔ RNA Processing
    (3, 4),   # DNA Replication ↔ Transcription
    (3, 6),   # DNA Replication ↔ DNA Repair
    
    # Gene expression control
    (4, 24),  # Transcription ↔ Gene Regulation
    (24, 18), # Gene Regulation ↔ Stress Response
    
    # Cell structure chain
    (8, 9),   # Membrane Biogenesis ↔ Lipid Synthesis
    (8, 21),  # Membrane Biogenesis ↔ Cell Wall
    (7, 21),  # Cell Division ↔ Cell Wall
    (10, 22), # Membrane Transport ↔ Secretion
    
    # Protein homeostasis
    (11, 12), # Protein Folding ↔ Protein Degradation
    (11, 23), # Protein Folding ↔ Quality Control
    (1, 11),  # Translation ↔ Protein Folding
    
    # Metabolism chain
    (14, 15), # Energy/ATP ↔ Glycolysis
    (15, 16), # Glycolysis ↔ Amino Acid Metabolism
    (13, 14), # Nucleotide Metabolism ↔ Energy/ATP
    (16, 17), # Amino Acid Metabolism ↔ Cofactor Synthesis
    
    # Support systems
    (17, 20), # Cofactor Synthesis ↔ Metal Homeostasis
    (19, 18), # Signal Transduction ↔ Stress Response
    
    # ncRNA and essential
    (2, 25),  # tRNA Synthetases ↔ ncRNA
    (26, 27), # Hypothetical Essential ↔ Other Essential
]


# =============================================================================
# GENE-TO-CLUSTER MAPPING FUNCTIONS
# =============================================================================

def assign_gene_to_cluster(
    gene_id: str,
    gene_name: str = "",
    category: str = "",
    function: str = ""
) -> int:
    """
    Assign a gene to one of the 28 functional clusters.
    
    Priority:
    1. Match gene name against patterns
    2. Match category/function keywords
    3. Default to cluster 27 (Other Essential)
    
    Args:
        gene_id: Gene identifier (e.g., "JCVISYN3A_0001")
        gene_name: Gene name (e.g., "dnaA", "rpsA")
        category: Category from BED file (e.g., "Ribosome Assembly")
        function: Function description
    
    Returns:
        Cluster ID (0-27)
    """
    # Use gene_name if available, else extract from gene_id
    name_to_check = gene_name.lower() if gene_name else gene_id.lower()
    category_lower = category.lower() if category else ""
    function_lower = function.lower() if function else ""
    
    # First pass: match gene name patterns
    for cluster_id, cluster_def in CLUSTER_DEFINITIONS.items():
        for pattern in cluster_def["patterns"]:
            if re.search(pattern, name_to_check, re.IGNORECASE):
                return cluster_id
    
    # Second pass: match category/function keywords
    combined_text = f"{category_lower} {function_lower}"
    for cluster_id, cluster_def in CLUSTER_DEFINITIONS.items():
        for cat_keyword in cluster_def["categories"]:
            if cat_keyword.lower() in combined_text:
                return cluster_id
    
    # Default: Other Essential
    return 27


def build_gene_cluster_mapping(gene_data: List[Dict]) -> Dict[str, int]:
    """
    Build mapping from all genes to their cluster IDs.
    
    Args:
        gene_data: List of gene dicts with 'gene_id', 'gene_name', 'category', 'function'
    
    Returns:
        Dict mapping gene_id -> cluster_id
    """
    mapping = {}
    for gene in gene_data:
        gene_id = gene.get('gene_id', '')
        gene_name = gene.get('gene_name', gene_id)
        category = gene.get('category', '')
        function = gene.get('function', '')
        
        cluster_id = assign_gene_to_cluster(gene_id, gene_name, category, function)
        mapping[gene_id] = cluster_id
    
    return mapping


def aggregate_cluster_levels(
    gene_data: List[Dict],
    gene_cluster_mapping: Dict[str, int],
    level_key: str = 'current_level'
) -> Dict[int, float]:
    """
    Aggregate transcription levels for each cluster.
    
    Args:
        gene_data: List of gene dicts with transcription data
        gene_cluster_mapping: Dict mapping gene_id -> cluster_id
        level_key: Key for transcription level in gene dict
    
    Returns:
        Dict mapping cluster_id -> normalized average level [0, 1]
    """
    cluster_levels = {i: [] for i in range(28)}
    
    for gene in gene_data:
        gene_id = gene.get('gene_id', '')
        level = gene.get(level_key, 0.0)
        
        # Normalize to [0, 1] if level is percentage (0-100)
        if level > 1.0:
            level = level / 100.0
        
        cluster_id = gene_cluster_mapping.get(gene_id, 27)
        cluster_levels[cluster_id].append(level)
    
    # Calculate average for each cluster
    result = {}
    for cluster_id, levels in cluster_levels.items():
        if levels:
            result[cluster_id] = sum(levels) / len(levels)
        else:
            result[cluster_id] = 0.0
    
    return result


def get_cluster_gene_counts(gene_cluster_mapping: Dict[str, int]) -> Dict[int, int]:
    """
    Count genes in each cluster.
    
    Args:
        gene_cluster_mapping: Dict mapping gene_id -> cluster_id
    
    Returns:
        Dict mapping cluster_id -> gene count
    """
    counts = {i: 0 for i in range(28)}
    for cluster_id in gene_cluster_mapping.values():
        counts[cluster_id] += 1
    return counts


def get_cluster_name(cluster_id: int) -> str:
    """Get human-readable name for a cluster."""
    if cluster_id in CLUSTER_DEFINITIONS:
        return CLUSTER_DEFINITIONS[cluster_id]["name"]
    return f"Cluster {cluster_id}"


def get_cluster_info(cluster_id: int) -> Dict:
    """Get full cluster definition."""
    return CLUSTER_DEFINITIONS.get(cluster_id, {
        "name": f"Unknown Cluster {cluster_id}",
        "patterns": [],
        "categories": [],
        "description": ""
    })
