"""
JCVI Genome Analyzer with AI Agents v2.0
Professional genome comparison with AI-powered promoter, transcription, and translation analysis
"""

__version__ = "2.0.0"
__author__ = "JCVI Bioinformatics"

from .data_structures import (
    Gene,
    PromoterAnnotation,
    BlastHit,
    PromoterPrediction,
    TranscriptionState,
    TranslationState,
    RegulationAnnotation,
    TranscriptionDynamics,
    TranslationDynamics
)

from .neural_networks import PromoterCNN, TranscriptionLSTM
from .ai_promoter_agent import PromoterAnalyzerAgent
from .ai_transcription_agent import TranscriptionDynamicsAgent
from .ai_translation_agent import AITranslationAgent
from .blast_manager import BlastManager
from .genome_manager import GenomeManager
from .workers import LoadGenomeWorker, BlastSearchWorker, AIAnalysisWorker

__all__ = [
    'Gene',
    'PromoterAnnotation', 
    'BlastHit',
    'PromoterPrediction',
    'TranscriptionState',
    'TranslationState',
    'RegulationAnnotation',
    'TranscriptionDynamics',
    'TranslationDynamics',
    'PromoterCNN',
    'TranscriptionLSTM',
    'PromoterAnalyzerAgent',
    'TranscriptionDynamicsAgent',
    'AITranslationAgent',
    'BlastManager',
    'GenomeManager',
    'LoadGenomeWorker',
    'BlastSearchWorker',
    'AIAnalysisWorker'
]
