"""
GEM Analysis Module for Beta-Cell Workload
==========================================

Genome-scale metabolic modeling integration with CWI workload analysis.

Modules:
- 01_download_gem_models: Download Recon3D and Human-GEM models
- 02_betacell_gem_analysis: Beta-cell specific GEM analysis
- 03_workload_metabolic_integration: Integrate CWI with metabolic modeling
"""

from pathlib import Path

MODULE_DIR = Path(__file__).parent
WORKLOAD_DIR = MODULE_DIR.parent.parent

__all__ = ['MODULE_DIR', 'WORKLOAD_DIR']
