"""
Beta-Cell Genome-Scale Metabolic Model Analysis
================================================

Performs constraint-based analysis of beta-cell metabolism:
1. Load and adapt human GEM for beta-cell context
2. Flux Balance Analysis (FBA) for metabolic states
3. Gene knockout analysis for workload genes
4. Flux variability analysis for key pathways

Author: Beta-Cell Workload Analysis Pipeline
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import warnings

warnings.filterwarnings('ignore')

# Paths
SCRIPT_DIR = Path(__file__).parent.resolve()
WORKLOAD_DIR = SCRIPT_DIR.parent.parent
DATA_DIR = WORKLOAD_DIR / "data" / "gem_models"
RESULTS_DIR = WORKLOAD_DIR / "results" / "gem_analysis"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Beta-cell specific objective reactions
BETACELL_OBJECTIVES = {
    "insulin_secretion": "DM_ins_LSEC",  # Insulin demand
    "atp_production": "ATPM",  # ATP maintenance
    "biomass": "biomass_reaction",  # Growth/maintenance
}

# Beta-cell metabolic constraints (relative to default)
BETACELL_CONSTRAINTS = {
    # High glucose uptake
    "EX_glc__D_e": (-20.0, 1000.0),  # Glucose exchange
    # Active pyruvate carboxylase (anaplerosis)
    "PC": (0.0, 1000.0),  # Pyruvate carboxylase
    # Low lactate production (unlike other cells)
    "EX_lac__L_e": (0.0, 5.0),  # Limit lactate export
    # Active mitochondrial metabolism
    "ATPS4m": (0.0, 1000.0),  # ATP synthase
}

# Workload-related genes for knockout analysis
WORKLOAD_GENES = {
    # MR-validated genes
    "PDX1": "transcription_factor",
    "SLC2A2": "glucose_transporter",
    "MAFA": "transcription_factor",
    # Metabolic capacity genes
    "GCK": "glucose_sensor",
    "INS": "insulin",
    "PCSK1": "proinsulin_processing",
    "PCSK2": "proinsulin_processing",
    # Stress response genes
    "DDIT3": "upr_terminal",
    "ATF4": "upr_adaptive",
    "XBP1": "upr_adaptive",
    # Mitochondrial genes
    "TFAM": "mitochondrial_biogenesis",
    "PPARGC1A": "mitochondrial_biogenesis",
}


class BetaCellGEM:
    """Beta-cell specific genome-scale metabolic model analysis."""

    def __init__(self, model_path: Optional[Path] = None):
        """
        Initialize beta-cell GEM analysis.

        Args:
            model_path: Path to SBML model file
        """
        self.model = None
        self.model_path = model_path
        self.original_bounds = {}

    def load_model(self, model_path: Optional[Path] = None) -> bool:
        """Load genome-scale metabolic model."""
        try:
            from cobra.io import read_sbml_model, load_model

            if model_path:
                self.model_path = model_path

            if self.model_path and self.model_path.exists():
                print(f"Loading model from: {self.model_path}")
                self.model = read_sbml_model(str(self.model_path))
            else:
                # Try to load bundled E. coli model for testing
                print("Loading test model (E. coli core)...")
                self.model = load_model("textbook")

            print(f"  Model ID: {self.model.id}")
            print(f"  Reactions: {len(self.model.reactions)}")
            print(f"  Metabolites: {len(self.model.metabolites)}")
            print(f"  Genes: {len(self.model.genes)}")

            # Store original bounds
            for rxn in self.model.reactions:
                self.original_bounds[rxn.id] = (rxn.lower_bound, rxn.upper_bound)

            return True

        except ImportError:
            print("ERROR: COBRApy not installed. Install with: pip install cobra")
            return False
        except Exception as e:
            print(f"ERROR loading model: {e}")
            return False

    def apply_betacell_constraints(self):
        """Apply beta-cell specific metabolic constraints."""
        if self.model is None:
            print("No model loaded")
            return

        print("\nApplying beta-cell specific constraints...")

        applied = []
        for rxn_id, bounds in BETACELL_CONSTRAINTS.items():
            try:
                rxn = self.model.reactions.get_by_id(rxn_id)
                rxn.bounds = bounds
                applied.append(rxn_id)
                print(f"  {rxn_id}: bounds set to {bounds}")
            except KeyError:
                # Reaction not in model
                pass

        print(f"  Applied constraints to {len(applied)} reactions")

    def run_fba(self, objective: Optional[str] = None) -> Dict:
        """
        Run Flux Balance Analysis.

        Args:
            objective: Objective reaction ID (uses model default if None)

        Returns:
            Dictionary with FBA results
        """
        if self.model is None:
            return {'status': 'error', 'message': 'No model loaded'}

        try:
            if objective:
                self.model.objective = objective

            solution = self.model.optimize()

            return {
                'status': solution.status,
                'objective_value': solution.objective_value,
                'fluxes': solution.fluxes.to_dict(),
                'objective_id': str(self.model.objective)
            }

        except Exception as e:
            return {'status': 'error', 'message': str(e)}

    def run_pfba(self) -> Dict:
        """
        Run parsimonious FBA (minimize total flux).

        Returns:
            Dictionary with pFBA results
        """
        if self.model is None:
            return {'status': 'error', 'message': 'No model loaded'}

        try:
            from cobra.flux_analysis import pfba

            solution = pfba(self.model)

            return {
                'status': solution.status,
                'objective_value': solution.objective_value,
                'fluxes': solution.fluxes.to_dict(),
                'total_flux': solution.fluxes.abs().sum()
            }

        except Exception as e:
            return {'status': 'error', 'message': str(e)}

    def run_fva(self, reactions: Optional[List[str]] = None,
                fraction_of_optimum: float = 0.9) -> pd.DataFrame:
        """
        Run Flux Variability Analysis.

        Args:
            reactions: List of reaction IDs (all if None)
            fraction_of_optimum: Fraction of optimal growth to maintain

        Returns:
            DataFrame with min/max flux for each reaction
        """
        if self.model is None:
            return pd.DataFrame()

        try:
            from cobra.flux_analysis import flux_variability_analysis

            fva_result = flux_variability_analysis(
                self.model,
                reaction_list=reactions,
                fraction_of_optimum=fraction_of_optimum
            )

            return fva_result

        except Exception as e:
            print(f"FVA error: {e}")
            return pd.DataFrame()

    def gene_knockout_analysis(self, genes: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Perform single gene knockout analysis.

        Args:
            genes: List of gene IDs to knockout (uses WORKLOAD_GENES if None)

        Returns:
            DataFrame with knockout effects
        """
        if self.model is None:
            return pd.DataFrame()

        try:
            from cobra.flux_analysis import single_gene_deletion

            if genes is None:
                genes = list(WORKLOAD_GENES.keys())

            # Filter to genes in model
            model_genes = {g.id for g in self.model.genes}
            valid_genes = [g for g in genes if g in model_genes]

            if not valid_genes:
                print("No valid genes found in model")
                # Try gene name matching
                for gene in genes:
                    matches = [g for g in model_genes if gene in g]
                    if matches:
                        valid_genes.extend(matches[:1])

            if not valid_genes:
                print("Using all model genes for demonstration")
                valid_genes = list(model_genes)[:20]

            print(f"Analyzing {len(valid_genes)} genes...")

            # Get wildtype growth
            wt_solution = self.model.optimize()
            wt_growth = wt_solution.objective_value

            results = []
            for gene_id in valid_genes:
                with self.model:
                    try:
                        gene = self.model.genes.get_by_id(gene_id)
                        gene.knock_out()
                        ko_solution = self.model.optimize()

                        results.append({
                            'gene': gene_id,
                            'wt_growth': wt_growth,
                            'ko_growth': ko_solution.objective_value,
                            'growth_ratio': ko_solution.objective_value / wt_growth if wt_growth > 0 else 0,
                            'status': ko_solution.status,
                            'essential': ko_solution.objective_value < 0.01 * wt_growth
                        })
                    except Exception as e:
                        results.append({
                            'gene': gene_id,
                            'wt_growth': wt_growth,
                            'ko_growth': np.nan,
                            'growth_ratio': np.nan,
                            'status': 'error',
                            'essential': False
                        })

            return pd.DataFrame(results)

        except Exception as e:
            print(f"Gene knockout error: {e}")
            return pd.DataFrame()

    def analyze_metabolic_state(self, state: str = "normal") -> Dict:
        """
        Analyze metabolic state corresponding to workload states.

        Args:
            state: One of 'resting', 'active', 'stressed', 'exhausted', 'failing'

        Returns:
            Dictionary with metabolic analysis results
        """
        if self.model is None:
            return {}

        # Define state-specific constraints
        state_constraints = {
            "resting": {
                "glucose_uptake": 5.0,  # Low glucose
                "description": "Low demand, balanced metabolism"
            },
            "active": {
                "glucose_uptake": 15.0,  # High glucose
                "description": "High insulin secretion demand"
            },
            "stressed": {
                "glucose_uptake": 20.0,  # Very high glucose
                "description": "UPR activated, increased demand"
            },
            "exhausted": {
                "glucose_uptake": 20.0,  # High glucose
                "atp_maintenance": 0.5,  # Reduced ATP capacity
                "description": "Metabolic dysfunction"
            },
            "failing": {
                "glucose_uptake": 10.0,  # Reduced glucose sensing
                "atp_maintenance": 0.3,  # Severely reduced capacity
                "description": "Beta-cell failure"
            }
        }

        if state not in state_constraints:
            state = "normal"
            constraints = {"glucose_uptake": 10.0, "description": "Normal metabolism"}
        else:
            constraints = state_constraints[state]

        print(f"\nAnalyzing metabolic state: {state}")
        print(f"  {constraints.get('description', '')}")

        with self.model:
            # Apply state-specific constraints
            # Try to find glucose exchange reaction
            glc_rxns = [r for r in self.model.reactions if 'glc' in r.id.lower() and 'EX_' in r.id]
            if glc_rxns:
                glc_rxn = glc_rxns[0]
                glc_rxn.lower_bound = -constraints["glucose_uptake"]

            # Run FBA
            solution = self.model.optimize()

            # Identify active pathways
            active_fluxes = solution.fluxes[solution.fluxes.abs() > 1e-6]

            return {
                'state': state,
                'growth_rate': solution.objective_value,
                'n_active_reactions': len(active_fluxes),
                'total_flux': solution.fluxes.abs().sum(),
                'top_fluxes': active_fluxes.abs().nlargest(20).to_dict()
            }

    def get_pathway_fluxes(self, pathway_name: str) -> pd.DataFrame:
        """Get fluxes for reactions in a specific pathway."""
        if self.model is None:
            return pd.DataFrame()

        # Run FBA first
        solution = self.model.optimize()

        # Find reactions in pathway (by subsystem)
        pathway_rxns = []
        for rxn in self.model.reactions:
            subsystem = getattr(rxn, 'subsystem', '')
            if pathway_name.lower() in subsystem.lower():
                pathway_rxns.append({
                    'reaction_id': rxn.id,
                    'reaction_name': rxn.name,
                    'subsystem': subsystem,
                    'flux': solution.fluxes[rxn.id],
                    'bounds': (rxn.lower_bound, rxn.upper_bound)
                })

        return pd.DataFrame(pathway_rxns)


def run_betacell_gem_analysis():
    """Run complete beta-cell GEM analysis pipeline."""
    print("=" * 60)
    print("BETA-CELL GENOME-SCALE METABOLIC ANALYSIS")
    print("=" * 60)

    # Initialize GEM analysis
    gem = BetaCellGEM()

    # Try to load human model, fall back to test model
    model_path = DATA_DIR / "Recon3D.xml"
    if not model_path.exists():
        model_path = DATA_DIR / "Human-GEM.xml"

    if model_path.exists():
        success = gem.load_model(model_path)
    else:
        print("\nNo human GEM found, using E. coli test model for demonstration")
        success = gem.load_model()

    if not success:
        print("Failed to load model")
        return None

    # Apply beta-cell constraints
    gem.apply_betacell_constraints()

    # Run FBA
    print("\n" + "-" * 40)
    print("FLUX BALANCE ANALYSIS")
    print("-" * 40)
    fba_result = gem.run_fba()
    print(f"  Status: {fba_result['status']}")
    print(f"  Objective value: {fba_result['objective_value']:.4f}")

    # Run pFBA
    print("\n" + "-" * 40)
    print("PARSIMONIOUS FBA")
    print("-" * 40)
    pfba_result = gem.run_pfba()
    print(f"  Status: {pfba_result['status']}")
    print(f"  Total flux: {pfba_result.get('total_flux', 0):.2f}")

    # Gene knockout analysis
    print("\n" + "-" * 40)
    print("GENE KNOCKOUT ANALYSIS")
    print("-" * 40)
    ko_results = gem.gene_knockout_analysis()
    if not ko_results.empty:
        print(ko_results.to_string(index=False))
        ko_results.to_csv(RESULTS_DIR / "gene_knockout_results.csv", index=False)

    # Analyze metabolic states
    print("\n" + "-" * 40)
    print("METABOLIC STATE ANALYSIS")
    print("-" * 40)
    state_results = []
    for state in ["resting", "active", "stressed", "exhausted", "failing"]:
        result = gem.analyze_metabolic_state(state)
        state_results.append(result)
        print(f"  {state}: growth={result['growth_rate']:.4f}, "
              f"active_rxns={result['n_active_reactions']}")

    state_df = pd.DataFrame(state_results)
    state_df.to_csv(RESULTS_DIR / "metabolic_states.csv", index=False)

    # Run FVA for key reactions
    print("\n" + "-" * 40)
    print("FLUX VARIABILITY ANALYSIS")
    print("-" * 40)
    fva_result = gem.run_fva(fraction_of_optimum=0.9)
    if not fva_result.empty:
        print(f"  Analyzed {len(fva_result)} reactions")
        # Save top variable reactions
        fva_result['range'] = fva_result['maximum'] - fva_result['minimum']
        top_variable = fva_result.nlargest(20, 'range')
        top_variable.to_csv(RESULTS_DIR / "fva_top_variable.csv")

    print("\n" + "=" * 60)
    print("GEM ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"\nResults saved to: {RESULTS_DIR}")

    return gem


if __name__ == "__main__":
    gem = run_betacell_gem_analysis()
