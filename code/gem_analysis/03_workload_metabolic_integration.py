"""
Workload-Metabolic Integration Analysis
=======================================

Integrates CWI workload scores with genome-scale metabolic modeling:
1. Maps workload states to metabolic phenotypes
2. Identifies metabolic bottlenecks in stressed cells
3. Predicts therapeutic targets based on metabolic analysis

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
DATA_DIR = WORKLOAD_DIR / "data"
CWI_DIR = WORKLOAD_DIR / "results" / "multiomics"
GEM_DIR = DATA_DIR / "gem_models"
RESULTS_DIR = WORKLOAD_DIR / "results" / "gem_analysis"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Workload state to metabolic phenotype mapping
WORKLOAD_METABOLIC_MAP = {
    "S1_Resting": {
        "glucose_uptake": 5.0,
        "insulin_demand": 0.2,
        "atp_requirement": 1.0,
        "upr_activity": 0.0,
        "description": "Low metabolic demand, efficient ATP production"
    },
    "S2_Active": {
        "glucose_uptake": 15.0,
        "insulin_demand": 0.8,
        "atp_requirement": 1.5,
        "upr_activity": 0.1,
        "description": "High glucose-stimulated insulin secretion"
    },
    "S3_Stressed": {
        "glucose_uptake": 18.0,
        "insulin_demand": 1.0,
        "atp_requirement": 1.8,
        "upr_activity": 0.5,
        "description": "UPR activated, ER stress response"
    },
    "S4_Exhausted": {
        "glucose_uptake": 12.0,
        "insulin_demand": 0.6,
        "atp_requirement": 1.2,
        "upr_activity": 0.8,
        "description": "Metabolic dysfunction, reduced capacity"
    },
    "S5_Failing": {
        "glucose_uptake": 8.0,
        "insulin_demand": 0.3,
        "atp_requirement": 0.8,
        "upr_activity": 1.0,
        "description": "Beta-cell failure, dedifferentiation"
    }
}

# Metabolic pathways affected by workload
WORKLOAD_PATHWAYS = {
    "glucose_sensing": {
        "genes": ["GCK", "SLC2A2", "G6PC2"],
        "reactions": ["HEX1", "GNK", "PGI"],
        "affected_by": ["capacity", "dediff"]
    },
    "glycolysis": {
        "genes": ["PFKFB2", "PFKL", "PKM", "ALDOA"],
        "reactions": ["PFK", "FBA", "PGK", "PYK"],
        "affected_by": ["demand"]
    },
    "tca_cycle": {
        "genes": ["CS", "IDH2", "OGDH", "MDH2", "PC"],
        "reactions": ["CS", "ICDHy", "AKGDH", "MDH", "PC"],
        "affected_by": ["capacity", "stress"]
    },
    "oxidative_phosphorylation": {
        "genes": ["NDUFS1", "SDHA", "UQCRC1", "COX4I1", "ATP5F1A"],
        "reactions": ["NADH2_u10m", "SUCD", "CYOR_u10m", "CYOOm", "ATPS4m"],
        "affected_by": ["capacity", "stress"]
    },
    "insulin_synthesis": {
        "genes": ["INS", "PCSK1", "PCSK2", "CPE"],
        "reactions": ["PROTS", "PROT_FOLD"],
        "affected_by": ["demand", "stress"]
    },
    "upr_metabolism": {
        "genes": ["XBP1", "ATF4", "DDIT3", "HSPA5"],
        "reactions": ["CHAP_FOLD", "ER_STRESS"],
        "affected_by": ["stress", "dediff"]
    }
}


class WorkloadMetabolicIntegrator:
    """Integrates CWI workload analysis with GEM metabolic modeling."""

    def __init__(self):
        """Initialize integrator."""
        self.cwi_data = None
        self.gem_model = None
        self.integration_results = {}

    def load_cwi_data(self, cwi_path: Optional[Path] = None) -> bool:
        """Load CWI workload scores."""
        if cwi_path is None:
            # Try default locations
            possible_paths = [
                CWI_DIR / "workload_scores.csv",
                WORKLOAD_DIR / "results" / "tables" / "workload_scores.csv",
                WORKLOAD_DIR / "results" / "multiomics" / "cwi_for_diablo.csv"
            ]
            for path in possible_paths:
                if path.exists():
                    cwi_path = path
                    break

        if cwi_path and cwi_path.exists():
            print(f"Loading CWI data from: {cwi_path}")
            self.cwi_data = pd.read_csv(cwi_path)

            # Standardize column names
            col_mapping = {
                'cwi_predicted': 'CWI',
                'workload_state': 'Workload_State',
                'Workload_State': 'Workload_State'
            }
            self.cwi_data.rename(columns=col_mapping, inplace=True)

            print(f"  Loaded {len(self.cwi_data)} cells")
            if 'Workload_State' in self.cwi_data.columns:
                print(f"  States: {self.cwi_data['Workload_State'].value_counts().to_dict()}")
            return True
        else:
            print("No CWI data found, generating simulated data...")
            self._generate_simulated_cwi()
            return True

    def _generate_simulated_cwi(self, n_cells: int = 500):
        """Generate simulated CWI data for testing."""
        np.random.seed(42)

        states = ["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]
        state_probs = [0.2, 0.3, 0.25, 0.15, 0.1]

        cell_states = np.random.choice(states, n_cells, p=state_probs)

        # Generate CWI based on state
        cwi_means = {"S1_Resting": 0.3, "S2_Active": 0.8, "S3_Stressed": 1.3,
                     "S4_Exhausted": 2.0, "S5_Failing": 3.5}

        cwi_values = [cwi_means[s] + np.random.normal(0, 0.2) for s in cell_states]
        cwi_values = np.clip(cwi_values, 0, 5)

        self.cwi_data = pd.DataFrame({
            'cell_id': [f'cell_{i}' for i in range(n_cells)],
            'CWI': cwi_values,
            'Workload_State': cell_states,
            'condition': np.random.choice(['ND', 'T2D'], n_cells, p=[0.5, 0.5])
        })

        print(f"  Generated {n_cells} simulated cells")

    def load_gem_model(self, model_path: Optional[Path] = None) -> bool:
        """Load genome-scale metabolic model."""
        try:
            from cobra.io import read_sbml_model, load_model

            if model_path and model_path.exists():
                print(f"Loading GEM from: {model_path}")
                self.gem_model = read_sbml_model(str(model_path))
            else:
                # Try default paths
                for filename in ["Recon3D.xml", "Human-GEM.xml"]:
                    path = GEM_DIR / filename
                    if path.exists():
                        print(f"Loading GEM from: {path}")
                        self.gem_model = read_sbml_model(str(path))
                        break

                if self.gem_model is None:
                    print("Using E. coli test model for demonstration")
                    self.gem_model = load_model("textbook")

            print(f"  Model: {self.gem_model.id}")
            print(f"  Reactions: {len(self.gem_model.reactions)}")
            return True

        except ImportError:
            print("COBRApy not installed")
            return False
        except Exception as e:
            print(f"Error loading model: {e}")
            return False

    def simulate_metabolic_state(self, state: str) -> Dict:
        """Simulate metabolism for a given workload state."""
        if self.gem_model is None:
            return {}

        phenotype = WORKLOAD_METABOLIC_MAP.get(state, WORKLOAD_METABOLIC_MAP["S2_Active"])

        with self.gem_model:
            # Apply state-specific constraints
            # Find glucose exchange
            glc_rxns = [r for r in self.gem_model.reactions
                        if 'glc' in r.id.lower() and ('EX_' in r.id or 'ex_' in r.id)]
            if glc_rxns:
                glc_rxns[0].lower_bound = -phenotype["glucose_uptake"]

            # Run FBA
            solution = self.gem_model.optimize()

            # Calculate pathway fluxes
            pathway_fluxes = {}
            for pathway, info in WORKLOAD_PATHWAYS.items():
                rxn_fluxes = []
                for rxn_id in info["reactions"]:
                    try:
                        flux = solution.fluxes.get(rxn_id, 0)
                        rxn_fluxes.append(abs(flux))
                    except:
                        pass
                pathway_fluxes[pathway] = np.mean(rxn_fluxes) if rxn_fluxes else 0

            return {
                'state': state,
                'growth_rate': solution.objective_value,
                'status': solution.status,
                'pathway_fluxes': pathway_fluxes,
                'phenotype': phenotype
            }

    def analyze_state_transitions(self) -> pd.DataFrame:
        """Analyze metabolic changes across workload states."""
        print("\nAnalyzing metabolic state transitions...")

        results = []
        for state in WORKLOAD_METABOLIC_MAP.keys():
            metabolic = self.simulate_metabolic_state(state)
            if metabolic:
                row = {
                    'state': state,
                    'growth_rate': metabolic['growth_rate'],
                    'glucose_uptake': WORKLOAD_METABOLIC_MAP[state]['glucose_uptake'],
                    'insulin_demand': WORKLOAD_METABOLIC_MAP[state]['insulin_demand'],
                    'upr_activity': WORKLOAD_METABOLIC_MAP[state]['upr_activity']
                }
                # Add pathway fluxes
                for pathway, flux in metabolic.get('pathway_fluxes', {}).items():
                    row[f'flux_{pathway}'] = flux
                results.append(row)

        return pd.DataFrame(results)

    def identify_metabolic_bottlenecks(self, state: str = "S4_Exhausted") -> pd.DataFrame:
        """Identify metabolic bottlenecks in stressed/exhausted cells."""
        if self.gem_model is None:
            return pd.DataFrame()

        print(f"\nIdentifying metabolic bottlenecks for {state}...")

        try:
            from cobra.flux_analysis import flux_variability_analysis

            # Simulate the state
            phenotype = WORKLOAD_METABOLIC_MAP.get(state, WORKLOAD_METABOLIC_MAP["S4_Exhausted"])

            with self.gem_model:
                # Apply constraints
                glc_rxns = [r for r in self.gem_model.reactions
                            if 'glc' in r.id.lower() and ('EX_' in r.id or 'ex_' in r.id)]
                if glc_rxns:
                    glc_rxns[0].lower_bound = -phenotype["glucose_uptake"]

                # Run FVA
                fva = flux_variability_analysis(
                    self.gem_model,
                    fraction_of_optimum=0.9
                )

            # Identify bottlenecks (reactions with narrow flux ranges)
            fva['range'] = fva['maximum'] - fva['minimum']
            fva['is_bottleneck'] = fva['range'] < 0.1

            # Get pathway assignments
            bottlenecks = []
            for rxn_id in fva[fva['is_bottleneck']].index[:50]:
                try:
                    rxn = self.gem_model.reactions.get_by_id(rxn_id)
                    bottlenecks.append({
                        'reaction_id': rxn_id,
                        'reaction_name': rxn.name,
                        'subsystem': getattr(rxn, 'subsystem', ''),
                        'flux_min': fva.loc[rxn_id, 'minimum'],
                        'flux_max': fva.loc[rxn_id, 'maximum'],
                        'flux_range': fva.loc[rxn_id, 'range']
                    })
                except:
                    pass

            return pd.DataFrame(bottlenecks)

        except Exception as e:
            print(f"Error identifying bottlenecks: {e}")
            return pd.DataFrame()

    def predict_therapeutic_targets(self) -> pd.DataFrame:
        """Predict therapeutic targets based on metabolic analysis."""
        if self.gem_model is None:
            return pd.DataFrame()

        print("\nPredicting therapeutic targets...")

        try:
            from cobra.flux_analysis import single_gene_deletion

            # Get genes affecting workload pathways
            target_genes = set()
            for pathway, info in WORKLOAD_PATHWAYS.items():
                target_genes.update(info["genes"])

            # Find matching genes in model
            model_genes = {g.id for g in self.gem_model.genes}
            valid_genes = []
            for gene in target_genes:
                if gene in model_genes:
                    valid_genes.append(gene)
                else:
                    # Try partial matching
                    matches = [g for g in model_genes if gene.lower() in g.lower()]
                    valid_genes.extend(matches[:1])

            if not valid_genes:
                valid_genes = list(model_genes)[:20]

            # Get wildtype growth
            wt_solution = self.gem_model.optimize()
            wt_growth = wt_solution.objective_value

            # Test gene knockouts
            results = []
            for gene_id in valid_genes[:30]:
                with self.gem_model:
                    try:
                        gene = self.gem_model.genes.get_by_id(gene_id)
                        gene.knock_out()
                        ko_solution = self.gem_model.optimize()

                        # Determine target potential
                        growth_ratio = ko_solution.objective_value / wt_growth if wt_growth > 0 else 0
                        is_essential = growth_ratio < 0.1
                        is_target = 0.3 < growth_ratio < 0.9

                        results.append({
                            'gene': gene_id,
                            'wt_growth': wt_growth,
                            'ko_growth': ko_solution.objective_value,
                            'growth_ratio': growth_ratio,
                            'essential': is_essential,
                            'therapeutic_potential': 'High' if is_target else ('Essential' if is_essential else 'Low'),
                            'pathway': self._get_gene_pathway(gene_id)
                        })
                    except Exception as e:
                        pass

            return pd.DataFrame(results)

        except Exception as e:
            print(f"Error predicting targets: {e}")
            return pd.DataFrame()

    def _get_gene_pathway(self, gene_id: str) -> str:
        """Get pathway for a gene."""
        for pathway, info in WORKLOAD_PATHWAYS.items():
            if gene_id in info["genes"] or any(gene_id in g for g in info["genes"]):
                return pathway
        return "other"

    def integrate_cwi_with_gem(self) -> pd.DataFrame:
        """Main integration: map CWI states to metabolic phenotypes."""
        if self.cwi_data is None:
            print("No CWI data loaded")
            return pd.DataFrame()

        print("\n" + "=" * 60)
        print("INTEGRATING CWI WITH GEM ANALYSIS")
        print("=" * 60)

        # Group cells by workload state
        state_summary = []

        for state in self.cwi_data['Workload_State'].unique():
            state_cells = self.cwi_data[self.cwi_data['Workload_State'] == state]
            n_cells = len(state_cells)
            mean_cwi = state_cells['CWI'].mean()

            # Get metabolic simulation
            metabolic = self.simulate_metabolic_state(state)

            row = {
                'workload_state': state,
                'n_cells': n_cells,
                'pct_cells': n_cells / len(self.cwi_data) * 100,
                'mean_cwi': mean_cwi,
                'metabolic_growth': metabolic.get('growth_rate', np.nan),
                'description': WORKLOAD_METABOLIC_MAP.get(state, {}).get('description', '')
            }

            # Add pathway fluxes
            for pathway, flux in metabolic.get('pathway_fluxes', {}).items():
                row[f'{pathway}_flux'] = flux

            state_summary.append(row)

        summary_df = pd.DataFrame(state_summary)

        # Sort by workload progression
        state_order = ["S1_Resting", "S2_Active", "S3_Stressed", "S4_Exhausted", "S5_Failing"]
        summary_df['state_order'] = summary_df['workload_state'].apply(
            lambda x: state_order.index(x) if x in state_order else 99
        )
        summary_df = summary_df.sort_values('state_order').drop('state_order', axis=1)

        return summary_df

    def generate_report(self) -> str:
        """Generate comprehensive integration report."""
        report = []
        report.append("=" * 60)
        report.append("WORKLOAD-METABOLIC INTEGRATION REPORT")
        report.append("=" * 60)

        if self.cwi_data is not None:
            report.append(f"\nCWI Data: {len(self.cwi_data)} cells")
            report.append(f"States: {self.cwi_data['Workload_State'].value_counts().to_dict()}")

        if self.gem_model is not None:
            report.append(f"\nGEM Model: {self.gem_model.id}")
            report.append(f"Reactions: {len(self.gem_model.reactions)}")
            report.append(f"Genes: {len(self.gem_model.genes)}")

        # Integration summary
        integration = self.integrate_cwi_with_gem()
        if not integration.empty:
            report.append("\n" + "-" * 40)
            report.append("STATE-METABOLIC MAPPING")
            report.append("-" * 40)
            report.append(integration[['workload_state', 'n_cells', 'mean_cwi',
                                        'metabolic_growth']].to_string(index=False))

        return "\n".join(report)


def main():
    """Main integration pipeline."""
    print("=" * 60)
    print("WORKLOAD-METABOLIC INTEGRATION ANALYSIS")
    print("=" * 60)

    # Initialize integrator
    integrator = WorkloadMetabolicIntegrator()

    # Load data
    integrator.load_cwi_data()
    integrator.load_gem_model()

    # Run integration analysis
    integration_results = integrator.integrate_cwi_with_gem()
    if not integration_results.empty:
        print("\n" + "-" * 40)
        print("WORKLOAD-METABOLIC INTEGRATION")
        print("-" * 40)
        print(integration_results.to_string(index=False))
        integration_results.to_csv(RESULTS_DIR / "workload_metabolic_integration.csv", index=False)

    # Analyze state transitions
    transitions = integrator.analyze_state_transitions()
    if not transitions.empty:
        print("\n" + "-" * 40)
        print("STATE TRANSITIONS")
        print("-" * 40)
        print(transitions[['state', 'growth_rate', 'glucose_uptake', 'upr_activity']].to_string(index=False))
        transitions.to_csv(RESULTS_DIR / "state_transitions.csv", index=False)

    # Identify bottlenecks
    bottlenecks = integrator.identify_metabolic_bottlenecks("S4_Exhausted")
    if not bottlenecks.empty:
        print("\n" + "-" * 40)
        print("METABOLIC BOTTLENECKS (Exhausted State)")
        print("-" * 40)
        print(bottlenecks.head(10).to_string(index=False))
        bottlenecks.to_csv(RESULTS_DIR / "metabolic_bottlenecks.csv", index=False)

    # Predict therapeutic targets
    targets = integrator.predict_therapeutic_targets()
    if not targets.empty:
        print("\n" + "-" * 40)
        print("THERAPEUTIC TARGETS")
        print("-" * 40)
        high_potential = targets[targets['therapeutic_potential'] == 'High']
        if not high_potential.empty:
            print(high_potential.to_string(index=False))
        targets.to_csv(RESULTS_DIR / "therapeutic_targets.csv", index=False)

    # Generate report
    report = integrator.generate_report()
    report_path = RESULTS_DIR / "integration_report.txt"
    with open(report_path, 'w') as f:
        f.write(report)

    print("\n" + "=" * 60)
    print("INTEGRATION COMPLETE")
    print("=" * 60)
    print(f"\nResults saved to: {RESULTS_DIR}")

    return integrator


if __name__ == "__main__":
    integrator = main()
