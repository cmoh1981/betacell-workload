#!/usr/bin/env python3
"""
run_pipeline.py - Run the complete deep learning workload pipeline

Usage:
    python run_pipeline.py              # Run all steps
    python run_pipeline.py --step 1     # Run only data preparation
    python run_pipeline.py --step 2     # Run only training
    python run_pipeline.py --step 3     # Run only validation
"""

import argparse
import subprocess
import sys
from pathlib import Path

SCRIPTS = {
    1: ("01_prepare_data.py", "Data Preparation"),
    2: ("03_train.py", "Model Training"),
    3: ("04_validate.py", "Model Validation")
}


def run_step(step: int):
    """Run a single pipeline step."""
    if step not in SCRIPTS:
        print(f"Invalid step: {step}. Valid steps: 1, 2, 3")
        return False

    script, name = SCRIPTS[step]
    script_path = Path(__file__).parent / script

    print("\n" + "=" * 70)
    print(f"STEP {step}: {name}")
    print(f"Running: {script}")
    print("=" * 70 + "\n")

    result = subprocess.run(
        [sys.executable, str(script_path)],
        cwd=script_path.parent
    )

    if result.returncode != 0:
        print(f"\nError in step {step}: {name}")
        return False

    print(f"\n✓ Step {step} complete: {name}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Run beta-cell workload deep learning pipeline"
    )
    parser.add_argument(
        "--step", type=int, choices=[1, 2, 3],
        help="Run specific step (1=prepare, 2=train, 3=validate)"
    )
    args = parser.parse_args()

    print("=" * 70)
    print("BETA-CELL WORKLOAD DEEP LEARNING PIPELINE")
    print("=" * 70)
    print("\nSteps:")
    for step, (script, name) in SCRIPTS.items():
        print(f"  {step}. {name} ({script})")

    if args.step:
        # Run single step
        success = run_step(args.step)
    else:
        # Run all steps
        print("\nRunning complete pipeline...")
        for step in [1, 2, 3]:
            success = run_step(step)
            if not success:
                print(f"\nPipeline stopped at step {step}")
                sys.exit(1)

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print("\nOutputs in: results/deep_learning/")
    print("  - training_data.h5ad")
    print("  - model_weights.pt")
    print("  - cwi_scores.csv")
    print("  - state_predictions.csv")
    print("  - gene_importance.csv")
    print("  - validation_report.json")
    print("  - figures/")


if __name__ == "__main__":
    main()
