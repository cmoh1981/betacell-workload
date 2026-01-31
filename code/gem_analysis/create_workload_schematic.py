"""
Beta-Cell Workload Hypothesis Schematic Generator
Creates a publication-quality figure showing diabetes progression
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Set up the figure
fig, ax = plt.subplots(1, 1, figsize=(16, 10))
ax.set_xlim(0, 16)
ax.set_ylim(0, 10)
ax.axis('off')

# Color scheme (colorblind-friendly)
colors = {
    'normal': '#2E8B57',      # Sea green
    'compensation': '#DAA520', # Goldenrod
    'stress': '#FF8C00',       # Dark orange
    'exhaustion': '#DC143C',   # Crimson
    'failure': '#8B0000',      # Dark red
    'arrow': '#4A4A4A',        # Dark gray
    'text': '#1A1A1A',         # Near black
    'light_bg': '#F5F5F5',     # Light gray background
}

# Title
ax.text(8, 9.5, 'Beta-Cell Workload Hypothesis: Diabetes Progression Model',
        fontsize=16, fontweight='bold', ha='center', va='center', color=colors['text'])

# Subtitle with demand/capacity relationship
ax.text(8, 9.0, 'Demand < Capacity  ───────────────────────────────────>  Demand >> Capacity',
        fontsize=10, ha='center', va='center', color=colors['text'], style='italic')

# Stage boxes parameters
box_width = 2.6
box_height = 5.5
box_y = 2.2
stage_x = [0.8, 3.6, 6.4, 9.2, 12.0]

# Stage data
stages = [
    {
        'name': 'NORMAL',
        'color': colors['normal'],
        'cwi': 'CWI < 0.3',
        'title_color': 'white',
        'markers': [
            ('GCK', 'up', 'Glucose sensor'),
            ('SLC2A2', 'up', 'GLUT2'),
            ('PDX1', 'up', 'Identity'),
            ('MAFA', 'up', 'Maturation'),
        ],
        'status': 'Low Workload\nHealthy Function',
        'detail': 'Insulin secretion\nmatches demand'
    },
    {
        'name': 'COMPENSATION',
        'color': colors['compensation'],
        'cwi': 'CWI 0.3-0.6',
        'title_color': 'black',
        'markers': [
            ('GCK', 'up', 'Increased'),
            ('INS', 'up', 'Hyperfunction'),
            ('HSPA5', 'up', 'ER expansion'),
        ],
        'status': 'Rising Workload\nMetabolic Strain',
        'detail': '43 metabolic\nbottlenecks emerge'
    },
    {
        'name': 'STRESS',
        'color': colors['stress'],
        'cwi': 'CWI 0.6-0.8',
        'title_color': 'white',
        'markers': [
            ('DDIT3', 'up', 'CHOP/Apoptosis'),
            ('XBP1', 'up', 'UPR active'),
            ('HSPA5', 'up', 'BiP saturated'),
            ('ATF4', 'up', 'Stress TF'),
        ],
        'status': 'ER Stress\nUPR Activated',
        'detail': 'Protein folding\noverload'
    },
    {
        'name': 'EXHAUSTION',
        'color': colors['exhaustion'],
        'cwi': 'CWI > 0.8',
        'title_color': 'white',
        'markers': [
            ('PDX1', 'down', 'Identity loss'),
            ('MAFA', 'down', 'Immaturity'),
            ('ALDH1A3', 'up', 'Dediff marker'),
            ('SOX9', 'up', 'Progenitor'),
        ],
        'status': 'Dedifferentiation\nIdentity Loss',
        'detail': 'Beta cells lose\nmature phenotype'
    },
    {
        'name': 'FAILURE',
        'color': colors['failure'],
        'cwi': 'Clinical T2D',
        'title_color': 'white',
        'markers': [
            ('Beta-cell', 'down', 'Mass < 50%'),
            ('Insulin', 'down', 'Secretion'),
            ('Glucose', 'up', 'Hyperglycemia'),
        ],
        'status': 'Clinical T2D\nBeta-Cell Failure',
        'detail': 'Fasting glucose\n> 126 mg/dL'
    },
]

# Draw stages
for i, (x, stage) in enumerate(zip(stage_x, stages)):
    # Main box
    box = FancyBboxPatch((x, box_y), box_width, box_height,
                         boxstyle="round,pad=0.05,rounding_size=0.2",
                         facecolor=stage['color'], edgecolor='black', linewidth=2,
                         alpha=0.9)
    ax.add_patch(box)

    # Stage name (header)
    ax.text(x + box_width/2, box_y + box_height - 0.4, stage['name'],
            fontsize=12, fontweight='bold', ha='center', va='center',
            color=stage['title_color'])

    # CWI indicator
    ax.text(x + box_width/2, box_y + box_height - 0.9, stage['cwi'],
            fontsize=9, ha='center', va='center', color=stage['title_color'],
            style='italic')

    # Status text
    ax.text(x + box_width/2, box_y + box_height - 1.6, stage['status'],
            fontsize=9, fontweight='bold', ha='center', va='center',
            color=stage['title_color'])

    # Molecular markers
    marker_y_start = box_y + box_height - 2.4
    for j, (gene, direction, desc) in enumerate(stage['markers']):
        y_pos = marker_y_start - j * 0.7

        # Arrow symbol
        arrow_symbol = '\u2191' if direction == 'up' else '\u2193'
        arrow_color = '#00AA00' if direction == 'up' else '#AA0000'

        # Gene name with arrow
        ax.text(x + 0.3, y_pos, f"{arrow_symbol} {gene}",
                fontsize=9, fontweight='bold', ha='left', va='center',
                color='white' if stage['title_color'] == 'white' else 'black')

        # Description
        ax.text(x + box_width - 0.2, y_pos, desc,
                fontsize=7, ha='right', va='center',
                color='white' if stage['title_color'] == 'white' else 'black',
                alpha=0.9)

    # Detail text at bottom
    ax.text(x + box_width/2, box_y + 0.5, stage['detail'],
            fontsize=8, ha='center', va='center',
            color=stage['title_color'], alpha=0.9)

# Draw connecting arrows
arrow_style = "Simple, tail_width=0.5, head_width=1.5, head_length=0.8"
for i in range(len(stage_x) - 1):
    arrow = FancyArrowPatch((stage_x[i] + box_width + 0.1, box_y + box_height/2),
                            (stage_x[i+1] - 0.1, box_y + box_height/2),
                            arrowstyle=arrow_style,
                            color=colors['arrow'], linewidth=2,
                            mutation_scale=15)
    ax.add_patch(arrow)

# Gradient bar at bottom (Workload Index)
gradient_y = 0.8
gradient_height = 0.4
for i in range(100):
    x_pos = 0.8 + i * 0.138
    # Create gradient from green to red
    r = min(1, i / 50)
    g = max(0, 1 - i / 50)
    b = 0
    rect = plt.Rectangle((x_pos, gradient_y), 0.14, gradient_height,
                          facecolor=(r, g, b), edgecolor='none')
    ax.add_patch(rect)

# Gradient bar border and labels
ax.plot([0.8, 14.6], [gradient_y, gradient_y], 'k-', linewidth=1)
ax.plot([0.8, 14.6], [gradient_y + gradient_height, gradient_y + gradient_height], 'k-', linewidth=1)
ax.plot([0.8, 0.8], [gradient_y, gradient_y + gradient_height], 'k-', linewidth=1)
ax.plot([14.6, 14.6], [gradient_y, gradient_y + gradient_height], 'k-', linewidth=1)

ax.text(7.7, gradient_y - 0.3, 'Composite Workload Index (CWI)',
        fontsize=11, fontweight='bold', ha='center', va='center', color=colors['text'])
ax.text(0.8, gradient_y + gradient_height + 0.15, 'Low', fontsize=9, ha='left', va='bottom')
ax.text(14.6, gradient_y + gradient_height + 0.15, 'High', fontsize=9, ha='right', va='bottom')

# Evidence sources box
evidence_box = FancyBboxPatch((14.8, 5.5), 1.0, 2.2,
                              boxstyle="round,pad=0.02,rounding_size=0.1",
                              facecolor=colors['light_bg'], edgecolor='gray', linewidth=1)
ax.add_patch(evidence_box)
ax.text(15.3, 7.5, 'Evidence', fontsize=8, fontweight='bold', ha='center', va='center')
ax.text(15.3, 7.1, 'Sources:', fontsize=8, fontweight='bold', ha='center', va='center')
evidence_items = ['MR Analysis', 'GEM Model', 'DisGeNET', 'LINCS']
for j, item in enumerate(evidence_items):
    ax.text(15.3, 6.6 - j*0.4, item, fontsize=7, ha='center', va='center')

# Key insight annotation
ax.annotate('Chronic metabolic demand\nexceeding capacity drives\nprogressive dysfunction',
            xy=(6.4, 2.2), xytext=(6.4, 1.5),
            fontsize=8, ha='center', va='top', style='italic',
            color=colors['text'])

plt.tight_layout()

# Save figure
output_path = r'C:\Users\ohbry\Aging and Metabolism Dropbox\Oh Chang-Myung\macbook\workload\figures\betacell_workload_model.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
print(f"Figure saved to: {output_path}")

# Also save as PDF for publication
pdf_path = output_path.replace('.png', '.pdf')
plt.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
print(f"PDF saved to: {pdf_path}")

plt.close()
print("Done!")
