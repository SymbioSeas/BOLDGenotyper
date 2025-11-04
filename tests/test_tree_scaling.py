#!/usr/bin/env python3
"""
Test automatic tree figure scaling based on number of tips.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from boldgenotyper import visualization, utils
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
import tempfile

# Setup logging
logger = utils.setup_logging(log_level="INFO")

print("=" * 70)
print("Testing Phylogenetic Tree Auto-Scaling")
print("=" * 70)

def create_test_tree(n_tips):
    """Create a simple test tree with n tips."""
    # Create terminal clades
    terminals = [Clade(name=f"tip_{i}") for i in range(n_tips)]

    # Create a simple star topology (all tips connect to root)
    root = Clade()
    root.clades = terminals

    tree = Tree(root=root)
    return tree

# Test with different numbers of tips
test_cases = [
    (10, "Small tree"),
    (30, "Medium tree"),
    (50, "Large tree"),
    (85, "Carcharhinus-like tree (85 consensus groups)")
]

print("\nExpected figure sizes for different tree sizes:")
print("-" * 70)

for n_tips, description in test_cases:
    tree = create_test_tree(n_tips)

    # Calculate what the figsize should be based on the logic
    height = max(10, min(50, n_tips * 0.3))
    width = 8 if n_tips <= 30 else min(14, 8 + (n_tips - 30) * 0.1)

    print(f"{description:40s} ({n_tips:3d} tips): {width:.1f}w × {height:.1f}h inches")

    # Write to temp file and test
    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        Phylo.write(tree, f, 'newick')
        tree_file = f.name

    output_file = Path(f"tests/output/test_tree_{n_tips}tips.png")
    output_file.parent.mkdir(exist_ok=True, parents=True)

    # Test the function - it should log the auto-calculated size
    visualization.plot_phylogenetic_tree(
        tree_file=tree_file,
        output_path=str(output_file),
        figsize=None,  # Trigger auto-scaling
        dpi=100  # Lower DPI for faster testing
    )

    # Clean up temp tree file
    Path(tree_file).unlink()

    print(f"  → Generated: {output_file}")

print("\n" + "=" * 70)
print("✓ Auto-scaling test completed")
print("=" * 70)
print("\nKey features:")
print("  • Height scales linearly: 0.3 inches per tip")
print("  • Minimum height: 10 inches")
print("  • Maximum height: 50 inches (caps at ~167 tips)")
print("  • Width: 8 inches for ≤30 tips, gradually increases to max 14 inches")
print("  • Can still override with explicit figsize parameter")
