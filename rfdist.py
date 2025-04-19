import sys
from Bio import Phylo
from io import StringIO

def get_splits(tree):
    terminals = set(clade.name for clade in tree.get_terminals())
    splits = set()
    for clade in tree.find_clades(order='preorder'):
        if clade.is_terminal():
            continue
        split = frozenset(t.name for t in clade.get_terminals())
        if 0 < len(split) < len(terminals):  # exclude trivial splits
            splits.add(split)
    return splits

def read_newick(file_path):
    with open(file_path, 'r') as f:
        newick = f.read().strip()
    return Phylo.read(StringIO(newick), "newick")

def compute_rf(tree1, tree2):
    splits1 = get_splits(tree1)
    splits2 = get_splits(tree2)
    shared = splits1 & splits2
    rf_distance = len(splits1 - shared) + len(splits2 - shared)
    max_splits = len(splits1 | splits2)
    return rf_distance, max_splits

def main(tree1_file, tree2_file):
    tree1 = read_newick(tree1_file)
    tree2 = read_newick(tree2_file)
    rf, max_rf = compute_rf(tree1, tree2)
    print(f"RF Distance: {rf}")
    print(f"Max Possible RF Distance: {max_rf}")
    print(f"Normalized RF Distance: {rf / max_rf:.4f}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rfdist.py <tree1.nwk> <tree2.nwk>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
