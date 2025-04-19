
import sys
import numpy as np
from collections import defaultdict

class NJTreeNode:
    def __init__(self, name):
        self.name = name
        self.left = None
        self.right = None
        self.branch_length = {}

    def to_newick(self):
        if self.left is None and self.right is None:
            return self.name
        left = self.left.to_newick()
        right = self.right.to_newick()
        l_len = self.branch_length[self.left]
        r_len = self.branch_length[self.right]
        return f"({left}:{l_len:.6f},{right}:{r_len:.6f})"

def read_phylip(filename):
    with open(filename) as f:
        lines = f.readlines()
        n = int(lines[0].strip())
        names = []
        matrix = []
        for line in lines[1:]:
            if not line.strip():
                continue
            parts = line.strip().split()
            names.append(parts[0])
            matrix.append(list(map(float, parts[1:])))
        D = np.array(matrix)
        return names, D

def write_newick(root, output_file):
    with open(output_file, 'w') as f:
        newick_str = root.to_newick() + ";\n"
        f.write(newick_str)

def neighbor_joining(names, D):
    print("⚙️ Starting neighbor-joining algorithm...")

    n = len(names)
    nodes = {i: NJTreeNode(names[i]) for i in range(n)}
    active = list(range(n))
    next_node_id = n

    while len(active) > 2:
        print(f"🔄 Iteration with {len(active)} active nodes...")
        m = len(active)
        r = {i: sum(D[i][j] for j in active if j != i) for i in active}

        Q = defaultdict(lambda: float('inf'))
        for i in active:
            for j in active:
                if i < j:
                    Q[(i, j)] = (m - 2) * D[i][j] - r[i] - r[j]

        i, j = min(Q, key=Q.get)
        print(f"🔗 Joining nodes: {i} and {j}")

        new_node = NJTreeNode(f"N{next_node_id}")
        next_node_id += 1

        delta = (r[i] - r[j]) / (m - 2)
        limb_i = 0.5 * D[i][j] + 0.5 * delta
        limb_j = 0.5 * D[i][j] - 0.5 * delta

        new_node.left = nodes[i]
        new_node.right = nodes[j]
        new_node.branch_length[nodes[i]] = limb_i
        new_node.branch_length[nodes[j]] = limb_j

        new_index = next_node_id - 1
        D = np.vstack([D, [0]*D.shape[1]])
        D = np.column_stack([D, [0]*D.shape[0]])
        for k in active:
            if k != i and k != j:
                D[new_index][k] = D[k][new_index] = (D[i][k] + D[j][k] - D[i][j]) / 2
        D[new_index][new_index] = 0.0

        active.remove(i)
        active.remove(j)
        active.append(new_index)
        nodes[new_index] = new_node

    i, j = active
    root = NJTreeNode("Root")
    root.left = nodes[i]
    root.right = nodes[j]
    root.branch_length[nodes[i]] = D[i][j] / 2
    root.branch_length[nodes[j]] = D[i][j] / 2
    return root

def main():
    if len(sys.argv) != 3:
        print("Usage: python nj_debug_final_fixed.py <input.phy> <output.nwk>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    print(f"📄 Reading: {input_file}")
    names, D = read_phylip(input_file)
    print(f"✅ Loaded {len(names)} sequences.")
    root = neighbor_joining(names, D)
    print("🌳 NJ tree constructed.")
    write_newick(root, output_file)
    print(f"✅ Newick tree saved to: {output_file}")

if __name__ == "__main__":
    main()
