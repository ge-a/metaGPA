import argparse
import matplotlib.pyplot as plt
import networkx as nx
import os
import numpy as np
from scipy.stats import ttest_1samp, t
from typing import List, Optional
from collections import deque
from Bio import Phylo
from io import StringIO


class TreeNode:
    def __init__(self, name: str = "", parent: Optional['TreeNode'] = None):
        self.name = name
        self.children: List['TreeNode'] = []
        self.parent: Optional['TreeNode'] = parent
        self.enrichment_value = self._parse_enrichment_value()
        self.children_enrichment_ratio: Optional[float] = None
        self.selected_leaves = 0
        self.unselected_leaves = 0

    def _parse_enrichment_value(self) -> int:
        if not self.name:
            return 0
        try:
            return int(self.name.split("_")[-1])
        except (IndexError, ValueError):
            return 0

    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def is_selected(self) -> bool:
        return self.name.startswith("S")

    def is_unselected(self) -> bool:
        return self.name.startswith("unselected")

    def add_child(self, child: 'TreeNode'):
        self.children.append(child)
        child.parent = self


def parse_newick_to_tree(newick_str: str) -> TreeNode:
    """Parses Newick tree into TreeNode structure using BioPython."""
    handle = StringIO(newick_str)
    tree = Phylo.read(handle, "newick")

    def convert(clade, parent=None):
        node = TreeNode(clade.name or "", parent)
        for child_clade in clade.clades:
            child_node = convert(child_clade, node)
            node.add_child(child_node)
        return node
    return convert(tree.clade)

def parse_tree(root: TreeNode):
    """Computes children_enrichment_ratio from leaves up to root."""

    def compute_node(node: TreeNode):
        if node.is_leaf():
            if node.is_selected():
                node.selected_leaves = 1
                node.unselected_leaves = 0
            elif node.is_unselected():
                node.selected_leaves = 0
                node.unselected_leaves = 1
            else:
                raise ValueError(f"Leaf node {node.name} does not start with S/unselected")
            node.children_enrichment_ratio = None
        else:
            selected = 0
            unselected = 0
            for child in node.children:
                selected += child.selected_leaves
                unselected += child.unselected_leaves
            node.selected_leaves = selected
            node.unselected_leaves = unselected
            if unselected > 0:
                node.children_enrichment_ratio = selected / unselected
            else:
                node.children_enrichment_ratio = -1
    # Bottom-up traversal
    def post_order_traversal(node: TreeNode):
        for child in node.children:
            post_order_traversal(child)
        compute_node(node)

    post_order_traversal(root)
    return root

def compute_significance_cutoff(root: TreeNode, alpha: float = 0.05) -> Optional[float]:
    """
    Computes the enrichment ratio cutoff required to be significantly greater than the root ratio at the given alpha.
    Excludes nodes with enrichment ratio <= 0 to avoid diluting signal.
    """
    enrichment_ratios = []

    def collect_ratios(node: TreeNode):
        if node is not root and node.children_enrichment_ratio is not None and node.children_enrichment_ratio > 0:
            enrichment_ratios.append(node.children_enrichment_ratio)
        for child in node.children:
            collect_ratios(child)

    collect_ratios(root)

    if not enrichment_ratios:
        print("No enrichment ratios found above zero.")
        return None

    sample = np.array(enrichment_ratios)
    sample_size = len(sample)
    sample_std = np.std(sample, ddof=1)
    root_ratio = root.children_enrichment_ratio or 0
    df = sample_size - 1
    t_crit = t.ppf(1 - alpha, df)

    cutoff_ratio = root_ratio + t_crit * (sample_std / np.sqrt(sample_size))
    return cutoff_ratio

def find_significant_nodes_above_cutoff(root: TreeNode, cutoff_ratio: float, min_depth: int = 3) -> List[TreeNode]:
    """
    Identifies nodes with enrichment ratios above the given cutoff and no ancestor also above the cutoff.
    """

    significant_nodes = []

    def get_max_child_depth(node: TreeNode) -> int:
        if node.is_leaf():
            return 0
        return 1 + max(get_max_child_depth(child) for child in node.children)

    def collect_significant_nodes(node: TreeNode):
        if (
            node is not root
            and node.children_enrichment_ratio is not None
            and node.children_enrichment_ratio > cutoff_ratio
            and get_max_child_depth(node) >= min_depth
        ):
            significant_nodes.append(node)
        for child in node.children:
            collect_significant_nodes(child)

    collect_significant_nodes(root)

    def has_significant_ancestor(node: TreeNode) -> bool:
        parent = node.parent
        while parent:
            if (
                parent.children_enrichment_ratio is not None
                and parent.children_enrichment_ratio > cutoff_ratio
                and get_max_child_depth(parent) >= min_depth
            ):
                return True
            parent = parent.parent
        return False

    top_significant_nodes = [n for n in significant_nodes if not has_significant_ancestor(n)]
    return top_significant_nodes

def find_leaves_under_high_ratio_nodes(root: TreeNode, threshold: float, min_depth: int = 3) -> List[TreeNode]:
    results = []

    def get_max_child_depth(node: TreeNode) -> int:
        if node.is_leaf():
            return 0
        return 1 + max(get_max_child_depth(child) for child in node.children)

    def collect_leaves(node: TreeNode) -> List[TreeNode]:
        """Collect all leaves under the given node."""
        if node.is_leaf():
            return [node]
        leaves = []
        for child in node.children:
            leaves.extend(collect_leaves(child))
        return leaves

    def dfs(node: TreeNode):
        if not node.is_leaf():
            if (
                node.children_enrichment_ratio is not None
                and node.children_enrichment_ratio > threshold
                and get_max_child_depth(node) >= min_depth
            ):
                leaves = collect_leaves(node)
                counts = {
                    "selected": sum(1 for leaf in leaves if leaf.is_selected()),
                    "unselected": sum(1 for leaf in leaves if leaf.is_unselected())
                }
                results.append((node, leaves, counts))
                return  # Do not recurse deeper
            for child in node.children:
                dfs(child)

    dfs(root)
    return results

def visualize_tree(root, save_path=None):
    G = nx.DiGraph()
    pos = {}
    labels = {}
    colors = []

    x_spacing = 10.0   
    y_spacing = 5

    def assign_positions(node, depth=0, x_counter=[0]):
        node_id = id(node)

        if node.is_leaf():
            x = x_counter[0]
            pos[node_id] = (x * x_spacing, -depth * y_spacing)
            x_counter[0] += 1

            # Leaf color by selection
            if node.is_selected():
                colors.append("blue")
            elif node.is_unselected():
                colors.append("gray")
        else:
            for child in node.children:
                assign_positions(child, depth + 1, x_counter)

            child_positions = [pos[id(child)][0] for child in node.children]
            center_x = sum(child_positions) / len(child_positions)
            pos[node_id] = (center_x, -depth * y_spacing)

            # Internal node color by enrichment
            ratio = node.children_enrichment_ratio
            if ratio == 0:
                colors.append("red")
            elif ratio is None:
                colors.append("gray")
            else:
                colors.append("green")

            labels[node_id] = f"{ratio:.2f}" if ratio is not None and ratio != 0 else ""

        G.add_node(node_id)
        for child in node.children:
            G.add_edge(node_id, id(child))

    assign_positions(root)

    # Plot
    plt.figure(figsize=(20, 12), constrained_layout=True)
    nx.draw(
        G,
        pos,
        with_labels=False,
        node_color=colors,
        node_size=300,
        width=1.0,
        arrows=False
    )
    nx.draw_networkx_labels(G, pos, labels, font_size=8, font_color="black")
    plt.title("Tree Visualization with Enrichment & Leaf Selection", fontsize=14)
    plt.axis('off')

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300)
        print(f"Saved visualization to: {save_path}")
    else:
        plt.show()
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Parse a Newick .tree file and compute enrichment ratios.")
    parser.add_argument("tree_file", type=str, help="Path to the .tree Newick file")
    parser.add_argument("--output", type=str, help="Path to save the output image (e.g. output/tree.png)")

    args = parser.parse_args()

    with open(args.tree_file, 'r') as f:
        newick_str = f.read()

    root = parse_newick_to_tree(newick_str)
    root = parse_tree(root)

    # Example output: show root stats
    print(f"Root selected leaves: {root.selected_leaves}")
    print(f"Root unselected leaves: {root.unselected_leaves}")
    print(f"Root enrichment ratio: {root.children_enrichment_ratio}")

     # Compute enrichment cutoff and significant nodes
    cutoff_ratio = compute_significance_cutoff(root, alpha=0.05)
    if cutoff_ratio is None:
        print("Unable to compute cutoff ratio.")
        return

    significant_nodes = find_significant_nodes_above_cutoff(root, cutoff_ratio)
    print(f"\nSignificance threshold (alpha=0.05): {cutoff_ratio:.4f}")
    print(f"Number of significant nodes: {len(significant_nodes)}")

    print("\nSignificant nodes with selected/unselected leaf counts:")

    significant_subtrees = find_leaves_under_high_ratio_nodes(root, cutoff_ratio)
    for i, (node, leaves, counts) in enumerate(significant_subtrees, start=1):
        print(f"{i}. Node: {node.name or '(unnamed)'}")
        print(f"   Enrichment ratio: {node.children_enrichment_ratio:.2f}")
        print(f"   Selected leaves: {counts['selected']}")
        print(f"   Unselected leaves: {counts['unselected']}")
        print(f"   Total leaves: {len(leaves)}")

    # Optional visualization
    visualize_tree(root, save_path=args.output)


if __name__ == "__main__":
    main()