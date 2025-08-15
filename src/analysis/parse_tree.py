import argparse
import matplotlib.pyplot as plt
import networkx as nx
import os
import numpy as np
from scipy.stats import ttest_1samp
from typing import List, Optional
from Bio import Phylo
from io import StringIO

# Create internal tree data structure for parsing Newick format trees
class TreeNode:
    def __init__(self, name: str = "", parent: Optional['TreeNode'] = None):
        self.name = name
        self.children: List['TreeNode'] = []
        self.parent: Optional['TreeNode'] = parent
        self.children_enrichment_perc: Optional[float] = None
        self.selected_leaves = 0
        self.unselected_leaves = 0
        self.p_value: Optional[float] = None

    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def is_selected(self) -> bool:
        return self.name.startswith("S")

    def is_unselected(self) -> bool:
        return self.name.startswith("unselected")

    def add_child(self, child: 'TreeNode'):
        self.children.append(child)
        child.parent = self

# Analyzer class to handle the parsing and analysis of TreeNode structures
class TreeAnalyzer:
    def __init__(self, newick_str: str):
        self.root = self.parse_newick_to_tree(newick_str)
        self.population_mean = 0.0
        self.ranked_nodes: List[TreeNode] = []

    # Parse Newick format string into a TreeNode structure
    def parse_newick_to_tree(self, newick_str: str) -> TreeNode:
        handle = StringIO(newick_str)
        tree = Phylo.read(handle, "newick")

        def convert(clade, parent=None):
            node = TreeNode(clade.name or "", parent)
            for child_clade in clade.clades:
                child_node = convert(child_clade, node)
                node.add_child(child_node)
            return node

        return convert(tree.clade)

    # Compute enrichment percentages for each node in the tree
    def compute_enrichment(self):
        def compute_node(node: TreeNode):
            if node.is_leaf():
                if node.is_selected():
                    node.selected_leaves = 1
                    node.unselected_leaves = 0
                elif node.is_unselected():
                    node.selected_leaves = 0
                    node.unselected_leaves = 1
                else:
                    raise ValueError(f"Leaf node {node.name} must start with S or unselected")
                node.children_enrichment_perc = None
            else:
                selected = 0
                unselected = 0
                for child in node.children:
                    selected += child.selected_leaves
                    unselected += child.unselected_leaves
                node.selected_leaves = selected
                node.unselected_leaves = unselected
                total = selected + unselected
                node.children_enrichment_perc = selected / total if total > 0 else 0

        def post_order_traversal(node: TreeNode):
            for child in node.children:
                post_order_traversal(child)
            compute_node(node)

        post_order_traversal(self.root)
        self.population_mean = self.root.children_enrichment_perc or 0.0

    # Assign p-values to each node based on selected and unselected leaves
    def assign_p_values(self):
        self.ranked_nodes.clear()

        def assign(node: TreeNode):
            total = node.selected_leaves + node.unselected_leaves
            if total > 1:
                sample = [1] * node.selected_leaves + [0] * node.unselected_leaves
                t_stat, p_two_tailed = ttest_1samp(sample, self.population_mean)
                node.p_value = p_two_tailed / 2 if t_stat > 0 else 1.0
            else:
                node.p_value = None

            if not node.is_leaf() and node.p_value is not None:
                self.ranked_nodes.append(node)

            for child in node.children:
                assign(child)

        assign(self.root)
        self.ranked_nodes.sort(key=lambda n: n.p_value)
    
    # Writes TT files with internal node statistics
    def write_tt(self, cutoff: float, output_path: str):
        total_selected = self.root.selected_leaves
        total_unselected = self.root.unselected_leaves

        rows = []

        def process_node(node: TreeNode):
            if node.is_leaf():
                return

            count_above = node.selected_leaves
            count_below = node.unselected_leaves

            frac_above = (count_above / total_selected) * 100 if total_selected else 0
            frac_below = (count_below / total_unselected) * 100 if total_unselected else 0

            node_id = node.name
            pval = node.p_value
            pval_str = f"{pval:.4e}" if pval is not None else "NA"

            # Collect data in a list (store pval as float or None for sorting)
            rows.append({
                "node_id": node_id,
                "pval": pval,
                "pval_str": pval_str,
                "frac_above": frac_above,
                "count_above": count_above,
                "frac_below": frac_below,
                "count_below": count_below,
            })

            for child in node.children:
                process_node(child)

        process_node(self.root)

        # Sort rows by pval, placing 'NA' (None) at the end
        rows.sort(key=lambda r: (float('inf') if r['pval'] is None else r['pval']))

        with open(output_path, "w") as out:
            out.write("Node_ID\tP_Value\t%AboveCutoff\tCountAbove\t%BelowCutoff\tCountBelow\n")
            for r in rows:
                out.write(f"{r['node_id']}\t{r['pval_str']}\t{r['frac_above']:.2f}\t{r['count_above']}\t{r['frac_below']:.2f}\t{r['count_below']}\n")

    # Print top nodes based on p-value and selected leaves
    def print_top_nodes(self, top_n=10, min_selected_frac=0.7, output_file=None) -> Optional[float]:
        total_selected = self.root.selected_leaves
        count = 0
        lowest_p = None

        def write(msg):
            if output_file:
                print(msg, file=output_file)

        write(f"\nTop {top_n} internal nodes by p-value:")
        if total_selected == 0:
            write("No selected leaves in the tree.")
            return None

        for node in self.ranked_nodes:
            node_selected_frac = node.selected_leaves / total_selected
            if node_selected_frac >= min_selected_frac:
                count += 1
                write(f"{count}. Node Name: {node.name} | p: {node.p_value:.2e} | enrichment: "
                      f"{node.children_enrichment_perc:.2f} | selected: {node.selected_leaves} | "
                      f"unselected: {node.unselected_leaves}")
                if lowest_p is None:
                    lowest_p = node.p_value
                if count >= top_n:
                    break

        return lowest_p

    def visualize_tree(self, save_path=None):
        G = nx.DiGraph()
        pos = {}
        labels = {}
        colors = []

        x_spacing = 15  
        y_spacing = 5

        def assign_positions(node, depth=0, x_counter=[0]):
            node_id = id(node)

            if node.is_leaf():
                x = x_counter[0]
                pos[node_id] = (x * x_spacing, -depth * y_spacing)
                x_counter[0] += 1

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

                perc = node.children_enrichment_perc
                if perc == 0:
                    colors.append("red")
                elif perc is None:
                    colors.append("gray")
                else:
                    colors.append("green")

                label = ""
                if perc is not None and perc != 0:
                    label += f"{perc:.2f}"
                if node.p_value is not None:
                    p_val_str = f"{node.p_value:.2e}"
                    if node.p_value < 0.05:
                        p_val_str += " *"
                        label += f"\n(p={p_val_str})"
                labels[node_id] = label

            G.add_node(node_id)
            for child in node.children:
                G.add_edge(node_id, id(child))

        assign_positions(self.root)

        plt.figure(figsize=(25, 15), constrained_layout=True)
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

# We only output a tree image if an output path is provided
def main():
    parser = argparse.ArgumentParser(description="Parse a Newick .tree file and compute enrichment ratios.")
    parser.add_argument("tree_file", type=str, help="Path to the .tree Newick file")
    parser.add_argument("--output", type=str, help="Path to save the output image (e.g. output/tree.png)")
    parser.add_argument("--write-tt", type=str, help="Write internal node stats to a .txt file")

    args = parser.parse_args()

    with open(args.tree_file, 'r') as f:
        newick_str = f.read()

    analyzer = TreeAnalyzer(newick_str)
    analyzer.compute_enrichment()
    analyzer.assign_p_values()

    if args.write_tt is not None:
        split_tt = args.write_tt.split(":")
        domain = split_tt[0]
        cutoff = split_tt[1]
        output_dir_path = split_tt[2]
        tt_output_path = output_dir_path + "/TT/" + domain + "_TT.txt"
        analyzer.write_tt(cutoff=cutoff, output_path=tt_output_path)
        print(f"Wrote internal node stats to: {tt_output_path}")

    if args.output:
        log_path = os.path.splitext(args.output)[0] + ".txt"
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        log_file = open(log_path, "w")
    else:
        log_file = None

    def log(msg):
        if log_file:
            print(msg, file=log_file)
    
    log(f"Root selected leaves: {analyzer.root.selected_leaves}")
    log(f"Root unselected leaves: {analyzer.root.unselected_leaves}")
    log(f"Root enrichment perc: {analyzer.root.children_enrichment_perc:.4f}")

    # Only prints the top nodes if an output file is provided
    lowest_p = analyzer.print_top_nodes(top_n=10, output_file=log_file)

    if args.output:
        analyzer.visualize_tree(save_path=args.output)
        log_file.close()

    if lowest_p is not None:
        print(f"{lowest_p:.6e}")
    else:
        print("None")


if __name__ == "__main__":
    main()