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
