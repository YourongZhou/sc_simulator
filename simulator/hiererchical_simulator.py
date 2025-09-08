import numpy as np
import pandas as pd
from scipy.stats import dirichlet
import xml.etree.ElementTree as ET
import re
import seaborn as sns

# -------------------------
# 解析树函数（dict 版本）
# -------------------------
def parse_tree(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    return parse_celltype(root.find("celltype"))

def parse_celltype(element, parent=None):
    cell = {
        "name": element.attrib["name"],
        "num_cells": int(element.attrib["num_cells"]),
        "gene_blocks": [],
        "parent": parent
    }
    
    # 解析 gene_block
    for gb in element.findall("gene_block"):
        block = {
            "num_genes": int(gb.attrib["num_genes"]),
            "express_strength": float(gb.attrib["express_strength"]),
            "background_strength": float(gb.attrib["background_strength"]),
            "type": gb.attrib.get("type", "normal")
        }
        cell["gene_blocks"].append(block)
    
    # 解析子 celltype
    cell["children"] = [parse_celltype(child, parent=cell) 
                        for child in element.findall("celltype")]
    return cell

# -------------------------
# 获取叶子节点（dict 版本）
# -------------------------
def get_leaf_nodes_dict(tree_dict):
    """递归获取所有叶子节点"""
    if not tree_dict.get("children"):
        return [tree_dict]
    leaves = []
    for child in tree_dict["children"]:
        leaves.extend(get_leaf_nodes_dict(child))
    return leaves


def generate_alpha_matrix(tree_dict):
    """
    生成 alpha_matrix，每一行对应一个 cell（cell_00001..），每列对应一个基因。
    - normal gene block: 在定义该 block 的节点的 descendant 叶子里使用 express_strength（cell-level 浮动），
      在其它叶子里使用 background_strength（cell-level 浮动）。
    - anticorrelation gene block: 将 num_genes 分为 pairs (half); 对每对 (a,b)：
        对每个叶子 L：
            if L is descendant(node): base = express_strength
            else: base = background_strength
            a_vals = linspace(base*0.9, base*1.1, n_cells_of_L)
            b_vals = 2*base - a_vals
        将 a_vals, b_vals 填入该叶子对应的那些行
    返回:
        alpha_df: pd.DataFrame, shape (total_cells, total_genes), index cell_00001...
    注意:
        如果 anticorrelation 的 num_genes 为奇数，最后一个基因将被作为 normal gene 处理（不会丢失）。
    """
    # --- helper: set parent pointers (if not already) ---
    def set_parent(node, parent=None):
        node["parent"] = parent
        for c in node.get("children", []):
            set_parent(c, node)
    set_parent(tree_dict)

    # --- 收集 nodes 顺序（preorder）和所有叶子 ---
    nodes = []
    def traverse(node):
        nodes.append(node)
        for c in node.get("children", []):
            traverse(c)
    traverse(tree_dict)

    leaf_nodes = get_leaf_nodes_dict(tree_dict)  # 保持你的已有函数

    # --- build all_genes list and gene->col map ---
    all_genes = []
    gene_info = {}  # gene_name -> metadata
    for node in nodes:
        for bi, gb in enumerate(node.get("gene_blocks", [])):
            gtype = gb.get("type", "normal")
            if gtype == "normal":
                for j in range(gb["num_genes"]):
                    gname = f"g{node['name']}_gb{bi}_{j}"
                    all_genes.append(gname)
                    gene_info[gname] = {
                        "node": node,
                        "type": "normal",
                        "express_strength": gb["express_strength"],
                        "background_strength": gb["background_strength"],
                        "block_idx": bi,
                        "gene_idx": j
                    }
            elif gtype == "anticorrelation":
                b_list = []
                half = gb["num_genes"] // 2
                # 如果是奇数，最后一个基因当 normal 处理
                for j in range(half):
                    ga = f"g{node['name']}_gb{bi}_{j}_a"
                    gbn = f"g{node['name']}_gb{bi}_{j}_b"
                    all_genes.append(ga)
                    # all_genes.append(gbn)
                    b_list.append(gbn)
                    gene_info[ga] = {
                        "node": node,
                        "type": "anticorrelation_a",
                        "express_strength": gb["express_strength"],
                        "background_strength": gb["background_strength"],
                        "block_idx": bi,
                        "pair_idx": j
                    }
                    gene_info[gbn] = {
                        "node": node,
                        "type": "anticorrelation_b",
                        "express_strength": gb["express_strength"],
                        "background_strength": gb["background_strength"],
                        "block_idx": bi,
                        "pair_idx": j
                    }
                all_genes = all_genes + b_list
                if gb["num_genes"] % 2 == 1:
                    # 额外的最后一个基因当 normal
                    j = half
                    gname = f"g{node['name']}_gb{bi}_{j}"
                    all_genes.append(gname)
                    gene_info[gname] = {
                        "node": node,
                        "type": "normal",
                        "express_strength": gb["express_strength"],
                        "background_strength": gb["background_strength"],
                        "block_idx": bi,
                        "gene_idx": j
                    }

    # col map
    col_index = {g:i for i,g in enumerate(all_genes)}

    # --- total rows & leaf -> row range mapping ---
    leaf_start = {}
    total_cells = 0
    for leaf in leaf_nodes:
        leaf_start[leaf["name"]] = total_cells
        total_cells += leaf["num_cells"]

    # prepare alpha matrix
    n_genes = len(all_genes)
    alpha = np.zeros((total_cells, n_genes), dtype=float)

    # helper to test ancestor relation (node is ancestor of leaf?)
    def is_ancestor(node, leaf):
        cur = leaf
        while cur is not None:
            if cur is node:
                return True
            cur = cur.get("parent")
        return False

    # --- Fill alpha by scanning nodes and their gene_blocks ---
    for node in nodes:
        for bi, gb in enumerate(node.get("gene_blocks", [])):
            gtype = gb.get("type", "normal")
            if gtype == "normal":
                # for each gene in block, fill each leaf's rows
                for j in range(gb["num_genes"]):
                    gname = f"g{node['name']}_gb{bi}_{j}"
                    col = col_index[gname]
                    for leaf in leaf_nodes:
                        start = leaf_start[leaf["name"]]
                        n = leaf["num_cells"]
                        if is_ancestor(node, leaf):
                            base = gb["express_strength"]
                        else:
                            base = gb["background_strength"]
                        vals = np.linspace(base * 0.9, base * 1.1, n)
                        alpha[start:start+n, col] = vals
            elif gtype == "anticorrelation":
                half = gb["num_genes"] // 2
                for j in range(half):
                    ga = f"g{node['name']}_gb{bi}_{j}_a"
                    gbn = f"g{node['name']}_gb{bi}_{j}_b"
                    col_a = col_index[ga]
                    col_b = col_index[gbn]
                    for leaf in leaf_nodes:
                        start = leaf_start[leaf["name"]]
                        n = leaf["num_cells"]
                        if is_ancestor(node, leaf):
                            base = gb["express_strength"]
                        else:
                            base = gb["background_strength"]
                        a_vals = np.linspace(base * 0.1, base * 2, n)
                        b_vals = 2 * base - a_vals
                        alpha[start:start+n, col_a] = a_vals
                        alpha[start:start+n, col_b] = b_vals
                # odd leftover -> treat as normal (if exists)
                if gb["num_genes"] % 2 == 1:
                    j = half
                    gname = f"g{node['name']}_gb{bi}_{j}"
                    col = col_index[gname]
                    for leaf in leaf_nodes:
                        start = leaf_start[leaf["name"]]
                        n = leaf["num_cells"]
                        if is_ancestor(node, leaf):
                            base = gb["express_strength"]
                        else:
                            base = gb["background_strength"]
                        vals = np.linspace(base * 0.1, base * 2, n)
                        alpha[start:start+n, col] = vals

    # --- cell ids in same order (leaf_nodes order) ---
    cell_ids = []
    gid = 1
    for leaf in leaf_nodes:
        for k in range(leaf["num_cells"]):
            cell_ids.append(f"ctype_{leaf['name']}_cell{gid+1:03d}")
            gid += 1

    alpha_df = pd.DataFrame(alpha, index=cell_ids, columns=all_genes)
    return alpha_df



def sample_cell_expressions(alpha_matrix):
    """
    根据每个 cell 的 alpha 向量逐行采样 Dirichlet，返回比例矩阵（cell x gene）
    """
    all_expressions = []
    for i in range(alpha_matrix.shape[0]):
        alpha_vec = np.maximum(alpha_matrix.iloc[i].values, 1e-10)
        expr = dirichlet.rvs(alpha_vec, size=1)[0]
        all_expressions.append(expr)

    result_df = pd.DataFrame(all_expressions,
                             columns=alpha_matrix.columns,
                             index=alpha_matrix.index)
    return result_df




# -------------------------
# 生成单细胞表达谱
# -------------------------
def generate_cell_expressions(alpha_matrix, cell_counts):
    all_expressions = []
    all_cell_types = []

    for cell_type in alpha_matrix.index:
        alpha_params = alpha_matrix.loc[cell_type].values
        alpha_params = np.maximum(alpha_params, 1e-10)
        n_cells = cell_counts[cell_type]
        expressions = dirichlet.rvs(alpha_params, size=n_cells)
        all_expressions.append(expressions)
        all_cell_types.extend([cell_type] * n_cells)

    all_expressions = np.vstack(all_expressions)
    result_df = pd.DataFrame(
        all_expressions,
        columns=alpha_matrix.columns,
        index=[f"{ctype}_cell{i+1}" for ctype, i in zip(all_cell_types, range(len(all_cell_types)))]
    )
    return result_df, all_cell_types

# -------------------------
# 根据总 UMI 生成整数计数矩阵
# -------------------------
def generate_final_counts(expression_matrix, total_expressions):
    final_counts = np.zeros_like(expression_matrix, dtype=np.int32)
    for i in range(len(total_expressions)):
        probabilities = expression_matrix.iloc[i].values
        probabilities = probabilities / probabilities.sum()
        counts = np.random.multinomial(n=int(total_expressions[i]), pvals=probabilities)
        final_counts[i] = counts
    df = pd.DataFrame(final_counts, index=expression_matrix.index, columns=expression_matrix.columns)
    return df
