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
def get_leaf_nodes_dict(celltype):
    if not celltype["children"]:
        return [celltype]
    leaves = []
    for child in celltype["children"]:
        leaves.extend(get_leaf_nodes_dict(child))
    return leaves

# -------------------------
# 生成 alpha 矩阵
# -------------------------
def generate_cell_gene_matrix(tree_dict):
    """
    生成 alpha 矩阵，每个叶子节点继承祖先 marker。
    支持任意数量 gene_blocks，express_strength 和 background_strength。
    """
    all_genes = []
    gene_info = {}

    # 给 dict 树添加 parent 属性
    def set_parent(node, parent=None):
        node["parent"] = parent
        for child in node["children"]:
            set_parent(child, node)

    set_parent(tree_dict)

    # 收集所有基因
    def collect_all_genes(node):
        for i, gb in enumerate(node.get("gene_blocks", [])):
            for j in range(gb["num_genes"]):
                gene_name = f"g{node['name']}_gb{i}_{j}"
                all_genes.append(gene_name)
                gene_info[gene_name] = {
                    "celltype": node["name"],
                    "express_strength": gb["express_strength"],
                    "background_strength": gb["background_strength"]
                }
        for child in node["children"]:
            collect_all_genes(child)

    collect_all_genes(tree_dict)

    # 获取叶子节点
    leaf_nodes = get_leaf_nodes_dict(tree_dict)
    data = []
    celltype_names = []

    for leaf in leaf_nodes:
        row = []
        inherited_markers = {}
        current = leaf
        while current:
            for i, gb in enumerate(current.get("gene_blocks", [])):
                num_genes = gb["num_genes"]
                gtype = gb.get("type", "normal")  # 默认为 normal
                
                # if gtype == "normal":
                # 随机微调 90%-110%
                strengths = np.linspace(gb["express_strength"]*0.1, gb["express_strength"]*2, gb["num_genes"])
                bg_strengths = np.linspace(gb["background_strength"]*0.1, gb["background_strength"]*2, gb["num_genes"])
                # strengths = gb["express_strength"] * np.random.uniform(0.9, 1.1, num_genes)
                # bg_strengths = gb["background_strength"] * np.random.uniform(0.9, 1.1, num_genes)

                for j in range(num_genes):
                    gene_name = f"g{current['name']}_gb{i}_{j}"
                    inherited_markers[gene_name] = {
                        "express_strength": strengths[j],
                        "background_strength": bg_strengths[j]
                    }

                # elif gtype == "anticorrelation":
                #     half = num_genes // 2
                #     # 前半部分普通 marker
                #     strengths_a = np.linspace(gb["express_strength"]*0.1, gb["express_strength"]*2, gb["num_genes"])
                #     bg_strengths_a = np.linspace(gb["background_strength"]*0.1, gb["background_strength"]*100, gb["num_genes"])
                #     # strengths_a = gb["express_strength"] * np.random.uniform(0.9, 1.1, half)
                #     # bg_strengths_a = gb["background_strength"] * np.random.uniform(0.9, 1.1, half)
                #     # 后半部分为负相关基因
                #     strengths_b = 2*gb["express_strength"] - strengths_a
                #     bg_strengths_b = 200*gb["background_strength"] - bg_strengths_a

                #     for j in range(half):
                #         gene_name_a = f"g{current['name']}_gb{i}_{j}_a"
                #         gene_name_b = f"g{current['name']}_gb{i}_{j}_b"
                #         inherited_markers[gene_name_a] = {
                #             "express_strength": strengths_a[j],
                #             "background_strength": bg_strengths_a[j]
                #         }
                #         inherited_markers[gene_name_b] = {
                #             "express_strength": strengths_b[j],
                #             "background_strength": bg_strengths_b[j]
                #         }
            current = current.get("parent")

        # 填充 alpha_matrix
        for g in all_genes:
            if g in inherited_markers:
                row.append(inherited_markers[g]["express_strength"])
            else:
                row.append(gene_info[g]["background_strength"] * np.random.uniform(0.9, 1.1))
        data.append(row)
        celltype_names.append(f"ctype_{leaf['name']}")

    df = pd.DataFrame(data, index=celltype_names, columns=all_genes)
    return df



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
