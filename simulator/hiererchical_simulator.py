import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import dirichlet
import xml.etree.ElementTree as ET
import re
import seaborn as sns
import matplotlib.pyplot as plt

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


def generate_data(xml_file='tree.xml', mu=10, sigma=0.3):
    print("Parsing xml tree file...")
    tree_dict = parse_tree(xml_file)

    print("Generating alpha matrix...")
    # alpha_matrix = generate_cell_gene_matrix(tree_dict) # 针对基因横向高低
    alpha_matrix = generate_alpha_matrix(tree_dict) # 针对纵向高低


    print("Sampling p's from Dirichlet distribution...")
    leaf_nodes = get_leaf_nodes_dict(tree_dict)
    cell_counts = {f'ctype_{leaf["name"]}': leaf["num_cells"] for leaf in leaf_nodes}

    # expression_matrix, cell_types = generate_cell_expressions(alpha_matrix, cell_counts) # 针对基因横向高低
    expression_matrix = sample_cell_expressions(alpha_matrix) # 针对纵向高低

    print("Sampling cell total UMIs from logNormal distribution...")
    total_expressions = []
    for leaf in leaf_nodes:
        n_cells = leaf["num_cells"]
        mean = leaf.get("mean", mu)
        sigma_val = leaf.get("sigma", sigma)
        leaf_expressions = np.random.lognormal(mean=mean, sigma=sigma_val, size=n_cells)
        total_expressions.extend(leaf_expressions)
    total_expressions = np.array(total_expressions)

    print("Sampling cells from multinomial distribution...")
    final_matrix = generate_final_counts(expression_matrix, total_expressions)
    return pd.DataFrame(final_matrix)

# def generate_data(xml_file='tree_10.xml', mu=10, sigma=0.3):
#     """
#     使用 Dirichlet-Multinomial 分布生成模拟单细胞表达数据。
    
#     流程说明：
#     1. 从 XML 文件解析细胞类型树。
#     2. 根据树结构生成 Dirichlet 分布的 alpha 矩阵。
#     3. 从 Dirichlet 分布中采样每个细胞类型的比例 (p)。
#     4. 根据采样的比例生成单细胞表达谱。
#     5. 从对数正态分布中采样每个细胞的总 UMI 数量。
#     6. 使用 Multinomial 分布生成最终的整数计数矩阵。
#     7. 将模拟矩阵保存为 CSV 文件。

#     参数
#     ----
#     xml_file : str, 可选
#         定义细胞类型树的 XML 文件路径（默认 'tree_10.xml'）。
#     mu : float, 可选
#         对数正态分布的均值，用于模拟每个细胞的总 UMI（默认 10）。
#     sigma : float, 可选
#         对数正态分布的标准差，用于模拟每个细胞的总 UMI（默认 0.3）。

#     返回
#     ----
#     pd.DataFrame
#         模拟生成的单细胞表达矩阵（行为细胞，列为基因）。
    
#     注意事项
#     ----
#     - 每个细胞类型的细胞数量由树的叶子节点确定。
#     - alpha 矩阵控制每个细胞类型的 Dirichlet 分布，影响基因表达比例。
#     - 每个细胞的总 UMI 独立从 log-normal 分布采样。
#     - Multinomial 分布生成整数计数，保证总 UMI 与采样一致。
#     - 输出文件保存路径为 './sim_data/csvs/'，文件名根据 XML 文件名生成。
#     """
    
#     # 根据xml生成树
#     print("Parsing xml tree file...")
#     root = parse_tree_from_xml(xml_file)
#     # 生成alpha矩阵
#     print("Generating alpha matrix...")
#     alpha_matrix = generate_cell_gene_matrix(root)
    
#     # 生成 Dirichlet的 p
#     # 设定每个细胞类型的细胞数量（从树的叶子节点获取）
#     print("Sampling p's from Dirichlet distribution...")
#     cell_counts = {f'ctype_{leaf.name}': leaf.val for leaf in get_leaf_nodes(root)}
    
#     # 生成单细胞表达谱
#     expression_matrix, cell_types = generate_cell_expressions(alpha_matrix, cell_counts)
    
#     # 生成logNormal总 UMI
#     print("Sampling cell total UMIs from logNormal distribution...")
#     total_expressions = []
#     for leaf in get_leaf_nodes(root):
#         n_cells = leaf.val
#         mean = getattr(leaf, "mean", mu)
#         sigma = getattr(leaf, "sigma", sigma)
        
#         leaf_expressions = np.random.lognormal(mean=mean, sigma=sigma, size=n_cells)
#         total_expressions.extend(leaf_expressions)

#     total_expressions = np.array(total_expressions)

    
#     # Multinomial生成最终矩阵
#     print("Sampling cells from multinomial distribution...")
#     final_matrix = generate_final_counts(expression_matrix, total_expressions)
    
#     return pd.DataFrame(final_matrix)


def get_data_and_visualize(tree_id, mu=10, sigma=0.3, 
                           tree_dir='../sim_data/trees', 
                           data_dir='../sim_data/adata'):
    """
    生成模拟单细胞表达数据，并绘制基因表达热图，同时保存 AnnData 对象。
    
    流程说明：
    1. 根据指定的 tree_id 从 XML 文件生成模拟单细胞数据（调用 generate_data）。
    2. 提取基因类别（gene category），用于热图列颜色映射。
    3. 提取细胞类型（cell type），用于热图行颜色映射。
    4. 为基因类别和细胞类型分配调色板。
    5. 绘制不聚类的单细胞表达热图（行为细胞，列为基因）。
    6. 在热图右侧添加细胞类型和基因类别图例。
    7. 将表达矩阵封装为 AnnData 对象，并保存为 h5ad 文件。
    
    参数
    ----
    tree_id : int
        树的编号，对应 XML 文件 tree_{tree_id}.xml。
    mu : float, 可选
        log-normal 分布的均值，用于模拟每个细胞总 UMI（默认 10）。
    sigma : float, 可选
        log-normal 分布的标准差，用于模拟每个细胞总 UMI（默认 0.3）。
    tree_dir : str, 可选
        存放 XML 树文件的目录（默认 './sim_data/trees'）。
    data_dir : str, 可选
        存放生成的 AnnData 对象的目录（默认 './sim_data/adata'）。
    
    返回
    ----
    data_df : pd.DataFrame
        模拟生成的单细胞表达矩阵（行为细胞，列为基因）。
    adata : AnnData
        对应的 AnnData 对象，包含表达矩阵和细胞类型信息 obs。
    
    注意事项
    ----
    - 生成的细胞名称格式为 "ctype_{cell_type}_cell{number}"。
    - 生成的基因名称以 "g{类别}_" 开头，如 "gA1_gene1"。
    - 热图行和列不进行聚类，仅用于可视化表达模式。
    """
    
    # 1. 生成模拟数据
    tree_file = f"{tree_dir}/tree_{tree_id}.xml"
    data_df = generate_data(tree_file, mu, sigma)

    # 2. 提取基因类别
    def extract_gene_category(gene_name):
        match = gene_name.split("_")
        if match:
            return match[0]
        return "Unknown"

    gene_categories = [extract_gene_category(gene) for gene in data_df.columns]
    gene_category_series = pd.Series(gene_categories, index=data_df.columns)

    # 3. 提取细胞类型
    def extract_cell_type(cell_name):
        match = re.match(r'ctype_(\w+)_cell\d+', cell_name)
        if match:
            return match.group(1)
        return None
    
    cell_types = [extract_cell_type(cell) for cell in data_df.index]
    if None in cell_types:
        raise ValueError("无法正确提取所有细胞类型，请检查索引格式")

    cell_type_series = pd.Series(cell_types, index=data_df.index)

    # 4. 分配颜色
    # 细胞类型颜色映射
    unique_cell_types = cell_type_series.unique()
    palette_cells = sns.color_palette("tab20", len(unique_cell_types))
    cell_type_colors = dict(zip(unique_cell_types, palette_cells))
    row_colors = cell_type_series.map(cell_type_colors)

    # 基因类别颜色映射
    unique_gene_categories = gene_category_series.unique()
    palette_genes = sns.color_palette("Set2", len(unique_gene_categories))
    gene_category_colors = dict(zip(unique_gene_categories, palette_genes))
    col_colors = gene_category_series.map(gene_category_colors)

    # 5. 绘制单细胞表达热图（不聚类）
    g = sns.clustermap(
        np.log1p(data_df),
        row_cluster=False,
        col_cluster=False,
        row_colors=row_colors,
        col_colors=col_colors,
        xticklabels=False,
        yticklabels=False,
        cbar_kws={'label': 'Expression'},
        figsize=(12, 8)
    )

    g.ax_heatmap.set_xlabel('Genes')
    g.ax_heatmap.set_ylabel('Cells')
    g.ax_heatmap.set_title("Gene Expression Heatmap")

    # 6. 添加图例
    # cell type legend
    handles_cell = [plt.Line2D([0], [0], marker='s', color=color, linestyle='') 
                    for color in cell_type_colors.values()]
    labels_cell = list(cell_type_colors.keys())

    # gene category legend
    handles_gene = [plt.Line2D([0], [0], marker='s', color=color, linestyle='') 
                    for color in gene_category_colors.values()]
    labels_gene = list(gene_category_colors.keys())

    # 在热图右侧添加两个 legend
    legend1 = g.ax_heatmap.legend(handles_cell, labels_cell, 
                                title="Cell Types", 
                                loc="upper left", bbox_to_anchor=(1.03, 1), 
                                borderaxespad=0)
    g.ax_heatmap.add_artist(legend1)  # 先画 cell types legend
    g.ax_heatmap.legend(handles_gene, labels_gene, 
                        title="Gene Categories", 
                        loc="lower left", bbox_to_anchor=(1.03, 0), 
                        borderaxespad=0)

    plt.show()

    # 7. 构建 AnnData 对象
    adata = sc.AnnData(X=data_df.values, 
                    obs=pd.DataFrame(index=data_df.index), 
                    var=pd.DataFrame(index=data_df.columns))

    # 添加细胞类型信息到 obs
    adata.obs['cell_type'] = [extract_cell_type(cell) for cell in adata.obs.index]

    # 保存 AnnData 对象
    print("Saving anndata to:", f"{data_dir}/adata_{tree_id}.h5ad")
    adata.write_h5ad(f"{data_dir}/adata_{tree_id}.h5ad")

    return data_df, adata