import re
from typing import Union, List
import numpy as np
import pandas as pd
from collections import deque
from scipy.stats import dirichlet
import matplotlib.pyplot as plt
import seaborn as sns
import xml.etree.ElementTree as ET

class TreeNode:
    def __init__(self, name, val=0):
        self.name = name
        self.val = val
        self.mean = 5
        self.sigma = 0.3
        self.children = []
        self.marker_genes = {}
        self.marker_count = 0

    def add_markers(self, count: int, tot_exp: float, 
                mean: float = 0, sigma: float = 0.1,
                high_ratio: float = 0.2, mid_ratio: float = 0.5, low_ratio: float = 0.3,
                high_scale: float = 1.25, low_scale: float = 0.8):
        """
        添加marker基因，分为高/中/低表达三类
        
        参数：
        count: marker基因数量
        tot_exp: marker总表达量
        high_ratio/mid_ratio/low_ratio: 不同类别的基因比例（和=1）
        high_scale: 高表达基因相对倍数
        low_scale: 低表达基因相对倍数
        """
        self.marker_count = count
        
        # 按比例分配基因数量
        n_high = int(count * high_ratio)
        n_mid  = int(count * mid_ratio)
        n_low  = count - n_high - n_mid
        
        # 生成高中低表达的基因，lognormal是为了取正数服从 alpha 设置
        raw_high = np.random.lognormal(mean=mean, sigma=sigma, size=n_high) * high_scale
        raw_mid = np.random.lognormal(mean=mean, sigma=sigma, size=n_mid)
        raw_low = np.random.lognormal(mean=mean, sigma=sigma, size=n_low) * low_scale
        
        raw_values = np.concatenate([raw_high, raw_mid, raw_low])
        # np.random.shuffle(raw_values) 
        
        # 归一化到 tot_exp
        normalized_values = raw_values / np.sum(raw_values) * tot_exp
        
        # 更新 marker_genes
        for i, expr in enumerate(normalized_values):
            gene_name = f'g{self.name}_{i+1}'
            self.marker_genes[gene_name] = expr


    def get_all_ancestor_genes(self):
        genes = {}
        current = self
        while current:
            genes.update(current.marker_genes)
            current = getattr(current, 'parent', None)  # 使用parent属性
        return genes
    
    def _get_root(self):
        """获取树的根节点"""
        current = self
        while hasattr(current, 'parent'):
            current = current.parent
        return current
    
    @staticmethod
    def _find_node_by_name(root, name):
        """在树中查找指定名称的节点"""
        if root.name == name:
            return root
        for child in root.children:
            result = TreeNode._find_node_by_name(child, name)
            if result:
                return result
        return None

def print_tree(node, level=0):
    """
    打印树结构
    """
    print("  " * level + f"{node.name}: 细胞数{node.val}，基因数{node.marker_count}，基因强度{node.marker_genes}")
    for child in node.children:
        print_tree(child, level + 1)

def get_leaf_nodes(root):
    """获取所有叶子节点"""
    leaves = []
    if not root.children:
        leaves.append(root)
    for child in root.children:
        leaves.extend(get_leaf_nodes(child))
    return leaves

def get_all_nodes_bfs(root):
    """广度优先遍历获取所有节点"""
    nodes = []
    queue = deque([root])
    while queue:
        node = queue.popleft()
        nodes.append(node)
        queue.extend(node.children)
    return nodes

def parse_tree_from_xml(file_path: str) -> TreeNode:
    """
    从XML文件解析生成树，自动进行类型转换
    
    参数：
    file_path: str, XML文件路径
    
    返回：
    TreeNode: 树的根节点
    """
    def convert_type(value: str, type_func):
        """转换字符串值到指定类型"""
        try:
            return type_func(value)
        except (ValueError, TypeError) as e:
            raise ValueError(f"无法转换值 '{value}' 到指定类型: {str(e)}")

    def _parse_node(element) -> TreeNode:
        name = element.get('name')
        num_cells = convert_type(element.get('num_cells'), int)
        num_genes = convert_type(element.get('num_genes'), int)
        tot_exp   = convert_type(element.get('tot_exp'), float)

        mean = float(element.get('mean', 10.0))   # 默认值可自定
        sigma = float(element.get('sigma', 0.3))

        high_ratio = float(element.get('high_ratio', 0.2))
        mid_ratio  = float(element.get('mid_ratio', 0.5))
        low_ratio  = float(element.get('low_ratio', 0.3))

        node = TreeNode(name, num_cells)
        node.add_markers(num_genes, tot_exp, 
                        high_ratio=high_ratio, mid_ratio=mid_ratio, low_ratio=low_ratio)
        
        node.mean = mean
        node.sigma = sigma
        
        for child_element in element.findall('celltype'):
            child_node = _parse_node(child_element)
            child_node.parent = node
            node.children.append(child_node)
            
        return node


    # 解析XML文件
    tree = ET.parse(file_path)
    root_element = tree.getroot()
    
    # 获取第一个celltype元素作为根节点
    first_celltype = root_element.find('celltype')
    if first_celltype is None:
        raise ValueError("XML文件中没有找到celltype节点")
        
    return _parse_node(first_celltype)
    
def generate_cell_gene_matrix(root):
    """生成细胞-基因表达矩阵"""
    # 获取所有叶子节点（细胞类型）
    leaf_nodes = get_leaf_nodes(root)
    
    # 按BFS顺序获取所有节点的所有基因
    all_genes = []
    for node in get_all_nodes_bfs(root):
        all_genes.extend([f'g{node.name}_{i+1}' for i in range(node.marker_count)])
    
    # 创建矩阵
    alpha_matrix = np.zeros((len(leaf_nodes), len(all_genes)))
    
    # 填充矩阵
    for i, cell in enumerate(leaf_nodes):
        # 获取该细胞类型表达的所有基因
        expressed_genes = cell.get_all_ancestor_genes()
        # 填充表达值
        for j, gene in enumerate(all_genes):
            if gene in expressed_genes:
                alpha_matrix[i, j] = expressed_genes[gene]
    
    # 创建DataFrame
    df = pd.DataFrame(
        alpha_matrix,
        index=[f'ctype_{cell.name}' for cell in leaf_nodes],
        columns=all_genes
    )
    
    return df

def generate_cell_expressions(alpha_matrix, cell_counts):
    """
    根据alpha矩阵和细胞数量生成单细胞基因表达谱
    
    参数：
    alpha_matrix: pandas DataFrame, 行是细胞类型，列是基因
    cell_counts: dict, 键是细胞类型名称，值是该类型的细胞数量
    
    返回：
    cell_expressions: pandas DataFrame, 每行是一个细胞的基因表达谱
    cell_types: list, 每个细胞对应的细胞类型
    """
    all_expressions = []  # 存储所有细胞的表达谱
    all_cell_types = []   # 存储每个细胞的类型标签
    
    # 对每个细胞类型进行处理
    for cell_type in alpha_matrix.index:
        # 获取该细胞类型的alpha参数
        alpha_params = alpha_matrix.loc[cell_type].values
        # 确保alpha参数都是正数
        alpha_params = np.maximum(alpha_params, 1e-10)
        
        # 获取该类型需要生成的细胞数量
        n_cells = cell_counts[cell_type]
        
        # 使用狄利克雷分布生成n_cells个表达谱
        expressions = dirichlet.rvs(alpha_params, size=n_cells)
        
        # 将生成的表达谱添加到结果中
        all_expressions.append(expressions)
        all_cell_types.extend([cell_type] * n_cells)
    
    # 将所有表达谱合并
    all_expressions = np.vstack(all_expressions)
    
    # 创建DataFrame
    result_df = pd.DataFrame(
        all_expressions,
        columns=alpha_matrix.columns,
        index=[f"{ctype}_cell{i+1}" for ctype, i in zip(all_cell_types, range(len(all_cell_types)))]
    )
    
    return result_df, all_cell_types

def generate_final_counts(expression_matrix, total_expressions):
    """
    根据表达比例矩阵和总UMI数生成最终的表达量矩阵
    
    参数：
    expression_matrix: pandas DataFrame, 行是细胞，列是基因的表达比例矩阵
    total_expressions: numpy array, 每个细胞的总UMI数
    
    返回：
    pandas DataFrame: 最终的表达量矩阵，保持与输入相同的索引和列名
    """
    if len(total_expressions) != expression_matrix.shape[0]:
        raise ValueError("总表达量数组长度与细胞数不匹配")
    
    # 初始化结果矩阵
    final_counts = np.zeros_like(expression_matrix, dtype=np.int32)
    
    # 对每个细胞生成表达量
    for i in range(len(total_expressions)):
        # 获取当前细胞的表达比例
        probabilities = expression_matrix.iloc[i].values
        
        # 确保概率和为1
        probabilities = probabilities / probabilities.sum()
        
        # 使用多项式分布生成具体表达量
        counts = np.random.multinomial(n=int(total_expressions[i]), 
                                     pvals=probabilities)
        
        final_counts[i] = counts
    
    # 转换为DataFrame，保持原有的索引和列名
    result_df = pd.DataFrame(final_counts,
                           index=expression_matrix.index,
                           columns=expression_matrix.columns)
    
    return result_df