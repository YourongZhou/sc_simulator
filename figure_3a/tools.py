import scanpy as sc
import pandas as pd
import numpy as np
from IPython.display import Image

from curr_method_test import *
import re
from sklearn.metrics import adjusted_mutual_info_score

import xml.etree.ElementTree as ET


def imbalance_entropy(adata, key="cell_type"):
    """
    计算基于 Shannon entropy 的类别不平衡指数 (归一化到 [0,1])
    
    参数:
        adata : AnnData
            输入的 AnnData 对象
        key : str
            obs 中的列名，例如 'cell_type'
    
    返回:
        float : 归一化 entropy (0=极度不平衡, 1=完全平衡)
    """
    counts = adata.obs[key].value_counts()
    print("Cell type counts:", counts)
    probs = counts / counts.sum()

    # Shannon entropy
    H = -(probs * np.log(probs)).sum()

    # 归一化
    H_norm = H #/ np.log(len(counts)) if len(counts) > 1 else 0.0
    return H_norm

def comparison_histogram(dubstepr_genes, anticorr_genes, mule_genes, all_markers, subtype_markers):
    # 所有出现过的基因全集（做混淆矩阵的背景空间）
    universe = set(dubstepr_genes) | set(anticorr_genes) | set(mule_genes) | set(all_markers) | set(subtype_markers)

    def compute_metrics(pred_genes, true_genes, universe):
        pred_set = set(pred_genes)
        true_set = set(true_genes)

        TP = len(pred_set & true_set)
        FP = len(pred_set - true_set)
        FN = len(true_set - pred_set)
        TN = len(universe - pred_set - true_set)

        precision = TP / (TP + FP) if (TP+FP) > 0 else 0
        sensitivity = TP / (TP + FN) if (TP+FN) > 0 else 0
        specificity = TN / (TN + FP) if (TN+FP) > 0 else 0
        accuracy = (TP + TN) / (TP+FP+FN+TN) if (TP+FP+FN+TN) > 0 else 0

        # AMI: 转成 0/1 标签向量
        pred_labels = [1 if g in pred_set else 0 for g in universe]
        true_labels = [1 if g in true_set else 0 for g in universe]
        ami = adjusted_mutual_info_score(true_labels, pred_labels)

        return ami, precision, sensitivity, specificity, accuracy

    # 三个方法 vs all_markers
    scores_all = {
        "DUBStepR": compute_metrics(dubstepr_genes, all_markers, universe),
        "Anti-correlation": compute_metrics(anticorr_genes, all_markers, universe),
        "MULE": compute_metrics(mule_genes, all_markers, universe),
    }

    # 三个方法 vs subtype_markers
    scores_sub = {
        "DUBStepR": compute_metrics(dubstepr_genes, subtype_markers, universe),
        "Anti-correlation": compute_metrics(anticorr_genes, subtype_markers, universe),
        "MULE": compute_metrics(mule_genes, subtype_markers, universe),
    }

    metrics_all = ["AMI", "Prediction", "Sensitivity", "Specificity", "Accuracy"]
    metrics_sub = ["AMI", "Precision", "Sensitivity", "Specificity", "Accuracy"]

    # -------- 绘图 ----------
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharey=True)

    bar_width = 0.18
    x_all = np.arange(len(metrics_all))
    x_sub = np.arange(len(metrics_sub))

    colors = {
        "DUBStepR": "tab:orange",
        "Anti-correlation": "tab:blue",
        "MULE": "tab:green",
        "All markers": "gray",
        "Subtype markers": "gray"
    }

    # all markers subplot
    for i, (method, scores) in enumerate(scores_all.items()):
        axes[0].bar(x_all + i*bar_width, scores, width=bar_width,
                    label=method, color=colors[method], alpha=0.8)

    axes[0].set_xticks(x_all + bar_width*len(scores_all)/2)
    axes[0].set_xticklabels(metrics_all)
    axes[0].set_title("all markers")

    # subtype markers subplot
    for i, (method, scores) in enumerate(scores_sub.items()):
        axes[1].bar(x_sub + i*bar_width, scores, width=bar_width,
                    label=method, color=colors[method], alpha=0.8)

    axes[1].set_xticks(x_sub + bar_width*len(scores_sub)/2)
    axes[1].set_xticklabels(metrics_sub)
    axes[1].set_title("subtype markers")

    axes[0].legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()
    plt.show()

    return scores_all, scores_sub

class TreeNode:
    def __init__(self, name, val=0):
        self.name = name
        self.val = val
        self.mean = 5
        self.sigma = 0.3
        self.children = []
        self.marker_genes = {}
        self.marker_count = 0
        self.background_ratio = 0

    def add_markers(self, count: int, tot_exp: float, 
                mean: float = 0, sigma: float = 0.1,
                high_ratio: float = 0.2, mid_ratio: float = 0.5, low_ratio: float = 0.3,
                high_scale: float = 1.25, low_scale: float = 0.8, background_ratio: float = 0.0):
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

        mu = float(element.get('mu', 10.0))   # 默认值可自定
        sigma = float(element.get('sigma', 0.3))

        high_ratio = float(element.get('high_ratio', 0.2))
        mid_ratio  = float(element.get('mid_ratio', 0.5))
        low_ratio  = float(element.get('low_ratio', 0.3))
        background_ratio = float(element.get('background_ratio', 0.01))

        high_scale = float(element.get('high_scale', 1.25))
        low_scale  = float(element.get('low_scale', 0.8))

        node = TreeNode(name, num_cells)
        node.add_markers(num_genes, tot_exp, 
                        high_ratio=high_ratio, mid_ratio=mid_ratio, low_ratio=low_ratio,
                        high_scale=high_scale, low_scale=low_scale)
        
        node.mu = mu
        node.sigma = sigma
        node.background_ratio = background_ratio

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