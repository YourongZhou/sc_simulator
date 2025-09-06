# 模拟单细胞基因表达数据生成器

本项目用于基于细胞系谱树（tree）结构，模拟符合生物学特性的单细胞基因表达矩阵。  
支持设定每个细胞类型的总基因表达量分布，以及高/中/低表达基因的比例，从而生成更接近真实生物学情况的数据。  

---

# 安装指南
```
 conda create --name sc_simulator --file requirements.txt
```
注：很多包不是必须的。未来需要改进。

## 功能特点

- **基于树结构的模拟**：输入树文件，自动构建细胞类型的谱系关系。  
- **多层次表达模式**：每个细胞类型可设定：
  - 总表达量 (`tot_exp`)  
  - 高表达基因比例 (`high_prop`)  
  - 中表达基因比例 (`mid_prop`)  
  - 低表达基因比例 (`low_prop`)  
- **随机化基因分布**：在高/中/低表达基因生成后，打乱顺序，避免固定模式。  
- **符合生物学的表达分布**：基因表达使用 Dirichlet-Multinomial + log-normal 分布模拟，常用于拟合真实单细胞 RNA-seq 的表达特征。  


## 树结构构建说明

树文件存放在 `./sim_data/trees/`，命名为 `tree_{id}.xml`。  

### 树格式
通过嵌套的 `<celltype>` 节点来表示父子关系，例如：  

```xml
<tree>
    <celltype name="Root" num_cells="240" num_genes="150" tot_exp="10000" high_prop="0.2" mid_prop="0.5" low_prop="0.3">
        <celltype name="R1" num_cells="150" num_genes="10" tot_exp="1000" high_prop="0.3" mid_prop="0.4" low_prop="0.3">
            <celltype name="A1" num_cells="70" num_genes="10" tot_exp="1000" high_prop="0.25" mid_prop="0.5" low_prop="0.25"/>
            <celltype name="A2" num_cells="80" num_genes="10" tot_exp="1000" high_prop="0.15" mid_prop="0.6" low_prop="0.25"/>
        </celltype>
        <celltype name="R2" num_cells="90" num_genes="10" tot_exp="1000" high_prop="0.1" mid_prop="0.7" low_prop="0.2"/>
    </celltype>
</tree>
```
通过缩进来表示层级关系：
没有缩进的行为根节点，跟节点的基因就是 housekeeping 基因。每多一级缩进（```\t```），代表进入子节点

# 许可

本项目遵循 MIT 协议，欢迎自由使用与修改。
