# 模拟单细胞基因表达数据生成器

本项目用于基于细胞系谱树（tree）结构，模拟符合生物学特性的单细胞基因表达矩阵。  
支持设定每个细胞类型的总基因表达量分布，以及高/中/低表达基因的比例，从而生成更接近真实生物学情况的数据。  

---

# 安装指南
```
 conda create --name sc_simulator --file requirements.txt
```
注：很多包不是必须的。未来需要改进。 

# 使用指南
## 树结构构建说明

树文件存放在 `./sim_data/trees/`，命名为 `tree_{id}.xml`。  

### 树格式
通过嵌套的 `<celltype>` 节点来表示父子关系，例如：  

```xml
<tree>
  <celltype name="Root" num_cells="100">
    <!-- 根节点的基因块 -->
    <gene_block num_genes="10" express_strength="30" background_strength="0.01" type="normal"/>
    <gene_block num_genes="6"  express_strength="20" background_strength="0.5"  type="anticorrelation"/>

    <!-- 大类 R1 -->
    <celltype name="R1" num_cells="60">
      <gene_block num_genes="8" express_strength="25" background_strength="0.2" type="normal"/>

      <!-- 子类 A1 -->
      <celltype name="A1" num_cells="30" mu="10" sigma="0.3">
        <gene_block num_genes="5" express_strength="20" background_strength="0.2" type="normal"/>
      </celltype>

      <!-- 子类 A2 -->
      <celltype name="A2" num_cells="30" mu="12" sigma="0.4">
        <gene_block num_genes="5" express_strength="18" background_strength="0.3" type="normal"/>
      </celltype>
    </celltype>

    <!-- 独立类 B -->
    <celltype name="B" num_cells="40" mu="11" sigma="0.35">
      <gene_block num_genes="7" express_strength="22" background_strength="0.25" type="normal"/>
    </celltype>
  </celltype>
</tree>
```

### 树结构说明
- tree：整个模拟的根。
- celltype：细胞类型节点，可以包含基因块 (gene_block) 和子细胞类型。
Root 定义housekeeping基因，也就是在所有细胞中都表达的基因。
R1 是一个大类，进一步细分为 A1 和 A2。
B 是独立的一类。

### celltype 节点参数
- name：细胞类型的名称，必须唯一。
- num_cells：该类型下要生成的细胞数。
- mu：表达分布的均值（可选，一般只在叶子节点使用）。
- sigma：表达分布的标准差（可选，一般只在叶子节点使用）。
示例：
```
<celltype name="A1" num_cells="30" mu="10" sigma="0.3">
```
表示子类 A1 包含 30 个细胞，细胞的总 UMI 服从logNormal(10 ,0.3)。

### gene_block 参数
- num_genes：基因数量。
- express_strength：marker 基因在本细胞类型中的表达强度。
- background_strength：marker 基因在非 marker 细胞类型中的背景表达强度。
- type：基因类型：
-- normal：普通 marker 基因，在本类型中高表达，在其他类型中低表达。
-- anticorrelation：负相关基因，按对生成：
在 marker 类型中：
基因 a ~ express_strength
基因 b = 2✖️express_strength - a；
在非 marker 类型中：
基因 a ~ background_strength
基因 b = 2✖️background_strength - a

示例1：
```
<gene_block num_genes="6" express_strength="20" background_strength="0.5" type="normal"/>
```
表示生成 6 个正常基因。

示例2：
```
<gene_block num_genes="6" express_strength="20" background_strength="0.5" type="anticorrelation"/>
```
表示生成 6 个基因，即 3 对负相关基因。



# 许可

本项目遵循 MIT 协议，欢迎自由使用与修改。
