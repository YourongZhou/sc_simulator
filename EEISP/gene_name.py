import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="filename prefix", type=str)
    parser.add_argument("file1", help="cdi_score_file", type=str)
    parser.add_argument("file2", help="eei_score_file", type=str)
    parser.add_argument("file3", help="expression_file (nonzero)", type=str)
    parser.add_argument("thre", help="threshold", type=str)

    args = parser.parse_args()
    
    # 1. 读取基因名称映射表 (file3)
    # 格式通常是: Index \t GeneName \t ...
    try:
        # 必须给列名，否则第一行可能丢失
        df_genes = pd.read_csv(args.file3, header=None, sep="\t", dtype=str)
        # 建立字典: { 'GeneIndex_String': 'GeneName' }
        # 假设第一列是索引，第二列是名字
        index_to_name = dict(zip(df_genes.iloc[:, 0], df_genes.iloc[:, 1]))
    except Exception as e:
        print(f"Error reading expression file: {e}")
        sys.exit(1)

    # 定义一个内部函数来处理 CDI 和 EEI 的逻辑，避免代码重复
    def process_score_file(input_path, output_path, mapping_dict):
        try:
            # 读取分数文件
            df_scores = pd.read_csv(input_path, header=None, sep="\t")
        except pd.errors.EmptyDataError:
            # 如果文件为空，直接创建一个空文件
            open(output_path, 'w').close()
            return 0
            
        written_count = 0
        seen_pairs = set() # 使用集合去重，彻底解决 IndexError
        
        with open(output_path, "w") as fout:
            # 遍历每一行
            for _, row in df_scores.iterrows():
                # 获取原始索引 (转为字符串以匹配字典键)
                idx1 = str(int(row[0])) 
                idx2 = str(int(row[1]))
                score = row[2]
                
                # 查找名字
                if idx1 in mapping_dict and idx2 in mapping_dict:
                    name1 = mapping_dict[idx1]
                    name2 = mapping_dict[idx2]
                    
                    # 排序以确保 (A, B) 和 (B, A) 被视为同一对
                    pair_key = tuple(sorted((name1, name2)))
                    
                    if pair_key not in seen_pairs:
                        # 写入结果: Name1 \t Name2 \t Score
                        fout.write(f"{name1}\t{name2}\t{score}\n")
                        seen_pairs.add(pair_key)
                        written_count += 1
        return written_count

    # 2. 处理 CDI
    data_file_cdi = f"{args.filename}_CDI_convert_data_thre{args.thre}.txt"
    count1 = process_score_file(args.file1, data_file_cdi, index_to_name)

    # 3. 处理 EEI
    data_file_eei = f"{args.filename}_EEI_convert_data_thre{args.thre}.txt"
    count2 = process_score_file(args.file2, data_file_eei, index_to_name)

    print("*************************************************************")
    print("Finish to convert to gene names! (Fixed Version)")
    print(f"The number of CDI gene pairs : {count1}")
    print(f"The number of EEI gene pairs : {count2}")
    print("*************************************************************")

if __name__ == '__main__':
    main()