import pandas as pd

def calculate_distances(input_file, output_file):
    # 读取Excel文件
    data = pd.read_excel(input_file)
    
    # 根据染色体和起始位置对数据进行排序
    sorted_data = data.sort_values(by=['Chr', 'Start'])
    
    # 计算同一染色体上连续条目之间的距离
    sorted_data['Distance'] = sorted_data.groupby('Chr')['Start'].diff().abs()
    
    # 将带有距离的数据保存到新的Excel文件
    sorted_data.to_excel(output_file, index=False)

# 输入文件路径
input_file_path = 'cp-cp-sp.xlsx'
# 输出文件路径
output_file_path = 'output_distance-sp.xlsx'

# 调用函数
calculate_distances(input_file_path, output_file_path)
