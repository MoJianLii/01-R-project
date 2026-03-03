import pandas as pd
from pathlib import Path
import numpy as np
import anndata
import time
import matplotlib.pyplot as plt
import os
import pandas as pd
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

download_base = Path('../../data/abc_atlas')
abc_cache = AbcProjectCache.from_cache_dir(download_base)

abc_cache.current_manifest

datasets = ['Zhuang-ABCA-1', 'Zhuang-ABCA-2', 'Zhuang-ABCA-3', 'Zhuang-ABCA-4']
example_section = {'Zhuang-ABCA-1': 'Zhuang-ABCA-1.057',
                   'Zhuang-ABCA-2': 'Zhuang-ABCA-2.037',
                   'Zhuang-ABCA-3': 'Zhuang-ABCA-3.010',
                   'Zhuang-ABCA-4': 'Zhuang-ABCA-4.002'}
base_dir = r"D:\GUOGUO_project\gsmap\R\analysefig"

cell = {}

# 只读 csv（如果目录里不止这 4 个，也不会误读）
csv_files = [
    f for f in os.listdir(base_dir)
    if f.endswith(".csv")
]

for i, f in enumerate(csv_files, start=1):
    # 构造统一 key 名
    key_name = f"Zhuang-ABCA-{i}"

    file_path = os.path.join(base_dir, f)

    # 读取 CSV
    df = pd.read_csv(file_path, dtype={"cell_label": str})
    df.set_index("cell_label", inplace=True)

    # 存入字典
    cell[key_name] = df

    # 可选：分组打印统计信息
    sdf = df.groupby("brain_section_label")
    print(
        key_name, ":",
        "Number of cells =", len(df),
        ", Number of sections =", sdf.ngroups
    )
cluster_details = abc_cache.get_metadata_dataframe(
    directory='WMB-taxonomy',
    file_name='cluster_to_cluster_annotation_membership_pivoted',
    keep_default_na=False
)
cluster_details.set_index('cluster_alias', inplace=True)

cluster_colors = abc_cache.get_metadata_dataframe(
    directory='WMB-taxonomy',
    file_name='cluster_to_cluster_annotation_membership_color',
)
cluster_colors.set_index('cluster_alias', inplace=True)
cell_extended = {}

for d in datasets :
    cell_extended[d] = cell[d].join(cluster_details, on='cluster_alias')
    cell_extended[d] = cell_extended[d].join(cluster_colors, on='cluster_alias')

def subplot_section(ax, xx, yy, cc = None, val = None, cmap = None) :
    
    if cmap is not None :
        ax.scatter(xx, yy, s=0.5, c=val, marker='.', cmap=cmap)
    elif cc is not None :
        ax.scatter(xx, yy, s=0.5, color=cc, marker='.')
    ax.set_ylim(11, 0)
    ax.set_xlim(0, 11)
    ax.axis('equal')
    ax.set_xticks([])
    ax.set_yticks([])
def plot_sections(cell_extended, example_section, cc = None, val = None, fig_width = 10, fig_height = 10, cmap = None) :
    
    fig, ax = plt.subplots(2, 2)
    fig.set_size_inches(fig_width, fig_height)
    
    for i, d in enumerate(cell_extended):
        
        pred = (cell_extended[d]['brain_section_label'] == example_section[d])
        section = cell_extended[d][pred] 
        
        if cmap is not None :
            subplot_section( ax.flat[i], section['x'], section['y'], val=section[val], cmap=cmap)
        elif cc is not None :
            subplot_section( ax.flat[i], section['x'], section['y'], section[cc])
            
        ax.flat[i].set_title(d)
        
    return fig, ax
ccf_dir = r"D:\GUOGUO_project\gsmap\R\analysefig\ccf"

# 列出所有 CSV 文件并排序（保证顺序和 key 一致）
ccf_files = sorted([f for f in os.listdir(ccf_dir) if f.endswith(".csv")])

# 创建字典存储 ccf 坐标
ccf_coordinates = {}

# 遍历 CSV，生成 key 名和读取 DataFrame
for i, f in enumerate(ccf_files, start=1):
    key_name = f"Zhuang-ABCA-{i}"  # 与之前 cell 字典 key 对齐

    file_path = os.path.join(ccf_dir, f)
    df = pd.read_csv(file_path, dtype={"cell_label": str})

    # 设置 cell_label 为 index
    df.set_index("cell_label", inplace=True)

    # 重命名列
    df.rename(
        columns={
            "x": "x_ccf",
            "y": "y_ccf",
            "z": "z_ccf"
        },
        inplace=True
    )

    # 存入字典
    ccf_coordinates[key_name] = df

    # 把 ccf 坐标 join 到 cell_extended
    if key_name in cell_extended:
        cell_extended[key_name] = cell_extended[key_name].join(df, how="inner")
    else:
        print(f"Warning: {key_name} not in cell_extended keys!")

# 打印信息确认
for k, v in ccf_coordinates.items():
    print(k, "number of cells =", len(v))
parcellation_annotation = abc_cache.get_metadata_dataframe(directory="Allen-CCF-2020",
                                                           file_name='parcellation_to_parcellation_term_membership_acronym')
parcellation_annotation.set_index('parcellation_index', inplace=True)
parcellation_annotation.columns = ['parcellation_%s'% x for x in  parcellation_annotation.columns]

parcellation_color = abc_cache.get_metadata_dataframe(directory="Allen-CCF-2020",
                                                      file_name='parcellation_to_parcellation_term_membership_color')
parcellation_color.set_index('parcellation_index', inplace=True)
parcellation_color.columns = ['parcellation_%s'% x for x in  parcellation_color.columns]

for d in datasets :
    cell_extended[d] = cell_extended[d].join(parcellation_annotation, on='parcellation_index')
    cell_extended[d] = cell_extended[d].join(parcellation_color, on='parcellation_index')   
    print(cell_extended[d].columns)
    print(parcellation_annotation.columns)
    print(parcellation_color.columns)


fig, ax = plot_sections(cell_extended, example_section, 'neurotransmitter_color')
res = fig.suptitle('Neurotransmitter Identity', fontsize=14)
plt.show()
fig, ax = plot_sections(cell_extended, example_section, 'subclass_color')
res = fig.suptitle('Cell Type Subclasses', fontsize=14)
plt.show()
fig, ax = plot_sections(cell_extended, example_section, 'parcellation_substructure_color')
res = fig.suptitle('Parcellation - substructure', fontsize=14)
plt.show()