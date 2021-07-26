
#conda create -n cpdb python=3.7
source activate cpdb
#pip install cellphonedb
cellphonedb method statistical_analysis test_meta.txt test_counts.txt --threads 40 --counts-data=gene_name

cellphonedb plot dot_plot
cellphonedb plot heatmap_plot test_meta.txt
