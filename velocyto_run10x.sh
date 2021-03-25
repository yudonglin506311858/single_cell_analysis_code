#export OPENBLAS_NUM_THREADS=1
conda activate velocity
velocyto run10x -m /data/yudonglin/velocity/mm10_rmsk.gtf /data/yudonglin/singlecell/ABFC20190816-04-分析结果/NoPro /data/yudonglin/velocity/refdata-gex-mm10-2020-A/genes/genes.gtf
