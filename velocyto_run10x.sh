export OPENBLAS_NUM_THREADS=1
conda activate velocity
velocyto run10x -m /data/yudonglin/velocity/mm10_rmsk.gtf /data/yudonglin/singlecell/ABFC20190816-04-分析结果/NoPro /data/yudonglin/velocity/refdata-gex-mm10-2020-A/genes/genes.gtf

export OPENBLAS_NUM_THREADS=1
conda activate velocity
nohup velocyto run10x -m /data/yudonglin/velocity/mm10_rmsk.gtf /data/yudonglin/software/cellranger-4.0.0/erythropoiesis_nonpro/ /data/yudonglin/velocity/refdata-gex-mm10-2020-A/genes/genes.gtf &


export OPENBLAS_NUM_THREADS=1
conda activate velocity
velocyto run10x -m /data/yudonglin/velocity/mm10_rmsk.gtf /data/yudonglin/software/cellranger-4.0.0/erythropoiesis_nonpro/ /data/yudonglin/velocity/refdata-gex-mm10-2020-A/genes/genes.gtf


export OPENBLAS_NUM_THREADS=1
conda activate velocity
velocyto run10x -m /data/yudonglin/velocity/mm10_rmsk.gtf /data/yudonglin/software/cellranger-4.0.0/nonpro1/ /data/yudonglin/velocity/refdata-gex-mm10-2020-A/genes/genes.gtf



export OPENBLAS_NUM_THREADS=1
conda activate velocity
velocyto run10x -m /data/yudonglin/velocity/mm10_rmsk.gtf /data/yudonglin/software/cellranger-4.0.0/pro/ /data/yudonglin/velocity/refdata-gex-mm10-2020-A/genes/genes.gtf


#两个样本
export OPENBLAS_NUM_THREADS=1
conda activate velocity
velocyto run10x -m /data/yudonglin/velocity/mm10_rmsk.gtf /data/yudonglin/software/cellranger-4.0.0/baso_merge/ /data/yudonglin/velocity/refdata-gex-mm10-2020-A/genes/genes.gtf
