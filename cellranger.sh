conda activate chipseq
perfetch SRR10127227 --max-size 100G
/data/yudonglin/miniconda3/envs/chipseq/bin/fastq-dump --split-files --gzip SRR10127227.sra
export PATH=/data/yudonglin/software/cellranger-4.0.0:$PATH
cellranger count --id=neutrophil \
                   --transcriptome=/data/yudonglin/velocity/refdata-gex-mm10-2020-A \
                   --fastqs=/data/yudonglin/ncbi/public/sra/ \
                   --sample=SRR10127227 \
                   --localcores=40 \
                   --localmem=64

#巨噬细胞
conda activate chipseq
prefetch -t ascp-a "/data/yudonglin/miniconda3/envs/chipseq/etc/asperaweb_id_dsa.openssh" SRR12493999 -O ~/example
mv SRR12493999_1.fastq.gz SRR12493999_S1_L001_R1_001.fastq.gz
mv SRR12493999_2.fastq.gz SRR12493999_S1_L001_R2_001.fastq.gz
export PATH=/data/yudonglin/software/cellranger-4.0.0:$PATH
cellranger count --id=macrophage \
                   --transcriptome=/data/yudonglin/velocity/refdata-gex-mm10-2020-A \
                   --fastqs=/data/yudonglin/example/SRR12493999/ \
                   --sample=SRR12493999 \
                   --localcores=40 \
                   --localmem=64


#两个样本的运行

libraries.csv
library_id,molecule_h5
baso,/data/yudonglin/software/cellranger-4.0.0/baso/outs/molecule_info.h5
baso_2,/data/yudonglin/software/cellranger-4.0.0/baso_2/outs/molecule_info.h5

export PATH=/data/yudonglin/software/cellranger-4.0.0:$PATH

cellranger aggr --id=baso_merge \
                  --csv=libraries.csv \
                  --normalize=mapped



