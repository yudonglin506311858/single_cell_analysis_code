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
export PATH=/data/yudonglin/software/cellranger-4.0.0:$PATH
cellranger count --id=macrophage \
                   --transcriptome=/data/yudonglin/velocity/refdata-gex-mm10-2020-A \
                   --fastqs=/data/yudonglin/example/SRR12493999/ \
                   --sample=SRR12493999 \
                   --localcores=40 \
                   --localmem=64




