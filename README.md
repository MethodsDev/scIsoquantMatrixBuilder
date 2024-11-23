# scIsoquantMatrixBuilder
build sparse matrices for genes, isoforms, and unambiguous read-assigned isoforms for sc-Isoquant outputs


## from Ref Annot Guided Isoquant

```
scIsoseqUtil.RefAnnotGuided.py --sample_id ${sample_id} \
                --bam ${sample_id}.aligned.sorted.bam \
                --transcript_model_reads ${sample_id}.transcript_model_reads.tsv.gz \
                --transcript_models_gtf ${sample_id}.transcript_models.gtf.gz
```

## from Ref quant only (no additional modeling) IsoQuant

```
scIsoseqUtil.RefOnly.py --sample_id ${sample_id} \
                        --ref_annot_gtf GRCh38.gencode.v39.annotation.gtf.gz \
                        --transcript_read_assignments ${sample_id}.read_assignments.tsv.gz 
```

