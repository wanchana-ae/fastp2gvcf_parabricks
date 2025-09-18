# fastp2gvcf_parabricks

## Description
Pipeline to Call SNP/Indels (GVCF file) by Parabricks

## Tutorials
[NVIDIA Parabricks Tutorials](https://docs.nvidia.com/clara/parabricks/latest/tutorials.html)

## Sample data
```bash
wget -O parabricks_sample.tar.gz "https://s3.amazonaws.com/parabricks.sample/parabricks_sample.tar.gz" 
```

## Metadata format for multi-samples
The metadata file for multiple samples should have the following format:

Sequence_ID,Sample_ID,/PATH_to_data/sample_1.fq.gz,/PATH_to_data/sample_2.fq.gz

- `Sequence_ID` → The ID of the sequencing run  
- `Sample_ID` → The name of the sample  
- `/PATH_to_data/sample_1.fq.gz` → Path to paired-end read 1  
- `/PATH_to_data/sample_2.fq.gz` → Path to paired-end read 2
