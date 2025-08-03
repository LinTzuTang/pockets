## This project supports automated batch processing of RNA structures via [RNA3DB](https://github.com/marcellszi/rna3db/tree/main/scripts) and pocket extraction using SLURM.  
Follow the steps below to process your data on an HPC cluster:

### 1. RNA Structure Processing

Submit the RNA preprocessing job with:

```bash
sbatch process_rna.sh
```

### 2. Pocket Extraction

```bash
sbatch extract_pockets.sh
```
