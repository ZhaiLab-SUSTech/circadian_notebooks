## Requirments
  1. Data Required
      - Genome file
        - Arabidopsis (Genome: TAIR10, GTF: Araport11)
      - Raw data
        10x data can be downloaded at [CNCB data base](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA021408)
  2. Softwares Required
      - CellRanger (Used for processing 10x raw data)
      - scanpy  (Used for downstream analysis of 10X data)
      - Seurat (Used for downstream analysis of 10X data)
      - scDblFinder (Used for removing putative doublets)
      - SCTransform (Used to implement 10X data normalization)
      - harmonypy (Used to implement 10X data integration)
      - stereopy (Used for downstream analysis of stereoseq data)
      - metacycle (Used to identify circadian genes)
      - AUCell (Used for calculating gene set expression score)
      - pyscenic (Used for calculating gene set expression score)
      - jpy_tools (A wrapper of single-cell analysis tools, which is available here [jpy_tools](https://github.com/liuzj039/jpy_tools/tree/master/jModule))
      - rpy2 (Used to implement invocation of R packages in python environment)
      - pygsea (Used to perform GO enrichmenth analysis)

## Main steps

### Preprocessing
1. Get 10X cell-gene matrix using Cell Ranger

### Analysis

These jupyter files contains the scripts needed for downstream analysis. Github often fails to preview large jupyter files, so you can preview these files using [nbviewer](https://nbviewer.org/github/ZhaiLab-SUSTech/circadian_notebooks/tree/master/). 

## Others
- The gene expression pattern can be explored at our [website](https://zhailab.bio.sustech.edu.cn/sc_circadian/)
  - If you found any bugs in our website, please reported [here](https://github.com/ZhaiLab-SUSTech/circadian_notebooks/issues/new)

## Notes
All notebooks share the same jupyter kernel.
