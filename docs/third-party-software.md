# Third-party software

The repository contains workflow source and small deterministic fixtures. It
does not redistribute the executables, reference databases, or container image
archives used by the workflow.

Runtime software is obtained through the pinned container references recorded
in `assets/reporting/tool_registry.json`. Each tool and container remains
subject to its upstream license and terms. The pipeline's MIT license applies
to this repository's workflow source; it does not replace upstream licenses.

The principal upstream projects are:

- [Nextflow](https://www.nextflow.io/)
- [MultiQC](https://multiqc.info/)
- [NanoPlot](https://github.com/wdecoster/NanoPlot)
- [pychopper](https://github.com/epi2me-labs/pychopper)
- [minimap2](https://github.com/lh3/minimap2) and
  [SAMtools](https://www.htslib.org/)
- [PacBio lima, Iso-Seq, and pbmm2](https://github.com/PacificBiosciences)
- [IsoQuant](https://github.com/ablab/IsoQuant)
- [SQANTI3](https://github.com/ConesaLab/SQANTI3)
- [FLAIR](https://github.com/BrooksLabUCSC/flair)
- [bambu](https://bioconductor.org/packages/bambu/)
- [edgeR](https://bioconductor.org/packages/edgeR/),
  [satuRn](https://bioconductor.org/packages/satuRn/),
  [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/IsoformSwitchAnalyzeR/),
  [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/), and
  [org.Hs.eg.db](https://bioconductor.org/packages/org.Hs.eg.db/)
- [SUPPA](https://github.com/comprna/SUPPA)
- [TransDecoder](https://github.com/TransDecoder/TransDecoder)

Users preparing redistribution or deployment should review the applicable
upstream licenses for their selected execution profile.
