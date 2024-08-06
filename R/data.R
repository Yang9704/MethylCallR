#' @title Illmunima manifest files
#' @description 
#' These files are the Illumina manifest file downloaded from 
#' the Illumina official website (\url{https://support.illumina.com/?tab=microarrays}). 
#' To reduce the file size, several columns have been removed.
#' The file header and information of control probes have also been deleted.
#' The information of control probes can be found in `?Control_probe_chain_file`.
#' 
#' These files are used for annotation of probes.
#' @format A data frame with 937691 rows and 28 columns (EPICv2)
#'
#' @usage data(IlluminaHumanMethylationEPIC.V2a1.hg38)
#' @keywords datasets
#' @name Illumina_manifest
#' @rdname Illumina_manifest
#' @aliases IlluminaHumanMethylationEPIC.V2a1.hg38, IlluminaHumanMethylationEPIC.V1b5.hg38, IlluminaHumanMethylation450K
#'
"IlluminaHumanMethylationEPIC.V2a1.hg38"
#' @name manifest
#' @rdname manifest
#' @format A data frame with 865918 rows and 27 columns (EPICv1)
#' @usage data(IlluminaHumanMethylationEPIC.V1b5.hg38)
"IlluminaHumanMethylationEPIC.V1b5.hg38"
#' @name manifest
#' @rdname manifest
#' @format A data frame with 485577 rows and 20 columns (450K)
#' @usage data(IlluminaHumanMethylation450K)
"IlluminaHumanMethylation450K"


#' @title SeSaMe manifest file
#' @description 
#' These files are the SeSaMe manifest file downloaded from 
#' the github \url{https://zwdzwd.github.io/InfiniumAnnotation}
#' Since Illumina manifest of 450K is not provide the genomic position based 
#' the hg38, this file is required for filtering mapping fail probes and 
#' generating address chain file for 450K.
#' @format A data frame with 486427 rows and 27 columns
#'
#' @references
#' Zhou, W. et al. (2017). Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Nucleic Acids Research. \url{https://doi.org/10.1093/nar/gkw967}
#'
#' @usage data(HM450.hg38.manifest)
#' @keywords datasets
#' @name SeSaMe_manifest
#' @rdname SeSaMe_manifest
#'
"HM450.hg38.manifest"


#' @title SNP related probe list
#' @description 
#' These files contain probe IDs where SNPs are located closely.
#' These files can download at the github pages of SeSaMe (\url{https://zwdzwd.github.io/InfiniumAnnotation}).
#' 
#' Each row represents a probe.
#' Documentation for each column can be found in "header description" of github pages of SeSaMe (\url{https://zwdzwd.github.io/InfiniumAnnotation}).
#'
#' These files are used for filtering probes.
#'
#' @format A data frame with 82740 rows and 36 columns (EPICv2)
#'
#' @references 
#' For 450K and EPICv1 - 
#' Zhou, W. et al. (2017). Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Nucleic Acids Research. \url{https://doi.org/10.1093/nar/gkw967}
#'
#' For EPICv2 - 
#' Kaur, D. et al. (2023). Comprehensive evaluation of the Infinium human MethylationEPIC v2 BeadChip. Epigenetics Communications. \url{https://doi.org/10.1186/s43682-023-00021-5}
#'
#' @usage data(SNPmask.EPICv2)
#' @keywords datasets
#' @name SNPmask
#' @rdname SNPmask
#' @aliases SNPmask.EPICv2, SNPmask.EPICv1, SNPmask.450K
#'
"SNPmask.EPICv2"
#' @name SNPmask
#' @rdname SNPmask
#' @usage data(SNPmask.EPICv1)
#' @format A data frame with 105555 rows and 36 columns (EPICv1)
"SNPmask.EPICv1"
#' @name SNPmask
#' @rdname SNPmask
#' @usage data(SNPmask.450K)
#' @format A data frame with 485577 rows and 66 columns (450K)
"SNPmask.450K"



#' @title Mapping failed probe list
#' @description 
#' These files are a subset of the Illumina manifest or 
#' SeSaMe manifest that include mapping failed probes.
#' 
#' `V2a1.hg38.mapping.fail` for EPICv2 (Illumina)
#' `V1b5.hg38.mapping.fail` for EPICv1 (Illumina)
#' `hg38.mapping.fail.450K` for 450K (SeSaMe)
#' 
#' Each row represents a probe.
#' Documentation for each column can be found in the Illumina product files (\url{https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html}) 
#' or github pages of SeSaMe (\url{https://zwdzwd.github.io/InfiniumAnnotation})
#'
#' These files are used for filtering probes.
#'
#' @format A data frame with 190 rows and 48 columns (EPICv2)
#'
#' @references 
#' For 450K - 
#' Zhou, W. et al. (2017). Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Nucleic Acids Research. \url{https://doi.org/10.1093/nar/gkw967}
#'
#' @usage data(V2a1.hg38.mapping.fail)
#' @keywords datasets
#' @name mapping_failed_probe
#' @rdname mapping_failed_probe
#' @aliases V2a1.hg38.mapping.fail, V1b5.hg38.mapping.fail, hg38.mapping.fail.450K
"V2a1.hg38.mapping.fail"
#' @name mapping_failed_probe
#' @rdname mapping_failed_probe
#' @format A data frame with 38 rows and 51 columns (EPICv1) 
#' @usage data(V1b5.hg38.mapping.fail)
"V1b5.hg38.mapping.fail"
#' @name mapping_failed_probe
#' @rdname mapping_failed_probe
#' @format A data frame with 32 rows and 27 columns (450K)
#' @usage data(hg38.mapping.fail.450K)
"hg38.mapping.fail.450K"


#' @title Modified EPICv2 manifest by TJ Peters
#' @description This file is modified version of EPICv2 manifest provided by Peters.
#' This file can be found in supplementary materials of Peters' study 
#' (doi : url{https://doi.org/10.1186/s12864-024-10027-5}).
#' 
#' This file is used for filtering probes and defining duplicated probe set in EPICv2.
#' @format A data frame with 936166 rows and 79 columns
#'
#' @references 
#' Peters, T.J. et al. (2024). Characterisation and reproducibility of the HumanMethylationEPIC v2.0 BeadChip for DNA methylation profiling. BMC Genomics. url{https://doi.org/10.1186/s12864-024-10027-5}
#'
#'
#'
#' @usage data(Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024)
#' @name Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024
#' @keywords datasets
"Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024"



#' @title Multi-Hit probe list in EPICv1
#' @description
#' This file contains the probe list that have multiple aligned sequences.
#' This file can be found in supplementary materials of Nordlund's study
#' (doi : url{https://doi.org/10.1186/gb-2013-14-9-r105}).
#' 
#' This file is used for filtering probes.
#' @format A data frame with 9341 rows and 37 columns   
#'
#' @references 
#' Nordlund, J. et al. (2013). Genome-wide signatures of differential DNA methylation in pediatric acute lymphoblastic leukemia. Genome Biology. url{https://doi.org/10.1186/gb-2013-14-9-r105}
#'
#' @usage data(EPIC.MultiHit.Nordlund_GenomeBiology_2013)
#' @name EPIC.MultiHit.Nordlund_GenomeBiology_2013
#' @keywords datasets
"EPIC.MultiHit.Nordlund_GenomeBiology_2013"


#' @title The list of control probes in Illumina methylation BeadChip
#' @description This file contains control probes in EPICv2, EPICv1, and 450K.
#' This file generated based on the Illumina manifests from 3 types of platform. 
#' 
#'
#' This file is used for calculating detection P-value and normalization step.
#' @format A data frame with 985 rows and 8 columns
#'
#'
#'
#'
#'
#'
#' @usage data(Control_probe_chain_file)
#' @name Control_probe_chain_file
#' @keywords datasets
"Control_probe_chain_file"


#' @title The preset data of duplicated probe to remove in EPICv2
#'
#' @description This file contains probe IDs of duplicated probe in EPICv2.
#' Using this file, the address chain file can be generated, which is a 
#' reference file to analysis EPICv2 data or convert the EPICv2 data 
#' into previous versions. 
#'
#' To customize this file, we recommend to use Illumina manifest file.
#' (custom <- IlluminaHumanMethylationEPIC.V2a1.hg38[removelist, ])
#' 
#' Please refer to this file to create customized data for the `MeCall.SetChainFile` function.
#' 
#' This file is used to address duplicated probe and generate `address_chain_file` for preprocessing of EPICv2.
#' @format A data frame with 6963 rows and 28 columns
#'
#'
#'
#'
#'
#'
#' @usage data(Duplicated.Probes.preset)
#' @name Duplicated.Probes.preset
#' @keywords datasets
"Duplicated.Probes.preset"


#' @title Flagged probe list
#' @description This file is a subset of the Illumina EPICv2 manifest that 
#'  includes flagged probes.
#' This file is included in the Infinium MethylationEPIC v2.0 Product Files
#' from Illumina official website (\url{https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html}). 
#' 
#' This file is used for filtering probes.
#'
#' @format A data frame with 50209 rows and 48 columns
#'
#' @usage data(FLAGprobes)
#' @name FLAGprobes
#' @keywords datasets
"FLAGprobes"
