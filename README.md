# GWAS-Imputation-Protocol-


## Requirements

**General**<br/>
-plink https://www.cog-genomics.org/plink2/

**Pre-Imputation-Processing**<br/>
-shapeit https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html
-checkVCF.py https://genome.sph.umich.edu/wiki/CheckVCF.py#Download
-modified HRC checking tool - HRC-1000G-check-bim-v4.2_AP_MODIFIED.pl [Mccarthy group](https://www.well.ox.ac.uk/~wrayner/tools/)

**Post-Imputation-Processing**<br/>
-R 3.5
-data.tables package
-optparse package
-HRC.r1-1.GRCh37.wgs.mac5.sites.tab from the [The Haplotype Reference Consortium](http://www.haplotype-reference-consortium.org/site)

## How to Run
**Create the following directory tree inside the Imputation Folder**<br/>
-/starting_plink_files
-/frequency_outputs
-/updated_plink_files
-/shapeit_phase_out
-/imputation/
-/imputation/pre-phased_vcfs

**For pre-Imputation-Processing**<br/>
Run: pipeline_pre_Imputation_Checking_Formatting.sh

Files inside imputation/pre-phased are now ready to be uploaded to Michigan Imputation Server https://imputationserver.readthedocs.io/en/latest/

**Michigan Imputation Server recommended parameters:**<br/>
Reference Panel: HRC r1.1 2016
Input Files: .vcf.gz and .vcf.tbi files inside the /imputation/pre-phased_vcfs directory
Phasing: ShapeIT v2.r790 (unphased)

**For post-Imputation-Processing**<br/>
Download your results (zip files) from Michigan Imputation Server inside the /imputation/ folder you created.
Output: A set of plink bfiles containing all the chromosomes merged. Only the variants with R2 quality greater than the threshold provided will be included in the file. The varians are now annotated to rsid numbers where possible. 
## 
