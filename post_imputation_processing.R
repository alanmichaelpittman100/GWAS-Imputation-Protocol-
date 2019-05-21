cat('#############################################
# Post Michigan Imputation Processing
    #############################################\n',sep='')
#pdir = dirname(sys.frame(1)$ofile)
source("dependencies.R")

option_list = list(
  make_option("--impute_dir", action="store", default=NA, type='character',
              help="Path to the imputation directory [required]"),
  make_option("--r2_threshold", action="store", default=NA, type='character',
              help="Imputation quality control threshold [required]"),
  make_option("--outNAME", action="store", default=NA, type='character',
              help="Name of the output plink file."),
  make_option("--zip_password", action="store", default=NA, type='character',
              help="Password to unzip the files downloaded from Michigan Imputation Server")
)

opt = parse_args(OptionParser(option_list=option_list))

impute_dir = "../GWAS_TUKraw/Imputation/tuk1M/imputation/"
impute_dir = opt$impute_dir
threshold = opt$r2_threshold
outNAME = opt$outNAME
outNAME = "tuk1M.biall.noDup.Rsq3"
passw = opt$zip_password

gsub(" ","",impute_dir)
gsub(" ","",threshold)
gsub(" ","",outNAME)

if(!is.na(passw)){
  gsub(" ","",passw)
}

prefix=""
chroms = paste0("chr",c(1:22))
suffix = "dose.vcf.gz"

files = paste0(impute_dir,prefix,chroms)
info_files_gz = paste0(impute_dir,prefix,chroms,".info.gz")
info_files = paste0(impute_dir,prefix,chroms,".info")

cat("##########PARAMETERS:###########\n")
cat(paste0("---Imputation Directory: ",impute_dir,"\n"))
cat(paste0("---Output file name: ",outNAME,"\n"))
cat(paste0("---R2 threshold: ",threshold,"\n"))
if(!is.na(passw)){
  cat(paste0("---ZIP password: ",passw,"\n"))
} else {cat(paste0("---ZIP password: No password was provided.",passw,"\n"))} 
cat("################################\n")

#Unzip 
cat('Decompressing ZIP files...\n',sep='')
system(paste0("for f in $(ls ",impute_dir,"*.zip -1) ; do unzip -P ",'"',passw,'"',' $f -d ',impute_dir,'; done'))

cat('Converting .vcf files to plink format and inserting RSIDs...\n',sep='')
opt
#hrc = data.frame(fread(hrc_ref))
#threshold=0.3


for(i in files){
  j=strsplit(i,"/chr")[[1]][2]
  
  cat('\n')
  cat('--------------------------------------------------------\n')
  cat(paste0('Processing ',i,'...\n'))
  #Decompress info file
  #system(paste0('gunzip -c ',i,'.info.gz > ',i,'.info'))
  
  #VCF to plink bfile format
  system(paste0(plink2,' --vcf ',i,'.dose.vcf.gz dosage=DS --threads 10 --exclude-if-info "R2<',threshold,'" --rm-dup exclude-mismatch --write-snplist --out ',i,'.TMPinit.biall.noDup.Rsq3'))
  
  cat('\n')
  cat(paste0('---INFO: Checking if ',i,'.TMPinit.biall.noDup.Rsq3.snplist was succesully created:\n'))
  stopifnot(file.exists(paste0(i,'.TMPinit.biall.noDup.Rsq3.snplist')))
  cat('-GOOD TO GO!\n')
  cat('\n')
  
  cat('---Convert variants IDs to rsid where possible\n')
  #Convert variants IDs to rsid where possible
  bim<-data.frame(fread(paste0(i,'.TMPinit.biall.noDup.Rsq3.snplist'),header = F))
  names(bim)<-"bim_ID"
  
  system(paste0("awk '$1 == \"",j,"\" { print $1",'"',"\t",'"',"$2",'"',"\t",'"',"$3 }' ",hrc_ref," > ", i, "_HRC_ref.txt"))
  
  
  HRC_ref<-data.frame(fread(paste0(i, "_HRC_ref.txt")))
  names(HRC_ref)<-paste0('HRC_ref_',names(HRC_ref))
  HRC_ref$HRC_ref_ID<-paste(HRC_ref[,1],HRC_ref[,2],sep=':')
  bim_HRC_ref<-merge(bim,HRC_ref, by.x='bim_ID', by.y='HRC_ref_ID')
  bim_HRC_ref = bim_HRC_ref[-which(duplicated(bim_HRC_ref$bim_ID)),]
  
  ID_update<-bim_HRC_ref[c('bim_ID','HRC_ref_V3')]
  ID_update<-ID_update[ID_update$HRC_ref_V3 != '.',]
  fwrite(ID_update, paste0(i,"_tempIDS",".txt"), col.names=F, sep=' ')
  
  system(paste0(plink2,' --vcf ',i,'.dose.vcf.gz dosage=DS --threads 10 --exclude-if-info "R2<',threshold,'" --rm-dup exclude-mismatch --update-name ',i,"_tempIDS.txt --export vcf vcf-dosage=GP --out ",i,'.VCFtmp.biall.noDup.Rsq3'))
  cat('\n')
  cat(paste0('---INFO: Checking if ',i,'.VCFtmp.biall.noDup.Rsq3.vcf was succesully created:\n'))
  stopifnot(file.exists(paste0(i,'.VCFtmp.biall.noDup.Rsq3.vcf')))
  cat('-GOOD TO GO!\n')
  cat('\n')
  
  #Remove Temporary Files
  system(paste0("rm ",i,"_HRC_ref.txt ",i,"_tempIDS.txt ",i,'.TMP*'))
  rm(bim,HRC_ref,bim_HRC_ref,ID_update)
  
  cat('--------------------------------------------------------\n')
  cat('\n')
  }


cat('---Merging per chromosome files...\n',sep='')

suffix = ".VCFtmp.biall.noDup.Rsq3.vcf"
bfiles = paste0(files,suffix)

cat(paste0('---INFO: Checking if all output vcf files were succesully created:\n'))
for(f in bfiles){
  stopifnot(file.exists(f))
}
cat('-GOOD TO GO!\n')
cat('\n')

bfiles_str = paste(bfiles,collapse=" ")
system(paste0("export PERL5LIB=",vcfperl,"; ",vcfconcat," ",bfiles_str," | gzip -c > ",impute_dir,outNAME,".vcf.gz"))