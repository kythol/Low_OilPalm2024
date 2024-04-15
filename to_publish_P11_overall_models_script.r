#!/usr/bin/env Rscript
options(echo=FALSE)
args <- commandArgs(trailingOnly =TRUE)

suppressMessages(library(rtracklayer))
suppressMessages(library(sqldf))
suppressMessages(library(dplyr))
suppressMessages(library(qdapRegex))
suppressMessages(library(intervals))
suppressMessages(library(plyr))
suppressMessages(library(data.table))

# Rscript to_publish_P11_overall_models_script.r input/testset_seqping.gff3 \
#     input/testset_mikado.gff3 \
#     input/testset_TSS.osc \
#     input/testset_mikado_proper_ATG.txt \
#     input/testset_gtf \
#     input/testset_maker_Liliopsida.txt \
#     input/testset_mikado_Liliopsida.txt \
#     input/testset_overlaps.txt

# prep output folder
ifelse(!dir.exists(file.path("tmp")), dir.create(file.path("tmp")), FALSE)
ifelse(!dir.exists(file.path("results")), dir.create(file.path("results")), FALSE)

# input
print("Load input data")
print(Sys.time())
seqpin<-readGFF(args[1], version=3)
mikado<-readGFF(args[2])
cage<-read.table(args[3],header=T)
mikado_proper<-read.table(args[4])
gtf_path=args[5]
maker_Liliopsida<-read.delim(args[6],sep="\t",header=T)
mikado_Liliopsida<-read.delim(args[7],sep="\t",header=T)
overlaps<-read.table(args[8])

print("Step 1. Maker and mikado predictions")
print(Sys.time())
# prepare data frame format for seqping annotation
seqpin_1<-as.data.frame(seqpin[, c("seqid","type","start","end","strand","ID","multiexonic","Parent")])
seqpin_1$seqid<-as.character(seqpin_1$seqid)
seqpin_1$type<-as.character(seqpin_1$type)
seqpin_1$start<-as.numeric(seqpin_1$start)
seqpin_1$end<-as.numeric(seqpin_1$end)
seqpin_1$strand<-as.character(seqpin_1$strand)
seqpin_1$ID<-as.character(seqpin_1$ID)
seqpin_1$multiexonic<-as.character(seqpin_1$multiexonic)
seqpin_1$Parent<-as.character(seqpin_1$Parent)
# get ATG coordinates for every model.
# Subset all CDS and for + strand obtain the min coordinate,
# and for - strand obtain the max coordinate and merge to the models file
unique_maker_isoforms<-unique(seqpin_1[seqpin_1["type"] == "mRNA",]$ID)
atg_list<-as.data.frame(matrix(0,ncol=2,nrow=0))
k=1
seqpin_1_CDS<-subset(seqpin_1, type == "CDS")
for( i in unique_maker_isoforms){
    a<-subset(seqpin_1_CDS, Parent == i)
    if( unique(a$strand) == "-"){
        c<-max(a[,c("start","end")])
    } else if(unique(a$strand) == "+"){
        c<-min(a[,c("start","end")])
    }
    yo<-data.frame(ID=i,ATG=c)
    atg_list<-rbind(atg_list,yo)
    # print(k)
    k=k+1
}
# the ATG is only located at the mRNA-type rows
maker<-left_join(seqpin_1, atg_list, by=c("ID"))

# repeat for mikado annotation.
# First change the data frame
mikado<-as.data.frame(mikado[, c("seqid","type","start","end","strand","Name","multiexonic","Parent")])
mikado$seqid<-as.character(mikado$seqid)
mikado$type<-as.character(mikado$type)
mikado$start<-as.numeric(mikado$start)
mikado$end<-as.numeric(mikado$end)
mikado$strand<-as.character(mikado$strand)
mikado$Name<-as.character(mikado$Name)
mikado$multiexonic<-as.character(mikado$multiexonic)
mikado$Parent<-as.character(mikado$Parent)
# get ATG coordinates for every model.
# Subset all CDS and for + strand obtain the min coordinate,
# and for - strand obtain the max coordinate and merge to the models file
unique_mikado_isoforms<-unique(mikado[mikado["type"] == "mRNA",]$Name)
atg_list_mikado<-as.data.frame(matrix(0,ncol=2,nrow=0))
k=1
mikado_CDS<-subset(mikado, type == "CDS")
for( i in unique_mikado_isoforms){
    a<-subset(mikado_CDS, Parent == i)
    if( unique(a$strand) == "-"){
        c<-max(a[,c("start","end")])
    } else if(unique(a$strand) == "+"){
        c<-min(a[,c("start","end")])
    }
    yo<-data.frame(Name=i,ATG=c)
    atg_list_mikado<-rbind(atg_list_mikado,yo)
    # print(k)
    k=k+1
}
mikado_1<-left_join(mikado, atg_list_mikado, by=c("Name"))

################################ CAGE #############################################
# add CAGE data. The inclusion criterion is that the CAGE tag is located within -3000 to -10 from the ATG for these gene isoform.
# One tag should be used per gene and if the tag is associated with a gene, it cannot be used again.
# First, we use maker models and then the remaining of the CAGE tags are used for the mikado models
# at this point we have separate annotations
# the models that do not have a CAGE tag are left in the data frame
# cage<-read.table("all_EG.merged.level2_clusters.osc",header=T)
# merge cage and models
print("Step 2. TSS information")
print(Sys.time())
cage_plus<-subset(cage, strand == "+")
cage_minus<-subset(cage, strand == "-")
maker_plus<-subset(maker, strand == "+")
maker_minus<-subset(maker, strand == "-")

maker_cage_plus<-sqldf("select * from maker_plus join cage_plus on maker_plus.seqid = cage_plus.chrom where cage_plus.pos between maker_plus.ATG - 3000 and maker_plus.ATG - 10")
maker_cage_plus$dist<-maker_cage_plus$pos - maker_cage_plus$ATG
maker_cage_minus<-sqldf("select * from maker_minus join cage_minus on maker_minus.seqid = cage_minus.chrom where cage_minus.pos between maker_minus.ATG + 10 and maker_minus.ATG + 3000")
maker_cage_minus$dist<- maker_cage_minus$ATG - maker_cage_minus$pos
maker_cage<-rbind(maker_cage_plus,maker_cage_minus)
# make a list of CAGE tags that are not associated with maker genes
a<-as.data.frame(table(maker_cage$id))
cage_id<-a$Var1
cage_list_left<-a[a["Freq"] == 0,]$Var1
maker_cage_list<-a[a["Freq"] != 0,]
# one CAGE tag can be associated with multiple isoforms of one gene,
# but it cannot be associated with several genes
# check this rule
filtered_maker_cage<-as.data.frame(matrix(0, ncol=18,nrow=0))
k=1
for( i in seq(1:dim(maker_cage_list)[1])){
    b<-maker_cage_list[i,]
    if(b$Freq == 1){
        # in case there is only one CAGE tag at all
        c<-maker_cage[maker_cage["id"] == paste0(b$Var1[[1]]),]
        filtered_maker_cage<-rbind(filtered_maker_cage,c)
        next
    } else {
        # in case there are >1 CAGE tags
        c<-maker_cage[maker_cage["id"] == paste0(b$Var1[[1]]),]
        if(length(unique(c$Parent)) == 1){
            # this tag corresponds to one gene name
            filtered_maker_cage<-rbind(filtered_maker_cage,c)
            next
        } else {
            # this tag corresponds to >1 gene name, choose the gene name (Parent) that is the closest to the ATG (min dist)
            d<-min(abs(c$dist))
            f<-c[c["dist"] == paste0("-",d),]$Parent
            c1<-c[c["Parent"] == paste0(f),]
            filtered_maker_cage<-rbind(filtered_maker_cage,c1)
        }
    }
    # print(k)
    k=k+1
}

# now check if every isoform has only one closest CAGE tag
filtered2_maker_cage<-as.data.frame(matrix(0, ncol=18,nrow=0))
k=1
for( i in unique(filtered_maker_cage$ID)){
    b<-filtered_maker_cage[filtered_maker_cage["ID"] == i,]
    if(dim(b)[1] == 1){
        # in case there is only one CAGE per isoform
        filtered2_maker_cage<-rbind(filtered2_maker_cage,b)
        next
    } else {
        # in case there are >1 CAGE tags per isoform
        if(length(unique(b$ID)) == 1){
            # all the CAGE tags are corresponding to the same isoform
            # pick the closest to ATG (dist)
            c<-b[b["dist"] == paste0("-",min(abs(b$dist))),]
            filtered2_maker_cage<-rbind(filtered2_maker_cage,c)
            next
        } else {
            # this tag corresponds to >1 gene name, choose the gene name (Parent) that is the closest to the ATG (min dist)
            # there were no such gene, I'm not sure what to do, let's see in the future
            # print(i)
        }
    }
    #print(k)
    #k=k+1
}
# final maker CAGE file
maker_CAGE<-data.frame(seqid=filtered2_maker_cage$seqid,start=filtered2_maker_cage$start,end=filtered2_maker_cage$end,strand=filtered2_maker_cage$strand, id=filtered2_maker_cage$ID,Parent=filtered2_maker_cage$Parent,ATG=filtered2_maker_cage$ATG,CAGE_id=filtered2_maker_cage$id,CAGE_TSS=filtered2_maker_cage$pos, CAGE_raw_counts=filtered2_maker_cage$raw.all_EG.merged.sorted_mapped_sorted, CAGE_norm_counts=filtered2_maker_cage$norm.all_EG.merged.sorted_mapped_sorted,dist_to_TSS=filtered2_maker_cage$dist)
# since we filtered CAGE-maker data frame, some CAGE tags were deleted as well.
# update a list of CAGE tags left after maker filtering
a<-as.data.frame(table(maker_CAGE$CAGE_id))
cage_id<-a$Var1
maker_cage_list1<-a[a["Freq"] != 0,]
cage_list_left1<-a[a["Freq"] == 0,]$Var1
print(paste0("the number of CAGE tags left increased from ",length(cage_list_left)," to ",length(cage_list_left1), " out of ",dim(cage)[1],". These CAGE tags will be used for mikado dataset mikado dataset merge with CAGE"))
# prepare CAGE dataset with only remaining CAGE tags
# update: check if there are any CAGE tags left for mikado, if not make a dummy file
if (length(cage_list_left1) == 0) {
    mikado_CAGE<-data.frame(seqid=NA,start=NA,end=NA,strand=NA, id=NA,Parent=NA,ATG=NA,CAGE_id=NA,CAGE_TSS=NA, CAGE_raw_counts=NA, CAGE_norm_counts=NA,dist_to_TSS=NA)
} else {
    cage_plus1<-cage_plus[cage_plus$id %in% cage_list_left1,]
    cage_minus1<-cage_minus[cage_minus$id %in% cage_list_left1,]

    mikado_plus<-subset(mikado_1, strand == "+")
    mikado_minus<-subset(mikado_1, strand == "-")

    mikado_cage_plus<-sqldf("select * from mikado_plus join cage_plus1 on mikado_plus.seqid = cage_plus1.chrom where cage_plus1.pos between mikado_plus.ATG - 3000 and mikado_plus.ATG - 10")
    mikado_cage_plus$dist<-mikado_cage_plus$pos - mikado_cage_plus$ATG
    mikado_cage_minus<-sqldf("select * from mikado_minus join cage_minus1 on mikado_minus.seqid = cage_minus1.chrom where cage_minus1.pos between mikado_minus.ATG + 10 and mikado_minus.ATG + 3000")
    mikado_cage_minus$dist<- mikado_cage_minus$ATG - mikado_cage_minus$pos
    mikado_cage<-rbind(mikado_cage_plus,mikado_cage_minus)
    # make a list of CAGE tags that were not associated with neither maker genes nor mikado genes
    a<-as.data.frame(table(cage$id))
    used_cage_ids<-unique(append(as.character(mikado_cage$id),as.character(maker_CAGE$CAGE_id)))
    print(paste0(length(used_cage_ids)," CAGE tags were used out of ",dim(cage)[1]))

    # one CAGE tag can be associated with multiple isoforms of one gene,
    # but it cannot be associated with several genes
    # check this rule
    a<-as.data.frame(table(mikado_cage$id))
    mikado_cage_list<-a[a["Freq"] != 0,]

    filtered_mikado_cage<-as.data.frame(matrix(0, ncol=18,nrow=0))
    for( i in seq(1:dim(mikado_cage_list)[1])){
        print(i)
        b<-mikado_cage_list[i,]
        c<-mikado_cage[mikado_cage["id"] == paste0(b$Var1[[1]]),]
        if(b$Freq == 1){
            # in case there is only one CAGE tag at all
            filtered_mikado_cage<-rbind(filtered_mikado_cage,c)
            next
        } else {
            # in case there are >1 CAGE tags
            if(length(unique(c$Parent)) == 1){
                # this tag corresponds to one gene name
                filtered_mikado_cage<-rbind(filtered_mikado_cage,c)
                next
            } else {
                # this tag corresponds to >1 gene name, choose the gene name (Parent) that is the closest to the ATG (min dist)
                d<-min(abs(c$dist))
                f<-c[c["dist"] == paste0("-",d),]$Parent
                c1<-c[c["Parent"] == paste0(f),]
                filtered_mikado_cage<-rbind(filtered_mikado_cage,c1)
            }
        }
    }


    # now check if every isoform has only one closest CAGE tag
    filtered2_mikado_cage<-as.data.frame(matrix(0, ncol=18,nrow=0))
    k=1
    for( i in unique(filtered_mikado_cage$Name)){
        # print(k)
        k=k+1
        b<-filtered_mikado_cage[filtered_mikado_cage["Name"] == i,]
        if(dim(b)[1] == 1){
            # in case there is only one CAGE per isoform
            filtered2_mikado_cage<-rbind(filtered2_mikado_cage,b)
            next
        } else {
            # in case there are >1 CAGE tags per isoform
            if(length(unique(b$Name)) == 1){
                # all the CAGE tags are corresponding to the same isoform
                # pick the closest to ATG (dist)
                c<-b[b["dist"] == paste0("-",min(abs(b$dist))),]
                filtered2_mikado_cage<-rbind(filtered2_mikado_cage,c)
                next
            } else {
                # this tag corresponds to >1 gene name, choose the gene name (Parent) that is the closest to the ATG (min dist)
                # there were no such genes so I'm not sure what to do, let's see in the future
                print(i)
            }
        }
    }
    # final mikado CAGE file
    mikado_CAGE<-data.frame(seqid=filtered2_mikado_cage$seqid,start=filtered2_mikado_cage$start,end=filtered2_mikado_cage$end,strand=filtered2_mikado_cage$strand, id=filtered2_mikado_cage$Name,Parent=filtered2_mikado_cage$Parent,ATG=filtered2_mikado_cage$ATG,CAGE_id=filtered2_mikado_cage$id,CAGE_TSS=filtered2_mikado_cage$pos, CAGE_raw_counts=filtered2_mikado_cage$raw.all_EG.merged.sorted_mapped_sorted, CAGE_norm_counts=filtered2_mikado_cage$norm.all_EG.merged.sorted_mapped_sorted,dist_to_TSS=filtered2_mikado_cage$dist)

}
# now let's go back to the full lists of models.
# sqldf only leaves rows that match the condition, thus we should merge the data back.
# it is possible to do full_join on this data
maker_tr<-subset(maker, type == "mRNA")
names(maker_tr)<-c("seqid","type","start","end","strand","id","multiexonic","Parent","ATG")
maker_with_cage<-full_join(maker_tr, maker_CAGE, by=c("seqid","start","end","strand","id","Parent","ATG"))
mikado_tr<-subset(mikado_1, type == "mRNA")
names(mikado_tr)<-c("seqid","type","start","end","strand","id","multiexonic","Parent","ATG")
mikado_with_cage<-full_join(mikado_tr, mikado_CAGE, by=c("seqid","start","end","strand","id","Parent","ATG"))
write.csv(maker_with_cage,paste0("tmp/","P11d_maker_with_cage_",Sys.Date(),".csv"))
write.csv(mikado_with_cage,paste0("tmp/","P11d_mikado_with_cage_",Sys.Date(),".csv"))
# in case of mikado, there were some models with questionable start/stop codons.
# Prepared a list of mikado models with proper start/stop codons
# that we will use to filter the models. Command: grep "^>" P11.mikado.CDS_J.fasta > P11_mikado_list_with_proper_atg_stop.txt
# mikado_proper<-read.table("P11_mikado_list_with_proper_atg_stop.txt")
mikado_proper<-data.frame(id=gsub(">","",mikado_proper$V1))
mikado_with_cage<-mikado_with_cage[mikado_with_cage$id %in% mikado_proper$id,]
################################ Models integration #############################################
# Models integration. We want to have a single list of models from both maker and mikado
# first - remove obviously identical models that are shared by two datasets (if they are identical,
# the CAGE tag should be in maker model not mikado model. If they are similar but not identical,
# e.g. the ATGs are different, there can be a different CAGE tag for mikado model as well
# overlapping models are not considered yet!

# check if the models from maker and mikado have identical start and end positions.
# In case they have, check if the CDS sequences from maker match
# CDS sequences from mikado. In case they do match perfectly, pick maker.
# take maker model one by one and look for a similar ones in mikado
print("Step 3. Model integration")
print(Sys.time())
# mikado1<-as.data.frame(mikado)
# mikado1$Parent<-as.character(mikado1$Parent)
together<-as.data.frame(matrix(0, ncol=16,nrow=0))
# mikado_delete -- the mikado models that are identical to maker models
mikado_delete<-as.data.frame(matrix(0, ncol=16,nrow=0))
# mikado_attention -- mikado models that have the same coordinates of mRNA
# but have different CDS coordinates (or number of CDS).
# for now we leave maker analogues of mikado models but we will have a look at them again
mikado_attention<-as.data.frame(matrix(0, ncol=17,nrow=0))
k=1
for( i in seq(1:dim(maker_with_cage)[1])){
    # print(k)
    k=k+1
    a<-maker_with_cage[i,]
    b<-mikado_with_cage[mikado_with_cage["seqid"] == a$seqid & mikado_with_cage["start"] == a$start & mikado_with_cage["end"] == a$end ,]
    if(dim(b)[1] == 0){
        together<-rbind(together, a)
        next
    } else if (dim(b)[1] == 1){
        # if there are two similar positions, let's compare their CDs and if they are the same we choose maker
        c<-seqpin_1[seqpin_1["Parent"] == a$id & seqpin_1["type"] == "CDS",]
        d<-mikado[mikado["Parent"] == b$id & mikado["type"] == "CDS",]
        c<-c[order(c$seqid,c$start,c$end),]
        d<-d[order(d$seqid,d$start,d$end),]
        rownames(c)<-NULL
        rownames(d)<-NULL
        if(identical(c[,c(1,3,4,5)],d[,c(1,3,4,5)]) == TRUE){
            # if all the CDS are identical, use maker data
            together<-rbind(together, a)
            mikado_delete<-rbind(mikado_delete,b)
        } else {
            together<-rbind(together, a)
            b$maker_id<-c(a$id)
            mikado_attention<-rbind(mikado_attention,b)
        }
    } else {
        # in case of several models with the same coordinates, check each separately
        for( j in unique(b$id)){
            g<-subset(b,id == j)
            c<-seqpin_1[seqpin_1["Parent"] == a$id & seqpin_1["type"] == "CDS",]
            d<-mikado[mikado["Parent"] == g$id & mikado["type"] == "CDS",]
            c<-c[order(c$seqid,c$start,c$end),]
            d<-d[order(d$seqid,d$start,d$end),]
            rownames(c)<-NULL
            rownames(d)<-NULL
            if(identical(c[,c(1,3,4,5)],d[,c(1,3,4,5)]) == TRUE){
                # if all the CDS are identical, use maker data
                mikado_delete<-rbind(mikado_delete,g)
            } else {
                # if they are not identical, mark the model with maker id and add to mikado_attention
                g$maker_id<-c(a$id)
                mikado_attention<-rbind(mikado_attention,g)
            }
            together<-rbind(together, a)
        }
    }
}
together<-together[!duplicated(together),]
# mikado_delete (even though some of them still have the CAGE tag) are to be deleted -- they are the same as in the maker
mikado_with_cage<-mikado_with_cage[!mikado_with_cage$id %in% mikado_delete$id,]
# filtering mikado_attention so that there are nrow = number of mikado models
# and if one mikado model matches several maker models they are listed in maker_id column
mikado_attention_filtered<-as.data.frame(matrix(0,ncol=15,nrow=0))
for(i in unique(mikado_attention$id)){
    a<-subset(mikado_attention,id==i)
    b<-a[,c(1:14)]
    b<-b[!duplicated(b),]
    if(dim(b)[1] == 1){
        b$maker_id<-paste(a$maker_id,collapse=",")
        mikado_attention_filtered<-rbind(mikado_attention_filtered,b)
    } else {
        # to see if there were any different models. There were none
        print(i)
    }
}
# merge data from mikado with mikado_attention to introduce maker_id column in mikado set
mikado_with_cage<-left_join(mikado_with_cage,mikado_attention_filtered,by=c('seqid','type','start','end','strand','id','multiexonic','Parent','ATG','CAGE_id','CAGE_TSS','CAGE_raw_counts','CAGE_norm_counts','dist_to_TSS'))
# introduce empty maker_id model in maker file
together$maker_id<-c(NA)
# merge maker data with remaining mikado models
together<-rbind(together,mikado_with_cage)
together<-together[order(together$seqid,together$strand,together$start,together$end),]
write.csv(together,paste0("tmp/","P11_preliminary_models_",Sys.Date(),".csv"))

################################ Expression #############################################
# Expression data integration. Right now I'm running expression data analysis on all datasets
# but using TPM count instead of coverage (and using transcript instead of exon in the code)
# all the GTF files are located at the folder
print("Step 4. RNAseq expression data")
print(Sys.time())
filenames <- list.files(path=gtf_path)
filenames<-data.frame(names=filenames)
filenames<-as.data.frame(filenames[grepl(".gtf",filenames[,1]),])
names <-apply(filenames, 1, function(x) strsplit(x, "[.]")[[1]][1])
for(i in names){
    print(i)
    filepath <- file.path(gtf_path,"/",paste(i,".stringtie.gtf",sep=""))
    assign(paste0("sample_",i), read.table(filepath,header = F,sep = "\t"))
    assign(paste0('sample_', i), transform(get(paste0('sample_', i)), sample =  rep(i)))
}
data<-do.call(rbind, mget(ls(pattern="sample_")))
data_transcript<-subset(data, V3 == "transcript")
data_exon<-subset(data, V3 == "exon")
data_exon$transcript_id<-unlist(rm_between(data_exon$V9,"transcript_id ",";",extract=T))
data_transcript$transcript_id<-unlist(rm_between(data_transcript$V9,"transcript_id ",";",extract=T))
data_transcript$TPM<-unlist(rm_between(data_transcript$V9,"TPM ",";",extract=T))
data_transcript1<-data_transcript[,c("transcript_id","TPM","sample")]

data_exon<-left_join(data_exon,data_transcript1,by=c("transcript_id","sample"))

models<-together
# correct for the name
names(mikado)[which(names(mikado) == "Name")] <- "id"
names(seqpin_1)[which(names(seqpin_1) == "ID")] <- "id"
cds_merged<-rbind(mikado, seqpin_1)
cds_merged<-cds_merged[cds_merged$id %in% models$id,]

# now find intervals of match
# first turn CDS into intervals data per scaffold
CDS_matches_total<-as.data.frame(matrix(0, ncol=23,nrow=0))
for( i in unique(cds_merged$seqid)){
    # print(i)
    a<-subset(cds_merged, seqid == i)
    # turn CDS into intervals
    b<-Intervals(a[,c("start","end")])
    # subset expression data for a scaffold
    a_exp<-subset(data_exon, V1 == i)
    if(dim(a_exp)[1] == 0){
        next
    }
    a$idd<-rownames(a)
    a_exp$iddd<-seq(1:dim(a_exp)[1])

    b_exp<-Intervals(a_exp[,c("V4","V5")])
    bb_overlap<-interval_overlap(b,b_exp)
    df<-ldply(bb_overlap,data.frame)
    if(dim(df)[1] == 0){
        next
    }
    names(df)<-c("idd","iddd")
    
    df<-inner_join(df,a,by=c("idd"))
    df<-inner_join(df,a_exp,by=c("iddd"))
    df$ratio<-c(NA)
    df$ratio<-ifelse(df$V4 <= df$start & df$V5 >= df$end,c(1),df$ratio)
    df$ratio<-ifelse(df$V4 > df$start & df$V5 >= df$end,c((df$end - df$V4)/(df$end - df$start)),df$ratio)
    df$ratio<-ifelse(df$V4 <= df$start & df$V5 < df$end,c((df$V5 - df$start)/(df$end - df$start)),df$ratio)
    df$ratio<-ifelse(df$V4 > df$start & df$V5 < df$end,c((df$V5 - df$V4)/(df$end - df$start)),df$ratio)
    CDS_matches_total<-rbind(CDS_matches_total,df)
}
# to calculate the ratio for an id rathrr than the CDS
aaaaa<-aggregate(CDS_matches_total$ratio, by=list(CDS_matches_total$id,CDS_matches_total$sample),mean)
names(aaaaa)<-c("id","sample","ratio")
CDS_matches_total$TPM<-as.numeric(CDS_matches_total$TPM)
aaaaa2<-aggregate(CDS_matches_total$TPM, by=list(CDS_matches_total$id,CDS_matches_total$sample),max)
names(aaaaa2)<-c("id","sample","TPM")
bb<-inner_join(aaaaa,aaaaa2,by=c("id","sample"))
#
mm<-models[,c("id","Parent")]
bb<-left_join(bb,mm,by=c("id"))
#
aaaaa_min<-aggregate(bb[bb["ratio"] > 0.8,]$TPM, by=list(bb[bb["ratio"] > 0.8,]$id),min)
aaaaa_median<-aggregate(bb[bb["ratio"] > 0.8,]$TPM, by=list(bb[bb["ratio"] > 0.8,]$id),median)
aaaaa_max<-aggregate(bb[bb["ratio"] > 0.8,]$TPM, by=list(bb[bb["ratio"] > 0.8,]$id),max)
bb$sample<-gsub("05WM","M05W",bb$sample)
bb$sample<-gsub("10WM","M10W",bb$sample)
bb$sample<-gsub("15WM","M15W",bb$sample)
bb$sample<-gsub("20WM","M20W",bb$sample)
bb$sample<-gsub("10WK","K10W",bb$sample)
bb$sample<-gsub("15WK","K15W",bb$sample)
bb$sample<-gsub("20WK","K20W",bb$sample)

for(i in unique(bb$sample)){
    assign(paste0("aa_",i), subset(bb, sample == i & ratio > 0.5))
    assign(paste0('aa_',i), get(paste0('aa_',i))[,c(1,3,4)])
}

names(aa_M05W)<-c("id","M05W_ratio","M05W_TPM")
names(aa_M10W)<-c("id","M10W_ratio","M10W_TPM")
names(aa_M15W)<-c("id","M15W_ratio","M15W_TPM")
names(aa_M20W)<-c("id","M20W_ratio","M20W_TPM")
names(aa_K10W)<-c("id","K10W_ratio","K10W_TPM")
names(aa_K15W)<-c("id","K15W_ratio","K15W_TPM")
names(aa_K20W)<-c("id","K20W_ratio","K20W_TPM")
names(aa_DRR053155)<-c("id","DRR053155_ratio","DRR053155_TPM")
names(aa_DRR053157)<-c("id","DRR053157_ratio","DRR053157_TPM")
names(aa_EC)<-c("id","EC_ratio","EC_TPM")
names(aa_NEC)<-c("id","NEC_ratio","NEC_TPM")
names(aa_EMB)<-c("id","EMB_ratio","EMB_TPM")
names(aa_GSD1)<-c("id","GSD1_ratio","GSD1_TPM")
names(aa_Plantlet)<-c("id","Plantlet_ratio","Plantlet_TPM")
names(aa_leaf_T0)<-c("id","leaf_T0_ratio","leaf_T0_TPM")
names(aa_leaf_T1)<-c("id","leaf_T1_ratio","leaf_T1_TPM")
names(aa_leaf_T2)<-c("id","leaf_T2_ratio","leaf_T2_TPM")
names(aa_PRJNA305816)<-c("id","PRJNA305816_ratio","PRJNA305816_TPM")
names(aa_SRR1612397)<-c("id","SRR1612397_ratio","SRR1612397_TPM")
names(aa_root_T1)<-c("id","root_T1_ratio","root_T1_TPM")
names(aa_root_T2S1)<-c("id","root_T2S1_ratio","root_T2S1_TPM")
names(aa_root_T2S3)<-c("id","root_T2S3_ratio","root_T2S3_TPM")
names(aa_trunk)<-c("id","trunk_ratio","trunk_TPM")

data<-data.frame(id=unique(bb$id))
for(i in unique(bb$sample)){
    data<- full_join(data, get(paste0('aa_',i)),by=c("id"))
}
names(aaaaa_min)<-c("id","TPM_min")
names(aaaaa_max)<-c("id","TPM_max")
names(aaaaa_median)<-c("id","TPM_median")
data<-full_join(aaaaa_max,data,by=c("id"))
data<-full_join(aaaaa_median,data,by=c("id"))
data<-full_join(aaaaa_min,data,by=c("id"))
expression<-data
types_of_tissues<-data.frame(tissue=names(expression)[c(6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)],group=c("mesocarp","kernel","mesocarp","kernel","mesocarp","kernel","mesocarp",rep("flowers",2),rep("embryos",3),rep("leaves",3),rep("embryos",2),"inflorescence",rep("root",3),"leaves","trunk"))
types_of_tissues$group <- as.factor(types_of_tissues$group)
# clean up TPMs with ratio < 0.8
for( i in gsub("_TPM","",types_of_tissues$tissue)){
    data[,grepl(paste0(i,"_TPM"),names(data))]<-ifelse(data[,grepl(paste0(i,"_ratio"),names(data))] < 0.8,c(NA),data[,grepl(paste0(i,"_TPM"),names(data))])
}
aa<-apply(data[,c(6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)],1,function(x) {
    a<-as.data.frame(x)
    a$tissue<-rownames(a)
    #print(a)
    a<-merge(a,types_of_tissues,by=c("tissue"))
    a<-a[!is.na(a$x),]
    b<-as.data.frame(table(a$group))
    #print(b)
    c<-data.frame(total=sum(b$Freq),mesocarp=b[b["Var1"] == "mesocarp",]$Freq,flowers=b[b["Var1"] == "flowers",]$Freq,inflorescence=b[b["Var1"] == "inflorescence",]$Freq,kernel=b[b["Var1"] == "kernel",]$Freq,embryos=b[b["Var1"] == "embryos",]$Freq,leaves=b[b["Var1"] == "leaves",]$Freq,trunk=b[b["Var1"] == "trunk",]$Freq,root=b[b["Var1"] == "root",]$Freq)
    c<-data.frame(total=sum(b$Freq),mesocarp=b[b["Var1"] == "mesocarp",]$Freq,flowers=b[b["Var1"] == "flowers",]$Freq,inflorescence=b[b["Var1"] == "inflorescence",]$Freq,kernel=b[b["Var1"] == "kernel",]$Freq,embryos=b[b["Var1"] == "embryos",]$Freq,leaves=b[b["Var1"] == "leaves",]$Freq,trunk=b[b["Var1"] == "trunk",]$Freq,root=b[b["Var1"] == "root",]$Freq)
    # print(c)
    return(c)
})
aa<-rbindlist(aa)
expression$Parent<-NULL
expression<-cbind(expression,aa)

write.csv(data,paste0("tmp/","expression_data_",Sys.Date(),".csv"))
# integrate with models file
together<-together[,1:15]
together<-left_join(together, expression, by=c("id"))
################################ BLAST #############################################
# upload BLAST results from Jane's files
# and rename them according to database
# Liliopsyda also has e-value for further classification
# maker_Liliopsida<-read.delim("maker.Lilopsida.isoforms-eval.txt",sep="\t",header=T)
print("Step 5. BLAST")
print(Sys.time())
names(maker_Liliopsida)<-c("id","title_Liliopsida","Liliopsida_eval")
# maker_Palms<-read.delim("maker.Palms.isoforms.txt",sep="\t",header=T)
maker_blast<-maker_Liliopsida
maker_blast$BLAST_Liliopsida<-ifelse(is.na(maker_blast$title_Liliopsida),c("."),c("+"))
# for the models file we will only leave +/- for all db,
# but also a title and e-value for Liliopsida
maker_blast<-maker_blast[,c("id","BLAST_Liliopsida","title_Liliopsida","Liliopsida_eval")]
maker_blast<-maker_blast[!duplicated(maker_blast),]
# upload BLAST results 
# and rename them according to database
# Liliopsida also has e-value for further classification
# mikado_Liliopsida<-read.delim("mikado.Lilopsida.isoforms-eval.txt",sep="\t",header=T)
names(mikado_Liliopsida)<-c("id","title_Liliopsida","Liliopsida_eval")
#merge all the the BLAST hits for maker
mikado_blast<-mikado_Liliopsida
mikado_blast$BLAST_Liliopsida<-ifelse(is.na(mikado_blast$title_Liliopsida),c("."),c("+"))
# for the models file we will only leave +/- for all db,
# but also a title and e-value for Liliopsida
mikado_blast<-mikado_blast[,c("id","BLAST_Liliopsida","title_Liliopsida","Liliopsida_eval")]
mikado_blast<-mikado_blast[!duplicated(mikado_blast),]
blast_set<-rbind(maker_blast,mikado_blast)
blast_set<-blast_set[!duplicated(blast_set),]
# check if there are multiple e-values and pick the minimum value
# first separate models that are unique from the repeated ones
a<-as.data.frame(table(blast_set$id))
blast_set_unique<-blast_set[blast_set$id %in% a[a["Freq"] == 1,]$Var1,]
blast_set_repeated<-blast_set[blast_set$id %in% a[a["Freq"] > 1,]$Var1,]
blast_set_filtered<-as.data.frame(matrix(0,ncol=7,nrow=0))
k=1
for( i in unique(blast_set_repeated$id)){
    # print(k)
    k=k+1
    a<-subset(blast_set_repeated, id == i)
    a$sum<-c(NA)
    if(dim(a)[1] == 1){
        blast_set_filtered<-rbind(blast_set_filtered,a)
        next
    } else {
        b<-a[,1:3]
        b<-b[!duplicated(b),]
        if(dim(b)[1] == 1){
            b<-a[a["Liliopsida_eval"] == min(a$Liliopsida_eval),]
            blast_set_filtered<-rbind(blast_set_filtered,b)
        } else if (dim(b)[1] > 1){
            c<-which(a$Liliopsida_eval == min(a$Liliopsida_eval))
            blast_set_filtered<-rbind(blast_set_filtered,a)
        }
    }
}
blast_set_filtered$sum<-NULL
blast_set_final<-rbind(blast_set_unique,blast_set_filtered)
# merge together with blast_set by id
together<-left_join(together, blast_set_final, by=c("id"))
together$title_Liliopsida<-as.character(together$title_Liliopsida)
together<-together[!duplicated(together),]
# finished finally?
write.csv(together,paste0("tmp/","P11_models_",Sys.Date(),".csv"))
# next is classification of the categories
#! add gene length according to CDS length (not ready)
# filter extra short sequences (<200 bp) (not ready)
########### classification of models #############################################
# models are classified according to columns blastLiliopsyda (+/.), CAGE (+/.) and EXP (+/.)
# prepare the columns
print("Step 6. Classification")
print(Sys.time())
together$EXP<-NULL # repeated so that I can define parameters together with cage and blast
together$EXP<-ifelse(together$TPM_max < 1, ".","+")
together$CAGE<-ifelse(together$CAGE_id == ".", ".","+")
together$blastLiliopsyda<-ifelse(together$BLAST_Liliopsida == ".", ".","+")
# now merge the +/.
together1<-together
together1[is.na(together1)]<-c(".")
together1$cat<-paste0(together1$blastLiliopsyda, together1$CAGE, together1$EXP)
# first assign all models 'no category'
together1$category<-rep("no category")
# assign category 5 for the models that have either expression or CAGE tag
together1$category<-ifelse( together1$CAGE == "+", "5A",together1$category)
together1$category<-ifelse( together1$EXP == "+", "5B",together1$category)
# assign category 4 for the models that have BLAST hit
together1$category<-ifelse( together1$blastLiliopsyda == "+", "4B",together1$category)
together1$category<-ifelse( together1$blastLiliopsyda == "+" & together1$EXP == "+", "4A",together1$category)
# now assign category 1 if CAGE, expression and BLAST are present
together1$category<-ifelse( together1$cat == "+++", "1",together1$category)
# now assign category 2 if CAGE and BLAST are present
together1$category<-ifelse( together1$cat == "++.", "2",together1$category)
# now assign category 3 if CAGE and expression are present
together1$category<-ifelse( together1$cat == ".++", "3",together1$category)
# differentiate between 1A (maker) and 1B (mikado)
together1$category<-ifelse(grepl("maker",together1$id) & together1$cat == "+++", "1A",together1$category)
together1$category<-ifelse(grepl("mikado",together1$id) & together1$cat == "+++", "1B",together1$category)
together1$category<-ifelse(grepl("maker",together1$id) & together1$cat == "++.", "2A",together1$category)
together1$category<-ifelse(grepl("mikado",together1$id) & together1$cat == "++.", "2B",together1$category)
# Category 1A – maker genes that have CAGE, expression and BLAST
# Category 1B - mikado genes that have CAGE, expression and BLAST
# Category 2A – maker genes that have CAGE and BLAST
# Category 2B - mikado genes that have CAGE and BLAST
# Category 3 – maker genes that have CAGE and expression
# Category 4 – maker genes that have BLAST (BLAST + Exp = 4A, BLAST only = 4B)
# Category 5 – maker genes that have either CAGE (5A) or expression (5B)
############ New Classes###########
# Defines new model classes
# Class 1: categories 1A and 1B (CAGE, BLAST & expression)
# Class 2: categories 2A and 2B and 3 (CAGE & BLAST/expression)
# Class 3: category 4A (BLAST & expression)
# Class 4: category 4B (only BLAST)
# Class 5: categories 5A and 5B
# Class 6: no category
together1$class<-c(NA)
together1$class<-ifelse(grepl("1",together1$category),c("1"),together1$class)
together1$class<-ifelse(grepl("2",together1$category),c("2"),together1$class)
together1$class<-ifelse(grepl("3",together1$category),c("2"),together1$class)
together1$class<-ifelse(grepl("4A",together1$category),c("3"),together1$class)
together1$class<-ifelse(grepl("4B",together1$category),c("4"),together1$class)
together1$class<-ifelse(grepl("5",together1$category),c("5"),together1$class)
together1$class<-ifelse(grepl("no category",together1$category),c("6"),together1$class)
together[is.na(together)]<-c(".")
write.csv(together1,paste0("results/","P11_models_categories_",Sys.Date(),".csv"))
#### overlaps
# overlaps<-read.table("OVERLAPS_P11.txt")
overlaps1<-overlaps[,c("V5","V12")]
names(overlaps1)<-c("id","overlapped_id")
overlaps1_filt<-as.data.frame(matrix(0,ncol=2,nrow=0))
k=0
for( i in unique(overlaps1$id)){
    # print(k)
    k=k+1
    a<-subset(overlaps1,id == i)
    if(dim(a)[1] == 1){
        overlaps1_filt<-rbind(overlaps1_filt,a)
    } else if (dim(a)[1] > 1){
        b<-data.frame(id=i,overlapped_id=paste(a$overlapped_id,collapse=","))
        overlaps1_filt<-rbind(overlaps1_filt,b)
    }
}
overlaps1_filt$id<-gsub("ID=","",overlaps1_filt$id)
overlaps1_filt$overlapped_id<-gsub("ID=","",overlaps1_filt$overlapped_id)
# merge with models
together1<-left_join(together1, overlaps1_filt,by=c("id"))
write.csv(together1,paste0("tmp/","P11_models_categories_",Sys.Date(),".csv"))
##################### final representatives ############################
print("Step 7. Representatives")
print(Sys.time())
models<-together1
cage_dist_sorgum <- subset(models, dist_to_TSS != '.')                   # subset only models with cage --> nrow(cage)=23884
cage_dist_sorgum$dist_to_TSS <- as.numeric(cage_dist_sorgum$dist_to_TSS)
cage_dist_sorgum$dist_to_TSS <- cage_dist_sorgum$dist_to_TSS * (-1)
theta = 96.08861    # sorgum
k = 1.353787        # sorgum
# calculate gamma distribution based on the parameters from Sorgum TSS data
cage_dist_sorgum$density <- dgamma(cage_dist_sorgum$dist_to_TSS, theta, k) 
cage_dist_sorgum <- as.data.frame(cage_dist_sorgum)
cage_dist_sorgum <- cage_dist_sorgum[, grepl('^(id)|dist_to_TSS|density', names(cage_dist_sorgum), perl=TRUE)]
# proceed with representatives
cage_sorgum<-cage_dist_sorgum[,c("id","density")]
names(cage_sorgum)[2]<-c("density_sorgum")
models<-left_join(models,cage_sorgum, by=c("id"))
# divide models that are one per gene or if there are multiple
# give categories a factor
# Maker Parent
seqpin_1<-as.data.frame(seqpin[, c("seqid","type","start","end","strand","ID","Name","Alias","multiexonic","Parent")])
seqpin_1$seqid<-as.character(seqpin_1$seqid)
seqpin_1$type<-as.character(seqpin_1$type)
seqpin_1$start<-as.numeric(seqpin_1$start)
seqpin_1$end<-as.numeric(seqpin_1$end)
seqpin_1$strand<-as.character(seqpin_1$strand)
seqpin_1$ID<-as.character(seqpin_1$ID)
seqpin_1$Name<-as.character(seqpin_1$Name)
seqpin_1$Alias<-as.character(seqpin_1$Alias)
seqpin_1$multiexonic<-as.character(seqpin_1$multiexonic)
seqpin_1$Parent<-as.character(seqpin_1$Parent)

master_parent<-seqpin_1[grepl("mikado",seqpin_1$Alias),c("ID","Name")]
names1 <-apply(master_parent, 1, function(x) strsplit(x["ID"], "[_]")[[1]][1])
names2 <-apply(master_parent, 1, function(x) strsplit(x["ID"], "[_]")[[1]][2])
names<-data.frame(id=master_parent$Name,Parent=c(paste0(names1,"_",names2)))
names<-names[names$id %in% models$id,]
for( i in 1:nrow(names)){
    # print(i)
    a<-names[i,]
    if(a$id %in% models$id){
        models[models["id"] == as.character(a$id),]$Parent<-as.character(a$Parent)
    } else {
        next
    }
}
# calculate CDS length
names(mikado)[6]<-"ID"
seqpin_1 <- seqpin_1[,names(mikado)]
for_cds_length<-rbind(seqpin_1,mikado)
for_cds_length<-subset(for_cds_length, type == "CDS")
for_cds_length$length<-for_cds_length$end - for_cds_length$start + 1
for_cds_length_agg<-aggregate(for_cds_length$length, by = list(for_cds_length$Parent), FUN= sum)
names(for_cds_length_agg)<-c("id","CDS_length")
models<-left_join(models, for_cds_length_agg, by=c("id"))
# resolve the issue with repeated Parents
models$Parent_1 <-apply(models, 1, function(x) strsplit(x["id"], "_R")[[1]][1])
models$Parent_1<-sub("\\.[^.]*$", "", models$Parent_1)
#
table(models$Parent == models$Parent_1)
models$check_parent<-c(NA)
models$check_parent<-ifelse(models$Parent == models$Parent_1, c("match"),c("attention"))
#
parent_db<-data.frame(maker_parent=as.character(models$Parent), mikado_parent=as.character(models$Parent_1))
parent_db$maker_parent<-as.character(parent_db$maker_parent)
parent_db$mikado_parent<-as.character(parent_db$mikado_parent)
parent_db<-parent_db[!duplicated(parent_db),]
parent_db<-parent_db[parent_db["maker_parent"] != parent_db["mikado_parent"],]
models$Parent_2<-models$Parent
for( i in parent_db$mikado_parent){
    a<-models[models["Parent_1"] == i,]
    if( length(unique(a$Parent)) == 1){
        next
    } else {
        a1<-parent_db[parent_db["mikado_parent"] == i,"maker_parent"]
        if( a1 %in% a$Parent){ #just a double-check
            models[models["Parent_1"] == i,]$Parent_2<-rep(a1)
        } else {
            print("ATTENTION")
        }
    }
}
# check again
table(models$Parent == models$Parent_2)
print("if the issue persists, sort through the incorrect Parent IDs manually")
models$Parent_old<-models$Parent
models$Parent<-models$Parent_2
#### representatives
# prepare for representatives: get only multi-isoform genes
a<-as.data.frame(table(models$Parent))
once<-a[a["Freq"] == 1,]$Var1
multiple_times<-a[a["Freq"] != 1,]$Var1
# if it is only one for gene, make it representative
models$representative<-c(NA)
models$representative<-ifelse(models$Parent %in% once, c("YES"),c("NO"))
models$category1<-factor(models$category,levels=c("1A","1B","2A","2B","3","4A","4B","5A","5B","no category"),ordered=T)
# if there are multiple models
multiple_models<-models[models["representative"] == "NO",]
representative_models<-as.data.frame(matrix(0, ncol=54,nrow=0))
what_the_cage<-as.data.frame(matrix(0, ncol=54,nrow=0))
for( i in unique(multiple_models$Parent)){
    a<-multiple_models[multiple_models["Parent"] == i,]
    #pick the best category
    b<-a[a["category"] == as.character(min(a$category1)),]
    if(dim(b)[1] == 1){
        b$representative<-c("YES")
        representative_models<-rbind(representative_models,b)
    } else {
        #for categories 1A and 1B
        if(any(b$category == "1A" | b$category == "1B")){
            if(all(is.na(b$density_sorgum))){
                what_the_cage<-rbind(what_the_cage,b)
                next
            }
            # first check the e-value
            b<-b[b["Liliopsida_eval"] == min(as.numeric(as.character(b$Liliopsida_eval)),na.rm=T),]
            b<-b[!is.na(b$id),]
            if (dim(b)[1] > 1){
                # next - density
                b<-b[b["density_sorgum"] == max(b$density_sorgum,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            if (dim(b)[1] > 1){
                    # next - number of tissues expressed
                b<-b[b["total"] == max(b$total,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            if (dim(b)[1] > 1){
                # next - max TPM
                b<-b[b["TPM_max"] == max(b$TPM_max,na.rm=T),]
            }
            if (dim(b)[1] > 1){
                # next - max TPM
                b<-b[b["CDS_length"] == max(b$CDS_length,na.rm=T),]
            }
            # if still there are more than 1 model, pick the first in the annotation
            if (dim(b)[1] > 1){
                if (any(grepl("maker_",b$id,fixed=T)) & any(grepl("mikado.",b$id,fixed=T))){
                    b<-b[grepl("maker_",b$id),]
                }
                if(all(grepl("maker_",b$id,fixed=T))){
                    yo<-sub('.*(?=.$)', '', b$id, perl=T)
                    yo<-yo[order(yo)]
                    b<-b[grepl(paste0("_R",yo[1]),b$id),]
                }
                if(all(grepl("mikado.",b$id,fixed=T))){
                    yo<-sub('.*\\.', '', b$id)
                    yo1<-gsub('(.*).\\w+', '\\1', b$id)
                    yo<-yo[order(as.numeric(yo))]
                    b<-b[grepl(paste0(yo1[1],".",yo[1]),b$id),]
                }
            }
            # if still there are more than 1 model, pick the first in the annotation
            if (dim(b)[1] > 1){
                if (any(grepl("maker_",b$id,fixed=T)) & any(grepl("mikado.",b$id,fixed=T))){
                    b<-b[grepl("maker_",b$id),]
                }
                if(all(grepl("maker_",b$id,fixed=T))){
                    yo<-sub('.*(?=.$)', '', b$id, perl=T)
                    yo<-yo[order(yo)]
                    b<-b[grepl(paste0("_R",yo[1]),b$id),]
                }
                if(all(grepl("mikado.",b$id,fixed=T))){
                    yo<-sub('.*\\.', '', b$id)
                    yo1<-gsub('(.*).\\w+', '\\1', b$id)
                    yo<-yo[order(as.numeric(yo))]
                    b<-b[grepl(paste0(yo1[1],".",yo[1]),b$id),]
                }
            }
            b$representative<-c("YES")
            representative_models<-rbind(representative_models,b)
        } else if (any(b$category == "2A" | b$category == "2B")){
            if(all(is.na(b$density_sorgum))){
                what_the_cage<-rbind(what_the_cage,b)
                next
            }
            b<-b[b["Liliopsida_eval"] == min(as.numeric(as.character(b$Liliopsida_eval)),na.rm=T),]
            b<-b[!is.na(b$id),]
            if (dim(b)[1] > 1){
                b<-b[b["density_sorgum"] == max(b$density_sorgum,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            if (dim(b)[1] > 1){
                # next - max TPM
                b<-b[b["CDS_length"] == max(b$CDS_length,na.rm=T),]
            }
            # if still there are more than 1 model, pick the first in the annotation
            if (dim(b)[1] > 1){
                if (any(grepl("maker_",b$id,fixed=T)) & any(grepl("mikado.",b$id,fixed=T))){
                    b<-b[grepl("maker_",b$id),]
                }
                if(all(grepl("maker_",b$id,fixed=T))){
                    yo<-sub('.*(?=.$)', '', b$id, perl=T)
                    yo<-yo[order(yo)]
                    b<-b[grepl(paste0("_R",yo[1]),b$id),]
                }
                if(all(grepl("mikado.",b$id,fixed=T))){
                    yo<-sub('.*\\.', '', b$id)
                    yo1<-gsub('(.*).\\w+', '\\1', b$id)
                    yo<-yo[order(as.numeric(yo))]
                    b<-b[grepl(paste0(yo1[1],".",yo[1]),b$id),]
                }
            }
            b$representative<-c("YES")
            representative_models<-rbind(representative_models,b)
        } else if (any(b$category == "3")){
            if(all(is.na(b$density_sorgum))){
                what_the_cage<-rbind(what_the_cage,b)
                next
            }
            if (dim(b)[1] > 1){
                b<-b[b["density_sorgum"] == max(b$density_sorgum,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            if (dim(b)[1] > 1){
                b<-b[b["total"] == max(b$total,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            if (dim(b)[1] > 1){
                b<-b[b["TPM_max"] == max(b$TPM_max,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            if (dim(b)[1] > 1){
                b<-b[b["CDS_length"] == max(b$CDS_length,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            # if still there are more than 1 model, pick the first in the annotation
            if (dim(b)[1] > 1){
                if (any(grepl("maker_",b$id,fixed=T)) & any(grepl("mikado.",b$id,fixed=T))){
                    b<-b[grepl("maker_",b$id),]
                }
                if(all(grepl("maker_",b$id,fixed=T))){
                    yo<-sub('.*(?=.$)', '', b$id, perl=T)
                    yo<-yo[order(yo)]
                    b<-b[grepl(paste0("_R",yo[1]),b$id),]
                }
                if(all(grepl("mikado.",b$id,fixed=T))){
                    yo<-sub('.*\\.', '', b$id)
                    yo1<-gsub('(.*).\\w+', '\\1', b$id)
                    yo<-yo[order(as.numeric(yo))]
                    b<-b[grepl(paste0(yo1[1],".",yo[1]),b$id),]
                }
            }
            b$representative<-c("YES")
            representative_models<-rbind(representative_models,b)
        } else if (any(b$category == "4A" | b$category == "4B")){
            if (dim(b)[1] > 1){
                b<-b[b["Liliopsida_eval"] == min(as.numeric(as.character(b$Liliopsida_eval)),na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            if (dim(b)[1] > 1){
                b<-b[b["total"] == max(b$total,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            if (dim(b)[1] > 1){
                b<-b[b["TPM_max"] == max(b$TPM_max,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            if (dim(b)[1] > 1){
                b<-b[b["CDS_length"] == max(b$CDS_length,na.rm=T),]
                b<-b[!is.na(b$id),]
            }

            if (dim(b)[1] > 1){
                if (any(grepl("maker_",b$id,fixed=T)) & any(grepl("mikado.",b$id,fixed=T))){
                    b<-b[grepl("maker_",b$id),]
                }
                if(all(grepl("maker_",b$id,fixed=T))){
                    yo<-sub('.*(?=.$)', '', b$id, perl=T)
                    yo<-yo[order(yo)]
                    b<-b[grepl(paste0("_R",yo[1]),b$id),]
                }
                if(all(grepl("mikado.",b$id,fixed=T))){
                    yo<-sub('.*\\.', '', b$id)
                    yo1<-gsub('(.*).\\w+', '\\1', b$id)
                    yo<-yo[order(as.numeric(yo))]
                    b<-b[grepl(paste0(yo1[1],".",yo[1]),b$id),]
                }
            }
            b$representative<-c("YES")
            representative_models<-rbind(representative_models,b)
        } else if (any(b$category == "5A" | b$category == "5B")){
            #if (dim(b)[1] > 1){
            if(all(b$category == "5A")){
                b<-b[b["density_sorgum"] == max(b$density_sorgum,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            #if (dim(b)[1] > 1){
            if(all(b$category == "5B")){
                b<-b[b["total"] == max(b$total,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            if (dim(b)[1] > 1){
                b<-b[b["CDS_length"] == max(b$CDS_length,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            # if still there are more than 1 model, pick the first in the annotation
            if (dim(b)[1] > 1){
                if (any(grepl("maker_",b$id,fixed=T)) & any(grepl("mikado.",b$id,fixed=T))){
                    b<-b[grepl("maker_",b$id),]
                }
                if(all(grepl("maker_",b$id,fixed=T))){
                    yo<-sub('.*(?=.$)', '', b$id, perl=T)
                    yo<-yo[order(yo)]
                    b<-b[grepl(paste0("_R",yo[1]),b$id),]
                }
                if(all(grepl("mikado.",b$id,fixed=T))){
                    yo<-sub('.*\\.', '', b$id)
                    yo1<-gsub('(.*).\\w+', '\\1', b$id)
                    yo<-yo[order(as.numeric(yo))]
                    b<-b[b["id"] == paste0(yo1[1],".",yo[1]),]
                }
            }
            b$representative<-c("YES")
            representative_models<-rbind(representative_models,b)
        } else if (any(b$category == "no category") & any(!is.na(b$total))){
            b<-b[b["total"] == max(b$total,na.rm=T),]
            b<-b[!is.na(b$id),]
            if (dim(b)[1] > 1){
                b<-b[b["CDS_length"] == max(b$CDS_length,na.rm=T),]
                b<-b[!is.na(b$id),]
            }
            # if still there are more than 1 model, pick the first in the annotation
            if (dim(b)[1] > 1){
                if (any(grepl("maker_",b$id,fixed=T)) & any(grepl("mikado.",b$id,fixed=T))){
                    b<-b[grepl("maker_",b$id),]
                }
                if(all(grepl("maker_",b$id,fixed=T))){
                    yo<-sub('.*(?=.$)', '', b$id, perl=T)
                    yo<-yo[order(yo)]
                    b<-b[grepl(paste0("_R",yo[1]),b$id),]
                }
                if(all(grepl("mikado.",b$id,fixed=T))){
                    yo<-sub('.*\\.', '', b$id)
                    yo1<-gsub('(.*).\\w+', '\\1', b$id)
                    yo<-yo[order(as.numeric(yo))]
                    b<-b[grepl(paste0(yo1[1],".",yo[1]),b$id),]
               }
            }
            b$representative<-c("YES")
            representative_models<-rbind(representative_models,b)
        } else if (any(b$category == "no category")){
            # if still there are more than 1 model, pick the first in the annotation
            if (dim(b)[1] > 1){
                if (any(grepl("maker_",b$id,fixed=T)) & any(grepl("mikado.",b$id,fixed=T))){
                    b<-b[grepl("maker_",b$id),]
                }
                if(all(grepl("maker_",b$id,fixed=T))){
                    yo<-sub('.*(?=.$)', '', b$id, perl=T)
                    yo<-yo[order(yo)]
                    b<-b[grepl(paste0("_R",yo[1]),b$id),]
                }
                if(all(grepl("mikado.",b$id,fixed=T))){
                    yo<-sub('.*\\.', '', b$id)
                    yo1<-gsub('(.*).\\w+', '\\1', b$id)
                    yo<-yo[order(as.numeric(yo))]
                    b<-b[grepl(paste0(yo1[1],".",yo[1]),b$id),]
                }
            }
            b$representative<-c("YES")
            representative_models<-rbind(representative_models,b)
        } else {
            print(i)
            break
        }
    }
}
multiple_models$representative<-ifelse(multiple_models$id %in% representative_models$id,c("YES"),c("NO"))
models$representative<-ifelse(models$id %in% representative_models$id, c("YES"),models$representative)
# check that there are no problems with CAGE data
dim(what_the_cage)[1] == 0
# save 
write.csv(models,paste0("results/","P11_models_representatives_",Sys.Date(),".csv"))
print("Finished")
print(Sys.time())
