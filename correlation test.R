######CpG######
{
value <- fread("CpG_percent_Mb.bed")
mut <- fread("mutation_counts_mb.tsv")
result <- merge(value, mut, by.x = "V4", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)



value <- fread("CpG_percent_500kb.bed")
mut <- fread("mutation_counts_500kb.tsv")
result <- merge(value, mut, by.x = "V4", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)



value <- fread("CpG_percent_100kb.bed")
mut <- fread("mutation_counts_100kb.tsv")
result <- merge(value, mut, by.x = "V4", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)



value <- fread("CpG_percent_10kb.bed")
mut <- fread("mutation_counts_10kb.tsv")
result <- merge(value, mut, by.x = "V4", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)

}


######GC#########
{
value <- fread("hg19_gcMb.bed")
mut <- fread("mutation_counts_mb.tsv")
result <- merge(value, mut, by.x = "V4", by.y = "V4")
result <- result[result$V5.x != ".",]
cor(as.numeric(result$V5.x), result$V5.y)
cor.test(as.numeric(result$V5.x), result$V5.y)


value <- fread("hg19_gc500kb.bed")
mut <- fread("mutation_counts_500kb.tsv")
result <- merge(value, mut, by.x = "V4", by.y = "V4")
result <- result[result$V5.x != ".",]
cor(as.numeric(result$V5.x), result$V5.y)
cor.test(as.numeric(result$V5.x), result$V5.y)



value <- fread("hg19_gc100kb.bed")
mut <- fread("mutation_counts_100kb.tsv")
result <- merge(value, mut, by.x = "V4", by.y = "V4")
result <- result[result$V5.x != ".",]
cor(as.numeric(result$V5.x), result$V5.y)
cor.test(as.numeric(result$V5.x), result$V5.y)



value <- fread("hg19_gc10kb.bed")
mut <- fread("mutation_counts_10kb.tsv")
result <- merge(value, mut, by.x = "V4", by.y = "V4")
result <- result[result$V5.x != ".",]
cor(as.numeric(result$V5.x), result$V5.y)
cor.test(as.numeric(result$V5.x), result$V5.y)


}


###conservation###
{
value <- fread("conservation_Mb.bed")
mut <- fread("mutation_counts_mb.tsv")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)


value <- fread("conservation_500kb.bed")
mut <- fread("mutation_counts_500kb.tsv")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)


value <- fread("conservation_100kb.bed")
mut <- fread("mutation_counts_100kb.tsv")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)


value <- fread("conservation_10kb.bed")
mut <- fread("mutation_counts_10kb.tsv")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)
}


###reptime###
{
	mut_num <- read.table("hg19_noncod_mut_numberMB.tsv", stringsAsFactors = F)
	mut_num <- mut_num[mut_num$V5 != ".",]
	mut_num$V5 <- as.numeric(mut_num$V5)
	
	reptime <- read.table("reptime_Mb.bed", stringsAsFactors = F)
	reptime <- reptime[reptime$V5 != ".",]
	reptime$V5 <- as.numeric(reptime$V5)
	result <- merge(reptime, mut_num, by.x = "V4", by.y = "V4")
	
	result <- result[result$V5 != 0, ]
	cor(result$V5, result$V9)
	cor.test(result$V5, result$V9)

}


###polII###
{
mut <- fread("mutation_counts_mb.tsv")
value <- fread("polII_mb.bed")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V3.x != 0,]
cor(result$V3.x, result$V5)
cor.test(result$V3.x, result$V5)



mut <- fread("mutation_counts_500kb.tsv")
value <- fread("polII_500kb.bed")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V3.x != 0,]
cor(result$V3.x, result$V5)
cor.test(result$V3.x, result$V5)



mut <- fread("mutation_counts_100kb.tsv")
value <- fread("polII_100kb.bed")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V3.x != 0,]
cor(result$V3.x, result$V5)
cor.test(result$V3.x, result$V5)



mut <- fread("mutation_counts_10kb.tsv")
value <- fread("polII_10kb.bed")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V3.x != 0,]
cor(result$V3.x, result$V5)
cor.test(result$V3.x, result$V5)



}


###mappability###
{
mut <- fread("mutation_counts_mb.tsv")
value <- fread("mean_mb_mappability.tsv")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)



mut <- fread("mutation_counts_500kb.tsv")
value <- fread("mean_500kb_mappability.tsv")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)



mut <- fread("mutation_counts_100kb.tsv")
value <- fread("mean_100kb_mappability.tsv")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)



mut <- fread("mutation_counts_10kb.tsv")
value <- fread("mean_10kb_mappability.tsv")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V8 != 0,]
cor(result$V8, result$V5.y)
cor.test(result$V8, result$V5.y)




}


###rec_rate###
{
	mut_num <- read.table("hg19_noncod_mut_numberMB.tsv", stringsAsFactors = F)
	mut_num <- mut_num[mut_num$V5 != ".",]
	mut_num$V5 <- as.numeric(mut_num$V5)
	
	rec_rata <- read.table("rec_rate.tsv", header = T)
	result <- merge(rec_rata, mut_num, by.x = "V4", by.y = "V4")
    
	cor()
}


###tfbs###
{
mut <- fread("mutation_counts_mb.tsv")
value <- fread("tfbs_mb.bed")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)



mut <- fread("mutation_counts_500kb.tsv")
value <- fread("tfbs_500kb.bed")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)



mut <- fread("mutation_counts_100kb.tsv")
value <- fread("tfbs_100kb.bed")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)


mut <- fread("mutation_counts_10kb.tsv")
value <- fread("tfbs_10kb.bed")
result <- merge(value, mut, by.x = "V1", by.y = "V4")
result <- result[result$V5.x != 0,]
cor(result$V5.x, result$V5.y)
cor.test(result$V5.x, result$V5.y)


	


}









