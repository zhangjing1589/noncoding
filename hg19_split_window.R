get_region_bed <- function(region_length)
{
  hg19_chr <- read.table("human_genome_chromosome_length.tsv", header = F)
  
  chr <- as.character(unique(hg19_chr$V1))
  chr_list <- vector("list", length = length(chr))
  names(chr_list)  <- chr
  
  for(i in chr)
  {
    end_pos <- hg19_chr[hg19_chr$V1 == i,]$V2
    df_end <- seq(0, end_pos, region_length)[-1]
    df_start <- df_end - region_length
    df_chr <- as.character(hg19_chr[hg19_chr$V1 == i,]$V1)
    df <- data.frame(df_chr, df_start, df_end)
    chr_list[[i]] <- df
  }
  
   result <- do.call("rbind", chr_list)
   result$name <- 1:dim(result)[1]
   write.table(result, "hg19_100kb.bed", sep = "\t", row.names = F, quote = F, col.names = F)

}  
