#script to visualize number of peak overlaps with ANTXR2 sweep region

#read in data
key1=read.csv('encode_key_1.txt',header=F)
colnames(key1)=c('id','tissue')
intersect=read.table('overlap_number.txt')
ids=read.table('encode_ids.txt')
intersect$id=ids$V1
joined_df <- merge(intersect, key1, by.x = "id", 
             by.y = "id", all.x = TRUE, all.y = FALSE)
counts=data.frame(joined_df$tissue,joined_df$V4)


#get stats for each tissue
stats=data.frame()
for(i in unique(joined_df$tissue)){
maximum=max(joined_df$V4[which(joined_df$tissue==i)])
minimum=min(joined_df$V4[which(joined_df$tissue==i)])
average=mean(joined_df$V4[which(joined_df$tissue==i)])
med=median(joined_df$V4[which(joined_df$tissue==i)])
stddev=sd(joined_df$V4[which(joined_df$tissue==i)])
df=data.frame(i,maximum,minimum,average,med,stddev)
colnames(df)=c('tissue','maximum','minimum','average','med','stddev')
stats=rbind(stats, df)
}


#plot
pdf('tissue_average.pdf')
barplot(stats$average,names.arg=stats$tissue)
dev.off()

#write out table of stats for each tissue
write.table(stats, row.names=F, quote=F, file="summary_tissues.txt",sep=',')
