
# Generate frame information at the 3' end. 

a=read.table("/Users/jwu/Downloads/5P_secondbatch_transcriptomebased_all_5end_stopcodon_RLchosed.txt")

sample <- unique(a$V4)

a[,5]=apply(a,1,function(x){b=as.numeric(x[1]);if(b%%3==0){return(0)}else if(b%%3==1){return(1)}else{return(2)}})

pdf("/Users/jwu/Downloads/5P_secondbatch_transcriptomebased_all_5end_stopcodon_RLchosed.pdf",height = 30,width=30)

for(i in 1:length(sample)){b=a[a$V4==sample[i],];p=ggplot(b)+geom_bar(aes(x=V1,y=V3+1,fill=factor(V5),width=0.8),stat = "identity")+facet_wrap(~V2,scales = "free",ncol=3)+coord_cartesian(xlim=c(-100,0))+scale_fill_manual(name="Class",values=color)+theme_bw()+ylab("Read counts")+xlab("Distance to stop codon")+theme(legend.position="non")+scale_x_continuous(breaks=seq(-100,42,by=6))+ggtitle(sample[i]);print(p)}

dev.off()
