library(Rmisc)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(cowplot)
library(metR)
library(reshape2)
library(vegan)
library(stringr)
library(sqldf)
library(iNEXT)
library(tibble)


# load data ---------------------------------------------------------------
dtraw<-read.csv("./nan2.csv",stringsAsFactors = FALSE)
dt<-dtraw %>%
  mutate(.,ba08=ifelse(D08>=0,pi*(D08/100/2)^2,0)) %>%
  mutate(.,ba13=ifelse(D13>=0,pi*(D13/100/2)^2,0)) %>%
  mutate(.,ba19=ifelse(D19>=0,pi*(D19/100/2)^2,0)) %>%
  group_by(.,tag) %>%
  mutate(.,sumba08=sum(ba08),sumba13=sum(ba13),sumba19=sum(ba19)) %>%
  ungroup(.) %>%
  mutate(.,x1y1=paste0(x1,y1)) %>%
  filter(.,D08!=0 & D08!=-1 & D08 !=-5)


# dtb0 --------------------------------------------------------------------
dtb0<-dt %>%
  filter(.,b==0) %>%
  mutate(.,new13=ifelse(D08==-4 & D13>0,1,0)) %>%
  mutate(.,new19=ifelse(D13==-4 & D19>0,1,0)) %>%
  mutate(.,live08=ifelse(D08==0 | D08==-1 | D08==-4 | D08==-5,0,1)) %>%
  mutate(.,live13=ifelse(D13==0 | D13==-1 | D13==-4 | D13==-5,0,1)) %>%
  mutate(.,live19=ifelse(S19==0 | S19==-1 | S19==-5,0,1))
plotnum<-dt %>%
  distinct(x1,y1)
plotha<-nrow(plotnum)/100

# 从羆计 ---------------------------------------------------------------------
dtstem<-dtb0 %>%
  summarise(.,stem08=sum(live08),stem13=sum(live13),stem19=sum(live19))

# 贺计 ---------------------------------------------------------------------
dtstemsp<-dtb0 %>%
  group_by(.,sp) %>%
  summarise(.,stem08=sum(live08),stem13=sum(live13),stem19=sum(live19)) %>%
  arrange(.,sp)
dtspeciesnum<-dtstemsp %>%
  summarise(.,sp08=sum(stem08>0),sp13=sum(stem13>0),sp19=sum(stem19>0))


# ネ妓┦计 -----------------------------------------------------------------
diver<-data.frame(t(dtstemsp),stringsAsFactors = F)
colnames(diver)<-diver[1,]
diveryear=c("2008","2013","2019")
divershannon<-c(indexshannon08=diversity(as.numeric(diver[-c(1,3:4),]),index = "shannon"),
                indexshannon13=diversity(as.numeric(diver[-c(1:2,4),]),index = "shannon"),
                indexshannon19=diversity(as.numeric(diver[-c(1:3),]),index = "shannon"))
diverj<-c(indexj08=diversity(as.numeric(diver[-c(1,3:4),])/log(specnumber(as.numeric(diver[-c(1,3:4),])))),
          indexj13=diversity(as.numeric(diver[-c(1:2,4),])/log(specnumber(as.numeric(diver[-c(1:2,4),])))),
          indexj19=diversity(as.numeric(diver[-c(1:3),])/log(specnumber(as.numeric(diver[-c(1:3),])))))
diverdf<-data.frame("だ"=diveryear,"indexshannon"=divershannon,"indexj"=diverj,row.names=1)

               
# sample site bootstrap --------------------------------------------------
site2019dt<-dtb0 %>%
  group_by(.,sp) %>%
  mutate(.,sumsp=sum(live19)) %>% 
  filter(.,sumsp>0) %>% 
  count(.,x1y1) %>% 
  count(.,x1y1) %>% 
  summarise(.,num=sum(n))
site2019df<-add_row(site2019dt,num=64,.before=1)
colnames(site2019df)<-NULL
site2019df<-as.numeric(as.matrix(site2019df[,-1]))

site2013dt<-dtb0 %>%
  group_by(.,sp) %>%
  mutate(.,sumsp=sum(live13)) %>%
  filter(.,sumsp>0) %>%
  count(.,x1y1) %>%
  count(.,x1y1) %>%
  summarise(.,num=sum(n))
site2013df<-add_row(site2013dt,num=64,.before = 1)
colnames(site2013df)<-NULL
site2013df<-as.numeric(as.matrix(site2013df[,-1]))

site2008dt<-dtb0 %>%
  group_by(.,sp) %>%
  mutate(.,sumsp=sum(live08)) %>%
  filter(.,sumsp>0) %>%
  count(.,x1y1) %>%
  count(.,x1y1) %>%
  summarise(.,num=sum(n))
site2008df<-add_row(site2008dt,num=64,.before = 1)
colnames(site2008df)<-NULL
site2008df<-as.numeric(as.matrix(site2008df[,-1]))

sitelist<-list("2008"=site2008df,"2013"=site2013df,"2019"=site2019df)

siteplot<-iNEXT(sitelist,q=0,datatype="incidence_freq")
ggiNEXT(siteplot,se=F)+
  theme_classic()+
  theme(legend.position = "bottom",plot.margin = unit(c(0,0.5,0,0.5),"cm"))
# ggsave("sample site.png",width = 15,height = 10,units = "cm")

# Hurlbert's Expected Number of Species -----------------------------------
library(benthos)

plotspnum2008<-dtb0 %>%
  dcast(formula = sp~x1y1,value.var = "live08",fun.aggregate = sum)

plotspnum2013<-dtb0 %>%
  dcast(formula = sp~x1y1,value.var = "live13",fun.aggregate = sum)

plotspnum2019<-dtb0 %>%
  dcast(formula = sp~x1y1,value.var = "live19",fun.aggregate = sum)

rqdf2008<-data.frame(s=numeric())
rqdf2013<-data.frame(s=numeric())
rqdf2019<-data.frame(s=numeric())

for (i in 1:64) {
  rqdf2008[i,1]<-hurlbert(plotspnum2008,taxon = sp,count = plotspnum2008[,i+1],n = 10)
}
for (i in 1:64) {
  rqdf2013[i,1]<-hurlbert(plotspnum2013,taxon = sp,count = plotspnum2013[,i+1],n = 10)
}
for (i in 1:64) {
  rqdf2019[i,1]<-hurlbert(plotspnum2019,taxon = sp,count = plotspnum2019[,i+1],n = 10)
}

rqdf2008 %<>% bind_cols(.,data.frame(year=rep("2008",64),stringsAsFactors = F))
rqdf2013 %<>% bind_cols(.,data.frame(year=rep("2013",64),stringsAsFactors = F))
rqdf2019 %<>% bind_cols(.,data.frame(year=rep("2019",64),stringsAsFactors = F))

rqdfall<-bind_rows(rqdf2008,rqdf2013,rqdf2019) %>%
  mutate_if(.,is.character,as.factor)

rqdfmean<-rqdfall %>% 
  group_by(.,year) %>%
  summarise(.,S=mean(s)) %>%
  mutate_if(.,is.factor,as.character)

rqlm<-lm(s~year,rqdfall)
rqout<-summary(rqlm)
rqstd<-data.frame(std=rqout$coefficients[,2],year=c("2008","2013","2019"),stringsAsFactors = F)
rownames(rqstd)<-NULL
rqdfall %<>% mutate_if(.,is.factor,as.character)
rqstd<-full_join(rqstd,rqdfmean,by="year")

ggplot(data = rqdfall,aes(x=year,y=s,group=1))+
  # geom_smooth(method = "lm",se=F,color="black")+
  geom_point(data = rqdfmean,aes(x=year,y=S),inherit.aes = F)+
  geom_errorbar(data = rqstd,aes(x=year,ymin=S-std,ymax=S+std),inherit.aes = F,width=0.2)+
  labs(y="Rarefaction index")+
  theme_classic()
# ggsave("Rarefaction index.png",width = 10,height = 10,units = "cm")


# Differences in the rarefaction index ------------------------------------
rqdfall_bycols<-bind_cols(rqdf2008,rqdf2013,rqdf2019) %>%
  select(.,c(1,3,5))
colnames(rqdfall_bycols)<-c("r2008","r2013","r2019")
rqdfall_diff<-rqdfall_bycols %>%
  mutate(.,diff1st=r2013-r2008) %>%
  mutate(.,diff2nd=r2019-r2013)
rqdiffdf<-data.frame(year=c("2008-2013","2013-2019"),
                     rdiffmean=c(mean(rqdfall_diff$diff1st),mean(rqdfall_diff$diff2nd)),
                     rdiffsd=c(sd(rqdfall_diff$diff1st),sd(rqdfall_diff$diff2nd)))
rqdiffci<-t(data.frame(CI(rqdfall_diff$diff1st),CI(rqdfall_diff$diff2nd))) %>%
  as.data.frame(.) %>%
  add_column(.,year=c("2008-2013","2013-2019"),.before = 1)
rownames(rqdiffci)<-NULL

ggplot(data = rqdiffdf,aes(x=year,y=rdiffmean))+
  geom_point()+
  geom_errorbar(data = rqdiffci,aes(x=year,ymin=lower,ymax=upper),inherit.aes = F,width=0.2)+
  geom_abline(slope = 0)+
  labs(y="Differences in the rarefaction index")+
  theme_classic()
# ggsave("fig3_Differences in the rarefaction index.png",width = 10,height = 10,units = "cm")


# The effect of mortality (ED) and recruitment (ER) on species diversity --------
dtdead<-dt %>%
  filter(.,b==0) %>%
  mutate(.,dead13=ifelse((D08!=0 & D08!=-1) & (D13==0 | D13==-1),1,0)) %>%
  mutate(.,dead19=ifelse((D13!=0 & D13!=-1) & (S19==0 | S19==-1),1,0))

deadspnum2013<-dtdead %>%
  dcast(formula = sp~x1y1,value.var = "dead13",fun.aggregate = sum)

deadspnum2019<-dtdead %>%
  dcast(formula = sp~x1y1,value.var = "dead19",fun.aggregate = sum)

newspnum2013<-dtb0 %>%
  dcast(formula = sp~x1y1,value.var = "new13",fun.aggregate = sum)

newspnum2019<-dtb0 %>%
  dcast(formula = sp~x1y1,value.var = "new19",fun.aggregate = sum)

splist<-deadspnum2013[,1]

ed2013dt<-plotspnum2008[,-1]-deadspnum2013[,-1]
ed2013dt<-add_column(ed2013dt,sp=splist,.before = 1)

ed2019dt<-plotspnum2013[,-1]-deadspnum2019[,-1]
ed2019dt<-add_column(ed2019dt,sp=splist,.before = 1)

ed2013df<-data.frame(s=numeric())
ed2019df<-data.frame(s=numeric())

for (i in 1:64) {
  ed2013df[i,1]<-hurlbert(ed2013dt,taxon = sp,count = ed2013dt[,i+1],n = 10)
}
for (i in 1:64) {
  ed2019df[i,1]<-hurlbert(ed2019dt,taxon = sp,count = ed2019dt[,i+1],n = 10)
}

ed2013df %<>% bind_cols(.,data.frame(year=rep("2013",64),stringsAsFactors = F))
ed2019df %<>% bind_cols(.,data.frame(year=rep("2019",64),stringsAsFactors = F))

edall<-bind_cols(rqdf2008,rqdf2013,ed2013df,ed2019df) %>%
  select(.,c(1,3,5,7))
colnames(edall)<-c("s2008","s2013","rd2013","rd2019")
edall %<>% mutate(.,diff2013=rd2013-s2008,diff2019=rd2019-s2013)
eddf<-data.frame(year=c("2008-2013","2013-2019"),
                 edmean=c(mean(edall$diff2013),mean(edall$diff2019)),
                 edsd=c(sd(edall$diff2013),sd(edall$diff2019)))
edci<-t(data.frame(CI(edall$diff2013),CI(edall$diff2019))) %>%
  as.data.frame(.) %>%
  add_column(.,year=c("2008-2013","2013-2019"),.before = 1)
rownames(edci)<-NULL

edplot<-ggplot(data = eddf,aes(x=year,y=edmean))+
  geom_point()+
  geom_errorbar(data = edci,aes(x=year,ymin=lower,ymax=upper),inherit.aes = F)+
  geom_abline(slope = 0)+
  labs(y="Effect of mortality on\n species diversity")+
  theme_classic()
edplot


er2013dt<-plotspnum2008[,-1]+newspnum2013[,-1]
er2013dt<-add_column(er2013dt,sp=splist,.before = 1)

er2019dt<-plotspnum2013[,-1]+newspnum2019[,-1]
er2019dt<-add_column(er2019dt,sp=splist,.before = 1)

er2013df<-data.frame(s=numeric())
er2019df<-data.frame(s=numeric())

for (i in 1:64) {
  er2013df[i,1]<-hurlbert(er2013dt,taxon = sp,count = er2013dt[,i+1],n = 10)
}
for (i in 1:64) {
  er2019df[i,1]<-hurlbert(er2019dt,taxon = sp,count = er2019dt[,i+1],n = 10)
}

er2013df %<>% bind_cols(.,data.frame(year=rep("2013",64),stringsAsFactors = F))
er2019df %<>% bind_cols(.,data.frame(year=rep("2019",64),stringsAsFactors = F))

erall<-bind_cols(rqdf2008,rqdf2013,er2013df,er2019df) %>%
  select(.,c(1,3,5,7))
colnames(erall)<-c("s2008","s2013","rr2013","rr2019")
erall %<>% mutate(.,diff2013=rr2013-s2008,diff2019=rr2019-s2013)
erdf<-data.frame(year=c("2008-2013","2013-2019"),
                 ermean=c(mean(erall$diff2013),mean(erall$diff2019)),
                 ersd=c(sd(erall$diff2013),sd(erall$diff2019)))
erci<-t(data.frame(CI(erall$diff2013),CI(erall$diff2019))) %>%
  as.data.frame(.) %>%
  add_column(.,year=c("2008-2013","2013-2019"),.before = 1)
rownames(erci)<-NULL

erplot<-ggplot(data = erdf,aes(x=year,y=ermean))+
  geom_point()+
  geom_errorbar(data = erci,aes(x=year,ymin=lower,ymax=upper),inherit.aes = F)+
  ylim(-0.05,0.15)+
  geom_abline(slope = 0)+
  labs(y="Effect of recruitment on\n species diversity")+
  theme_classic()
erplot

bindplot<-plot_grid(edplot,erplot,ncol = 2,labels = "AUTO")
bindplot
# ggsave("fig4_ed and er plot.png",width = 20,height = 10,units = "cm")
