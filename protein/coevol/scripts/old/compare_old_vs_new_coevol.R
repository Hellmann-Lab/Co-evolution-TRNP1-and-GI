#old vs new Coevol results

coevol_estimates_new<-readRDS("/data/share/htp/TRNP1/paper_data/protein/coevol/results/for_figures/coevol_3phenos_31sp_summarized.rds")

coevol_estimates_old<-readRDS("/data/share/htp/TRNP1/paper_data/data_tables/old/coevol_30sp_ready4plotting.rds")


names(coevol_estimates_new)<-paste0("new.",names(coevol_estimates_new))
names(coevol_estimates_old)<-paste0("old.",names(coevol_estimates_old))

coevol_estimates_comb<-full_join(coevol_estimates_new,coevol_estimates_old, by=c("new.species"="old.species"))

ggplot(coevol_estimates_comb, aes(x=new.omega,y=old.omega))+scale_x_log10()+scale_y_log10()+geom_point()+
  geom_text(aes(label=new.species))

ggplot(coevol_estimates_comb, aes(x=new.BodyM,y=old.BodyM))+scale_x_log10()+scale_y_log10()+geom_point()+
  geom_text(aes(label=new.species))
ggplot(coevol_estimates_comb, aes(x=new.EQ,y=old.EQ))+scale_x_log10()+scale_y_log10()+geom_point()+
  geom_text(aes(label=new.species))
