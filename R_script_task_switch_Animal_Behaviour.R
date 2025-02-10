#Script for analysis and data visualisation for the "Fungal infection intensity drives accelerated age polyethism without alteration in immune response in ants" in Animal Behaviour

#Cox mixed models
library(coxme)
task=read.csv('taskswitch.csv', sep=',')
head(task)

#cox mixed model 
cox.task.mm <- coxme(Surv(switch, status) ~ lvl_inf  + rickia + (1|nested/mother), data = task)
cox1=summary(cox.task.mm)
cox1

#task switch - plot - the probability of becoming forager regarding infection intensity
library(survival)
library(ggplot2)
library(dplyr)
library(survminer)

task=read.csv('taskswitch.csv', sep=',')
task$infcat=cut(task$lvl_inf, c(0,1,25,100), right = FALSE)

task$rickia=factor(task$rickia)
task$infcat=factor(task$infcat)

fit.task.cat <- survfit(Surv(age, status) ~ infcat, data=task)
p2=ggsurvplot(fit.task.cat,  conf.int=TRUE , legend.labs=c("uninfected", "infected (<25 thalli)", 'higly infected (>25 thalli)'),
              legend.title="", ylab='probability of staying a intranidal worker', xlab='time [days]',legend = c(0.23, 0.21),palette=c('#619cff','#F8766D','darkred'))

p2_df <- p2$plot$data
sum_df_cat <-summary(fit.task.cat, times = p2_df$time) 

columns_cat <- lapply(c(2:10,15:16), function(x) sum_df_cat[x])
tbl_cat <- do.call(data.frame, columns_cat)

custom_tbl_cat <- data.frame(time=tbl_cat$time, 
                             n.risk=tbl_cat$n.risk, 
                             n.event=tbl_cat$n.event,
                             cum.inc=1-tbl_cat$surv, 
                             std.err=1-tbl_cat$std.err,
                             upper=1-tbl_cat$lower,
                             lower=1-tbl_cat$upper,
                             strata=tbl_cat$strata)

custom_tbl_cat$strata=factor(custom_tbl_cat$strata, levels=c('infcat=[0,1)', 'infcat=[1,25)', 'infcat=[25,100)'), labels=c('uninfected',"infected (<25 thalli)", 'higly infected (>25 thalli)'))

task%>%group_by(infcat,status)%>%count()

#probability of foraging infection level (categories)
ggplot(custom_tbl_cat, aes(x=time,y=cum.inc,color=strata))+
  geom_line(linewidth=1.5)+
  geom_ribbon(data = custom_tbl_cat, aes(x = time, ymin = lower, ymax = upper, fill = strata), alpha = 0.05, color=NA)+
  labs(x='Worker age [days]',y='Proportion of foraging workers',fill='Infection status:', color="Infection status:")+
  theme_classic()+
  scale_color_manual(labels=c('Uninfected \n n = 62',"Infected (<25 thalli) \n n = 35", 'Higly Infected (>25 thalli) \n n = 14'),values=c('blue','red','darkred'))+
  scale_fill_manual(labels=c('Uninfected \n n = 62',"Infected (<25 thalli) \n n = 35", 'Higly Infected (>25 thalli) \n n = 14'),values=c('darkblue','red','darkred'))+
  scale_x_continuous(limits=c(0, 74), breaks = c(0,10,20,30,40,50,60,70))+
  scale_y_continuous(limits=c (0,1), labels = scales::percent)+
  #theme(legend.position = c(0.2,0.8))+
  theme(legend.position = 'top')

# Immune functions (PO an PPO) and workers size 
library(lme4)
library(ggplot2)
library(car)
library(DHARMa)
library(MuMIn)

ppo=read.csv("ppo.csv", sep=',')

ppo$colony=factor(ppo$colony)
ppo$nestid=factor(ppo$nestid)
ppo$rickia=factor(ppo$rickia)

ggplot(ppo, aes(x=head_len*head_wid, y = age))+
  geom_point()+
  geom_smooth(method = 'lm')


#head surface
ppo$head=ppo$head_len*ppo$head_wid


data.size=subset(ppo, head!="NA")

#infection effect on worker size
model.size <- lmer(head ~ rickia + (1|colony), data = data.size)
vif(model.size)
testDispersion(model.size)

model.size <- lmer(head ~ rickia + (1|colony), data = data.size, REML = FALSE)
model.reduced = lmer(head ~ (1|colony), data = data.size, REML = FALSE)

anova(model.size, model.reduced)


ggplot(ppo, aes(x=rickia, y=head, fill = rickia))+
  geom_boxplot()

kruskal.test(head~rickia, data=data.size)



ppo$size.index=ppo$head/916521 # mean head area
ppo$po.rel=ppo$po/ppo$size.index #standardising PO
ppo$ppo.rel=ppo$ppo/ppo$size.index #standardising PPO


#testing for correlation between PO and PPO

cor.test(ppo$po.rel, ppo$ppo.rel, method=c("pearson"))
cor.test(ppo$po.rel, ppo$ppo.rel, method=c("spearman"))


#Standardised PO (active phenoloxidase)
data.po=subset(ppo, po.rel!="NA")

model.po <- lmer(po.rel ~ rickia + lvl_inf + age  + task + (1|nestid), data = data.po)

model.po1 <- lm(po.rel ~ rickia + lvl_inf + age + task , data = data.po)
anova(model.po, model.po1)

vif(model.po1)
testDispersion(model.po1)

options(na.action = "na.fail")
po.dred=dredge(model.po1)

summary(model.avg(po.dred, subset = delta <= 2))


#Standardised PPO (total phenoloxidase)
data.ppo=subset(ppo, ppo.rel!="NA")

model.ppo <- lmer(ppo.rel ~ rickia + lvl_inf + age + task + (1|colony/nestid), data = data.ppo)
isSingular(model.ppo)

model.ppo1 <- lm(ppo.rel ~ rickia + lvl_inf + age  + task, data = data.ppo)

anova(model.ppo, model.ppo1)

testDispersion(model.ppo1)
vif(model.ppo1)

options(na.action = "na.fail")
ppo.dred=dredge(model.ppo1)
summary(model.avg(ppo.dred, subset = delta <= 2))



#figure phenoloxidase

ppo_scaled=read.csv("ppo_scaled.csv", sep = '\t')

summary(ppo_scaled)
y2_max <- 0.006021 #PPO max value
y1_max <- 0.00148 #PO max value

ppo_scaled$type=factor(ppo_s$type, levels = c('po', 'ppo'))
ggplot(ppo_scaled, aes(x = age, color = type))+
  geom_point(aes(y = scaled))+
  geom_smooth(aes(y=scaled), method ='lm')+
  scale_y_continuous(
    name = "Standardised PO (Vmax)",
    sec.axis = sec_axis(~ . * (y2_max / y1_max), name = "Standardised PPO (Vmax)"))+
  scale_color_manual(
    name = "Phenoloxidase",
    values = c("po" = "#F8766D", "ppo" = "#00BFC4"),
    labels = c("active", "total")
  )+scale_x_continuous(name = 'Worker age [days]', breaks = c(10,20,30,40,50,60, 70))+
  theme_classic() +
  theme(
    axis.line.y.right = element_line(color = "#00BFC4", linewidth = 0.7),
    axis.line.y.left = element_line(color = "#F8766D", linewidth = 0.7),
    axis.ticks.y.right = element_line(color = "#00BFC4"),
    axis.ticks.y.left = element_line(color = "#F8766D"),
    legend.position = "top"
  )