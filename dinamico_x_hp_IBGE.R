library(BayesMortalityPlus)
library(readxl)
library(MortalityLaws)
library(gridExtra)
library(tidyr)
library(dplyr)
options(scipen=999)

data_tot <- read_xls("IBGE 2020/ambos_os_sexos.xls")
#Ex_tot <- as.numeric(data_tot$...4[-c(1:4,45:60,101:111)])
#Dx_tot <- as.numeric(data_tot$...3[-c(1:4,45:60,101:111)])
qx_tot <- as.numeric(data_tot$...2[-c(1:4,45:60,101:111)])/1000
x <- 0:79

projs <- readODS::read_ods("IBGE 2020/projecoes_2018_populacao_idade_simples_2010_2060_20201209.ods",
                           sheet = "BR")

Ex_tot <- projs[197:276,12]   
Dx_tot <- Ex_tot * qx_tot

qx2 <- 1 - exp(-(Dx_tot/Ex_tot))
y <- log(qx2)

###ajustes pop total
fit1 <- hp(x, Ex_tot, Dx_tot, model="binomial")
fit2 <- dlm(y, delta=0.9)
#modelo HP usa thin=10, pois estamos interessados em amostras nao autocorrelacionadas

plotly:::ggplotly(plot(list(fit1,fit2),labels = c("HP","DLM")))
###dados já são graduados.
###dados projetados da exposição pra 2020, obitos obtidos a partir da relacao Dx = Ex*qx


##cadeias hp
plot(fit1)
plot_chain(fit1)   #ok
plot_chain(fit1, type = "acf")    #ok
summary(fit1)   ## 20.4, ok

saveRDS(fit1, "IBGE 2020/hp_79_IBGE_BR.rds")

#dlm
plot(fit2)
summary(fit2)

plot_chain(fit2, param = "sigma2")  #ok
gridExtra::grid.arrange(plot_chain(dlm_tot, param="sigma2"), 
                        plot_chain(dlm_tot, param="sigma2", type="acf"),
                        ncol = 2)
plot_chain(fit2, param = c("mu[0]","mu[6]","mu[10]",
                           "mu[30]","mu[60]","mu[79]"))  #idades criticas ok (0,6,10,40,79

saveRDS(fit2, "IBGE 2020/dlm_79_IBGE_BR.rds")



################ TENTATIVA 2
hp_tot <- readRDS("IBGE 2020/hp_79_IBGE_BR.rds")
dlm_tot <- readRDS("IBGE 2020/dlm_79_IBGE_BR.rds")

#### R EQM
dist_euc(qx_tot, fitted(hp_tot)$qx_fitted)
sqrt(sum((log(fitted(hp_tot)$qx_fitted) - log(qx_tot))^2)/length(qx_tot)) #dist. euclid. HP   #0.7465 #melhor 
dist_euc(qx_tot, fitted(dlm_tot)$qx_fitted)
sqrt(sum((log(fitted(dlm_tot)$qx_fitted) - log(qx_tot))^2)/length(qx_tot)) #dist. euclid. DLM #1.06  

#### distancia euclid. 50+ (inflexoes)
dist_euc(qx_tot[40:79], fitted(hp_tot)$qx_fitted[40:79])
sqrt(sum((log(fitted(hp_tot)$qx_fitted[40:79]) - log(qx_tot[40:79]))^2)/length(qx_tot[0:79])) #0.2525
dist_euc(qx_tot[40:79], fitted(dlm_tot)$qx_fitted[40:79])
sqrt(sum((log(fitted(dlm_tot)$qx_fitted[40:79]) - log(qx_tot[40:79]))^2)/length(qx_tot[0:79])) #0.1748 #melhor

#### distancia euclid. 0:40 (idades baixas)
dist_euc(qx_tot[0:40], fitted(hp_tot)$qx_fitted[0:40])
sqrt(sum((log(fitted(hp_tot)$qx_fitted[0:40]) - log(qx_tot[0:40]))^2)/length(qx_tot[0:40])) #0.697  #melhor
dist_euc(qx_tot[0:40], fitted(dlm_tot)$qx_fitted[0:40])
sqrt(sum((log(fitted(dlm_tot)$qx_fitted[0:40]) - log(qx_tot[0:40]))^2)/length(qx_tot[0:40])) #1.045 

#### HP SE COMPORTOU MELHOR NAS IDADES BAIXAS (0:40 ANOS), ENQUANTO QUE O MODELO DLM
#SE SAI MELHOR CAPTURANDO AS INFLEXOES E AS IDADES AVANÇADAS. NO GERAL, O MODELO HP
#SE AJUSTOU MELHOR AOS DADOS DO USA, COM MENOR DIST. EUCLID. NUM GERAL



######## IBGE MASC x FEM:
data_fem <- read_xls("IBGE 2020/mulheres.xls")
#Ex_fem <- as.numeric(data_fem$...4[-c(1:4,45:60,101:111)])
#Dx_fem <- as.numeric(data_fem$...3[-c(1:4,45:60,101:111)])
qx_fem <- as.numeric(data_fem$...2[-c(1:4,45:60,101:111)])/1000

Ex_fem <- projs[102:181,12]
Dx_fem <- Ex_fem * qx_fem

optf <- MortalityLaw(x, Dx_fem, Ex_fem, law="HP3")
coeff <- optf$coefficients

hp_ibge_fem <- hp(x, Ex_fem, Dx_fem, model = "binomial",
                  m = c(NA,NA,NA,NA,NA,17,NA,NA),
                  v = c(NA,NA,NA,NA,NA,1,NA,NA))
###priori F = 17 adquirida visualmente do ajuste

plot(hp_ibge_fem)
summary(hp_ibge_fem)

plot_chain(hp_ibge_fem, type="trace")
ggsave("ibge_fem_mcmc.png", device = png(width = 700, height = 500))
plot_chain(hp_ibge_fem, type="acf")
ggsave("ibge_fem_acf.png", device = png(width = 700, height = 500))
saveRDS(hp_ibge_fem, "IBGE 2020/hp_ibge_fem.rds")



data_hom <- read_xls("IBGE 2020/homens.xls")
#Ex_hom <- as.numeric(data_hom$...4[-c(1:4,45:60,101:111)])
#Dx_hom <- as.numeric(data_hom$...3[-c(1:4,45:60,101:111)])
qx_hom <- as.numeric(data_hom$...2[-c(1:4,45:60,101:111)])/1000

Ex_hom <- projs[6:85,12]
Dx_hom <- Ex_hom * qx_hom


hp_ibge_hom <- hp(x, Ex_hom, Dx_hom, model = "binomial")
summary(hp_ibge_hom)

plot(hp_ibge_hom)
plot_chain(hp_ibge_hom, type="trace")
ggsave("ibge_masc_mcmc.png", device = png(width = 700, height = 500))
plot_chain(hp_ibge_hom, type="acf")
ggsave("ibge_masc_acf.png", device = png(width = 700, height = 500))
saveRDS(hp_ibge_hom, "IBGE 2020/hp_ibge_hom.rds")



fits <- list(hp_ibge_hom,
             hp_ibge_fem)
saveRDS(fits, "IBGE 2020/ajustes_hp_IBGE_param.rds")
fits <- readRDS("IBGE 2020/ajustes_hp_IBGE_param.rds")
hp_ibge_hom <- fits[[1]]
hp_ibge_fem <- fits[[2]]

plot(list(hp_ibge_hom, hp_ibge_fem, hp_tot), labels = c("Masculino", "Feminino", "Total"),
     colors = c("slateblue","darkred","seagreen"), plotIC = FALSE)  ##ok (testar cores mais bonitas kk)
ggsave("ibge_mascxfem.png", device = png(width = 700, height = 500))


plot(fits)

params <- matrix(nrow=2,ncol=8)
colnames(params) <- c("A","B","C","D","E","F","G","H")
row.names(params) <- c("Masc", "Fem")
params_inf <- matrix(nrow=2,ncol=8)
colnames(params_inf) <- c("A","B","C","D","E","F","G","H")
row.names(params_inf) <- c("Masc", "Fem")
params_sup <- matrix(nrow=2,ncol=8)
colnames(params_sup) <- c("A","B","C","D","E","F","G","H")
row.names(params_sup) <- c("Masc", "Fem")


for(i in 1:2){
  for(j in 1:8){
    params[i,j] <- summary(fits[[i]])[j,"mean"]
    params_inf[i,j] <- summary(fits[[i]])[j, "2.5%"]
    params_sup[i,j] <- summary(fits[[i]])[j, "97.5%"]
  }
}

write.table(params,"IBGE 2020/params_IBGE.txt")
write.table(params_inf,"IBGE 2020/params_inf_IBGE.txt")
write.table(params_sup,"IBGE 2020/params_sup_IBGE.txt")


params <- data.frame(params)
params_inf <- data.frame(params_inf)
params_sup <- data.frame(params_sup)
###CRIAR TABELA DOS PARAMETROS
###ATUALIZAR O RELATORIO 2 COM ESSES NOVOS RESULTADOS

###DIFERENTES DELTAS:
fit2_1 <- dlm(y, delta=0.85)
fit2_2 <- dlm(y, delta=0.8)

fits <- list(dlm_tot, fit2_1, fit2_2)
saveRDS(fits, "IBGE 2020/dlm_deltas_IBGE_BR.rds")

fits <- readRDS("IBGE 2020/dlm_deltas_IBGE_BR.rds")
plot(fits, labels = c("0.9", "0.85", "0.8"))
ggsave("delta_IC.png", device = png(width = 700, height = 500))

plot(fits, labels = c("0.9", "0.85", "0.8"), plotIC = FALSE)
ggsave("delta_noIC.png", device = png(width = 700, height = 500))

hp_tot <- readRDS("IBGE 2020/hp_79_IBGE_BR.rds")

dist_euc <- function(x1,x2){
  sqrt(sum((log(x1,base=10) - log(x2, base=10))^2)/length(x2))
}

##vis ajustes:
plot(list(hp_tot,fits[[2]],fits[[3]]), labels=c("HP","DLM 0.85","DLM 0.8"))
ggsave("deltahp_IC.png", device = png(width = 700, height = 500))

plot(list(hp_tot,fits[[2]],fits[[3]]), labels=c("HP","DLM 0.85","DLM 0.8"), plotIC = FALSE)
ggsave("deltahp_noIC.png", device = png(width = 700, height = 500))

fit2_1 <- fits[[2]]
fit2_2 <- fits[[3]]

###total:
dist_euc(fitted(hp_tot)$qx_fitted, qx_tot)   ##0.324
dist_euc(fitted(fit2_1)$qx_fitted, qx_tot)   ##0.288
dist_euc(fitted(fit2_2)$qx_fitted, qx_tot)   ##0.178

###idades baixas
dist_euc(fitted(hp_tot)$qx_fitted[0:40], qx_tot[0:40])  ##0.30
dist_euc(fitted(fit2_1)$qx_fitted[0:40], qx_tot[0:40])  ##0.28
dist_euc(fitted(fit2_2)$qx_fitted[0:40], qx_tot[0:40])  ##0.17

dist_euc(fitted(hp_tot)$qx_fitted[40:79], qx_tot[40:79])  ##0.30
dist_euc(fitted(fit2_1)$qx_fitted[40:79], qx_tot[40:79])  ##0.28
dist_euc(fitted(fit2_2)$qx_fitted[40:79], qx_tot[40:79])

#####diminuir o fator de desconto ajudou muito o DLM, superando o modelo HP nas distancias euclideanas
library(ggplot2)

params <- read.table("IBGE 2020/params_IBGE.txt")
params_inf <- read.table("IBGE 2020/params_inf_IBGE.txt")
params_sup <- read.table("IBGE 2020/params_sup_IBGE.txt")

par(mfrow=c(3,3))
names <- c("A","B","C","D","E","F","G","H")
data <- data.frame(matrix(nrow=16,ncol=5))
for(i in 1:8){
  data[(i*2)-1,] <- c(params[names[i]][1,],params_inf[names[i]][1,],
                      params_sup[names[i]][1,],"Masculino", names[i])
  data[i*2,] <- c(params[names[i]][2,],params_inf[names[i]][2,],
                  params_sup[names[i]][2,],"Feminino", names[i])
}
colnames(data) <- c("point", "IC_inf", "IC_sup", "group", "param")
print(ggplot(data.frame(data)) + 
  geom_pointrange(aes(x = group, y = point, ymin = IC_inf, ymax = IC_sup)) + 
  labs( x = "", y ="", title = "Parâmetros") + facet_wrap(vars(param), scales = "free"))
ggsave("ibge_mascxfem_params.png", device = png(width = 700, height = 500))

plot(list(hp_tot,dlm_tot), labels = c("HP", "DLM"))
plot(hp_tot)
ggsave("hptot_IC.png", device = png(width = 700, height = 500))
library(ggplot2)
plot(list(hp_tot,dlm_tot), labels = c("HP", "DLM"), plotIC = FALSE)
ggsave("hpxdlm_noIC.png", device = png(width = 700, height = 500))


gridExtra::grid.arrange(plot(fits[[1]], labels = "DLM"), 
                        plot(hp_tot, labels = "HP"),
                        ncol = 2)
ggsave("dlmdeltas_IC.png", device = png(width = 900, height = 500))


#################### EXTRAPOLACAO
pred_dlm_100 <- predict(fits[[3]], h = 21)
closed <- hp_close(hp_tot, "hp", max_age = 100)

p <- plot(closed, plotIC = F, plotData = F, labels = "HP", age = 0:100) + 
  geom_line(data = fitted(fits[[3]]), aes(x = age, y = qx_fitted, col = "DLM")) +
  geom_line(data = pred_dlm_100, aes(x = Ages, y = qx_fitted, col = "DLM Predito", linetype = "longdash")) +
  ggplot2::scale_colour_manual(name = NULL, values = scales::hue_pal()(3), labels = c("DLM", "DLM predito", "HP")) +
  ggplot2::scale_linetype_manual(name = NULL, values = c("solid", "longdash"))+
  ggplot2::guides(fill = "none",
                  linetype = "none") 

g <- plot(closed, plotIC = F, plotData = F, labels = "HP", age = 70:100) + 
  geom_line(data = fitted(fits[[3]])[71:80,], aes(x = age, y = qx_fitted, col = "DLM")) +
  geom_line(data = pred_dlm_100, aes(x = Ages, y = qx_fitted, col = "DLM Predito", linetype = "longdash")) +
  ggplot2::scale_colour_manual(name = NULL, values = scales::hue_pal()(3), labels = c("DLM", "DLM predito", "HP")) +
  ggplot2::scale_linetype_manual(name = NULL, values = c("solid", "longdash"))+
  ggplot2::guides(fill = "none",
                  linetype = "none") + 
  ggplot2::scale_y_continuous(trans = 'log10', breaks = 10^-seq(3,0),
                              limits = 10^-c(3,0), labels = scales::comma) 

d <- gridExtra::grid.arrange(plot(closed, plotIC = F, plotData = F, labels = "HP", age = 0:100) + 
                          geom_line(data = fitted(fits[[3]]), aes(x = age, y = qx_fitted, col = "DLM")) +
                          geom_line(data = pred_dlm_100, aes(x = Ages, y = qx_fitted, col = "DLM Predito")) +
                          ggplot2::scale_colour_manual(name = NULL, values = scales::hue_pal()(3), labels = c("DLM", "DLM predito", "HP")) +
                          ggplot2::guides(fill = "none",
                                          linetype = "none"), 
                        plot(closed, plotIC = F, plotData = F, labels = "HP", age = 70:100) + 
                          geom_line(data = fitted(fits[[3]])[71:80,], aes(x = age, y = qx_fitted, col = "DLM")) +
                          geom_line(data = pred_dlm_100, aes(x = Ages, y = qx_fitted, col = "DLM Predito", linetype = "dashed")) +
                          ggplot2::scale_colour_manual(name = NULL, values = scales::hue_pal()(3), labels = c("DLM", "DLM predito", "HP")) +
                          ggplot2::guides(fill = "none",
                                          linetype = "none"),
                        ncol = 2)  ### trocar ylim da g?

d <- gridExtra::grid.arrange(p,
                             g,
                             ncol = 2)  ### trocar ylim da g?

plot(close_hp_110, labels = "fitted", plotIC = F) +
  geom_line(data = fitted(close_hp_110)[71:111,], aes(x = age, y = qx_fitted, col = "closed")) +
  scale_colour_manual(name = NULL, values = c("data" = "grey3",
                                              "fitted" = "red3",
                                              "closed" = "green4"))
ggsave("fechamento.png", device = png(width = 500, height = 360))

ggsave("predicoes_dois.png", d, device = png(width = 700, height = 500))
