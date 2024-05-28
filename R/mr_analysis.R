
mr_analysis <- function(mr_data){
#接下来就可以进行MR分析了，在这里作者定义了5种方法，包括固定效应和随机效应模型
cat("Begin to analyze MR, it may take a long time......")
res<-mr(mr_data)
res
write.xlsx(res,"analysis results/mr_analysis_result.xlsx")

#生成了结果，结论是一样的，精神病对前臂骨质疏松没关联。

#进行了一个（MR-PRESSO）检验，这个也是多水平效应检验，P值应该要大于0.05
cat("Begin the MR-PRESSO test, it will take a long time......")
pres <- mr_presso(BetaOutcome="beta.outcome",
          BetaExposure ="beta.exposure",
          SdOutcome ="se.outcome",
          SdExposure = "se.exposure",
          OUTLIERtest =TRUE,DISTORTIONtest = TRUE,
          data =mr_data, NbDistribution = 1000,SignifThreshold = 0.05)

pres

write.xlsx(pres,"analysis results/mr_presso.xlsx")

# 查看离群值
Outlier <- pres$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`

if(!is.null(Outlier) & is.numeric(Outlier) )
  mr_data <- mr_data[-Outlier,] # 剔除离群值61和70

#异质性检验，没有异质性的。
mr_hete <- mr_heterogeneity(mr_data,
                 method_list=c("mr_egger_regression", "mr_ivw"))
mr_hete
write.xlsx(mr_hete,"analysis results/mr_heterogeneity.xlsx")

#多水平校验，这里是没有多水平效应的
pleio<- mr_pleiotropy_test(mr_data)
pleio

write.xlsx(pleio,"analysis results/mr_pleiotropy_test.xlsx")

#生成OR和可信区间
OR<-generate_odds_ratios(res)
write.xlsx(OR,"analysis results/OR value.xlsx")

#Leave-one-out analysis是指逐步剔除SNP后观察剩余的稳定性，理想的是剔除后变化不大
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6)
par(mar = c(0,4,2,0))
png(filename = "analysis results/1.mr_leaveoneout.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

single<- mr_leaveoneout(mr_data)
p1<- mr_leaveoneout_plot(single)
print(p1)

dev.off()

#散点图
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(filename = "analysis results/2.mr_scatter_plot.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

p2 <- mr_scatter_plot(res,mr_data)
print(p2)

dev.off()


#绘制森林图
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(filename = "analysis results/3.mr_forest_plot.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

res_single<- mr_singlesnp(mr_data)
p3 <- mr_forest_plot(res_single)
print(p3)

dev.off()

#绘制漏斗图，主要是看蓝线周围的散点是否对称
while (!is.null(dev.list()))  dev.off()#关闭Plots
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(filename = "analysis results/4.funnel_plot.png",
    width=800, height=600,units = "px",res = 100)
par(mai = c(1, 1, 1, 1)) # 边距 c(bottom, left, top, right)

p4 <- mr_funnel_plot(res_single)
print(p4)

dev.off()#关闭Plots

show_result <- list(result=p4)

return(show_result)

}



