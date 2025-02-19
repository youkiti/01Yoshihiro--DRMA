# ------------------------------------------
# 1. 初期設定
# ------------------------------------------
# 環境変数を設定
Sys.setenv(LANGUAGE="en_US.UTF-8")
# メモリの中のデータをクリア
rm(list=ls())   
# コンソールをクリア
cat("\014")   
# 既存のグラフィックデバイスをオフにする
if(!is.null(dev.list())) dev.off()

# ------------------------------------------
# 2. パッケージのロード
# ------------------------------------------
library(dosresmeta)  # メタ分析のためのパッケージ
library(rms)        # 回帰モデリングツール
library(dplyr)      # データ操作ツール
library(openxlsx)   # Excelファイルの読み込みツール

# ------------------------------------------
# 3. データのロード
# ------------------------------------------
# Excelファイルからデータを読み込む
df <- openxlsx::read.xlsx("20231209_srws-psg.xlsx")

#
print(df)

# ------------------------------------------
# 4. 関数の定義
# ------------------------------------------

# 用量反応曲線をプロットする関数
plot_dosres <- function(mod, ylab, ylim){
  dosex_bin <- data.frame(dose= seq(0, 3, length.out = 100))
  xref_bin <- 0
  with(predict(mod, dosex_bin, xref_bin, exp = TRUE),{
    plot (get("rcs(dose, knots)dose"), pred, type = "l",
          ylim = ylim, 
          ylab = ylab, 
          xlab = "Brexpiprazole (mg)", 
          cex.lab  = 1, 
          cex.axis = 1,
          log = "y", 
          bty = "l", 
          las = 1,
    ) 
    abline(h=1,lty=3) 
    matlines(get("rcs(dose, knots)dose"), cbind(ci.ub, ci.lb), col = 1, lty = "dashed")})
}

# VPCをプロットする関数
plot_vpc <- function(data, mod, se){
  plot(data$dose[!is.na(data[,se])], vpc(mod), xlab = "BRE(mg)")
  lines(lowess(data$dose[!is.na(data[,se])], vpc(mod)))
}

# 制御イベント率を考慮して予測を行う関数
pred_md <- function(mod, blodds){
  newdata <- data.frame(dose = c(0,1,2,3))
  pred_md <- predict(mod, newdata = newdata, xref = 0, expo = TRUE)
  pred_md %>%
    mutate(N.lb = 100*((blodds)*ci.lb)/(1+(blodds*ci.lb)))->pred_md
  pred_md %>%
    mutate(N = 100*((blodds)*pred)/(1+(blodds*pred)))->pred_md
  pred_md %>%
    mutate(N.ub = 100*((blodds)*ci.ub)/(1+(blodds*ci.ub)))->pred_md
  return(round(pred_md, 2))
}

# 指定された効果の量での投与量を推定する関数
doseEff <- function(p, dose, Ep, trunc = FALSE){
  max <- max(Ep)
  EDmax <- dose[which.min(abs(Ep - max))]
  if (trunc == TRUE){
    if (EDmax == max(dose)) return(data.frame(p = NA, ED = NA, Ep = NA))      
  }
  ED <- apply(matrix(p), 1, function(x)
    dose[which.min(abs(Ep[dose < EDmax] - x * max(Ep)))])
  return(data.frame(p, ED, Ep = p * max(Ep)))
}

# ------------------------------------------
# 5. メインの分析: 2値アウトカムデータ
# ------------------------------------------
knots<-c(1,2,3)
#こちらで別途、knotsを四分位で指定することも可能
#knots<-quantile(df$dose, probs = c(0.25, 0.5, 0.75))

mod1 <- dosresmeta(formula = E_logor ~ rcs(dose, knots), 
                   type = 'cc', 
                   id = studyID, 
                   se = E_se, 
                   cases = N_responders_arm, 
                   n = N_arm, 
                   data = df[!(is.na(df$N_responders_arm)),], 
                   method="ml",
                   proc="1stage") 
summary(mod1)

plot_dosres(mod1, ylab="Response (OR)", ylim=c(0.75,2))
#ファイルとして保存
png(filename = "dichotomous_dose_response_curve.png", width = 1920, height = 1080, res = 300)
plot_dosres(mod1, ylab="Response (OR)", ylim=c(0.75,2))
dev.off()


plot_vpc(df, mod1, "E_se")
#ファイルとして保存
png(filename = "dichotomous_vpc_plot.png", width = 1920, height = 1080, res = 300)
plot_vpc(df, mod1, "E_se")
dev.off()


predicted_outcome <- pred_md(mod1, blodds=18.3/(100-18.3))
print(predicted_outcome)

p <- c(.5, .95, 1)
newdata <- data.frame(dose = seq(0, 3, length.out = 100))
edp <- with(predict(mod1, newdata), doseEff(p = p, dose = newdata$dose, Ep = pred, trunc = FALSE))
edp %>%
  mutate(ExpEp = exp(Ep))->edp
print(round(edp, 4))

# ------------------------------------------
# 6. メインの分析: 連続データ
# ------------------------------------------
df_cont <- subset(df, !is.na(MADRS))
mod2 <- dosresmeta(formula= MADRS ~ rcs(dose, knots),
                   id=studyID, 
                   sd=MADRS_SD,
                   n=MADRS_N,
                   covariance="smd",
                   data=df_cont,
                   proc="1stage") 
summary(mod2)

newdata<-data.frame(dose=seq(0,3,length.out = 100)) 
with(predict(mod2,newdata,order=TRUE),{ 
  plot(get("rcs(dose, knots)dose"),pred,type="l", ylim=c(-0.6,0),xlab="brexpiprazole(mg)",ylab="SMD") 
  matlines(get("rcs(dose, knots)dose"), cbind(ci.ub, ci.lb), col = 1, lty = "dashed")})

#ファイルとして保存
png(filename = "continuous_dose_response_curve.png", width = 1920, height = 1080, res = 300)
with(predict(mod2, newdata, order=TRUE), {
  plot(get("rcs(dose, knots)dose"), pred, type="l", ylim=c(-0.6,0), xlab="Brexpiprazole (mg)", ylab="SMD")
  matlines(get("rcs(dose, knots)dose"), cbind(ci.ub, ci.lb), col = 1, lty = "dashed")
})
dev.off()


