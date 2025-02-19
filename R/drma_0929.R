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
library(dosresmeta) # メタ分析のためのパッケージ
library(rms)       # 回帰モデリングツール
library(dplyr)     # データ操作ツール
library(openxlsx)  # Excelファイルの読み込みツール
library(gtsummary)

# ------------------------------------------
# 3. データのロード
# ------------------------------------------
# Excelファイルからデータを読み込む
df <- openxlsx::read.xlsx("ds_transfusion_infection.xlsx")

# データフレームの列名を確認
print(names(df))

df %>% tbl_summary()

# ------------------------------------------
# 4. データの準備
# ------------------------------------------
# 列名を確認し、必要に応じて変更
# 例: もし 'dose' という列が 'Distribution_in_group' の代わりに使われている場合
if("dose" %in% names(df) && !("Distribution_in_group" %in% names(df))) {
  df <- df %>% rename(Distribution_in_group = dose)
}

# データに合うように修正
df$N_arm <- df$n_group
df$N_responders_arm <- df$n_infection

# odds ratioと標準誤差を計算 ここでlogをとっているので、変数名はlogor
df$E_logor <- log((df$n_infection / (df$n_group - df$n_infection))) 
df$E_se <- sqrt(1/df$n_infection + 1/(df$n_group - df$n_infection)) 

# 新しいユニークなID列を作成 <- ここの意図は？
df$unique_id <- paste(df$Study_ID, df$Distribution_in_group, sep = "_")

df %>% select(Study_ID, Distribution_in_group, E_logor, E_se, N_responders_arm, N_arm) %>% head()

df %>% select(Study_ID, Distribution_in_group, E_logor, E_se, N_responders_arm, N_arm) %>% tbl_summary()

df <- df %>% select(Study_ID, Distribution_in_group, E_logor, E_se, N_responders_arm, N_arm)

df$Study_ID <- as.factor(df$Study_ID)


# ------------------------------------------
# 5. 関数の定義
# ------------------------------------------

# 用量反応曲線をプロットする関数 (修正済み)
plot_dosres <- function(mod, ylab, ylim) {
  dosex_bin <- data.frame(Distribution_in_group = seq(0, 6, length.out = 100))
  xref_bin <- 0
  with(predict(mod, dosex_bin, xref_bin, exp = TRUE), {
    plot(
      get("rcs(Distribution_in_group, knots)Distribution_in_group"),
      pred,
      type = "l",
      ylim = ylim,
      ylab = ylab,
      xlab = "Distribution in Group",
      cex.lab = 1,
      cex.axis = 1,
      log = "y",
      bty = "l",
      las = 1
    )
    abline(h = 1, lty = 3)
    matlines(get("rcs(Distribution_in_group, knots)Distribution_in_group"), cbind(ci.ub, ci.lb), col = 1, lty = "dashed")
  })
}

# VPCをプロットする関数 (修正済み)
plot_vpc <- function(data, mod, se){
  plot(data$Distribution_in_group[!is.na(data[,se])], vpc(mod), xlab = "Distribution in Group")
  lines(lowess(data$Distribution_in_group[!is.na(data[,se])], vpc(mod)))
}

# 制御イベント率を考慮して予測を行う関数 (修正済み)
pred_md <- function(mod, blodds){
  newdata <- data.frame(Distribution_in_group = c(1.0, 2.6, 4.3, 5.6))
  pred_md <- predict(mod, newdata = newdata, xref = 0, expo = TRUE)
  pred_md %>%
    mutate(N.lb = 100*((blodds)*ci.lb)/(1+(blodds*ci.lb)))->pred_md
  pred_md %>%
    mutate(N = 100*((blodds)*pred)/(1+(blodds*pred)))->pred_md
  pred_md %>%
    mutate(N.ub = 100*((blodds)*ci.ub)/(1+(blodds*ci.ub)))->pred_md
  return(round(pred_md, 2))
}

# 指定された効果の量での投与量を推定する関数 (修正済み)
doseEff <- function(p, Distribution_in_group, Ep, trunc = FALSE){
  max <- max(Ep)
  EDmax <- Distribution_in_group[which.min(abs(Ep - max))]
  if (trunc == TRUE){
    if (EDmax == max(Distribution_in_group)) return(data.frame(p = NA, ED = NA, Ep = NA))    
  }
  ED <- apply(matrix(p), 1, function(x)
    Distribution_in_group[which.min(abs(Ep[Distribution_in_group < EDmax] - x * max(Ep)))])
  return(data.frame(p, ED, Ep = p * max(Ep)))
}

# ------------------------------------------
# 6. メインの分析: 2値アウトカムデータ
# ------------------------------------------
# ノットをデータの範囲内に設定
knots <- quantile(df$Distribution_in_group, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

#dfの各変数の型
str(df)


head(df)

# dosresmeta 関数の実行
model <- dosresmeta(formula = E_logor ~ rcs(Distribution_in_group, knots = knots),
                    id = Study_ID,
                    se = E_se,
                    type = "ir",
                    cases = N_responders_arm,
                    n = N_arm,
                    data = df)

summary(mod1)

################
# 線形モデル
# メタ回帰モデルのフィッティング（自然スプライン）
res_spline <- rma(yi = E_logor, sei = E_se, mods = ~ ns(Distribution_in_group, df = 3), data = df_filtered, method = "REML")

# モデルの概要を表示
summary(res_spline)

library(ggplot2)
library(splines)  # For natural splines (ns)
library(metafor)  # For meta-analysis functions

# 予測用のデータ範囲を生成
new_data <- data.frame(Distribution_in_group = seq(min(df_filtered$Distribution_in_group), 
                                                   max(df_filtered$Distribution_in_group), 
                                                   length.out = 100))

# 予測値を取得
predictions <- predict(res_spline, newmods = ns(new_data$Distribution_in_group, df = 3))

# 予測結果をデータフレームに追加
new_data$pred <- predictions$pred

# 元データの散布図とスプラインフィットをggplotで描画
ggplot() +
  # 元データの散布図
  geom_point(data = df_filtered, aes(x = Distribution_in_group, y = E_logor), color = "gray", size = 3) +
  
  # スプラインフィットの予測線
  geom_line(data = new_data, aes(x = Distribution_in_group, y = pred), color = "blue", size = 1) +
  
  # ラベルとタイトルの追加
  labs(title = "Spline Fit for Meta-Regression",
       x = "Distribution_in_group",
       y = "Predicted E_logor") +
  
  # ggplotのテーマを調整
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))



##################





# 用量反応曲線をプロット
plot_dosres(mod1, ylab = "Infection Risk (OR)", ylim = c(0.5, 3)) 

# ファイルとして保存
png(
  filename = "dichotomous_dose_response_curve.png",
  width = 1920,
  height = 1080,
  res = 300
)
plot_dosres(mod1, ylab = "Infection Risk (OR)", ylim = c(0.5, 3))
dev.off()

# VPCをプロット (データに合わせて修正)
plot(df$Distribution_in_group, vpc(mod1), xlab = "Distribution in Group")
lines(lowess(df$Distribution_in_group, vpc(mod1)))

# ファイルとして保存
png(
  filename = "dichotomous_vpc_plot.png",
  width = 1920,
  height = 1080,
  res = 300
)
plot(df$Distribution_in_group, vpc(mod1), xlab = "Distribution in Group")
lines(lowess(df$Distribution_in_group, vpc(mod1)))
dev.off()

# doseEff を使う部分 (修正済み)
p <- c(.5, .95, 1)
newdata <- data.frame(Distribution_in_group = seq(0, 6, length.out = 100))
edp <- with(predict(mod1, newdata), doseEff(p = p, Distribution_in_group = newdata$Distribution_in_group, Ep = pred, trunc = FALSE))
edp %>%
  mutate(ExpEp = exp(Ep))->edp
print(round(edp, 4))