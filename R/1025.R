# ------------------------------------------
# 1. 初期設定
# ------------------------------------------
Sys.setenv(LANGUAGE="en_US.UTF-8")
rm(list=ls())   
cat("\014")   
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
# 3. データのロードと準備
# ------------------------------------------
# Excelファイルからデータを読み込む
df <- openxlsx::read.xlsx("ds_transfusion_infection_001.xlsx")

df

# オッズ比 (E_or)、ログオッズ比 (E_logor)、標準誤差 (E_se) を計算
# df <- df %>%
#   mutate(
#     E_or = (n_infection / (n_group - n_infection)),   # オッズ比の計算
#     E_logor = log(E_or),                              # ログオッズ比の計算
#     E_se = sqrt(1/n_infection + 1/(n_group - n_infection)) # 標準誤差の計算
#  )

# # 基準行 (Referent Dose) をデータフレームに追加
# referent_row <- data.frame(
#   Study_ID = "Referent",              # 適切なIDを入力
#   Distribution_in_group = 0,          # 基準用量
#   n_group = 1,                        # グループサイズ（1に設定）
#   n_infection = 0,                    # 感染者数（0に設定）
#   E_or = 1,                           # オッズ比は1（対照群）
#   E_logor = 0,                        # ログオッズ比は0（対照群）
#   E_se = 0                            # 標準誤差は0（対照群）
# )
# 
# # 既存のデータフレームと結合
# df <- rbind(referent_row, df)


# データを試験ごとにグループ化し、各試験内で最小用量を基準とする処理

df <- df %>%
  # まず、各試験でDistributionの異なる値の数を計算
  group_by(Study_ID) %>%
  mutate(dist_unique_count = n_distinct(Distribution_in_group)) %>%
  ungroup() %>%
  # 異なる用量が2つ以上ある試験のみを残す
  filter(dist_unique_count >= 2) %>%
  # 不要になった列を削除
  select(-dist_unique_count) %>%
  # ここから相対値の計算
  group_by(Study_ID) %>%
  mutate(
    # 各群のオッズを計算
    odds = (n_infection) / (n_group - n_infection),
    # 最小Distribution_in_groupの行のオッズをリファレンスとして使用
    ref_odds = odds[which.min(Distribution_in_group)],
    # オッズ比を計算
    E_or = odds / ref_odds,
    # log(OR)を計算
    E_logor = log(E_or),
    # 標準誤差の計算 (1/a + 1/b + 1/c + 1/d の平方根)
    E_se = sqrt(1/n_infection + 1/(n_group-n_infection) + 
                  1/n_infection[which.min(Distribution_in_group)] + 
                  1/(n_group-n_infection)[which.min(Distribution_in_group)])
  ) %>%
  ungroup() %>%
  # reference行の値を1, 0, NAに設定
  mutate(
    E_or = if_else(Distribution_in_group == ave(Distribution_in_group, Study_ID, FUN = min), 1, E_or),
    E_logor = if_else(Distribution_in_group == ave(Distribution_in_group, Study_ID, FUN = min), 0, E_logor),
    E_se = if_else(Distribution_in_group == ave(Distribution_in_group, Study_ID, FUN = min), NA_real_, E_se)
  ) %>%
  # 作業用の中間列を削除
  select(-odds, -ref_odds)


df

# ------------------------------------------
# 4. メインの分析: 2値アウトカムデータ
# ------------------------------------------
# ノットをデータの範囲内に設定
knots <- quantile(df$Distribution_in_group, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

# dosresmeta モデルの構築
mod1 <- dosresmeta(
  formula = E_logor ~ rcs(Distribution_in_group, knots),
  type = 'cc',
  id = Study_ID,
  se = E_se,
  cases = n_infection,
  n = n_group,
  data = df[!is.na(df$n_infection), ],
  method = "ml",
  proc = "1stage"
)
summary(mod1)

# 用量反応曲線をプロット
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

plot_dosres(mod1, ylab = "Infection Risk (OR)", ylim = c(0.5, 3))

# VPCのプロット
plot_vpc <- function(data, mod, se){
  plot(data$Distribution_in_group[!is.na(data[,se])], vpc(mod), xlab = "Distribution in Group")
  lines(lowess(data$Distribution_in_group[!is.na(data[,se])], vpc(mod)))
}

plot_vpc(df, mod1, "E_se")

# 指定された効果の量での投与量を推定する関数を修正
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


# doseEff を使って推定
p <- c(.5, .95, 1)
newdata <- data.frame(Distribution_in_group = seq(0, 6, length.out = 100))
edp <- with(predict(mod1, newdata), doseEff(p = p, Distribution_in_group = newdata$Distribution_in_group, Ep = pred, trunc = FALSE))
edp %>%
  mutate(ExpEp = exp(Ep))->edp
print(round(edp, 4))
