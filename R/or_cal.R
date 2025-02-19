# 必要なパッケージを読み込み
library(dosresmeta)
library(rms)
library(dplyr)
library(openxlsx)
library(gtsummary)

# Excelファイルからデータを読み込む
df <- openxlsx::read.xlsx("ds_transfusion_infection_001.xlsx")

# オッズ比 (E_or)、ログオッズ比 (E_logor)、標準誤差 (E_se) を計算
df <- df %>%
  mutate(
    E_or = (n_infection / (n_group - n_infection)),   # オッズ比の計算
    E_logor = log(E_or),                              # ログオッズ比の計算
    E_se = sqrt(1/n_infection + 1/(n_group - n_infection)) # 標準誤差の計算
  )

# 計算結果を新しいExcelファイルに書き出し
output_file <- "calculated_odds_logor_se.xlsx"
openxlsx::write.xlsx(df, file = output_file)

# 保存が完了したことを表示
print(paste("Excelファイルが作成されました:", output_file))
