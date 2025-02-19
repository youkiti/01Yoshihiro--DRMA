library(tidyverse)

# データの読み込み
ds0121 <- read.csv("ds0121.csv")

# 変換対象の列名のプレフィックスを指定
prefixes <- c("Distribution_in_group", "n_group", "n_infection", "n_SepsisSepticshock", "n_Bacteremia", "n_AnyInfection", "n_OtherInfection")

# データの変換とNAの削除、Study_idの保持
ds0121_long <- ds0121 %>%
  pivot_longer(
    cols = starts_with(prefixes),
    names_to = c(".value", "Group"),
    names_pattern = "(.*)(\\d)"
  ) %>%
  filter(!is.na(Distribution_in_group)) %>%
  select(Study_id, everything()) %>%
  select(-Group)


# infection_types を定義
infection_types <- c("n_infection", "n_SepsisSepticshock", "n_Bacteremia", "n_AnyInfection", "n_OtherInfection")

# ループでデータセットを作成
for (infection_type in infection_types) {
  # infection_type に対応するデータセットを作成
  subset_data <- ds0121_long %>%
    select(Study_id, Distribution_in_group, n_group, all_of(infection_type)) %>%
    filter(!is.na(!!sym(infection_type))) %>%
    rename(n_infection = all_of(infection_type)) %>%
  mutate(
    E_or = (n_infection / (n_group - n_infection)),   # オッズ比の計算
    E_logor = log(E_or),                              # ログオッズ比の計算
    E_se = sqrt(1/n_infection + 1/(n_group - n_infection)) # 標準誤差の計算
  )
  
  # データセット名を作成
  dataset_name <- paste0("dataset_", infection_type)
  
  # データセットを環境に割り当て
  assign(dataset_name, subset_data)
  
  # CSVファイルとして保存 (オプション)
  file_name <- paste0(dataset_name, ".csv")
  write.csv(subset_data, file_name, row.names = FALSE)
}

# 作成されたデータセットの確認 (例)
head(dataset_n_infection)