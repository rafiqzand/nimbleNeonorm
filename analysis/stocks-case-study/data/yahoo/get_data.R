dir.create("analysis/stocks-case-study/data/yahoo", recursive = TRUE, showWarnings = FALSE)

for (tk in c("ADRO", "PTBA")) {
  px <- quantmod::getSymbols(paste0(tk, ".JK"), src = "yahoo",
                             from = "2010-01-01", to = "2019-12-31",
                             auto.assign = FALSE)
  px <- data.frame(date = zoo::index(px), close = as.numeric(quantmod::Cl(px)))
  write.csv(px, file.path("analysis/stocks-case-study/data/yahoo",
                          paste0(tk, ".csv")), row.names = FALSE)
  cat(tk, ":", nrow(px), "rows\n")
}
