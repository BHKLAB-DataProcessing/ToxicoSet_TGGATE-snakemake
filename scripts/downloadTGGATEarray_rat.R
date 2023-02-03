require(downloader)
library(curl)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
download_dir <- args[1]

curl_download(
  "https://orcestradata.blob.core.windows.net/toxico/TGGATES_rat_CEL.tar.gz", 
  destfile =file.path(download_dir, "TGGATES_rat_CEL.tar.gz")
)
# untar(file.path(download_dir, "TGGATES_rat_CEL.tar.gz"), exdir = download_dir)

# file.remove(file.path(download_dir, "TGGATES_rat_CEL.tar.gz"))
