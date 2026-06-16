# .onLoad <- function(libname, pkgname){
#   swapStars()
# }

.onAttach <- function(...) {
  s1 <- " ____  _                                   ____     \n"
  s2 <- "/ ___|| | ___   _ ___  ___ __ _ _ __   ___|  _  \\    \n"
  s3 <- "\\ __ \\  |/ / | | / __|/ __/ _` | '_ \\ / _ \\ |_) | \n"
  s4 <- " __)  |   <| |_| \\ _ \\ (_| (_| | |_) |  __/  _ <  \n"
  s5 <- paste0("|____/|_|\\_\\ __  |___/\\ __\\ _,_| .__/ \\ __|_| \\ |  v", packageVersion("skyscapeR"), "\n")
  s6 <- "             __/ |             | |               \n"
  s7 <- "            |___/              |_|               \n"
  packageStartupMessage(paste("\n\nWelcome to\n", s1, s2, s3, s4, s5, s6, s7, sep=""))
  swapStars()
}
