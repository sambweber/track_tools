
# The adjust.duplicateTimes function taken from package trip

adjust_duplicate_times <- function (time) {
  dups <- duplicated(time)
  if (any(dups)) {
    time[dups] <- time[dups] + 1
    time <- Recall(time)
  }
  time
}
