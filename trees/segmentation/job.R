## execute segmentation in parallel with batches of stands
source("trees/segmentation/setup.R")
jobStartTime = Sys.time()
setwd(file.path(getwd(), "../.."))

workerNumber = 1

chunkSize = 100 # stands
standIDsToSegment = elliottStands$STD_ID[seq(1:100) + chunkSize * (workerNumber - 1)]
standsSegmented = 0
for (standID in standIDsToSegment)
{
  standStartTime = Sys.time()
  process_stand(standID)
  standsSegmented = standsSegmented + 1
  cat(sprintf("%d: %.1f s total (%d of %d).", standID, standsSegmented, length(standIDsToSegment), difftime(Sys.time(), standStartTime, units = "secs")), fill = TRUE)
}
