#' IsoVision Data
#'
#' @param path.bed The path to the .bed file specifying the feature track.
#' @param path.out The path to the .out file from the CDD with conserved domain annotations.
#' @param path.rank The path to the SHaMsPipe .tab file specifying isoform ranks.
#' @param path.fasta The path to the .fasta file with ORF information in the headers.
#'
#' @return Data ready for plot with IsoVision
isoviz_data <- function(path.bed, path.out, path.rank, path.fasta) {
  require(data.table)
  require(dplyr)

  bed <- fread(path.bed, header = F)
  out <- fread(path.out, header = F)
  rank <- fread(path.rank)
  fasta <- fread(path.fasta, header = F) %>% .[seq(1, nrow(.[]), by = 2)]

  if(bed[, unique(V6)] %>% length == 2) {
    message('.bed file has isoforms on both strands. This is not supported!')
    return(NULL)
  }

  # Normalize names
  bed[, V4 := sapply(strsplit(V4, '_'), '[[', 1)]
  out[, V1 := sapply(strsplit(sapply(strsplit(sapply(strsplit(V1, '>'), '[[', 2), ' '), '[[', 1), '_'), '[[', 1)]
  rank[, isoformName := sapply(strsplit(isoformName, '_'), '[[', 1)]

  fasta[, isoformName := sapply(strsplit(sapply(strsplit(sapply(strsplit(V1, '>'), '[[', 2), ' '), '[[', 1), '_'), '[[', 1)]
  fasta[, c('orf', 'V1') := list(sapply(strsplit(sapply(strsplit(V1, ':'), tail, 1), '\\('), '[[', 1), NULL)]

  starts <- lapply(strsplit(bed[, V12], ','), as.integer)
  sizes <- lapply(strsplit(bed[, V11], ','), as.integer)
  sizes.sum <- lapply(sizes, cumsum)

  blocks <- rbindlist(lapply(mapply(cbind,
                                    setNames(lapply(sizes.sum, function(x) c(0, head(x, -1))), bed[, V4]),
                                    sizes.sum,
                                    starts,
                                    mapply('+', sizes, starts)), as.data.table), idcol = 'isoformName') %>%
    setnames(paste0('V', 1:4), c('start', 'end', 'cStart', 'cEnd')) %>%
    merge(rank[, .(isoformName, total, ratio_total_reads)], by = 'isoformName', all.x = T) %>%
    merge(fasta, by = 'isoformName', all.x = T) %>% setorder(-total)

  for(i in 1:nrow(blocks)) {
    blocks[i, ORF := as.integer(sapply(strsplit(orf, '-'), tail, 1)) >= start & end >= as.integer(sapply(strsplit(orf, '-'), head, 1))]
  }

  blocks[, orf := NULL]

  out[, c('V4', 'V5') := list(3 * V4, 3 * V5)]

  strand.neg <- bed[1, V6] == '-'
  blocks[, out[V2 == 'specific', unique(V9)] := list(F)]
  for(i in which(out[, V2] == 'specific')) {
    #if(strand.neg) {
    #  eend <- tail(blocks[isoformName == out[1, V1], end], 1)
    #  blocks[isoformName == out[i, V1] & out[i, V5] <= eend - start & eend - end >= out[i, V4], eval(quote(out[i, V9])) := 1]
    #} else
      blocks[isoformName == out[i, V1] & out[i, V5] >= start & end >= out[i, V4], eval(quote(out[i, V9])) := T]
  }

  blocks
}
