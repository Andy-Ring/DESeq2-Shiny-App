
#Convert ID's Function
convertIDs <- function(counts, session = NULL) {
  tryCatch({
    rownames(counts) <- sub("\\..*", "", rownames(counts))
    anno <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = rownames(counts),
                                  columns = c("ENSEMBL", "SYMBOL"),
                                  keytype = "ENSEMBL")
    counts <- cbind(ENSEMBL = rownames(counts), counts)
    counts1 <- dplyr::left_join(as.data.frame(counts), anno, by = "ENSEMBL")
    counts1$SYMBOL <- ifelse(is.na(counts1$SYMBOL), counts1$ENSEMBL, counts1$SYMBOL)
    counts1 <- counts1[!duplicated(counts1$SYMBOL), ]
    rownames(counts1) <- counts1$SYMBOL
    counts1 <- counts1[, !(names(counts1) %in% c("ENSEMBL", "SYMBOL"))]
    return(as.matrix(counts1))
  }, error = function(e) {
    if (!is.null(session)) {
      showNotification("ENSEMBL IDs not detected, skipping conversion!", type = "warning", duration = 5)
    } else {
      message("ENSEMBL IDs not detected, skipping conversion!")
    }
    return(counts)
  })
}