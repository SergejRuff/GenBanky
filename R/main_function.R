#' Read GenBank Files with Single or Multiple Records
#'
#' @description
#' Reads GenBank format files (.gb/.gbk) containing either single or multiple sequence records.
#' Modified from the original geneviewer implementation to handle multi-record files and
#' store all data in structured S4 objects.
#'
#' @details
#' This enhanced version handles both use cases:
#' \itemize{
#'   \item Single files containing one sequence record
#'   \item Files containing multiple sequence records (separated by //)
#'   \item Directories containing multiple GenBank files
#' }
#'
#' Key improvements over original geneviewer::read_gbk():
#' \itemize{
#'   \item Full parsing of all GenBank sections (FEATURES, ORIGIN, REFERENCES)
#'   \item Structured S4 object output maintaining original file organization
#'   \item Automatic handling of duplicate accession numbers across files
#'   \item Preservation of multi-record file structure when processing directories
#' }
#'
#' @param path Path to either:
#' \itemize{
#'   \item A directory containing GenBank files (.gb/.gbk)
#'   \item A single GenBank file (with one or multiple records)
#' }
#'
#' @return Returns a \code{GenBankyObj} S4 object containing:
#' \itemize{
#'   \item \code{sequences}: Named list of \code{GBSequence} objects (named by accession)
#'   \item Attribute \code{changed_accessions}: Named vector showing accession modifications
#' }
#'
#' @examples
#' \dontrun{
#'
#' example_gbk_path <- system.file("extdata", "flavivirus.gb", package = "GenBanky")
#' imported_data <- GenBanky(example_gbk)
#' }
#'
#'
#' @author
#' Sergej Ruff
#'
#' @export
GenBanky <- function(path) {
  if (dir.exists(path)) {
    files <- list.files(path, pattern = "\\.gbk$|\\.gb$", full.names = TRUE)
    if (length(files) == 0) stop("No GenBank files found in directory")
    records <- unlist(lapply(files, function(f) split_records(readLines(f))), recursive = FALSE)
  } else if (file.exists(path)) {
    records <- split_records(readLines(path))
  } else {
    stop("Path does not exist")
  }

  sequences <- lapply(records, process_record)
  accessions <- sapply(sequences, function(x) x@ACCESSION)


  changed_acc <- character(0)
  if(any(duplicated(accessions))) {
    dupes <- accessions[duplicated(accessions) | duplicated(accessions, fromLast = TRUE)]
    unique_dupes <- unique(dupes)


    msg <- paste("Found duplicate accessions:\n",
                 paste("-", unique_dupes, collapse = "\n"),
                 "\nAppending numbers to make unique:")


    unique_accessions <- make.unique(accessions)
    changes <- unique_accessions != accessions
    changed_acc <- setNames(unique_accessions[changes], accessions[changes])


    change_msgs <- sapply(seq_along(changed_acc), function(i) {
      sprintf("Changed '%s' to '%s'", names(changed_acc)[i], changed_acc[i])
    })

    warning(paste(msg, paste(change_msgs, collapse = "\n"), sep = "\n"))
    accessions <- unique_accessions
  }

  names(sequences) <- accessions
  gb_collection <- new("GenBankyObj", sequences = sequences)


  attr(gb_collection, "changed_accessions") <- changed_acc
  return(gb_collection)
}
