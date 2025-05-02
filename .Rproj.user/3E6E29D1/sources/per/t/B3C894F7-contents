#' @title Split GenBank records
#' @description Internal function to split multi-record GenBank files
#' @keywords internal
split_records <- function(lines) {
  ends <- which(lines == "//")
  starts <- c(1, ends[-length(ends)] + 1)
  lapply(seq_along(ends), function(i) {
    if(starts[i] > ends[i]) character(0) else lines[starts[i]:ends[i]]
  })
}


#' @title Process GenBank record
#' @description Internal parser for individual GenBank records
#' @keywords internal
process_record <- function(record_lines) {
  methods::new("GBSequence",
      LOCUS = parse_section(record_lines, "LOCUS"),
      DEFINITION = parse_section(record_lines, "DEFINITION"),
      ACCESSION = parse_section(record_lines, "ACCESSION"),
      VERSION = parse_section(record_lines, "VERSION"),
      DBLINK = parse_section(record_lines, "DBLINK"),
      KEYWORDS = parse_section(record_lines, "KEYWORDS"),
      SOURCE = parse_source(record_lines),
      ORGANISM = parse_organism(record_lines),
      REFERENCE = parse_references(record_lines),
      COMMENT = parse_comment(record_lines),
      FEATURES = parse_features(record_lines),
      ORIGIN = parse_origin(record_lines)
  )
}


#' @title Parse GenBank section
#' @description Internal section parser for GenBank files
#' @keywords internal
parse_section <- function(lines, section) {
  section_lines <- grep(paste0("^\\s{0,2}", section, "\\s"), lines, value = TRUE)
  if(length(section_lines) > 0) {
    gsub("\\s+", " ", sub(paste0("^\\s*", section, "\\s+"), "", section_lines[1]))
  } else ""
}


#' @title Parse GenBank source
#' @description Internal source parser for GenBank files
#' @keywords internal
parse_source <- function(lines) {
  source_line <- grep("^SOURCE", lines, value = TRUE)
  if(length(source_line) > 0) {
    sub("^SOURCE\\s+", "", source_line[1])
  } else ""
}


#' @title Parse GenBank organism
#' @description Internal organism parser for GenBank files
#' @keywords internal
parse_organism <- function(lines) {
  org_line <- grep("^\\s+ORGANISM", lines, value = TRUE)
  if(length(org_line) > 0) {
    organism <- sub("^\\s+ORGANISM\\s+", "", org_line[1])
    taxonomy <- paste(grep("^\\s+\\S+;", lines, value = TRUE), collapse = " ")
    paste(organism, taxonomy, sep = "\n")
  } else ""
}


#' @title Parse GenBank reference
#' @description Internal reference parser for GenBank files
#' @keywords internal
parse_references <- function(lines) {
  ref_starts <- grep("^REFERENCE", lines)
  lapply(ref_starts, function(start) {
    end <- grep("^\\s{0,2}[A-Z]", lines[-(1:start)], invert = TRUE)[1] + start - 1
    if(is.na(end)) end <- length(lines)
    ref_lines <- lines[start:end]

    methods::new("GBReference",
        number = sub("^REFERENCE\\s+(\\d+).*", "\\1", ref_lines[1]),
        authors = paste(grep("^\\s+AUTHORS", ref_lines, value = TRUE) |>
                          sub("^\\s+AUTHORS\\s+", "", x = _), collapse = " "),
        title = paste(grep("^\\s+TITLE", ref_lines, value = TRUE) |>
                        sub("^\\s+TITLE\\s+", "", x = _), collapse = " "),
        journal = paste(grep("^\\s+JOURNAL", ref_lines, value = TRUE) |>
                          sub("^\\s+JOURNAL\\s+", "", x = _), collapse = " "),
        pubmed = paste(grep("^\\s+PUBMED", ref_lines, value = TRUE) |>
                         sub("^\\s+PUBMED\\s+", "", x = _), collapse = " "))
  })
}




#' @title Parse GenBank features
#' @description Internal feature parser for GenBank files
#' @keywords internal
parse_features <- function(lines) {
  feat_lines <- get_section_lines(lines, "FEATURES")
  features <- list()
  current_feat <- NULL

  for(line in feat_lines) {
    if(grepl("^\\s{5}\\w+", line)) {
      if(!is.null(current_feat)) {
        features[[current_feat$type]] <- c(features[[current_feat$type]],
                                           methods::new("GBFeature",
                                               type = current_feat$type,
                                               location = current_feat$location,
                                               qualifiers = current_feat$qualifiers))
      }
      current_feat <- list(
        type = sub("^\\s{5}(\\w+).*", "\\1", line),
        location = sub("^\\s{5}\\w+\\s+(.*)", "\\1", line),
        qualifiers = list()
      )
    } else if(grepl("^\\s{21}/", line)) {
      qual <- sub("^\\s+/", "", line)
      key <- sub("=.*", "", qual)
      value <- sub('.*?"(.*?)"', "\\1", qual)
      current_feat$qualifiers[[key]] <- c(current_feat$qualifiers[[key]], value)
    }
  }
  if(!is.null(current_feat)) {
    features[[current_feat$type]] <- c(features[[current_feat$type]],
                                       methods::new("GBFeature",
                                           type = current_feat$type,
                                           location = current_feat$location,
                                           qualifiers = current_feat$qualifiers))
  }
  features
}


#' @title parse_origin
#' @description Internal helper for origin parsing
#' @keywords internal
parse_origin <- function(lines) {
  origin_lines <- get_section_lines(lines, "ORIGIN")[-1]
  paste0(gsub("[^a-zA-Z]", "", unlist(strsplit(origin_lines, "\\s+"))),
         collapse = "")
}


#' @title parse_comment
#' @description Internal helper for comment parsing
#' @keywords internal
parse_comment <- function(lines) {
  comment_lines <- get_section_lines(lines, "COMMENT")
  paste(sub("^COMMENT\\s+", "", comment_lines), collapse = "\n")
}


#' @title Get section lines
#' @description Internal helper for section extraction
#' @keywords internal
get_section_lines <- function(lines, section) {
  start <- grep(paste0("^\\s{0,2}", section, "\\b"), lines)
  if(length(start) == 0) return(character(0))
  next_section <- grep("^\\s{0,2}[A-Z]", lines[-(1:start[1])])[1]
  end <- if(!is.na(next_section)) start[1] + next_section - 2 else length(lines)
  lines[start[1]:end]
}
