#' GenBank Feature Class
#'
#' @description
#' S4 class representing a GenBank feature with type, location, and qualifiers
#'
#' @slot type Character specifying feature type (e.g., "CDS", "gene", "mRNA")
#' @slot location Character describing genomic location (e.g., "1..500")
#' @slot qualifiers List of key-value pairs with feature attributes
#'
#' @author Sergej Ruff
#' @keywords internal
#' @aliases GBFeature-class
setClass("GBFeature",
         slots = list(
           type = "character",
           location = "character",
           qualifiers = "list"
         ))


#' GenBank Reference Class
#'
#' @description
#' S4 class representing a GenBank reference entry
#'
#' @slot number Character reference index
#' @slot authors Character vector of authors
#' @slot title Character reference title
#' @slot journal Character journal information
#' @slot pubmed Character PubMed ID
#'
#' @author Sergej Ruff
#' @keywords internal
#' @aliases GBReference-class
setClass("GBReference",
         slots = list(
           number = "character",
           authors = "character",
           title = "character",
           journal = "character",
           pubmed = "character"
         ))


#' GenBank Sequence Class
#'
#' @description
#' S4 class representing a complete GenBank record
#'
#' @slot LOCUS Character locus information
#' @slot DEFINITION Character sequence definition
#' @slot ACCESSION Character accession number
#' @slot VERSION Character version information
#' @slot DBLINK Character database links
#' @slot KEYWORDS Character keywords
#' @slot SOURCE Character source organism
#' @slot ORGANISM Character organism taxonomy
#' @slot REFERENCE List of GBReference objects
#' @slot COMMENT Character comments
#' @slot FEATURES List of GBFeature objects
#' @slot ORIGIN Character DNA sequence
#'
#' @author Sergej Ruff
#' @keywords internal
#' @aliases GBSequence-class
setClass("GBSequence",
         slots = list(
           LOCUS = "character",
           DEFINITION = "character",
           ACCESSION = "character",
           VERSION = "character",
           DBLINK = "character",
           KEYWORDS = "character",
           SOURCE = "character",
           ORGANISM = "character",
           REFERENCE = "list",
           COMMENT = "character",
           FEATURES = "list",
           ORIGIN = "character"
         ))


#' GenBank Collection Class
#'
#' @description
#' S4 class containing multiple GenBank records
#'
#' @slot sequences List of GBSequence objects, named by accession numbers
#'
#' @details
#' Contains attribute "changed_accessions" showing any modified accession names
#' when duplicates were present. Access with \code{attr(object, "changed_accessions")}
#'
#' @author Sergej Ruff (modified from geneviewer package)
#' @keywords internal
#' @aliases GenBankyObj-class
setClass("GenBankyObj",
         slots = list(
           sequences = "list"
         ))
