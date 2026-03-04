# R/fct_channels.R
# Helpers for reading / harmonising flow cytometry channel metadata
# (NAME/DESC mapping + scatter channel detection)

#' Extract parameter table (NAME/DESC) from a flowFrame
#'
#' @param ff A `flowCore::flowFrame`
#' @return A data.frame with columns `name` and `desc`
#' @noRd

# Read parameter table from fcs file
param_table <- function(ff) {
  pd <- flowCore::pData(flowCore::parameters(ff))
  data.frame(
    name = as.character(pd$name),
    desc = {x <- as.character(pd$desc); x[is.na(x)] <- ""; x},
    stringsAsFactors = FALSE
  )
}
#' Build consistent NAME <-> label maps across multiple flowFrames
#'
#' Uses the union of channel NAMEs across frames. For each NAME, chooses the
#' most common DESC (case-insensitive). If DESC is missing, falls back to NAME.
#' Duplicate display labels are disambiguated by appending ` [NAME]`.
#'
#' @param frames List of `flowCore::flowFrame`
#' @return A list with `union_names`, `name_to_label`, and `label_to_name`
#' @noRd
build_name_label_maps <- function(frames) {
  tabs <- lapply(frames, param_table)
  union_names <- Reduce(union, lapply(tabs, `[[`, "name"))

  # collect candidate descs per name across files
  desc_map <- lapply(union_names, function(nm) {
    v <- unlist(lapply(tabs, function(tb) tb$desc[tb$name == nm]), use.names = FALSE)
    v[!is.na(v) & nzchar(v)]
  })
  names(desc_map) <- union_names

  # majority label (case-insensitive), fallback to name
  pick_desc <- function(v) {
    if (!length(v)) return(NA_character_)
    key <- tolower(v)
    v[ match(names(sort(table(key), decreasing = TRUE))[1], key) ]
  }
  chosen  <- vapply(desc_map, pick_desc, character(1))
  display <- ifelse(is.na(chosen) | !nzchar(chosen), union_names, chosen)

  # disambiguate duplicate labels by appending [NAME] (i.e. if raw file contains labels for fluorescent channels)
  dups <- duplicated(display) | duplicated(display, fromLast = TRUE)
  display[dups] <- paste0(display[dups], " [", union_names[dups], "]")

  list(
    union_names   = union_names,
    name_to_label = setNames(display, union_names),
    label_to_name = setNames(union_names, display)
  )
}

#' Identify Forward Scatter (FSC) channel names
#'
#' Searches a vector of channel names and returns the first match corresponding
#' to either Forward Scatter Area (FSC-A) or Forward Scatter Height (FSC-H).
#' Common separator variations such as `"FSC-A"`, `"FSC.A"`, `"FSC_A"`,
#' `"FSC-Area"`, and `"FSC.Height"` are supported.
#'
#' If no matching channel is found, `NA_character_` is returned.
#'
#' @param cn Character vector of channel names (e.g., `colnames(flowFrame)`).
#' @param which Character scalar specifying which FSC measurement to return:
#'   `"A"` for area (FSC-A) or `"H"` for height (FSC-H).
#'
#' @return Character scalar giving the matched FSC channel name, or
#'   `NA_character_` if no suitable channel is found.
#' @noRd
# find forward scatter area and height channels
pick_fsc <- function(cn, which = c("A","H")) {
  which <- match.arg(which)
  pat <- if (which == "A") "^FSC[._-]?(A|Area)$" else "^FSC[._-]?(H|Height)$"
  hit <- cn[grepl(pat, cn, ignore.case = TRUE)]
  if (length(hit)) hit[1] else NA_character_
}
