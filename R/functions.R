
#' @export

# functions

# pre-process functions
preprocess_2 <- function(ff){
  ff <- flowCore::compensate(ff, flowCore::keyword(ff)$SPILL)
  ff <- flowCore::transform(ff,
                            flowCore::estimateLogicle(ff,
                                                      colnames(flowCore::keyword(ff)$SPILL)))
  ff
}

channel_select <- function(ff){
  library(data.table)
  ff = ff[,c((flowCore::colnames(ff) %like% ".A") | flowCore::colnames(ff) %like% "Time") | flowCore::colnames(ff) %like% "FSC.H" | flowCore::colnames(ff) %in% marker_panel$Marker]

  ff
}
# Functions for .fsc data importing and mapping to metadata

# Below sections are to allow for cbind.fill -- package origin had depreciated
cbind_fill<-function(...,fill=NULL)
{
  inputs<-list(...)
  inputs<-lapply(inputs,vert)
  maxlength<-max(unlist(lapply(inputs,len)))
  bufferedInputs<-lapply(inputs,buffer,length.out=maxlength,fill,preserveClass=FALSE)
  return(Reduce(cbind.data.frame,bufferedInputs))
}

vert<-function(object)
{
  #result<-as.data.frame(cbind(as.matrix(object)))
  if(is.list(object))
    object<-cbind(object)
  object<-data.frame(object)

  return(object)
}

len <- function(data)
{
  result<-ifelse(is.null(nrow(data)),length(data),nrow(data))
  return(result)
}

buffer<-function(x,length.out=len(x),fill=NULL,preserveClass=TRUE)
{
  xclass<-class(x)
  input<-lapply(vert(x),unlist)
  results<-as.data.frame(lapply(input,rep,length.out=length.out))
  if(length.out>len(x) && !is.null(fill))
  {
    results<-t(results)
    results[(length(unlist(x))+1):length(unlist(results))]<-fill
    results<-t(results)
  }
  if(preserveClass)
    results<-as2(results,xclass)
  return(results)
}

# Now process each with the preprocess function which applies the set spillover file and normalises with biexponential logicle function
preprocess <- function(ff) {
  ff <- compensate(ff, ff@description$SPILL)
  ff <- transform(ff, flowCore::transformList(colnames(ff@description$SPILL), flowCore::logicleTransform()))
  ff@exprs <- scale(ff@exprs)
  return(ff)
}

# function to convert wellIDs to match those attached to flowjo objects
checkWellNames = function(wellNames) {
  #' Check and convert well names to the appropriate format
  #' E.g. A1 -> A01
  #'
  #' @param wellNames String vector (e.g. c("A1", "A2", "B3", "B4"))
  #'
  #' @return String vector
  #'
  o = wellNames
  rows = substr(wellNames, 1, 1)
  stopifnot(all(rows %in% toupper(letters)[1:8]))
  columns = as.integer(substr(wellNames, 2, 10))
  stopifnot(all(columns >= 1 & columns <= 12))
  columns = as.character(columns)
  columns = sapply(columns, function(x) {
    if (nchar(x) == 1) {
      return(paste0("0", x))
    } else {
      return(x)
    }
  })
  return(paste0(rows, columns))
}

# Function to extract labels from fcs file paired with workspace --
workspace_cell_labels <- function(flowset, workspace, cell_types=cell_types_fixed, ws_group){
  groups = fj_ws_get_sample_groups(open_flowjo_xml(workspace))#, execute = FALSE)
  group_pos = groups[match(ws_group, groups$groupName), 2] + 1 # returns the group ID for the matched group +1 as the numbering begins at 0
  gating_manual = GetFlowJoLabels(flowset,
                                  workspace,
                                  group = group_pos,
                                  cellTypes = cell_types,
                                  get_data = TRUE)
  manual_labels = do.call("c",lapply(gating_manual,function(x){x$manual}))
  rm(gating_manual)
  return(manual_labels)
}


# wiht this the condition is a on points that deviate from FSC.H ~ FSC.A linearly, plus removing debrit (v. small vals in both)

# Function to remove data 5% from the line
remove_outliers <- function(x, y, p = 0.05, slope, intercept) {
  # Get the predicted values
  predicted <- slope * x + intercept

  # Calculate the residuals
  residuals <- y - predicted

  # Calculate the threshold for outliers
  threshold <- quantile(abs(residuals), 1 - p)

  # Identify the outliers
  outliers <- which(abs(residuals) > threshold)

  # Return the new data without outliers
  return(outliers)
}
