#' @keywords internal
close.txtProgressBar <- function (con, ...) 
{
  con$kill()
  invisible(NULL)
}