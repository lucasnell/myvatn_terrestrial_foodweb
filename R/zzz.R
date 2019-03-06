
#' Set ggplot2 theme on loading of package
#'
#' @noRd
#'
#' @importFrom ggplot2 theme_set
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 %+replace%
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 margin
#' @importFrom ggplot2 element_text
#'
.onLoad <- function(libname, pkgname) {

    theme_set(theme_classic() %+replace%
                  theme(panel.grid = element_blank(),
                        strip.background = element_blank(),
                        legend.margin = margin(0,0,0,0),
                        strip.text = element_text(size = 12,
                                                  margin = margin(b = 1, t = 0, 0, 0),
                                                  face = "bold"),
                        legend.text = element_text(size = 10),
                        axis.text = element_text(size = 10, color = "black"),
                        axis.title.y = element_text(angle = 90,
                                                    margin = margin(0,15,0,0)),
                        axis.title.x = element_text(margin = margin(15,0,0,0))))

}
