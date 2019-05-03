
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
                                                  margin = margin(b = 1, t = 0, 0, 0)),
                        legend.text = element_text(size = 10),
                        axis.text = element_text(size = 10, color = "black"),
                        axis.title.y = element_text(angle = 90,
                                                    margin = margin(0,15,0,0)),
                        axis.title.x = element_text(margin = margin(15,0,0,0))))

}



#' Color palette.
#'
#' @export
#' @noRd
#'
color_pal <- function() {
    # < order of colors: green, red, purple, pink, light green, yellow >
    # rows correspond to `RColorBrewer::brewer.pal(6, "Dark2")`
    rgb_mat <- rbind(c(27,158,119), c(217,95,2), c(117,112,179),
                     c(231,41,138), c(102,166,30), c(230,171,2))
    # Switching detritus with plant and detritivore with herbivore
    # Order is now "detritivore", "herbivore", "predator", "soil", "detritus", "plant"
    rgb_mat <- rgb_mat[c(2, 5, 3, 6, 4, 1),]
    # Multipliers for each color, < 1 makes it darker
    rgb_mults <- c(0.5, 1.3, 1.1,
                   1.05, 0.9, 1)
    cp <- apply(rgb_mat * matrix(rgb_mults, 6, 3), 1,
                function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
    return(cp)
}
