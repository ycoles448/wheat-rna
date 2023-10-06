# Override for plot_venn
# Need labels with backgrounds
plot_venn <- function (x, show_intersect, set_color, set_size, label, label_geom,
    label_alpha, label_color, label_size, label_percent_digit,
    label_txtWidth, edge_lty, edge_size, ...)
{
    venn <- Venn(x)
    data <- process_data(venn)
    p <- ggplot()
    region.params <- list(data = data@region, mapping = aes_string(fill = "count"))
    edge.params <- list(data = data@setEdge, mapping = aes_string(color = "id"),
        show.legend = FALSE)
    if (utils::packageVersion("ggplot2") >= "3.4.0") {
        edge.params$linetype <- edge_lty
        edge.params$linewidth <- edge_size
    }
    else {
        edge.params$lty <- edge_lty
        edge.params$size <- edge_size
    }
    text.params <- list(data = data@setLabel, mapping = aes_string(label = "name"),
        size = set_size, inherit.aes = FALSE)
    region.layer <- do.call("geom_sf", region.params)
    edge.layer <- do.call("geom_sf", edge.params)
    text.layer <- do.call("geom_sf_label", text.params)
    p <- p + region.layer + edge.layer + text.layer + theme_void()
    if (label != "none" & show_intersect == FALSE) {
        region_label <- data@region %>% dplyr::filter(.data$component ==
            "region") %>% dplyr::mutate(percent = paste(round(.data$count *
            100/sum(.data$count), digits = label_percent_digit),
            "%", sep = "")) %>% dplyr::mutate(both = paste(.data$count,
            paste0("(", .data$percent, ")"), sep = "\n"))
        if (label_geom == "label") {
            p <- p + geom_sf_label(aes_string(label = label),
                data = region_label, alpha = label_alpha, color = label_color,
                size = label_size, lineheight = 0.85, label.size = NA)
        }
        if (label_geom == "text") {
            p <- p + geom_sf_text(aes_string(label = label),
                data = region_label, alpha = label_alpha, color = label_color,
                size = label_size, lineheight = 0.85)
        }
    }
    if (show_intersect == TRUE) {
        items <- data@region %>% dplyr::rowwise() %>% dplyr::mutate(text = yulab.utils::str_wrap(paste0(.data$item,
            collapse = " "), width = label_txtWidth)) %>% sf::st_as_sf()
        label_coord = sf::st_centroid(items$geometry) %>% sf::st_coordinates()
        p <- ggplot(items) + geom_sf(aes_string(fill = "count")) +
            geom_sf_text(aes_string(label = "name"), data = data@setLabel,
                inherit.aes = F) + geom_text(aes_string(label = "count",
            text = "text"), x = label_coord[, 1], y = label_coord[,
            2], show.legend = FALSE) + theme_void()
        ax <- list(showline = FALSE)
        p <- plotly::ggplotly(p, tooltip = c("text")) %>% plotly::layout(xaxis = ax,
            yaxis = ax)
    }
    p
}

# Override for ggVennDiagram
function (x, category.names = names(x), show_intersect = FALSE,
    set_color = "black", set_size = NA, label = c("both", "count",
        "percent", "none"), label_alpha = 0.5, label_geom = c("label",
        "text"), label_color = "black", label_size = NA, label_percent_digit = 0,
    label_txtWidth = 40, edge_lty = "solid", edge_size = 1, ...)
{
    if (!is.list(x)) {
        stop(simpleError("ggVennDiagram() requires at least a list."))
    }
    names(x) <- category.names
    dimension <- length(x)
    label <- match.arg(label)
    ## label_geom <- match.arg(label_geom)
    if (dimension <= 7) {
        plot_venn(x, show_intersect = show_intersect, set_color = set_color,
            set_size = set_size, label = label, label_alpha = label_alpha,
            label_geom = label_geom, label_color = label_color,
            label_size = label_size, label_percent_digit = label_percent_digit,
            label_txtWidth = label_txtWidth, edge_lty = edge_lty,
            edge_size = edge_size, ...)
    }
    else {
        stop("Only support 2-7 dimension Venn diagram.")
    }
}
