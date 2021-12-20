theme_jk <- function (base_size = 11, base_family = "", base_line_size = base_size/22, 
          base_rect_size = base_size/22) 
{
  theme_bw(base_size = base_size, base_family = base_family, 
             base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(
      axis.text.x = element_text(angle = 90,vjust = 0.5),
      strip.background = element_rect(fill = "white"),
      panel.background  = element_blank(),
      legend.key = element_rect(fill="transparent", colour=NA)
      
      
    )
}
