

# MAss_rt_filterfunction_working hopefully --------------------------------
# input_data = x
mass_rt_filter_20200309 <- function(input_data) {
  
  x1 <- input_data %>% arrange(D_m_z, C_rt_Start)
  x2 <- x1 %>% pull(D_m_z)
  x3 <- x1 %>% pull(C_rt_Start)
  x4 <- x1 %>% pull(C_rt_End)
  y = 1
  # y = 30
  dsammelliste <- tibble()
  while (y <= (length(x2))) {
    # selection for masses:
    c <- c()
    masss <- x2[y]
    while (abs(x2[y] - masss) <= delta_mz2 & y <= (length(x2))) {
      c[y] <- y
      y = y + 1
    }
    y = y - 1
    c <- na.omit(c)
    # y;c
    # if more than one entry per similar mass, then check RT overlaps... by increasing a rt window:
    if (length(c) > 1) {
      x5 <-
        x1 %>% slice(c) %>% arrange(C_rt_Start) # contains entries with similar mass.
      
      z = 1
      repeat {
        # while (z<nrow(x5)) {
        w <- z
        x6 <- x5 %>% slice(z)
        start <- x6$C_rt_Start
        end <- x6$C_rt_End
        b <- x5 %>% slice(z+1) %>% pull(C_rt_Start)
        x7 <- x5 %>% slice(z+1) %>% pull(C_rt_End)
        while (all((b >= start & b <= end), z <= nrow(x5))) {
          end <- max(end, x7)
          z = z + 1
          b <- x5 %>% slice(z) %>% pull(C_rt_Start)
          x7 <- x5 %>% slice(z) %>% pull(C_rt_End)
        }
        if (length(w:z)>1) {
          z = z - 1
          sammelliste1 <-
            x5 %>% slice(w:z) %>%  mutate(E_Group = c[z],
                                          E_Start = start,
                                          E_End = end) # need to report this part as group:
          dsammelliste <-
            bind_rows(dsammelliste, sammelliste1) # pro massenfenster wenn duplikate vorhanden.
          z = z + 1  
        } else{z = z + 1}
        
        # print(z)
        if (z >= nrow(x5)) {
          break
        }
      }
      y = y + 1
    } else{ y = y + 1}
    # print(y)
  }
  dsammelliste
}
# 2nd part doing the data thing.