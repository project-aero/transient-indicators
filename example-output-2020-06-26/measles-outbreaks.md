---
title: "Application of transient indicators to measles outbreak data"
author: "Eamon O'Dea"
date: "6/22/2020"
output: html_document
---



## Data loading

Outbreak data from this paper: 
Hall V, Banerjee E, Kenyon C, et al. Measles Outbreak — Minnesota April–May 2017. MMWR Morb Mortal Wkly Rep 2017;66:713–717. DOI: http://dx.doi.org/10.15585/mmwr.mm6627a1external

and the paper here:
https://www.cdc.gov/mmwr/preview/mmwrhtml/mm6406a5.htm?s_cid=mm6406a5_w



```r
outs <- list()

outs$mn2017 <- 
tribble(~month, ~day, ~onsets,
        3, 30, 1,
        4, 8, 1,
        4, 10, 1,
        4, 11, 4,
        4, 12, 1,
        4, 13, 1,
        4, 18, 1,
        4, 19, 1,
        4, 20, 5,
        4, 21, 2,
        4, 22, 3,
        4, 23, 1,
        4, 24, 3,
        4, 25, 1,
        4, 26, 2,
        4, 27, 1, 
        4, 28, 1,
        4, 29, 1,
        4, 30, 2,
        5, 1, 1,
        5, 2, 4,
        5, 3, 2,
        5, 4, 1,
        5, 5, 1, 
        5, 6, 1,
        5, 7, 4,
        5, 8, 1,
        5, 9, 2,
        5, 10, 3,
        5, 12, 1, 
        5, 13, 1,
        5, 14, 1,
        5, 15, 1,
        5, 16, 2,
        5, 17, 3, 
        5, 19, 1,
        5, 26, 1, 
        5, 27, 1) %>%
  mutate(year = 2017)

outs$ca2014 <- 
  tribble(~month, ~day, ~year, ~onsets,
          12, 27, 2014, 1,
          12, 28, 2014, 1,
          12, 29, 2014, 2,
          12, 30, 2014, 2,
          12, 31, 2014, 6,
          1, 1, 2015, 5,
          1, 2, 2015, 8,
          1, 3, 2015, 4,
          1, 4, 2015, 4,
          1, 5, 2015, 2,
          1, 6, 2015, 2,
          1, 7, 2015, 1,
          1, 8, 2015, 4,
          1, 9, 2015, 0,
          1, 10, 2015, 2,
          1, 11, 2015, 2,
          1, 12, 2015, 6,
          1, 13, 2015, 4,
          1, 14, 2015, 3,
          1, 15, 2015, 6,
          1, 16, 2015, 5,
          1, 17, 2015, 0,
          1, 18, 2015, 3,
          1, 19, 2015, 6,
          1, 20, 2015, 3,
          1, 21, 2015, 1,
          1, 22, 2015, 3,
          1, 23, 2015, 2, 
          1, 24, 2015, 1,
          1, 25, 2015, 2,
          1, 26, 2015, 5,
          1, 27, 2015, 2,
          1, 28, 2015, 2,
          2, 1, 2015, 3,
          2, 2, 2015, 1,
          2, 3, 2015, 3,
          2, 4, 2015, 1,
          2, 7, 2015, 1,
          2, 8, 2015, 1) 

odf <- bind_rows(outs, .id="outbreak")
```

## Analysis


```r
linear_model <- function(df){
  data <- df[1:which.max(df$totcases), ]
  stats::lm(totcases ~ obtime, data = data)
}

react_stat <- function(mod){
  coef(mod)["obtime"]
}

max_amp <- function(df){
  max(df$totcases)
}

rdf <- odf %>% 
  group_by(outbreak) %>% 
  mutate(date = lubridate::make_date(year = year, month = month, day = day),
         totcases = cumsum(onsets),
         obtime = (date - min(date)) / lubridate::ddays(1)) %>%
  nest() %>%
  mutate(lmod = purrr::map(data, linear_model)) %>%
  mutate(stat = purrr::map_dbl(lmod, react_stat)) %>%
  mutate(maxamp = purrr::map_dbl(data, max_amp))

rdf
```

```
## # A tibble: 2 x 5
## # Groups:   outbreak [2]
##   outbreak           data lmod    stat maxamp
##   <chr>    <list<df[,7]>> <list> <dbl>  <dbl>
## 1 mn2017         [38 × 7] <lm>    1.43     65
## 2 ca2014         [39 × 7] <lm>    2.76    110
```

The California outbreak has greater reactivity and maximum amplification.

## Figures for publication


```r
theme_set(new = theme_minimal())

pmeasles <-
  rdf %>% unnest(data) %>% ggplot(aes(x = obtime, y = totcases, color = outbreak)) +
  geom_step() +
  scale_color_manual(
    values = c("#999999", "#E69F00"),
    breaks = c("ca2014", "mn2017"),
    labels = c("California December 2014", "Minnesota April 2017"),
    name = "Outbreak:"
  ) +
  theme(legend.position = "top") +
  labs(x = "Days since first case", y = "Number of measles cases") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

pmeasles
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

```r
ggsave("emergence_measles_analysis.pdf", pmeasles, width = 90, units = "mm")
```

```
## Saving 90 x 178 mm image
```

