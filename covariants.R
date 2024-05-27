library(tidyverse)
library(rjson)

voc0 <- fromJSON(file="data/perCountryData.json")
voc_us <- voc$regions[[1]]$distributions[[1]]$distribution

tot_seq <- tibble(
    date=ymd(map_chr(voc_us, "week")),
    n=map_int(voc_us, "total_sequences")
)

ggplot(tot_seq, aes(date, n)) +
    geom_point()
