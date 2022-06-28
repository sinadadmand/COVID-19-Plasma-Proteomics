library(ggplot2)
library(dplyr)


library(readxl)
data2 <- gene_overlap_analysis_results <- read_excel("C:/Users/sinad/Dropbox/aKoc/aProject/Covid/COMBAT validation part/gene overlap analysis results.xlsx")
data2 <- filter(data2, Overlaps_with_COMBAT!= 'Total without steiner')
data2 <- filter(data2, Overlaps_with_COMBAT!= 'community5')
data2 <- filter(data2, Overlaps_with_COMBAT!= 'Total with all first neighbors')


ggplot(data2, aes(x=reorder(Overlaps_with_COMBAT, -adj_p_value), y=x, size = -1*log(adj_p_value,10), color = overlap_percentage*100)) +
geom_point(alpha=0.7) +
scale_size(range = c(4, 25)) + scale_colour_gradient(low = "#76c8c8", high = "red2")



