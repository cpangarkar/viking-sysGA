
library(ggrepel)
assay = assay_list[2]

cdf = df |> filter(assay=="BUN" & species=="feline")

p = 
  cdf |> 
  ggplot(aes(x=predicate, y=vital)) + 
  geom_point() + 
  geom_abline(slope=1, intercept = 0) +
  geom_label(aes(label=disc))+
  scale_x_log10() + 
  scale_y_log10() + 
  #geom_label_repel(aes(label=disc))+
  labs(title = cdf$assay[1])

windows()
ggplot2:::print.ggplot(p)