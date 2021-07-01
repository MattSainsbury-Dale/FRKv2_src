library("ggplot2")
library("ggpubr")
library("dplyr")
FRK_df <- read.csv("./intermediates/Am_FRK.csv")
georob_df <- read.csv("./intermediates/Am_georob.csv")
# PaC2011_pred <- read.csv("./data/Am_Paul_Cressie_2011_block_results.csv",
                         # head = TRUE)

## Combine into one df
combined_df <- FRK_df %>% 
  dplyr::select(c("variable", "value", "area_sqrt", "Framework", "Scheme")) %>% 
  # rbind(PaC2011, georob_df) %>% 
  rbind(georob_df) %>% 
  mutate(fwk_sch = paste(Framework, Scheme, sep = ": ")) %>%
  ## alter the labels to change the facet_wrap titles:
  mutate(variable = factor(
    variable, 
    # labels = c("'Block prediction /' ~ 10^3 ~ 'counts' ~ min^-1", 
    #            "'RMSPE from block prediction /' ~ 10^3 ~ 'counts' ~ min^-1")
    labels = c("'Block prediction'", "'RMSPE from block prediction'")
    ))

combined_df$Framework <- as.character(combined_df$Framework)
combined_df$Scheme <- as.character(combined_df$Scheme)


combined_df$Framework[combined_df$Framework == "georob: optimal block prediction"] <- "georob"


ggsave( 
  ggplot(data = combined_df %>% 
           subset(!grepl("georob: permanence", Framework)) #%>%
           # subset(!grepl("Paul", Framework))
         ,
         aes(x = area_sqrt, lty = Framework, colour = Scheme, group = fwk_sch)) +
    geom_line(aes(y = value)) +
    facet_wrap(~variable, scales = "free", labeller = label_parsed) + 
    labs(x = "Block size (m)", y = "", 
         colour = "Blocking Scheme") +
    theme_bw() + scale_y_continuous(labels = scales::scientific) + 
    theme(text = element_text(size = 20), 
          strip.text = element_text(size = 20)),
  filename = "Am_comparison.png", device = "png", width = 13.6, height = 4.5,
  path = "./img/"
)

# labs(y = as.expression(bquote("Block prediction /" ~ 10^3 ~ "counts" ~ min^-1)))
# labs(y = as.expression(bquote("RMSPE from block prediction /" ~ 10^3 ~ "counts" ~ min^-1)))

## Size of largest block used by Paul and Cressie (2011)
# sqrt(((2.85575 - 2.8532) * (2.19835 - 2.196) ) * 10^5 * 10^5)
