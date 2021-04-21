# flashLIT search term results

# flashLIT: fast, largely automated, systematic handling of literature.
# workingconservation@gmail.com

# ==================================================================

# quick analysis of the numerical results of alternative search terms
# aimed at narrowing down what sets we would like to look at in more detail

res <- read_csv("data/SearchTermResults.csv") 

res <- res %>% 
  mutate(Terms = factor(Terms, levels = unique(res$Terms)[c(1,2,3,4)])) %>% 
  mutate(Domain = factor(Domain, levels = unique(res$Domain)[c(1,2)])) %>% 
  mutate(Version = factor(Version, levels = unique(res$Version)[c(1,2,3,4)]))


ggplot(res) + 
  geom_col(aes(y=N, x=Version, fill=Domain), position = "dodge") +
  geom_text(aes(y=N, x=Version, label = N, group=Domain), position = position_dodge(width = .8), angle = 90, hjust = -0.1, vjust = 0.5) +
  facet_grid(~Terms) +
  scale_y_continuous(expand = expansion(mult = c(0, .2))) +
  scale_fill_viridis_d() +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))


# use lm to summarise differences:
mod <- lm(N ~ Terms * Domain * Version, data = res)
summary(mod)
