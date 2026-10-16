library(ggplot2)

res <- read.csv("replicates.csv")

res$parameter <- factor(res$parameter, levels = c("m","a","t_peak","eta","sigma"),
   labels = c("m (cases/day)", "a", "t_peak (day)", "eta", "sigma (days)"))

truth <- data.frame(parameter = levels(res$parameter), true = c(0.1, 0.8, 200, 0.6, 20))

p <- ggplot(res, aes(x = "", y = estimate)) +
  geom_hline(data = truth, aes(yintercept = true), linetype = "dashed") +
  geom_jitter(width = 0.15, alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "maximum-likelihood estimate") + theme_minimal() +
  theme(axis.text.x = element_blank())

ggsave("replicates_r.png", p, width = 10, height = 3, dpi = 150)
