pacman::p_load(rstanarm, broom, dplyr, ggplot2, tidyr)
data("radon", package = "rstanarm")

fit_radon_1 <- stan_glm(log_radon ~ 1,
                        data = radon, chains = 1)
fit_radon_2 <- stan_glm(log_radon ~ -1 + county,
                        data = radon, chains = 1)
fit_radon_3 <- stan_glmer(log_radon ~ (1 + floor | county),
                          data = radon, chains = 1)


alpha_1 <- tidyMCMC(fit_radon_1, conf.int = TRUE) %>%
  filter(term == "(Intercept)") %>%
  # add county
  mutate(county = list(unique(as.character(radon[["county"]])))) %>%
  unnest(county) %>%
  select(-term) %>%
  mutate(model = "Complete")

alpha_2 <- tidyMCMC(fit_radon_2, conf.int = TRUE) %>%
  filter(str_detect(term, "^county")) %>%
  mutate(county = str_replace(term, "^county", "")) %>%
  select(-term) %>%
  mutate(model = "No")

alphas <- as.matrix(fit_radon_3, regex_pars = "^b\\[")
alpha_mean <- as.matrix(fit_radon_3, pars = "(Intercept)")
alpha_3 <- sweep(alphas, 1, alpha_mean, FUN = "+") %>%
  as_tibble() %>%
  gather(term, value) %>%
  group_by(term) %>%
  summarise(estimate = mean(value), conf.low = quantile(value, 0.025),
            conf.high = quantile(value, 0.975)) %>%
  ungroup() %>%
  mutate(county = str_match(term, "county:(.*)\\]")[ , 2]) %>%
  select(-term) %>%
  mutate(model = "Partial")

all_models <-
  bind_rows(alpha_1, alpha_2, alpha_3) %>%
  # reorder county by size
  mutate(county = fct_reorder(county, estimate, mean))

ggplot(all_models, aes(x = county, y = estimate, ymin = conf.low, ymax = conf.high,
                       color = model)) +
  geom_pointrange(position = position_dodge(width = 1)) +
  coord_flip()

