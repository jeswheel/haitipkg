library(tidyverse)
library(pomp)

path_epi <- system.file("extdata", "iffits_epi232_nprofs3_nmif50_nparticles100_passTmtrue_seed20190415.rds", package = 'haitipkg')
path_end <- system.file("extdata", "iffits_end232_nmif50_nparticles100_passTmfalse_seed4172019.rds", package = 'haitipkg')

IF_epi <- readRDS(path_epi)
IF_end <- readRDS(path_end)

h1LeeStartsEpi <- IF_epi %>%
  filter(loglik > -2250, beta1 < 20, nu >= 0.9) %>%
  rename(tau = theta, E_0 = E.0, I_0 = I.0, S_0 = S.0,
         A_0 = A.0, R_0 = R.0) %>%
  mutate(sig_sq = 0, pop_0 = 10911819) %>%
  select(-model, -loglik.se, -nfail.min, -nfail.max, -etime)

h1LeeStartsEnd <- IF_end %>%
  filter(loglik > -2400) %>%
  rename(tau = theta, E_0 = E.0, I_0 = I.0, S_0 = S.0,
          A_0 = A.0, R_0 = R.0, pop_0 = N0, incid_0 = incid.0) %>%
  mutate(sig_sq = 0) %>%
  select(-model, -loglik.se, -nfail.min, -nfail.max, -etime)

# GGally::ggpairs(h1LeeStartsEpi, columns = c("rho", "tau", 'beta1', 'nu', 'loglik'))
# GGally::ggpairs(h1LeeStartsEnd, columns = c("rho", "tau", 'beta1', 'nu', 'loglik'))

usethis::use_data(
  h1LeeStartsEpi,
  h1LeeStartsEnd,
  overwrite = TRUE
)


