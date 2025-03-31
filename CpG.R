library(ggplot2)
library(moments)
library(methods)
library(rstan)
library(bestNormalize)
library(dplyr)
library(tidyr)
library(rstudioapi)
library(ratematrix)
library(bayesplot)
library(loo)


data <- read.csv("Smoker_Epigenetic_df.csv")
data <- na.omit(data)

AgeDistribution <- function() {
  ggplot(data, aes(x = Age)) +
    geom_histogram(bins = 30, fill = 4, color = "black") +
    labs(title = "Age Distribution", x = "Age", y = "Count") +
    theme_minimal()
}

AgeDensityBySmokingStatus <- function() {
  ggplot(data, aes(x = Age, fill = Smoking.Status)) +
  geom_density(alpha = 0.5) +
  labs(title = "Age Density by Smoking Status", x = "Age", y = "Density") +
  theme_minimal()
}

MethylationBetaValuecg00050873 <- function() {
  ggplot(data, aes(x = cg00050873, fill = Smoking.Status)) +
  geom_density(alpha = 0.5) +
  labs(title = "Methylation Distribution (cg00050873) by Smoking Status",
       x = "Methylation Beta Value", y = "Density") +
  theme_minimal()
}
MethylationBetaValuecg00212031SS <- function() {
  ggplot(data, aes(x = cg00212031, fill = Smoking.Status)) +
    geom_density(alpha = 0.5) +
    labs(title = "Methylation Distribution (cg00212031) by Smoking Status",
         x = "Methylation Beta Value", y = "Density") +
    theme_minimal()
}

AgeVScg00050873Methylationcg00212031G <- function() {
  ggplot(data, aes(x = Age, y = cg00050873, color = Gender)) +
  geom_point(alpha = 0.7, size = 3) +
  labs(title = "Age vs. cg00050873 Methylation",
       x = "Age",
       y = "Methylation Value (cg00050873)",
       color = "Gender") +
  theme_minimal()
}

BoxPlotOfcg00050873MethylationbySmokingStatus <- function() {
  ggplot(data, aes(x = Smoking.Status, y = cg00050873, fill = Smoking.Status)) +
  geom_boxplot() +
  labs(title = "Box Plot of cg00050873 Methylation by Smoking Status",
       x = "Smoking Status",
       y = "Methylation Value (cg00050873)") +
  theme_minimal()
}


data_long <- pivot_longer(
  data,
  cols = starts_with("cg"),
  names_to = "CpG_site",
  values_to = "MethylationValue"
)

AgeVSMethylationValuesforDifferentCpGSites <- function() {
  ggplot(data_long, aes(x = Age, y = MethylationValue, color = Gender)) +
  geom_point(alpha = 0.4) +
  facet_wrap(~ CpG_site, scales = "free_y") +
  labs(title = "Age vs. Methylation Values for Different CpG Sites",
       x = "Age", y = "Methylation Value") +
  theme_minimal()
}

data_females <- data_long %>% filter(Gender == " f")
data_males   <- data_long %>% filter(Gender == " m")



FemaleMethylationBySmokingStatusforDifferentCpGSites <- function() {
  ggplot(data_females, aes(x = MethylationValue, fill = Smoking.Status)) +
    geom_density(alpha = 0.7) +
    facet_wrap(~ CpG_site, scales = "free") +
    labs(title = "Female Methylation Distribution by Smoking Status for Different CpG Sites",
         x = "Methylation Beta Value", y = "Density") +
    theme_minimal()
}
MaleMethylationBySmokingStatusforDifferentCpGSites <- function() {
  ggplot(data_males, aes(x = MethylationValue, fill = Smoking.Status)) +
    geom_density(alpha = 0.7) +
    facet_wrap(~ CpG_site, scales = "free") +
    labs(title = "Male Methylation Distribution by Smoking Status for Different CpG Sites",
         x = "Methylation Beta Value", y = "Density") +
    theme_minimal()
}


FemaleAgeVsMethylationBySmokingStatusForCpGSites <- function() {
  ggplot(data_females, aes(x = Age, y = MethylationValue, color = Smoking.Status)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "loess", se = FALSE) +
    facet_wrap(~ CpG_site, scales = "free_y") +
    labs(title = "Female: Age vs Methylation by Smoking Status for Different CpG Sites",
         x = "Age",
         y = "Methylation Beta Value") +
    theme_minimal()
}

MaleAgeVsMethylationBySmokingStatusForCpGSites <- function() {
  ggplot(data_males, aes(x = Age, y = MethylationValue, color = Smoking.Status)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "loess", se = FALSE) +
    facet_wrap(~ CpG_site, scales = "free_y") +
    labs(title = "Male: Age vs Methylation by Smoking Status for Different CpG Sites",
         x = "Age",
         y = "Methylation Beta Value") +
    theme_minimal()
}




data <- data %>% mutate(smoker = ifelse(Smoking.Status == "current", 1, 0))
data <- data %>% mutate(male = ifelse(Gender == " m", 1, 0))
data <- data %>% mutate(Age_std = scale(Age))
X <- as.matrix(data %>% select(Age_std, male, smoker))
meth_cols <- setdiff(grep("^cg", names(data), value = TRUE), c("cg00050873", "cg00212031", "cg00455876", "cg02011394", "cg02050847", "cg02494853", "cg02842889", "cg02839557", "cg03244189", "cg03706273"))
Y <- as.matrix(data[, meth_cols])

N <- nrow(X)
K <- ncol(X)
J <- ncol(Y)

stan_data <- list(
  N = N,
  K = K,
  J = J,
  X = X,
  Y = Y
)




fit <- stan(file = "CpGStanCheck.stan", data = stan_data, iter = 2000,
             chains = 4, warmup = 500, seed = 123, control = list(max_treedepth = 15))


selected_params <- c("beta0[1]", "beta0[2]", "beta0[3]",
                     "beta[1,1]", "beta[2,1]", "beta[3,1]",
                     "sigma[1]", "sigma[2]", "sigma[3]")

traceplot(fit, pars = selected_params)



checkConvergence <- function(model) {
  library(rstan)
  rhat_values <- summary(model)$summary[, "Rhat"]
  if (any(rhat_values > 1.01, na.rm = TRUE)) {
    warning("Some parameters have Rhat > 1.01, indicating potential convergence issues.")
  } else {
    message("All Rhat values are <= 1.01. Convergence looks good!")
  }
  sampler_params <- get_sampler_params(model, inc_warmup = FALSE)
  num_divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))

  if (num_divergences > 0) {
    warning(paste("There are", num_divergences, "divergent transitions. Consider reparameterizing or increasing adapt_delta."))
  } else {
    message("No divergent transitions detected.")
  }

  ess_values <- summary(model)$summary[, "n_eff"]
  if (any(ess_values < 100, na.rm = TRUE)) {
    warning("Some parameters have low effective sample sizes (<100). Consider increasing iterations.")
  } else {
    message("Effective sample sizes are reasonable.")
  }
}

checkConvergence(fit)


rhat_values <- summary(fit)$summary[, "Rhat"]

mcmc_rhat(rhat_values)

log_lik <- extract_log_lik(fit, merge_chains = FALSE)

loo_result <- loo(log_lik)

print(loo_result)
plot(loo_result)


y <- data$cg00213748
y_rep <- as.matrix(rstan::extract(fit, pars = "y_rep")$y_rep)
num_samples <- length(rstan::extract(fit, pars = "lp__")$lp__)
print(num_samples)
num_samples <- 6000
num_obs <- 10*length(y)
y_rep <- matrix(y_rep, nrow = num_samples, ncol = num_obs, byrow = TRUE)




y = data$cg02004872
y_rep = as.matrix(rstan::extract(fit, pars = "y_rep")$y_rep)

yrep_reshaped <- array(y_rep, dim = c(6000, 621, 10))
y_rep_matrix <- yrep_reshaped[1:100, , 4]
ppc_dens_overlay(y = y, yrep = y_rep_matrix)






