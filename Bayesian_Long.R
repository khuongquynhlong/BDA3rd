library(nimble)
library(coda)
library(tidybayes)
library(tidyverse)

#---------- Read data

df <- readRDS("Data/GSPS.RData")


glimpse(df)

sapply(df, function(x){sum(is.na(x))})


#----------  Set up theme
#===============================================================================
mytheme <- function(...) {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14,color = "grey10",  face = "bold", hjust = 0.5),
      plot.subtitle = element_text(face = "italic", color = "gray10", size = 14),
      plot.caption = element_text(face = "italic", size = 14, color = "gray10"),
      axis.line = element_line(linetype = "solid"),
      axis.text.x = element_text(color = "gray10", size = 14),
      axis.text.y = element_text(color = "gray10", size = 14),
      # axis.ticks = element_blank(),
      axis.title.x = element_text(color = "gray10", size = 14),
      axis.title.y = element_text(color = "gray10", size = 14),
      panel.grid.minor = element_blank(),
      # panel.grid.major = element_blank(),
      plot.background = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 14, face = "bold"),
      legend.direction = "horizontal",
      legend.position = "top",
      legend.background = element_rect(fill = NA, color = NA),
      legend.text = element_text(size = 14),
      legend.key.width = unit(4, "line"),
      strip.text = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = NA, color = NA)
    )
}


#---------- Data exploration
#===============================================================================
df %>% filter(handdum == 0) %>% select(handper)
df %>% filter(handdum == 1) %>% select(handper)

df %>% select(age, hsat, handper, hhninc, educ, docvis, hospvis) %>%
  gather(key = "var", value = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(fill = "navy", color = "white") +
  facet_wrap(~var, ncol = 3, scales = "free") +
  mytheme()


par(mfrow = c(1, 2))
hist(df$hhninc, prob = T)
lines(density(df$hhninc))
hist(log(df$hhninc), probability = T)
lines(density(log(df$hhninc)))
par(mfrow = c(1, 1))

table(df$year)
hist(df$educ)

table(df$public)

#---------- Standardized variables
#===============================================================================
std <- function(x){(x - mean(x))/sd(x)}

#---------- Question 1
#===============================================================================
# Set up the data given
model_data <- list(working     = df$working,
                   female      = df$female,
                   year        = std(df$year),
                   age         = std(df$age),
                   hsat        = std(df$hsat),
                   handdum     = df$handdum,
                   handper     = std(df$handper),
                   hhninc      = std(df$hhninc),
                   hhkids      = df$hhkids,
                   educ        = std(df$educ),
                   married     = df$married,
                   docvis      = std(df$docvis),
                   hospvis     = std(df$hospvis),
                   public      = df$public)
model_constant <- list(N = nrow(df))

# Initial value
initial_values <- list(beta = rep(0, 12))

# specify MCMC details
params <- c("beta")

# Vague prior
#-------------------------------------------------------------------------------
# Write nimble model, logistic regression
logit_mod1_code <- nimbleCode(
  {
    # Data model
    for (i in 1:N) {
      logit(p[i]) <- beta[1] + beta[2]*female[i] + beta[3]*year[i] + beta[4]*age[i] + 
        beta[5]*hsat[i] + beta[6]*handper[i] + beta[7]*hhninc[i] + beta[8]*hhkids[i] +
        beta[9]*educ[i] + beta[10]*married[i] + beta[11]*docvis[i] + beta[12]*hospvis[i]
      working[i] ~ dbern(p[i])
    }
    # Prior
    for (i in 1:12){
      beta[i]~ dnorm(0.0,1.0E-4)
    }
})


logit_mod1 <- nimbleMCMC(code = logit_mod1_code,
                        data = model_data,
                        constants = model_constant,
                        inits = initial_values,
                        monitors = params,
                        niter = 5000,
                        nburnin = 1000,
                        nchains = 2,
                        thin = 1,
                        samplesAsCodaMCMC = TRUE,
                        WAIC = TRUE)

# Convert  into mcmc.list 

model1_mcmc <- as.mcmc.list(logit_mod1$samples)

gelman.diag(model1_mcmc)
gelman.plot(model1_mcmc)

# Produce general summary of obtained MCMC sampling

plot(model1_mcmc)
summary(model1_mcmc)
densplot(model1_mcmc)

# Extract MCMC
model1_mcmc %>% gather_draws(`beta[1]`) %>% median_hdi() %>% 
  select(".value",".lower",".upper")


# WAIC
logit_mod1$WAIC





#---------- Question 2
#===============================================================================

# Initial value for GLMM
initial_re <- list(beta = rep(0, 12), mu0 = 0, tau0 = 1, b0 = rnorm(nrow(df)))

# specify MCMC details
params_re <- c("beta", "mu0", "tau0")

# Vague prior
#-------------------------------------------------------------------------------
# Write nimble model, logistic regression + random intercept
logit_RE_code <- nimbleCode(
  {
    # Data model
    for (i in 1:N) {
      logit(p[i]) <- beta[1] + beta[2]*female[i] + beta[3]*year[i] + beta[4]*age[i] + 
        beta[5]*hsat[i] + beta[6]*handper[i] + beta[7]*hhninc[i] + beta[8]*hhkids[i] +
        beta[9]*educ[i] + beta[10]*married[i] + beta[11]*docvis[i] + beta[12]*hospvis[i] +
        b0[i]
      working[i] ~ dbern(p[i])
      b0[i] ~ dnorm(mu0, var = tau0)
    }
    # Prior
    for (i in 1:12){
      beta[i]~ dnorm(0.0, 1.0E-4)
    }
    mu0 ~ dnorm(0, 1.0E-4)
    tau0 ~ dinvgamma(2, 1)
  })


logit_RE <- nimbleMCMC(code = logit_RE_code,
                         data = model_data,
                         constants = model_constant,
                         inits = initial_re,
                         monitors = params_re,
                         niter = 5000,
                         nburnin = 1000,
                         nchains = 2,
                         thin = 1,
                         samplesAsCodaMCMC = TRUE,
                         WAIC = TRUE)

# Convert  into mcmc.list 

logit_RE_mcmc <- as.mcmc.list(logit_RE$samples)

gelman.diag(logit_RE_mcmc)
gelman.plot(logit_RE_mcmc)

# Produce general summary of obtained MCMC sampling

plot(logit_RE_mcmc)
summary(logit_RE_mcmc)
densplot(logit_RE_mcmc)















