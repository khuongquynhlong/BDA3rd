library(nimble)
library(coda)
library(tidybayes)
library(tidyverse)

#---------- Read data

df <- readRDS("GSPS.RData")


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


#---------- Question 1
#===============================================================================
# Set up the data given

model_data <- list(working = df$working,
                   female = df$female,
                   year = df$year,
                   age = df$age,
                   hsat = df$hsat,
                   handdum = df$handdum,
                   handper = df$handper,
                   hhninc = df$hhninc,
                   hhkids = df$hhkids,
                   educ = df$educ,
                   married = df$married,
                   docvis = df$docvis,
                   hospvis = df$hospvis,
                   public = df$public)
model_constant <- list(N = nrow(df))


# Initial value
initial_values <- list(beta = c(0, 0))

# specify MCMC details
params <- c("beta")


# Vague prior
#-------------------------------------------------------------------------------
# Write nimble model, logistic regression
logit_mod1 <- nimbleCode(
  {
    # Data model
    for (i in 1:N) {
      logit(p[i]) <- beta[1] + beta[2]*female[i]
      working[i] ~ dbern(p[i])
    }
    # Prior
    for (i in 1:2){
      beta[i]~ dnorm(0.0,1.0E-4)
    }
})


logit_mod1 <- nimbleMCMC(code = logit_mod1,
                        data = model_data,
                        constants = model_constant,
                        inits = initial_values,
                        monitors = params,
                        niter = 5000,
                        nburnin = 1000,
                        nchains = 3,
                        thin = 1,
                        samplesAsCodaMCMC = TRUE,
                        WAIC = TRUE)



# Convert  into mcmc.list 

model1_mcmc <- as.mcmc.list(logit_mod1$samples)

# Produce general summary of obtained MCMC sampling

plot(model1_mcmc)
summary(model1_mcmc)

# WAIC
logit_mod1$WAIC




















