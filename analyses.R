library(gee) # For fitting GEE models.
library(jmv) # For t-tests.
library(ltm) # For Cronbach's alpha.
library(emmeans) # For hypothesis testing using estimated marginal means.
library(effectsize) # For eta squared.
library(rcompanion) # For epsilon squared.
library(boot) # For bootstrapping pseudo-R^2 confidence intervals.
library(dplyr) # For selecting columns.
library(tidyr) # For pivoting tables.



# [LOAD THE DATA]

# Load study data from experiment 1 in long format.
sl1 <- read.csv("exp1-long.tsv", sep = "\t")
sl1$group <- relevel(as.factor(sl1$group), ref="OG") # Set outgroup as the reference level.
sl1$subject <- as.factor(sl1$subject) # Convert subject ID to a factor.
sl1$chain <- as.factor(sl1$chain) # Convert chain ID to a factor.

# Load study data from experiment 1 in wide format.
sw1 <- read.csv("exp1-wide.tsv", sep = "\t")
sw1$subject <- as.factor(sw1$subject) # Convert subject ID to a factor.
sw1$chain <- as.factor(sw1$chain) # Convert chain ID to a factor.

# Load study data from experiment 2 in long format.
sl2 <- read.csv("exp2-long.tsv", sep = "\t")
sl2$group <- relevel(as.factor(sl2$group), ref="OG") # Set outgroup as the reference level.
sl2$subject <- as.factor(sl2$subject) # Convert subject ID to a factor.
sl2$models <- as.factor(sl2$models) # Convert the number of observed individuals (behavioral models) to factor.

# Load study data from experiment 2 in wide format.
sw2 <- read.csv("exp2-wide.tsv", sep = "\t")
sw2$subject <- as.factor(sw2$subject) # Convert subject ID to a factor.
sw2$models <- as.factor(sw2$models) # Convert the number of observed individuals (behavioral models) to factor.

# Load study data from experiment 2 for computing Cronbach's alpha.
cron <- read.csv("exp2-cronbach.tsv", sep = "\t")

# Fixed correlational structure, where responses are perfectly negatively correlated. Tiny offsets are required for the gee package to run.
fixed <- matrix(c(1, -1 + 10^-10, -1 + 10^-10, 1), 2, 2)



# [HELPER FUNCTIONS]

# Convert Wald F-statistics to Wald X^2 statistics.
F_to_chisq.gee <- function(tests, digits=3, eps=0.001) {
  output <- data.frame(matrix(ncol=4, nrow=0, dimnames=list(NULL, c("model term", "df", "chisq", "p.value"))))
  for(i in 1:nrow(tests)) {
    chisq <- format(round(as.numeric(tests[["F.ratio"]][i]) * as.numeric(tests[["df1"]][i]), digits=digits), nsmall=digits)
    output[nrow(output) + 1,] <- list(tests[["model term"]][i], tests[["df1"]][i], chisq, format.pval(tests[["p.value"]][i], digits=digits, eps=eps))
  }
  return(output)
}

# Calculate the McKelvey-Zavoina pseudo-R^2.
r2.mz <- function(model, data) {
  y <- predict(model, type="link")
  y.hat <- mean(y)
  term <- sum((y - y.hat)^2)
  return(term / (term + nrow(data) * (pi^2)/3))
}

# Calculate the partial McKelvey-Zavoina pseudo-R^2.
r2.mz.partial <- function(model, reduced, data) {
  m <- r2.mz(model, data)
  r <- r2.mz(reduced, data)
  return((m - r) / (1 - r))
}

# Bootstrapping function to get CIs for the McKelvey-Zavoina pseudo-R^2 and partial R^2.
boot.r2.mz <- function(data, indices, formula, dv) {
  gc()
  d <- data[indices,]
  d$subject <- as.factor(1:nrow(d))
  d$group <- NULL
  d <- pivot_longer(d, cols=c(paste(dv, "IG", sep=""), paste(dv, "OG", sep="")), names_to="group", values_to=dv)
  if(is.list(formula)) {
    model <- gee(formula[[1]], data=d, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject, maxiter=2500)
    reduced <- gee(formula[[2]], data=d, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject, maxiter=2500)
    return(r2.mz.partial(model, reduced, d))
  }
  else {
    model <- gee(formula, data=d, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject, maxiter=2500)
    return(r2.mz(model, d))
  }
}

# Convert R to R^2.
r_to_r2 <- function(corr, digits=2) {
  r <- corr$estimate[[1]]
  lower <- ifelse(sign(r) == sign(corr$conf.int[[1]]), corr$conf.int[[1]], 0)
  upper <- ifelse(sign(r) == sign(corr$conf.int[[2]]), corr$conf.int[[2]], 0)
  r2 <- format(round(r^2, digits=digits), nsmall=digits)
  LCL <- format(round(min(lower^2, upper^2), digits=digits), nsmall=digits)
  UCL <- format(round(max(lower^2, upper^2), digits=digits), nsmall=digits)
  level <- paste(attr(corr$conf.int,"conf.level") * 100, "%", sep="")
  return(data.frame(matrix(c(r2, LCL, UCL, level), ncol=4, nrow=1, dimnames=list(NULL, c("r2", "LCL", "UCL", "conf.level")))))
}

# Post hoc tests to determine which level(s) of an x*2 table differ.
chisq.posthoc <- function(x, y, adjust="bonferroni") {
  output <- data.frame(matrix(ncol=7, nrow=0, dimnames=list(NULL, c("model term", "df", "chisq", "p.value", "V", "LCL", "UCL"))))
  df <- data.frame("x" = as.factor(x), "y" = as.factor(y))
  for(i in 1:nlevels(df$x)) {
    df$bin <- ifelse(df$x == levels(df$x)[[i]], 1, 0)
    t <- chisq.test(df$bin, df$y)
    es <- effectsize(t)
    output[nrow(output) + 1,] <- list(levels(df$x)[[i]], t$parameter[[1]], t$statistic[[1]], t$p.value, es$Cramers_v, es$CI_low, es$CI_high)
  }
  output$p.value <- p.adjust(output$p.value, method=adjust)
  return(output)
}



# [EXPERIMENT 1]

# [INGROUP COPYING BIAS]

m0 <- gee(cbind(agree, agreeTotal - agree) ~ group, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
emmeans(m0, ~ group, type="response")
pairs(.Last.value, reverse=TRUE)
confint(.Last.value)

# Proportions.
prop.table(table(sw1$agreeBias))


# [INGROUP ATTENTIONAL BIAS]

m0 <- gee(cbind(obs, obsTotal - obs) ~ group, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
emmeans(m0, ~ group, type="response")
pairs(.Last.value, reverse=TRUE)
confint(.Last.value)

# Proportions.
prop.table(table(sw1$obsBias))


# [INGROUP RELIABILITY RATINGS]

# Proportions.
pre <- table(sw1$pre)
post <- table(sw1$post)
prop.table(pre)
prop.table(post)

# Most explicitly reject the notion that the ingroup is more reliable.
binom.test(pre[["equal"]] + pre[["less"]], sum(pre))
binom.test(post[["equal"]] + post[["less"]], sum(post))

# Of those who believe one group is more reliable, the ingroup has a transient advantage.
binom.test(pre[["more"]], pre[["more"]] + pre[["less"]])
binom.test(post[["more"]], post[["more"]] + post[["less"]])


# [CULTURAL DIVERGENCE]

# Amount of cultural divergence.
mean(sw1[["hammingIG"]])
mean(sw1[["hammingOG"]])
ttestPS(data=sw1, pairs = list(list(i1="hammingIG", i2="hammingOG")), meanDiff=TRUE, ci=TRUE, effectSize=TRUE, ciES=TRUE)

# Rate of cultural divergence.
m0 <- lm(hammingDelta ~ generation, data=sw1)
summary(m0)
r_to_r2(cor.test(sw1$hammingDelta, sw1$generation, conf.level=.9))



# [EXPERIMENT 2]

# [CRONBACH'S ALPHA]

cronbach.alpha(select(cron, c("warm", "friendly")), CI=TRUE, B=5000)
cronbach.alpha(select(cron, c("competent", "capable")), CI=TRUE, B=5000)


# [INGROUP COPYING BIAS BY CONDITION]

# Omnibus test.
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*models, data=subset(sl2, agreeTotal > 0), family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
pairs(pairs(emmeans(m0, ~ group | models, type="response"), reverse=TRUE), by=NULL, adjust="bonferroni")
r2 <- boot(data=subset(sw2, agreeTotal > 0), statistic=boot.r2.mz, R=1000, formula=list(cbind(agree, agreeTotal - agree) ~ group*models, cbind(agree, agreeTotal - agree) ~ group), dv="agree", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Ingroup copying bias by condition.
emmeans(m0, ~ group | models, type="response")
pairs(.Last.value, reverse=TRUE)
confint(.Last.value)


# [INGROUP PERCEPUTUAL BIASES]

# Pre-game bias in perceived warmth and competence.
ttestPS(data=sw2, pairs = list(list(i1="preWarmIG", i2="preWarmOG")), meanDiff=TRUE, ci=TRUE, effectSize=TRUE, ciES=TRUE)
ttestPS(data=sw2, pairs = list(list(i1="preCompIG", i2="preCompOG")), meanDiff=TRUE, ci=TRUE, effectSize=TRUE, ciES=TRUE)

# Post-game bias in perceived warmth and competence.
ttestPS(data=sw2, pairs = list(list(i1="postWarmIG", i2="postWarmOG")), meanDiff=TRUE, ci=TRUE, effectSize=TRUE, ciES=TRUE)
ttestPS(data=sw2, pairs = list(list(i1="postCompIG", i2="postCompOG")), meanDiff=TRUE, ci=TRUE, effectSize=TRUE, ciES=TRUE)


# [CULTURAL DIVERGENCE]

# Amount of cultural divergence by condition.
m0 <- aov(hammingDelta ~ models, data=sw2)
summary(m0)
effectsize(m0)

# Rate of cultural divergence by condition.
m0 <- aov(hammingDelta ~ models*generation, data=sw2)
summary(m0)
effectsize(m0)

# Amount of cultural divergence.
mean(sw2[["hammingIG"]])
mean(sw2[["hammingOG"]])
ttestPS(data=sw2, pairs = list(list(i1="hammingIG", i2="hammingOG")), meanDiff=TRUE, ci=TRUE, effectSize=TRUE, ciES=TRUE)

# Rate of cultural divergence.
m0 <- lm(hammingDelta ~ generation, data=sw2)
summary(m0)
r_to_r2(cor.test(sw2$hammingDelta, sw2$generation, conf.level=.9))



# [SUPPLEMENTAL ANALYSES OF COPYING]

# [DIRECT VS. INDIRECT BIAS]

# Ingroup copying bias after controlling for attention.
m0 <- gee(cbind(agree, agreeMax - agree) ~ group, data=subset(sl1, agreeMax > 0), family=binomial, corstr="exchangeable", scale.fix=TRUE, id=subject)
m0 <- gee(cbind(agree, agreeMax - agree) ~ group, data=subset(sl2, agreeMax > 0), family=binomial, corstr="exchangeable", scale.fix=TRUE, id=subject)
pairs(emmeans(m0, ~ group, type="response"), reverse=TRUE)
confint(.Last.value)

# Ingroup copying bias after controlling for attentional and perceptual biases, experiment 1.
sl1$preBias <- ifelse(sl1$pre == "more" | sl1$pre == "less", "Y", "N")
sl1$postBias <-ifelse(sl1$post == "more" | sl1$post == "less", "Y", "N")
m0 <- gee(cbind(agree, agreeMax - agree) ~ group*preBias + group*postBias, data=subset(sl1, agreeMax > 0), family=binomial, corstr="exchangeable", scale.fix=TRUE, id=subject)
pairs(emmeans(m0, ~ group | preBias + postBias, type="response"), reverse=TRUE)
confint(.Last.value)

# Ingroup copying bias after controlling for attentional and perceptual biases, experiment 2.
m0 <- gee(cbind(agree, agreeMax - agree) ~ group*preWarmDelta + group*preCompDelta + group*postWarmDelta + group*postCompDelta, data=subset(sl2, agreeMax > 0), family=binomial, corstr="exchangeable", scale.fix=TRUE, id=subject)
pairs(ref_grid(m0, at=list(preCompDelta=0, postCompDelta=0, preWarmDelta=0, postWarmDelta=0), type="response"), reverse=TRUE)
confint(.Last.value)


# [RELIABILITY AND COPYING]

# Pre-game ratings omnibus test.
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*pre, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(agree, agreeTotal - agree) ~ group*pre, cbind(agree, agreeTotal - agree) ~ group), dv="agree", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Ingroup copying bias by pre-game rating.
emmeans(m0, ~ group | pre, type="response")
pairs(.Last.value, reverse=TRUE)
confint(.Last.value)

# Post-game ratings omnibus test.
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*post, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
pairs(pairs(emmeans(m0, ~ group | post, type="response"), reverse=TRUE), by=NULL, adjust="bonferroni")
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(agree, agreeTotal - agree) ~ group*post, cbind(agree, agreeTotal - agree) ~ group), dv="agree", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Ingroup copying bias by post-game rating.
emmeans(m0, ~ group | post, type="response")
pairs(.Last.value, reverse=TRUE)
confint(.Last.value)


# [WARMTH, COMPETENCE, AND COPYING]

# Is the ingroup copying bias affected by pre-game ratings?
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*preWarmDelta + group*preCompDelta, data=subset(sl2, agreeTotal > 0), family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
pairs(emtrends(m0, ~ group, var="preWarmDelta"), reverse=TRUE)
pairs(emtrends(m0, ~ group, var="preCompDelta"), reverse=TRUE)
confint(.Last.value)
c(exp(.Last.value$asymp.LCL), exp(.Last.value$estimate), exp(.Last.value$asymp.UCL))

# Pre-game ratings as categorical measures.
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*preWarmBias + group*preCompBias, data=subset(sl2, agreeTotal > 0), family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject, maxiter=1000)
pairs(emmeans(m0, ~ group | preWarmBias, type="response"), reverse=TRUE)
pairs(emmeans(m0, ~ group | preCompBias, type="response"), reverse=TRUE)
confint(.Last.value)

# Is the ingroup copying bias affected by post-game ratings?
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*postWarmDelta + group*postCompDelta, data=subset(sl2, agreeTotal > 0), family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
pairs(emtrends(m0, ~ group, var="postWarmDelta"), reverse=TRUE)
pairs(emtrends(m0, ~ group, var="postCompDelta"), reverse=TRUE)
confint(.Last.value)
c(exp(.Last.value$asymp.LCL), exp(.Last.value$estimate), exp(.Last.value$asymp.UCL))

# Post-game ratings as categorical measures.
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*postWarmBias + group*postCompBias, data=subset(sl2, agreeTotal > 0), family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
pairs(emmeans(m0, ~ group | postWarmBias, type="response"), reverse=TRUE)
pairs(emmeans(m0, ~ group | postCompBias, type="response"), reverse=TRUE)
confint(.Last.value)



# [SUPPLEMENTAL EXAMINATION OF DATA QUALITY, PREREGISTRATION, AND MISCELLANEOUS RESULTS]

# [DATA QUALITY]

# Country of origin.
table(sw1$IP_Hub_Country_Code)
prop.table(.Last.value)

# VPS.
table(sw1$IP_Hub_VPS)
prop.table(.Last.value)

# Flagged based on country of origin or VPS use.
table(sw1$IP_Hub_recommend_block)
prop.table(.Last.value)

# Ingroup copying bias by whether participant is flagged.
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*IP_Hub_recommend_block, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(agree, agreeTotal - agree) ~ group*IP_Hub_recommend_block, cbind(agree, agreeTotal - agree) ~ group), dv="agree", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Ingroup attentional bias by whether participant is flagged.
m0 <- gee(cbind(obs, obsTotal - obs) ~ group*IP_Hub_recommend_block, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(obs, obsTotal - obs) ~ group*IP_Hub_recommend_block, cbind(obs, obsTotal - obs) ~ group), dv="obs", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Pre-game ratings by whether participant is flagged.
chisq.test(sw1$pre, sw1$IP_Hub_recommend_block)
effectsize(.Last.value)

# Post-game ratings by whether participant is flagged.
chisq.test(sw1$post, sw1$IP_Hub_recommend_block)
effectsize(.Last.value)
chisq.posthoc(sw1$post, sw1$IP_Hub_recommend_block, adjust="bonferroni")

# Amount of cultural divergence by whether participant is flagged.
m0 <- aov(hammingDelta ~ IP_Hub_recommend_block, data=sw1)
summary(m0)
effectsize(m0)

# Rate of cultural divergence by whether participant is flagged.
m0 <- aov(hammingDelta ~ generation*IP_Hub_recommend_block, data=sw1)
summary(m0)
effectsize(m0)

# Effect of post-game ratings on copying by whether participant is flagged. If bootstrapping creates rank-deficient models, try wrapping the bootstrapping function in a try-catch statement.
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*post*IP_Hub_recommend_block, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(agree, agreeTotal - agree) ~ group*post*IP_Hub_recommend_block, cbind(agree, agreeTotal - agree) ~ group*post), dv="agree", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Effect of post-game ratings on attention by whether participant is flagged. If bootstrapping creates rank-deficient models, try wrapping the bootstrapping function in a try-catch statement.
m0 <- gee(cbind(obs, obsTotal - obs) ~ group*post*IP_Hub_recommend_block, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(obs, obsTotal - obs) ~ group*post*IP_Hub_recommend_block, cbind(obs, obsTotal - obs) ~ group*post), dv="obs", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Do flagged individuals show the same qualitative patterns in their post-game ratings?
pre <- table(subset(sw1, IP_Hub_recommend_block==1)$pre)
post <- table(subset(sw1, IP_Hub_recommend_block==1)$post)

# Most explicitly reject the notion that the ingroup is more reliable.
binom.test(post[["equal"]] + post[["less"]], sum(post))

# Of those who believe one group is more reliable, post-game ratings are evenly split.
binom.test(post[["more"]], post[["more"]] + post[["less"]])


# [COPYING AND ATTENTION]

# Correlation between copying and attention.
cor.test(sw1$agreeIG, sw1$obsIG)
cor.test(sw2$agreeIG, sw2$obsIG)

# Participants who prefer to observe the outgroup show an outgroup copying bias
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group, data=subset(sl1, obsBias == "OG"), family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group, data=subset(sl2, obsBias == "OG"), family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
pairs(emmeans(m0, ~ group, type="response"), reverse=TRUE)
confint(.Last.value)


# [RELIABILITY AND ATTENTION]

# Pre-game ratings omnibus test.
m0 <- gee(cbind(obs, obsTotal - obs) ~ group*pre, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(obs, obsTotal - obs) ~ group*pre, cbind(obs, obsTotal - obs) ~ group), dv="obs", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Ingroup attentional bias by pre-game rating.
emmeans(m0, ~ group | pre, type="response")
pairs(.Last.value, reverse=TRUE)
confint(.Last.value)

# Post-game ratings omnibus test.
m0 <- gee(cbind(obs, obsTotal - obs) ~ group*post, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
pairs(pairs(emmeans(m0, ~ group | post, type="response"), reverse=TRUE), by=NULL, adjust="bonferroni")
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(obs, obsTotal - obs) ~ group*post, cbind(obs, obsTotal - obs) ~ group), dv="obs", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Ingroup attentional bias by post-game rating.
emmeans(m0, ~ group | post, type="response")
pairs(.Last.value, reverse=TRUE)
confint(.Last.value)


# [ATTENTION BY CONDITION]

# Omnibus test.
m0 <- gee(cbind(obs, obsTotal - obs) ~ group*models, data=sl2, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
pairs(pairs(emmeans(m0, ~ group | models, type="response"), reverse=TRUE), by=NULL, adjust="bonferroni")
r2 <- boot(data=sw2, statistic=boot.r2.mz, R=1000, formula=list(cbind(obs, obsTotal - obs) ~ group*models, cbind(obs, obsTotal - obs) ~ group), dv="obs", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Ingroup attentional bias by condition.
emmeans(m0, ~ group | models, type="response")
pairs(.Last.value, reverse=TRUE)
confint(.Last.value)


# [CHANGES IN PERCEIVED WARMTH AND COMPETENCE]

# Decrease in both biases over the course of the game.
ttestOneS(data=sw2, vars=vars(warmDelta, compDelta), effectSize=TRUE, ciES=TRUE)

# The change in warmth ratings does not vary with condition.
anovaNP(formula=warmDelta ~ models, data=sw2)
epsilonSquared(sw2$warmDelta, sw2$models, ci=TRUE, conf=.9, R=2000, type="bca")

# The change in competence ratings does not vary with condition.
anovaNP(formula=compDelta ~ models, data=sw2)
epsilonSquared(sw2$compDelta, sw2$models, ci=TRUE, conf=.9, R=2000, type="bca")


# [RATINGS AND FINAL SCORES]

# Pre-game reliability ratings.
m0 <- aov(finalScore ~ pre, data=sw1)
summary(m0)
effectsize(m0)
pairs(emmeans(m0, ~ pre, type="response"), adjust="bonferroni")

# Post-game reliability ratings.
m0 <- aov(finalScore ~ post, data=sw1)
summary(m0)
effectsize(m0)
pairs(emmeans(m0, ~ post, type="response"), adjust="bonferroni")

# Pre-game competence ratings.
m0 <- lm(finalScore ~ preCompIG + preCompOG, data=sw2)
summary(m0)
eta_squared(m0)

# Post-game competence ratings.
m0 <- lm(finalScore ~ postCompIG + postCompOG, data=sw2)
summary(m0)
eta_squared(m0)


# [DYNAMICS]

# Ingroup copying bias over generations, experiment 1.
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*generation, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(agree, agreeTotal - agree) ~ group*generation, cbind(agree, agreeTotal - agree) ~ group), dv="agree", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Ingroup copying bias over generations, experiment 2.
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*generation, data=subset(sl2, agreeTotal > 0), family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=subset(sw2, agreeTotal > 0), statistic=boot.r2.mz, R=1000, formula=list(cbind(agree, agreeTotal - agree) ~ group*generation, cbind(agree, agreeTotal - agree) ~ group), dv="agree", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Ingroup attentional bias over generations, experiment 1.
m0 <- gee(cbind(obs, obsTotal - obs) ~ group*generation, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(obs, obsTotal - obs) ~ group*generation, cbind(obs, obsTotal - obs) ~ group), dv="obs", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Ingroup attentional bias over generations, experiment 2.
m0 <- gee(cbind(obs, obsTotal - obs) ~ group*generation, data=sl2, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=sw2, statistic=boot.r2.mz, R=1000, formula=list(cbind(obs, obsTotal - obs) ~ group*generation, cbind(obs, obsTotal - obs) ~ group), dv="obs", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")


# [CHAINS]

# Ingroup copying bias by chain.
m0 <- gee(cbind(agree, agreeTotal - agree) ~ group*chain, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(agree, agreeTotal - agree) ~ group*chain, cbind(agree, agreeTotal - agree) ~ group), dv="agree", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Ingroup attentional bias by chain.
m0 <- gee(cbind(obs, obsTotal - obs) ~ group*chain, data=sl1, family=binomial, corstr="fixed", R=fixed, scale.fix=TRUE, id=subject, maxiter=500)
F_to_chisq.gee(joint_tests(m0))
r2 <- boot(data=sw1, statistic=boot.r2.mz, R=1000, formula=list(cbind(obs, obsTotal - obs) ~ group*chain, cbind(obs, obsTotal - obs) ~ group), dv="obs", parallel="multicore", ncpus=4)
print(r2)
boot.ci(r2, conf=.9, type="bca")

# Amount of cultural divergence by chain.
m0 <- aov(hammingDelta ~ chain, data=sw1)
summary(m0)
eta_squared(m0)

# Rate of cultural divergence by chain.
m0 <- aov(hammingDelta ~ generation*chain, data=sw1)
summary(m0)
eta_squared(m0)


# [INDIVIDUAL LEARNING BY CONDITION]

m0 <- aov(il ~ models, data=sw2)
summary(m0)
effectsize(m0)
