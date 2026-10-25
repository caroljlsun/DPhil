# Use packages belonging to this R installation
# with the help of chatgpt

.libPaths(.Library, include.site = FALSE)

library(glmmTMB)
library(emmeans)

cat("R installation:", R.home(), "\n")
print(sessionInfo())

# Verify numerical calculations before fitting
A <- matrix(1, 106, 106)
diag(A) <- 2
C <- chol(A)
stopifnot(max(abs(crossprod(C) - A)) < 1e-10)

df <- readRDS("analysis_df.rds")

# Preserve the intended categorical predictors
for (nm in c("knockdown", "plate", "row_f", "col_f")) {
  df[[nm]] <- factor(df[[nm]])
}
df$knockdown <- relevel(df$knockdown, ref = "AAVS1")

m_clean <- glmmTMB(
  log(spot_area) ~ knockdown + plate + row_f + col_f,
  family = t_family(link = "identity"),
  data = df,
  control = glmmTMBControl(
    parallel = 1,
    optCtrl = list(iter.max = 3000, eval.max = 5000)
  )
)

V <- as.matrix(vcov(m_clean)$cond)

print(list(
  convergence = m_clean$fit$convergence,
  message = m_clean$fit$message,
  pdHess = m_clean$sdr$pdHess,
  variance_range = range(diag(V)),
  Student_t_df = family_params(m_clean)
))

# Stop if uncertainty estimates remain invalid
stopifnot(
  m_clean$fit$convergence == 0,
  isTRUE(m_clean$sdr$pdHess),
  all(is.finite(V))
)
invisible(chol((V + t(V)) / 2))

emm <- emmeans(
  m_clean, ~ knockdown,
  data = df,
  nuisance = c("plate", "row_f", "col_f"),
  wt.nuis = "equal"
)

cmp <- contrast(emm, "trt.vs.ctrl", ref = "AAVS1")

v <- diag(vcov(cmp))
print(range(v))
stopifnot(all(is.finite(v)), all(v > 0))

set.seed(123)
results <- as.data.frame(summary(
  cmp,
  type = "response",
  infer = c(TRUE, TRUE),
  adjust = "mvt"
))

write.csv(results, "knockdown_vs_AAVS1.csv", row.names = FALSE)
saveRDS(m_clean, "model_clean.rds")
cat("Analysis completed successfully.\n")