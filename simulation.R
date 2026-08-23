rm(list = ls())

library(dplyr)

# Monte Carlo comparison of GCV, RR, MLE, PRESS, Min Sol'n, and Min Data
# Based on ill-conditioned linear model with Laplace transform design

n = 21
p = 10

# choose 21 positive sample points s_i (equally spaced)
s = seq(0.1
        , 1.5, length.out = n)

# construct tridiagonal matrix for Gauss–Laguerre quadrature
alpha = numeric(p)
offdiag = numeric(p)
M = matrix(0, p, p)

for (i in 1:p){
  alpha[i] = 2*i - 1
  offdiag[i] = -i
  M[i,i] = alpha[i]
  if (i < p){
    M[i, i+1] = offdiag[i]
    M[i+1, i] = offdiag[i]
  }
}

# eigen decomposition to obtain quadrature nodes and weights
eig = eigen(M)

# nodes t_k (eigenvalues)
tk = eig$values

# weights w_k from first component of normalized eigenvectors
wk = eig$vectors[1, ]^2

# sort nodes and weights in increasing order
ord = order(tk)
tk = tk[ord]
wk = wk[ord]

# construct design matrix X using Laplace transform discretization
# X_{ik} = w_k * exp(t_k * (1 - s_i))
X = matrix(0, nrow = n, ncol = p)
for (i in 1:n){
  for (k in 1:p){
    X[i,k] = wk[k] * exp(tk[k] * (1 - s[i]))
  }
}

# quick verification: dimension, condition number, numerical range
dim(X)
kappa(X)
range(X)

# generate smooth decreasing coefficient vector beta
beta_true = exp(-0.5 * (0:(p-1)))

# compute initial signal g = X beta
g0 = X %*% beta_true

# rescale beta so that ||X beta||^2 matches the value reported in the paper
scale_factor = sqrt(370.84 / sum(g0^2))
beta_true = scale_factor * beta_true

# final signal
g = X %*% beta_true

# verify signal magnitude
sum(g^2)

# sigma squared choices
sigmasquared = c(1e-8, 1e-6, 1e-4, 1e-2)

# lambda choices
lambda_choices = 10^seq(-6, 2, length.out = 50)

# Identity matrix matching the column dimension of the input matrix
identity_matrix <- function(X){
  return(diag(ncol(X)))
}

# Estimate beta using ridge
ridge_beta <- function(lambda, X, y){
  n = nrow(X)
  beta = solve(t(X)%*%X + n*lambda*identity_matrix(X))%*%t(X)%*%y
  return(beta)
}

# Ridge hat matrix A(lambda)
A_lambda <- function(lambda, X){
  n = nrow(X)
  return(X%*%solve(t(X)%*%X + n*lambda*identity_matrix(X))%*%t(X))
}

# Trace function
tr <- function(M){
  return(sum(diag(M)))
}

# squared norm
norm_squared <- function(x){
  return(sum(x^2))
}
# GCV Function
gcv_value <- function(lambda, X, y){
  A = A_lambda(lambda, X)
  M = identity_matrix(A)-A
  n = nrow(X)
  return((1/n)*norm_squared(M%*%y)/(tr(M)/n)^2)
}

# Compute fitted signal y_hat = X beta_hat
y_estimate <- function(X, bhat){
  return(X%*%bhat)
}

# Domain Loss
domain_loss <- function(btrue, bhat){
  return(norm_squared(btrue-bhat))
}
# Range Loss
range_loss <- function(gtrue,yhat){
  n = length(gtrue)
  return((1/n)*norm_squared(gtrue-yhat))
}

# Residual Vector Function
residual_vector <- function(A,y){
  return((identity_matrix(A) - A)%*%y)
}

# Leverages Function
leverages <- function(A){
  h = numeric(nrow(A))
  for(j in 1:nrow(A)){
    h[j] = A[j,j]
  }
  return(h)
}

# PRESS-adjusted residuals
r_PRESS <- function(r, h){
  r_P = numeric(length(r))
  for(i in 1:length(r)){
    r_P[i] = r[i]/(1-h[i])
  }
  return(r_P)
}

# PRESS Value Function
press_value <- function(lambda, X, y){
  A = A_lambda(lambda,X)
  r = residual_vector(A,y)
  h = leverages(A)
  r_P = r_PRESS(r,h)
  n = nrow(X)
  return((1/n)*norm_squared(r_P))
}

mle_value <- function(lambda, X, y){
  A = A_lambda(lambda,X) 
  M = identity_matrix(A)-A 
  n = nrow(X) 
  return(as.numeric((1/n)*(t(y)%*%M%*%y)/((det(M))^(1/n)))) 
}

# Estimate sigma^2 using OLS residuals for the RR criterion
sigma_squared_estimate <- function(X, y){
  n = nrow(X)
  p = ncol(X)
  M = X%*%solve(t(X)%*%X)%*%t(X)
  return((1/(n-p))*norm_squared((identity_matrix(M)-M)%*%y))
}

# Range Risk (RR)
rr_value <- function(lambda, X, y){
  n = nrow(X)
  A = A_lambda(lambda,X)
  M = identity_matrix(A) - A
  sigma2_estimate = sigma_squared_estimate(X,y)
  return((1/n)*norm_squared(M%*%y)-(sigma2_estimate*tr(M)/n)+sigma2_estimate)
}

# Evaluate dataset
evaluate_dataset <- function(y,X,g,beta_true,lambda_choices){
  gcv_list = numeric(length(lambda_choices))
  rr_list = numeric(length(lambda_choices))
  mle_list = numeric(length(lambda_choices))
  press_list = numeric(length(lambda_choices))
  domain_loss_list = numeric(length(lambda_choices))
  range_loss_list = numeric(length(lambda_choices))

  for (i in 1:length(lambda_choices)){
    lambda = lambda_choices[i]

    beta_hat = ridge_beta(lambda, X, y)
    y_hat = y_estimate(X,beta_hat)

    gcv_list[i] = gcv_value(lambda, X, y)
    rr_list[i] = rr_value(lambda, X, y)
    mle_list[i] = mle_value(lambda, X, y)
    press_list[i] = press_value(lambda, X, y)
    domain_loss_list[i] = domain_loss(beta_true, beta_hat)
    range_loss_list[i] = range_loss(g,y_hat)
  }

  idx_gcv = which.min(gcv_list)
  idx_rr = which.min(rr_list)
  idx_mle = which.min(mle_list)
  idx_press = which.min(press_list)
  idx_min_soln = which.min(domain_loss_list)
  idx_min_data = which.min(range_loss_list)

  lambda_gcv = lambda_choices[idx_gcv]
  lambda_rr = lambda_choices[idx_rr]
  lambda_mle = lambda_choices[idx_mle]
  lambda_press = lambda_choices[idx_press]
  lambda_min_soln = lambda_choices[idx_min_soln]
  lambda_min_data = lambda_choices[idx_min_data]

  min_domain_loss = domain_loss_list[idx_min_soln]
  min_range_loss = range_loss_list[idx_min_data]

  ID_GCV = domain_loss_list[idx_gcv]/min_domain_loss
  IR_GCV = range_loss_list[idx_gcv]/min_range_loss

  ID_RR = domain_loss_list[idx_rr]/min_domain_loss
  IR_RR = range_loss_list[idx_rr]/min_range_loss

  ID_MLE = domain_loss_list[idx_mle]/min_domain_loss
  IR_MLE = range_loss_list[idx_mle]/min_range_loss

  ID_PRESS = domain_loss_list[idx_press]/min_domain_loss
  IR_PRESS = range_loss_list[idx_press]/min_range_loss

  ID_min_soln = domain_loss_list[idx_min_soln]/min_domain_loss
  IR_min_soln = range_loss_list[idx_min_soln]/min_range_loss

  ID_min_data = domain_loss_list[idx_min_data]/min_domain_loss
  IR_min_data = range_loss_list[idx_min_data]/min_range_loss


  return(c(lambda_gcv, lambda_rr,lambda_mle, lambda_press,lambda_min_soln, lambda_min_data,
           ID_GCV, IR_GCV,
           ID_RR, IR_RR,
           ID_MLE, IR_MLE,
           ID_PRESS, IR_PRESS,
           ID_min_soln, IR_min_soln,
           ID_min_data, IR_min_data))
}

# single-run results across all four sigma^2 values
results = matrix(nrow = 4, ncol = 18)
set.seed(2026)
for (i in 1:length(sigmasquared)){
  eps = rnorm(n, 0, sqrt(sigmasquared[i]))
  y = g + eps
  results[i,] = evaluate_dataset(y,X,g,beta_true, lambda_choices)
}
colnames(results) <- c(
  "lambda_gcv",
  "lambda_rr",
  "lambda_mle",
  "lambda_press",
  "lambda_min_soln",
  "lambda_min_data",
  "ID_GCV",
  "IR_GCV",
  "ID_RR",
  "IR_RR",
  "ID_MLE",
  "IR_MLE",
  "ID_PRESS",
  "IR_PRESS",
  "ID_min_soln",
  "IR_min_soln",
  "ID_min_data",
  "IR_min_data"
)
rownames(results) <- paste0("sigma2=", sigmasquared)

# monte carlo simulation
monte_carlo <- function(B,sigmasquared,X,g,beta_true,lambda_choices){
  n = nrow(X)
  results = matrix(nrow = B*length(sigmasquared), ncol = 18)
  for(j in 1:B){
    set.seed(2026 + j)
    for (i in 1:length(sigmasquared)){
      eps = rnorm(n, 0, sqrt(sigmasquared[i]))
      y = g + eps
      results[length(sigmasquared)*(j - 1) + i,] = evaluate_dataset(y,X,g,beta_true, lambda_choices)
    }
  }
  return(results)
}



# 100 Monte Carlo Simulations

# Monte Carlo results across B replications
B = 100 # simulations
mc_results = monte_carlo(B,sigmasquared,X,g,beta_true, lambda_choices)
# Naming columns and rows
colnames(mc_results) <- c(
  "lambda_gcv",
  "lambda_rr",
  "lambda_mle",
  "lambda_press",
  "lambda_min_soln",
  "lambda_min_data",
  "ID_GCV",
  "IR_GCV",
  "ID_RR",
  "IR_RR",
  "ID_MLE",
  "IR_MLE",
  "ID_PRESS",
  "IR_PRESS",
  "ID_min_soln",
  "IR_min_soln",
  "ID_min_data",
  "IR_min_data"
)

rep_id = rep(1:B, each = length(sigmasquared))
sigma_id = rep(sigmasquared, times = B)

# Put it into a dataframe
mc_results = as.data.frame(mc_results)

mc_results$replication = rep_id
mc_results$sigma2 = sigma_id

# Keep only replication, sigma2, and inefficiency measures
mc_results = mc_results[, c(
  "replication", "sigma2",
  "ID_GCV", "IR_GCV",
  "ID_RR", "IR_RR",
  "ID_MLE", "IR_MLE",
  "ID_PRESS", "IR_PRESS",
  "ID_min_soln", "IR_min_soln",
  "ID_min_data", "IR_min_data"
)]

# Summary table for mean of 100 simulations

summary_table <- mc_results %>%
  group_by(sigma2) %>%
  summarise(
    ID_GCV = mean(ID_GCV), IR_GCV = mean(IR_GCV),
    ID_RR = mean(ID_RR), IR_RR = mean(IR_RR),
    ID_MLE = mean(ID_MLE), IR_MLE = mean(IR_MLE),
    ID_PRESS = mean(ID_PRESS), IR_PRESS = mean(IR_PRESS),
    ID_min_soln = mean(ID_min_soln), IR_min_soln = mean(IR_min_soln),
    ID_min_data = mean(ID_min_data), IR_min_data = mean(IR_min_data)
  )

# Create one compact summary table from summary_table
make_summary_compact <- function(summary_table) {
  summary_compact <- data.frame(
    sigma2 = character(),
    Method = character(),
    ID = numeric(),
    IR = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(summary_table)) {
    row <- summary_table[i, ]

    temp <- data.frame(
      sigma2 = rep(as.character(row$sigma2), 6),
      Method = c("GCV", "RR", "MLE", "PRESS", "Min Sol'n", "Min Data"),
      ID = c(
        row$ID_GCV,
        row$ID_RR,
        row$ID_MLE,
        row$ID_PRESS,
        row$ID_min_soln,
        row$ID_min_data
      ),
      IR = c(
        row$IR_GCV,
        row$IR_RR,
        row$IR_MLE,
        row$IR_PRESS,
        row$IR_min_soln,
        row$IR_min_data
      ),
      stringsAsFactors = FALSE
    )

    summary_compact <- rbind(summary_compact, temp)
  }

  summary_compact$ID <- round(summary_compact$ID, 2)
  summary_compact$IR <- round(summary_compact$IR, 2)

  return(summary_compact)
}

# Build the compact summary table
summary_compact <- make_summary_compact(summary_table)

# View it
summary_compact

# Optional: save it
# write.csv(summary_compact, "summary_compact.csv", row.names = FALSE)

make_summary_latex_wide <- function(summary_compact) {
  # Keep methods in desired order
  method_order <- c("GCV", "RR", "MLE", "PRESS", "Min Sol'n", "Min Data")
  
  # Keep sigma order
  sigma_order <- c("1e-08", "1e-06", "1e-04", "0.01")
  
  # Build wide table manually
  wide_tab <- data.frame(Method = method_order, stringsAsFactors = FALSE)
  
  for (sigma_val in sigma_order) {
    temp <- summary_compact[summary_compact$sigma2 == sigma_val, c("Method", "ID", "IR")]
    temp <- temp[match(method_order, temp$Method), ]
    
    wide_tab[[paste0("ID_", sigma_val)]] <- temp$ID
    wide_tab[[paste0("IR_", sigma_val)]] <- temp$IR
  }
  
  cat("\\begin{table}[htbp]\n")
  cat("\\centering\n")
  cat("\\caption{Mean inefficiency measures ($I_D$, $I_R$) over 100 Monte Carlo replications}\n")
  cat("\\begin{tabular}{lcccccccc}\n")
  cat("\\toprule\n")
  cat(" & \\multicolumn{2}{c}{$\\sigma^2 = 10^{-8}$}")
  cat(" & \\multicolumn{2}{c}{$\\sigma^2 = 10^{-6}$}")
  cat(" & \\multicolumn{2}{c}{$\\sigma^2 = 10^{-4}$}")
  cat(" & \\multicolumn{2}{c}{$\\sigma^2 = 10^{-2}$} \\\\\n")
  cat("Method & $I_D$ & $I_R$ & $I_D$ & $I_R$ & $I_D$ & $I_R$ & $I_D$ & $I_R$ \\\\\n")
  cat("\\midrule\n")
  
  for (i in 1:nrow(wide_tab)) {
    row <- wide_tab[i, ]
    cat(
      row$Method, " & ",
      formatC(row[[2]], format = "f", digits = 2), " & ",
      formatC(row[[3]], format = "f", digits = 2), " & ",
      formatC(row[[4]], format = "f", digits = 2), " & ",
      formatC(row[[5]], format = "f", digits = 2), " & ",
      formatC(row[[6]], format = "f", digits = 2), " & ",
      formatC(row[[7]], format = "f", digits = 2), " & ",
      formatC(row[[8]], format = "f", digits = 2), " & ",
      formatC(row[[9]], format = "f", digits = 2), " \\\\\n",
      sep = ""
    )
  }
  
  cat("\\bottomrule\n")
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
}

# Print Out
make_summary_latex_wide(summary_compact)



# Generate plots for report

# install.packages(c("dplyr", "tidyr", "ggplot2", "patchwork", "scales", "cowplot"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)


# 1. Setup

mc_results$sigma2 <- as.numeric(mc_results$sigma2)

method_levels <- c("GCV", "MLE", "PRESS", "RR")

method_colors <- c(
  "GCV"   = "#1B9E77",
  "MLE"   = "#D95F02",
  "PRESS" = "#7570B3",
  "RR"    = "#E7298A"
)

sigma_breaks <- c(1e-8, 1e-6, 1e-4, 1e-2)

# 2. Prepare data

# Mean ID
mean_id_plot_data <- mc_results %>%
  group_by(sigma2) %>%
  summarise(
    GCV   = mean(ID_GCV),
    MLE   = mean(ID_MLE),
    PRESS = mean(ID_PRESS),
    RR    = mean(ID_RR),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -sigma2,
    names_to = "Method",
    values_to = "Mean_ID"
  )

# Mean IR
mean_ir_plot_data <- mc_results %>%
  group_by(sigma2) %>%
  summarise(
    GCV   = mean(IR_GCV),
    MLE   = mean(IR_MLE),
    PRESS = mean(IR_PRESS),
    RR    = mean(IR_RR),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -sigma2,
    names_to = "Method",
    values_to = "Mean_IR"
  )

# Boxplot ID (sigma^2 = 1e-4)
box_id_data <- mc_results %>%
  filter(sigma2 == 1e-4) %>%
  select(replication, ID_GCV, ID_MLE, ID_PRESS, ID_RR) %>%
  pivot_longer(
    cols = -replication,
    names_to = "Method",
    values_to = "ID"
  )

# Boxplot IR (sigma^2 = 1e-4)
box_ir_data <- mc_results %>%
  filter(sigma2 == 1e-4) %>%
  select(replication, IR_GCV, IR_MLE, IR_PRESS, IR_RR) %>%
  pivot_longer(
    cols = -replication,
    names_to = "Method",
    values_to = "IR"
  )


# 3. Clean method names and order
mean_id_plot_data$Method <- factor(mean_id_plot_data$Method, levels = method_levels)
mean_ir_plot_data$Method <- factor(mean_ir_plot_data$Method, levels = method_levels)

box_id_data$Method <- factor(
  gsub("^ID_", "", box_id_data$Method),
  levels = method_levels
)

box_ir_data$Method <- factor(
  gsub("^IR_", "", box_ir_data$Method),
  levels = method_levels
)


# 4. Theme

common_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(size = 24, face = "plain", hjust = 0),
    axis.title = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 22)
  )

# 5. Plots

sigma_breaks <- c(1e-8, 1e-6, 1e-4, 1e-2)
sigma_labels <- c(expression(10^-8), expression(10^-6), expression(10^-4), expression(10^-2))

# (a) Mean ID
p1 <- ggplot(mean_id_plot_data, aes(x = sigma2, y = Mean_ID, color = Method)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  scale_x_log10(
    breaks = sigma_breaks,
    labels = sigma_labels,
    minor_breaks = NULL
  ) +
  scale_color_manual(values = method_colors) +
  labs(
    title = "(a) Mean Domain Inefficiency",
    x = expression(sigma^2),
    y = expression(I[D]),
    color = "Method"
  ) +
  common_theme

# (b) Mean IR
p2 <- ggplot(mean_ir_plot_data, aes(x = sigma2, y = Mean_IR, color = Method)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  scale_x_log10(
    breaks = sigma_breaks,
    labels = sigma_labels,
    minor_breaks = NULL
  ) +
  scale_color_manual(values = method_colors) +
  labs(
    title = "(b) Mean Range Inefficiency",
    x = expression(sigma^2),
    y = expression(I[R]),
    color = "Method"
  ) +
  common_theme

# (c) Box ID
p3 <- ggplot(box_id_data, aes(x = Method, y = ID, fill = Method)) +
  geom_boxplot(alpha = 0.75) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = expression("(c) Distribution of " * I[D] * " at " * sigma^2 == 10^{-4}),
    x = "Method",
    y = expression(I[D])
  ) +
  common_theme +
  theme(legend.position = "none")

# (d) Box IR
p4 <- ggplot(box_ir_data, aes(x = Method, y = IR, fill = Method)) +
  geom_boxplot(alpha = 0.75) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = expression("(d) Distribution of " * I[R] * " at " * sigma^2 == 10^{-4}),
    x = "Method",
    y = expression(I[R])
  ) +
  common_theme +
  theme(legend.position = "none")

# 6. Combine 2x2
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")

legend_plot <- ggplot(
  mean_id_plot_data,
  aes(x = sigma2, y = Mean_ID, color = Method)
) +
  geom_line(linewidth = 2) +
  geom_point(size = 4) +
  scale_color_manual(values = method_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.8, "cm")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 2.5,
        size = 4
      ),
      nrow = 1
    )
  )
library(grid)
library(cowplot)

legend <- cowplot::get_legend(legend_plot)
library(patchwork)

# 2x2 plot
combined_plot <- (p1 | p2) / (p3 | p4)

# Extract legend
legend <- cowplot::get_legend(legend_plot)

# Wrap both pieces as elements, then stack
final_plot <- patchwork::wrap_elements(combined_plot) /
  patchwork::wrap_elements(legend) +
  plot_layout(heights = c(10, 1))

# Show final
final_plot

# Save final
ggsave(
  "combined_monte_carlo_figure.png",
  final_plot,
  width = 12,
  height = 9,
  dpi = 500
)
