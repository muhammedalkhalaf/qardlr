# Script to generate example dataset for qardl package
# Run this script to regenerate data/qardl_sim.rda

set.seed(42)

n <- 200
burnin <- 100
n_total <- n + burnin

# True parameters
phi1 <- 0.4
phi2 <- 0.2
gamma1 <- 0.5
gamma2 <- 0.3
sigma_u <- 1
sigma_x <- 0.5

# Generate covariates (random walks)
x1 <- cumsum(rnorm(n_total, sd = sigma_x))
x2 <- cumsum(rnorm(n_total, sd = sigma_x))

# Generate dependent variable
y <- numeric(n_total)
u <- rnorm(n_total, sd = sigma_u)

# Initialize
y[1] <- u[1]
y[2] <- u[2]

# Generate QARDL(2,2) process
for (t in 3:n_total) {
  y[t] <- phi1 * y[t-1] + phi2 * y[t-2] +
          gamma1 * x1[t] + gamma2 * x2[t] + u[t]
}

# Remove burn-in
y <- y[(burnin + 1):n_total]
x1 <- x1[(burnin + 1):n_total]
x2 <- x2[(burnin + 1):n_total]

# Create data frame
qardl_sim <- data.frame(
  y = y,
  x1 = x1,
  x2 = x2
)

# Save to data/
save(qardl_sim, file = "data/qardl_sim.rda", compress = "xz")

cat("Dataset qardl_sim created with", nrow(qardl_sim), "observations\n")
cat("True parameters:\n")
cat("  phi1 =", phi1, ", phi2 =", phi2, "\n")
cat("  gamma1 =", gamma1, ", gamma2 =", gamma2, "\n")
cat("  Long-run beta1 =", gamma1 / (1 - phi1 - phi2), "\n")
cat("  Long-run beta2 =", gamma2 / (1 - phi1 - phi2), "\n")
