# ============================================================================
# MATRIX EXPONENTIAL AND LIE GROUPS - R IMPLEMENTATION
# MATH 7370 Linear Algebra and Matrix Analysis
# ============================================================================


# Vector of required packages
#packages <- c("expm", "matlib", "pracma", "ggplot2", "plotly", "gridExtra")

# Install any that are missing
#install.packages(setdiff(packages, rownames(installed.packages())))

# Load them
#lapply(packages, library, character.only = TRUE)

# Required libraries
library(expm)      # Matrix exponential
library(matlib)    # Matrix utilities
library(pracma)    # Practical math functions
library(ggplot2)   # Plotting
library(plotly)    # 3D plotting
library(gridExtra) # Multiple plots

# ============================================================================
# 1. MATRIX EXPONENTIAL ALGORITHMS
# ============================================================================

# Method 1: Direct power series
matrix_exp_series <- function(X, terms = 20) {
  n <- nrow(X)
  result <- diag(n)
  term <- diag(n)
  
  for (m in 1:terms) {
    term <- (term %*% X) / m
    result <- result + term
  }
  
  return(result)
}

# Method 2: Eigenvalue decomposition
matrix_exp_eigen <- function(X) {
  eigen_decomp <- eigen(X)
  V <- eigen_decomp$vectors
  D <- diag(exp(eigen_decomp$values))
  
  result <- V %*% D %*% solve(V)
  return(Re(result))  # Return real part if complex
}

# Method 3: Compare all methods
compare_methods <- function(X) {
  cat("=== Matrix Exponential Methods Comparison ===\n\n")
  
  # Method 1: Direct series
  t1 <- system.time({
    exp1 <- matrix_exp_series(X)
  })
  
  # Method 2: Eigenvalue
  t2 <- system.time({
    exp2 <- matrix_exp_eigen(X)
  })
  
  # Method 3: expm package (Padé approximation)
  t3 <- system.time({
    exp3 <- expm(X)
  })
  
  cat(sprintf("Series method: %.6f seconds\n", t1[3]))
  cat(sprintf("Eigen method:  %.6f seconds\n", t2[3]))
  cat(sprintf("expm package:  %.6f seconds\n", t3[3]))
  cat(sprintf("\nMax difference (Series vs expm): %.6e\n", max(abs(exp1 - exp3))))
  cat(sprintf("Max difference (Eigen vs expm):  %.6e\n", max(abs(exp2 - exp3))))
  
  return(list(series = exp1, eigen = exp2, expm = exp3))
}

# ============================================================================
# 2. VERIFY THEOREM 16.15 PROPERTIES
# ============================================================================

verify_theorem_16_15 <- function() {
  cat("\n=== Theorem 16.15 Properties Verification ===\n\n")
  
  # Property 1: exp(0) = I
  cat("Property 1: exp(0) = I\n")
  zero_mat <- matrix(0, 3, 3)
  cat(sprintf("Verified: %s\n\n", all.equal(expm(zero_mat), diag(3))))
  
  # Property 2: exp(X^T) = (exp(X))^T
  set.seed(123)
  X <- matrix(rnorm(9), 3, 3)
  cat("Property 2: exp(X^T) = (exp(X))^T\n")
  cat(sprintf("Verified: %s\n\n", 
              isTRUE(all.equal(expm(t(X)), t(expm(X))))))
  
  # Property 3: Adjoint (conjugate transpose)
  cat("Property 3: exp(X*) = (exp(X))*\n")
  cat(sprintf("Verified: %s\n\n", 
              isTRUE(all.equal(expm(Conj(t(X))), Conj(t(expm(X)))))))
  
  # Property 4: Conjugation A·exp(X)·A^(-1) = exp(A·X·A^(-1))
  A <- matrix(rnorm(9), 3, 3)
  A_inv <- solve(A)
  cat("Property 4: A·exp(X)·A^(-1) = exp(A·X·A^(-1))\n")
  lhs <- A %*% expm(X) %*% A_inv
  rhs <- expm(A %*% X %*% A_inv)
  cat(sprintf("Verified: %s\n\n", isTRUE(all.equal(lhs, rhs, tolerance = 1e-10))))
  
  # Property 5: det(exp(X)) = exp(tr(X))
  cat("Property 5: det(exp(X)) = exp(tr(X))\n")
  det_expX <- det(expm(X))
  exp_trX <- exp(sum(diag(X)))
  cat(sprintf("det(exp(X)) = %.6f\n", det_expX))
  cat(sprintf("exp(tr(X))  = %.6f\n", exp_trX))
  cat(sprintf("Verified: %s\n\n", isTRUE(all.equal(det_expX, exp_trX))))
  
  # Property 6: Commuting matrices
  X_diag <- diag(c(1, 2, 3))
  Y_diag <- diag(c(4, 5, 6))
  cat("Property 6: exp(X+Y) = exp(X)·exp(Y) when [X,Y]=0\n")
  cat(sprintf("Verified: %s\n\n", 
              isTRUE(all.equal(expm(X_diag + Y_diag), 
                               expm(X_diag) %*% expm(Y_diag)))))
  
  # Property 7: Lie product formula
  X2 <- matrix(rnorm(4), 2, 2)
  Y2 <- matrix(rnorm(4), 2, 2)
  m <- 100
  cat("Property 7: Lie product formula (m=100)\n")
  lie_product <- Reduce(`%*%`, replicate(m, expm(X2/m) %*% expm(Y2/m), 
                                         simplify = FALSE))
  direct <- expm(X2 + Y2)
  error <- norm(lie_product - direct, "F")
  cat(sprintf("Error: %.6e\n", error))
}

# ============================================================================
# 3. BAKER-CAMPBELL-HAUSDORFF FORMULA
# ============================================================================

# Lie bracket (commutator)
commutator <- function(X, Y) {
  return(X %*% Y - Y %*% X)
}

# BCH expansion
bch_expansion <- function(X, Y, order = 4) {
  Z <- X + Y
  
  if (order >= 2) {
    Z <- Z + 0.5 * commutator(X, Y)
  }
  
  if (order >= 3) {
    Z <- Z + (1/12) * commutator(X, commutator(X, Y))
    Z <- Z - (1/12) * commutator(Y, commutator(X, Y))
  }
  
  if (order >= 4) {
    Z <- Z - (1/24) * commutator(Y, commutator(X, commutator(X, Y)))
  }
  
  return(Z)
}

# Test BCH formula
test_bch_formula <- function() {
  cat("\n=== Baker-Campbell-Hausdorff Formula ===\n\n")
  
  # Use small matrices for convergence
  scale <- 0.1
  set.seed(456)
  X <- scale * matrix(rnorm(9), 3, 3)
  Y <- scale * matrix(rnorm(9), 3, 3)
  
  # Direct computation
  direct <- expm(X) %*% expm(Y)
  
  cat(sprintf("Scale factor: %.2f\n", scale))
  cat(sprintf("||[X,Y]||: %.6e\n\n", norm(commutator(X, Y), "F")))
  
  errors <- numeric(4)
  for (order in 1:4) {
    Z <- bch_expansion(X, Y, order)
    approx <- expm(Z)
    error <- norm(direct - approx, "F")
    errors[order] <- error
    cat(sprintf("Order %d: Error = %.6e\n", order, error))
  }
  
  # Error vs scale
  cat("\nError vs Scale:\n")
  scales <- c(0.01, 0.05, 0.1, 0.2, 0.5)
  for (s in scales) {
    X_s <- s * X / scale
    Y_s <- s * Y / scale
    Z4 <- bch_expansion(X_s, Y_s, order = 4)
    err <- norm(expm(X_s) %*% expm(Y_s) - expm(Z4), "F")
    cat(sprintf("Scale %.2f: Error = %.6e\n", s, err))
  }
  
  return(errors)
}

# ============================================================================
# 4. SO(3) - ROTATION GROUP
# ============================================================================

# Convert vector to skew-symmetric matrix
skew_symmetric <- function(v) {
  matrix(c(0, v[3], -v[2],
           -v[3], 0, v[1],
           v[2], -v[1], 0), 3, 3, byrow = TRUE)
}

# Verify SO(3) properties
verify_so3_properties <- function() {
  cat("\n=== SO(3) Properties ===\n\n")
  
  set.seed(789)
  v <- rnorm(3)
  X <- skew_symmetric(v)
  R <- expm(X)
  
  cat(sprintf("R^T·R = I: %s\n", 
              isTRUE(all.equal(t(R) %*% R, diag(3), tolerance = 1e-10))))
  cat(sprintf("det(R) = 1: %s\n", 
              isTRUE(all.equal(det(R), 1, tolerance = 1e-10))))
  cat(sprintf("X is skew-symmetric: %s\n", 
              isTRUE(all.equal(X, -t(X), tolerance = 1e-10))))
}

# Visualize SO(3) rotation
visualize_so3_rotation <- function() {
  # Rotation axis
  axis <- c(1, 1, 1) / sqrt(3)
  X <- skew_symmetric(axis)
  
  # Initial vector
  v0 <- c(1, 0, 0)
  
  # Generate rotation path
  t_values <- seq(0, 2*pi, length.out = 100)
  path <- t(sapply(t_values, function(t) {
    as.vector(expm(t * X) %*% v0)
  }))
  
  # Create sphere data
  theta <- seq(0, 2*pi, length.out = 50)
  phi <- seq(0, pi, length.out = 50)
  sphere_grid <- expand.grid(theta = theta, phi = phi)
  
  sphere_x <- sin(sphere_grid$phi) * cos(sphere_grid$theta)
  sphere_y <- sin(sphere_grid$phi) * sin(sphere_grid$theta)
  sphere_z <- cos(sphere_grid$phi)
  
  # Create 3D plot
  fig <- plot_ly() %>%
    add_surface(
      x = matrix(sphere_x, 50, 50),
      y = matrix(sphere_y, 50, 50),
      z = matrix(sphere_z, 50, 50),
      opacity = 0.2,
      colorscale = list(c(0, 1), c("lightblue", "lightblue")),
      showscale = FALSE,
      name = "Unit Sphere"
    ) %>%
    add_trace(
      type = "scatter3d",
      mode = "lines",
      x = path[,1], y = path[,2], z = path[,3],
      line = list(color = "red", width = 5),
      name = "Rotation Path"
    ) %>%
    add_trace(
      type = "scatter3d",
      mode = "markers",
      x = v0[1], y = v0[2], z = v0[3],
      marker = list(size = 8, color = "green"),
      name = "Start"
    ) %>%
    add_trace(
      type = "scatter3d",
      mode = "markers",
      x = path[100,1], y = path[100,2], z = path[100,3],
      marker = list(size = 8, color = "red"),
      name = "End"
    ) %>%
    layout(
      title = sprintf("SO(3) Rotation: exp(tX)v, axis=[%.2f,%.2f,%.2f]", 
                      axis[1], axis[2], axis[3]),
      scene = list(
        xaxis = list(title = "X"),
        yaxis = list(title = "Y"),
        zaxis = list(title = "Z"),
        aspectmode = "cube"
      )
    )
  
  return(fig)
}

# ============================================================================
# 5. SU(2) - SPECIAL UNITARY GROUP
# ============================================================================

# Pauli matrices
pauli_matrices <- function() {
  sigma_x <- matrix(c(0, 1, 1, 0), 2, 2)
  sigma_y <- matrix(c(0, -1i, 1i, 0), 2, 2)
  sigma_z <- matrix(c(1, 0, 0, -1), 2, 2)
  
  return(list(x = sigma_x, y = sigma_y, z = sigma_z))
}

# Verify SU(2) double cover
verify_su2_double_cover <- function() {
  cat("\n=== SU(2) Double Cover Verification ===\n\n")
  
  sigma <- pauli_matrices()
  
  # X ∈ su(2): traceless, anti-Hermitian
  X <- 1i * sigma$x / 2
  
  cat(sprintf("X is traceless: %s\n", 
              isTRUE(all.equal(sum(diag(X)), 0+0i, tolerance = 1e-10))))
  cat(sprintf("X is anti-Hermitian: %s\n\n", 
              isTRUE(all.equal(X, -Conj(t(X)), tolerance = 1e-10))))
  
  # Compute exp(2πX)
  exp_2pi <- expm(2*pi*X)
  
  cat("exp(2π·X):\n")
  print(exp_2pi)
  cat(sprintf("\nexp(2π·X) = -I: %s\n", 
              isTRUE(all.equal(exp_2pi, -diag(2), tolerance = 1e-10))))
  
  # Show exp(4πX) = I
  exp_4pi <- expm(4*pi*X)
  cat(sprintf("exp(4π·X) = I: %s\n", 
              isTRUE(all.equal(exp_4pi, diag(2), tolerance = 1e-10))))
  
  # Plot determinant evolution
  t_values <- seq(0, 4*pi, length.out = 100)
  dets <- sapply(t_values, function(t) det(expm(t*X)))
  
  df <- data.frame(
    t_over_pi = t_values / pi,
    real_part = Re(dets),
    imag_part = Im(dets)
  )
  
  p <- ggplot(df) +
    geom_line(aes(x = t_over_pi, y = real_part, color = "Real part"), size = 1) +
    geom_line(aes(x = t_over_pi, y = imag_part, color = "Imag part"), size = 1) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.5, color = "red") +
    geom_vline(xintercept = c(2, 4), linetype = "dashed", alpha = 0.5, 
               color = c("green", "blue")) +
    labs(title = "SU(2): det(exp(tX)) along path",
         x = "t/π", y = "det(exp(tX))",
         color = "") +
    theme_minimal() +
    theme(legend.position = "top")
  
  return(p)
}

# ============================================================================
# 6. QUANTUM ROTATION OPERATORS
# ============================================================================

# Quantum rotation operator
rotation_operator <- function(theta, n) {
  # Normalize n
  n <- n / sqrt(sum(n^2))
  
  sigma <- pauli_matrices()
  
  # J·n = (σ·n)/2
  J_dot_n <- (n[1]*sigma$x + n[2]*sigma$y + n[3]*sigma$z) / 2
  
  # U = exp(-i·θ·J·n)
  U <- expm(-1i * theta * J_dot_n)
  
  return(U)
}

# Verify quantum rotation properties
verify_rotation_properties <- function() {
  cat("\n=== Quantum Rotation Properties ===\n\n")
  
  theta <- pi / 4
  n <- c(1, 0, 0)
  
  U <- rotation_operator(theta, n)
  
  cat(sprintf("U is unitary (U†U = I): %s\n", 
              isTRUE(all.equal(Conj(t(U)) %*% U, diag(2), tolerance = 1e-10))))
  cat(sprintf("det(U) = 1: %s\n", 
              isTRUE(all.equal(det(U), 1+0i, tolerance = 1e-10))))
  
  # Verify composition
  U1 <- rotation_operator(pi/3, n)
  U2 <- rotation_operator(pi/6, n)
  U_combined <- rotation_operator(pi/2, n)
  
  cat("\nRotation composition:\n")
  cat(sprintf("U(π/3)·U(π/6) = U(π/2): %s\n", 
              isTRUE(all.equal(U1 %*% U2, U_combined, tolerance = 1e-10))))
  
  # Verify action on spin states
  cat("\nAction on spin states:\n")
  spin_up <- matrix(c(1, 0), 2, 1) + 0i
  spin_down <- matrix(c(0, 1), 2, 1) + 0i
  
  U_flip <- rotation_operator(pi, c(1, 0, 0))
  result <- U_flip %*% spin_up
  expected <- 1i * spin_down
  cat(sprintf("|↑⟩ → i|↓⟩: %s\n", 
              isTRUE(all.equal(result, expected, tolerance = 1e-10))))
}

# ============================================================================
# 7. LIE BRACKET VERIFICATION
# ============================================================================

verify_lie_bracket_properties <- function() {
  cat("\n=== Lie Bracket Properties ===\n\n")
  
  set.seed(321)
  X <- matrix(rnorm(9), 3, 3)
  Y <- matrix(rnorm(9), 3, 3)
  Z <- matrix(rnorm(9), 3, 3)
  
  # Antisymmetry
  antisym <- isTRUE(all.equal(commutator(X, Y), -commutator(Y, X)))
  cat(sprintf("Antisymmetry [X,Y] = -[Y,X]: %s\n", antisym))
  
  # Jacobi identity
  term1 <- commutator(X, commutator(Y, Z))
  term2 <- commutator(Y, commutator(Z, X))
  term3 <- commutator(Z, commutator(X, Y))
  jacobi <- isTRUE(all.equal(term1 + term2 + term3, 
                             matrix(0, 3, 3), tolerance = 1e-10))
  cat(sprintf("Jacobi identity: %s\n", jacobi))
  
  # SO(3) closure
  cat("\nSO(3) Lie algebra closure:\n")
  v1 <- rnorm(3)
  v2 <- rnorm(3)
  X_so3 <- skew_symmetric(v1)
  Y_so3 <- skew_symmetric(v2)
  bracket <- commutator(X_so3, Y_so3)
  is_skew <- isTRUE(all.equal(bracket, -t(bracket), tolerance = 1e-10))
  cat(sprintf("[X,Y] is skew-symmetric: %s\n", is_skew))
  
  # Cross product correspondence
  expected <- skew_symmetric(pracma::cross(v1, v2))
  cat(sprintf("[skew(v1),skew(v2)] = skew(v1×v2): %s\n", 
              isTRUE(all.equal(bracket, expected, tolerance = 1e-10))))
}

# ============================================================================
# 8. COMPLETE TEST SUITE
# ============================================================================

run_complete_test_suite <- function() {
  cat(rep("=", 60), "\n", sep = "")
  cat("MATRIX EXPONENTIAL & LIE GROUPS - COMPLETE TEST SUITE\n")
  cat(rep("=", 60), "\n", sep = "")
  
  # Test 1: Algorithm comparison
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("TEST 1: Matrix Exponential Algorithms\n")
  cat(rep("=", 60), "\n", sep = "")
  set.seed(100)
  X_test <- matrix(rnorm(16), 4, 4)
  compare_methods(X_test)
  
  # Test 2: Theorem 16.15
  verify_theorem_16_15()
  
  # Test 3: BCH Formula
  test_bch_formula()
  
  # Test 4: SO(3)
  verify_so3_properties()
  
  # Test 5: SU(2)
  p_su2 <- verify_su2_double_cover()
  
  # Test 6: Quantum rotations
  verify_rotation_properties()
  
  # Test 7: Lie brackets
  verify_lie_bracket_properties()
  
  # Generate visualization
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("GENERATING VISUALIZATIONS...\n")
  cat(rep("=", 60), "\n", sep = "")
  
  fig_so3 <- visualize_so3_rotation()
  
  cat("\n✓ All tests completed!\n")
  cat("✓ Use print(p_su2) to view SU(2) plot\n")
  cat("✓ Use print(fig_so3) to view SO(3) 3D rotation\n")
  
  return(list(su2_plot = p_su2, so3_plot = fig_so3))
}

# ============================================================================
# RUN ALL TESTS
# ============================================================================

# Execute the complete test suite
results <- run_complete_test_suite()

# Display plots
print(results$su2_plot)
print(results$so3_plot)
