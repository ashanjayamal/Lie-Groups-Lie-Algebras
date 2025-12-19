# ============================================================================
# TEST RUNNER AND VERIFICATION SCRIPT
# This script runs all tests and checks answers systematically
# MATH 7370 - Matrix Exponential and Lie Groups Project
# ============================================================================

# Clear workspace
rm(list = ls())

# Load required libraries
required_packages <- c("expm", "matlib", "pracma", "ggplot2", "plotly", "gridExtra")

cat("Checking and installing required packages...\n")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
cat("✓ All packages loaded successfully!\n\n")

# Source the main implementation file
cat("Loading main implementation...\n")
if (file.exists("matrix_exponential.R")) {
  source("matrix_exponential.R")
  cat("✓ Main file loaded successfully!\n\n")
} else {
  stop("ERROR: matrix_exponential.R not found! Please ensure it's in the same directory.")
}

# ============================================================================
# TESTING FRAMEWORK
# ============================================================================

# Counter for test results
test_counter <- list(passed = 0, failed = 0, total = 0)

# Test assertion function
assert_test <- function(condition, test_name, details = "") {
  test_counter$total <<- test_counter$total + 1
  
  if (isTRUE(condition)) {
    test_counter$passed <<- test_counter$passed + 1
    cat(sprintf("✓ PASS: %s\n", test_name))
    if (details != "") cat(sprintf("  Details: %s\n", details))
    return(TRUE)
  } else {
    test_counter$failed <<- test_counter$failed + 1
    cat(sprintf("✗ FAIL: %s\n", test_name))
    if (details != "") cat(sprintf("  Details: %s\n", details))
    return(FALSE)
  }
}

# Print test summary
print_test_summary <- function() {
  cat("\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("TEST SUMMARY\n")
  cat(rep("=", 60), "\n", sep = "")
  cat(sprintf("Total tests:  %d\n", test_counter$total))
  cat(sprintf("Passed:       %d (%.1f%%)\n", 
              test_counter$passed, 
              100 * test_counter$passed / test_counter$total))
  cat(sprintf("Failed:       %d (%.1f%%)\n", 
              test_counter$failed, 
              100 * test_counter$failed / test_counter$total))
  cat(rep("=", 60), "\n", sep = "")
  
  if (test_counter$failed == 0) {
    cat("\n🎉 ALL TESTS PASSED! 🎉\n\n")
  } else {
    cat("\n⚠️  SOME TESTS FAILED - Review output above ⚠️\n\n")
  }
}

# ============================================================================
# TEST SUITE 1: MATRIX EXPONENTIAL ALGORITHMS
# ============================================================================

test_matrix_exponential_algorithms <- function() {
  cat("\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("TEST SUITE 1: Matrix Exponential Algorithms\n")
  cat(rep("=", 60), "\n", sep = "")
  
  set.seed(12345)
  X <- matrix(rnorm(9), 3, 3)
  
  # Test 1.1: Series converges to expm
  exp_series <- matrix_exp_series(X, terms = 30)
  exp_exact <- expm(X)
  error_series <- max(abs(exp_series - exp_exact))
  assert_test(error_series < 1e-6, 
              "Series method converges to expm",
              sprintf("Error: %.2e", error_series))
  
  # Test 1.2: Eigenvalue method matches expm
  exp_eigen <- matrix_exp_eigen(X)
  error_eigen <- max(abs(exp_eigen - exp_exact))
  assert_test(error_eigen < 1e-6,
              "Eigenvalue method matches expm",
              sprintf("Error: %.2e", error_eigen))
  
  # Test 1.3: exp(0) = I
  zero_mat <- matrix(0, 3, 3)
  assert_test(all.equal(expm(zero_mat), diag(3)),
              "exp(0) = Identity matrix")
  
  # Test 1.4: Series method with diagonal matrix (known answer)
  D <- diag(c(1, 2, 3))
  exp_D_series <- matrix_exp_series(D, terms = 20)
  exp_D_exact <- diag(exp(c(1, 2, 3)))
  assert_test(max(abs(exp_D_series - exp_D_exact)) < 1e-10,
              "Series method correct for diagonal matrix",
              sprintf("||exp(D) - diag(e^λ)||: %.2e", 
                      max(abs(exp_D_series - exp_D_exact))))
  
  cat("\n")
}

# ============================================================================
# TEST SUITE 2: THEOREM 16.15 PROPERTIES
# ============================================================================

test_theorem_properties <- function() {
  cat(rep("=", 60), "\n", sep = "")
  cat("TEST SUITE 2: Theorem 16.15 Properties\n")
  cat(rep("=", 60), "\n", sep = "")
  
  set.seed(23456)
  X <- matrix(rnorm(9), 3, 3)
  
  # Test 2.1: exp(X^T) = (exp(X))^T
  assert_test(isTRUE(all.equal(expm(t(X)), t(expm(X)), tolerance = 1e-10)),
              "Property: exp(X^T) = (exp(X))^T")
  
  # Test 2.2: exp(X*) = (exp(X))*
  X_complex <- X + 1i * matrix(rnorm(9), 3, 3)
  assert_test(isTRUE(all.equal(expm(Conj(t(X_complex))), 
                               Conj(t(expm(X_complex))), 
                               tolerance = 1e-10)),
              "Property: exp(X*) = (exp(X))*")
  
  # Test 2.3: Conjugation property
  A <- matrix(rnorm(9), 3, 3)
  A_inv <- solve(A)
  lhs <- A %*% expm(X) %*% A_inv
  rhs <- expm(A %*% X %*% A_inv)
  assert_test(max(abs(lhs - rhs)) < 1e-9,
              "Property: A·exp(X)·A^(-1) = exp(A·X·A^(-1))",
              sprintf("Error: %.2e", max(abs(lhs - rhs))))
  
  # Test 2.4: det(exp(X)) = exp(tr(X))
  det_expX <- det(expm(X))
  exp_trX <- exp(sum(diag(X)))
  assert_test(abs(det_expX - exp_trX) < 1e-10,
              "Property: det(exp(X)) = exp(tr(X))",
              sprintf("det(exp(X))=%.6f, exp(tr(X))=%.6f", det_expX, exp_trX))
  
  # Test 2.5: Commuting matrices
  D1 <- diag(c(1, 2, 3))
  D2 <- diag(c(4, 5, 6))
  assert_test(max(abs(expm(D1 + D2) - expm(D1) %*% expm(D2))) < 1e-10,
              "Property: exp(X+Y)=exp(X)·exp(Y) when [X,Y]=0")
  
  # Test 2.6: Non-commuting matrices (should NOT satisfy product rule)
  X_nc <- matrix(c(1, 2, 0, 3), 2, 2)
  Y_nc <- matrix(c(0, 1, 1, 0), 2, 2)
  error_nc <- max(abs(expm(X_nc + Y_nc) - expm(X_nc) %*% expm(Y_nc)))
  assert_test(error_nc > 1e-5,
              "Non-commuting matrices: exp(X+Y) ≠ exp(X)·exp(Y)",
              sprintf("Error: %.2e (should be large)", error_nc))
  
  # Test 2.7: Lie product formula convergence
  m_values <- c(10, 50, 100)
  errors <- numeric(length(m_values))
  for (i in seq_along(m_values)) {
    m <- m_values[i]
    lie_prod <- Reduce(`%*%`, replicate(m, expm(X_nc/m) %*% expm(Y_nc/m), 
                                        simplify = FALSE))
    errors[i] <- norm(lie_prod - expm(X_nc + Y_nc), "F")
  }
  # Error should decrease as m increases
  assert_test(all(diff(errors) < 0),
              "Lie product formula: convergence as m increases",
              sprintf("Errors: m=10: %.2e, m=50: %.2e, m=100: %.2e", 
                      errors[1], errors[2], errors[3]))
  
  cat("\n")
}

# ============================================================================
# TEST SUITE 3: BAKER-CAMPBELL-HAUSDORFF FORMULA
# ============================================================================

test_bch_formula <- function() {
  cat(rep("=", 60), "\n", sep = "")
  cat("TEST SUITE 3: Baker-Campbell-Hausdorff Formula\n")
  cat(rep("=", 60), "\n", sep = "")
  
  # Test 3.1: Commutator definition
  set.seed(34567)
  X <- matrix(rnorm(4), 2, 2)
  Y <- matrix(rnorm(4), 2, 2)
  comm <- commutator(X, Y)
  expected_comm <- X %*% Y - Y %*% X
  assert_test(max(abs(comm - expected_comm)) < 1e-14,
              "Commutator [X,Y] = XY - YX")
  
  # Test 3.2: Commutator antisymmetry
  assert_test(max(abs(commutator(X, Y) + commutator(Y, X))) < 1e-14,
              "Commutator antisymmetry: [X,Y] = -[Y,X]")
  
  # Test 3.3: BCH for commuting matrices (reduces to X+Y)
  D1 <- diag(c(1, 2))
  D2 <- diag(c(3, 4))
  Z_bch <- bch_expansion(D1, D2, order = 4)
  assert_test(max(abs(Z_bch - (D1 + D2))) < 1e-14,
              "BCH reduces to X+Y for commuting matrices")
  
  # Test 3.4: BCH convergence with scale
  scale <- 0.05
  X_small <- scale * X
  Y_small <- scale * Y
  
  errors_by_order <- numeric(4)
  for (order in 1:4) {
    Z <- bch_expansion(X_small, Y_small, order)
    approx <- expm(Z)
    direct <- expm(X_small) %*% expm(Y_small)
    errors_by_order[order] <- norm(approx - direct, "F")
  }
  
  # Higher order should give smaller error
  assert_test(errors_by_order[4] < errors_by_order[2],
              "BCH: Higher order gives better approximation",
              sprintf("Order 2: %.2e, Order 4: %.2e", 
                      errors_by_order[2], errors_by_order[4]))
  
  # Test 3.5: BCH second-order term
  # For small X, Y: exp(X)exp(Y) ≈ exp(X + Y + [X,Y]/2)
  Z2 <- bch_expansion(X_small, Y_small, order = 2)
  expected_Z2 <- X_small + Y_small + 0.5 * commutator(X_small, Y_small)
  assert_test(max(abs(Z2 - expected_Z2)) < 1e-14,
              "BCH order 2 formula correct")
  
  cat("\n")
}

# ============================================================================
# TEST SUITE 4: SO(3) ROTATION GROUP
# ============================================================================

test_so3_group <- function() {
  cat(rep("=", 60), "\n", sep = "")
  cat("TEST SUITE 4: SO(3) Rotation Group\n")
  cat(rep("=", 60), "\n", sep = "")
  
  # Test 4.1: Skew-symmetric construction
  set.seed(45678)
  v <- rnorm(3)
  X <- skew_symmetric(v)
  assert_test(max(abs(X + t(X))) < 1e-14,
              "skew_symmetric(v) is skew-symmetric")
  
  # Test 4.2: exp(skew) ∈ SO(3) - orthogonality
  R <- expm(X)
  assert_test(max(abs(t(R) %*% R - diag(3))) < 1e-10,
              "exp(skew) is orthogonal: R^T·R = I",
              sprintf("||R^T·R - I||: %.2e", max(abs(t(R) %*% R - diag(3)))))
  
  # Test 4.3: exp(skew) ∈ SO(3) - determinant = 1
  det_R <- det(R)
  assert_test(abs(det_R - 1) < 1e-10,
              "exp(skew) has det = 1",
              sprintf("det(R) = %.10f", det_R))
  
  # Test 4.4: SO(3) preserves vector length
  v_test <- rnorm(3)
  v_rotated <- as.vector(R %*% v_test)
  len_original <- sqrt(sum(v_test^2))
  len_rotated <- sqrt(sum(v_rotated^2))
  assert_test(abs(len_original - len_rotated) < 1e-10,
              "SO(3) rotation preserves vector length",
              sprintf("Original: %.6f, Rotated: %.6f", len_original, len_rotated))
  
  # Test 4.5: Commutator of skew-symmetric is skew-symmetric
  v1 <- rnorm(3)
  v2 <- rnorm(3)
  X1 <- skew_symmetric(v1)
  X2 <- skew_symmetric(v2)
  comm <- commutator(X1, X2)
  assert_test(max(abs(comm + t(comm))) < 1e-14,
              "so(3) is closed under commutator")
  
  # Test 4.6: Commutator corresponds to cross product
  expected <- skew_symmetric(pracma::cross(v1, v2))
  assert_test(max(abs(comm - expected)) < 1e-14,
              "[skew(v1), skew(v2)] = skew(v1 × v2)")
  
  # Test 4.7: One-parameter subgroup property
  t1 <- 0.5
  t2 <- 0.3
  R1 <- expm(t1 * X)
  R2 <- expm(t2 * X)
  R_sum <- expm((t1 + t2) * X)
  assert_test(max(abs(R1 %*% R2 - R_sum)) < 1e-10,
              "One-parameter subgroup: exp(t1·X)·exp(t2·X) = exp((t1+t2)·X)")
  
  cat("\n")
}

# ============================================================================
# TEST SUITE 5: SU(2) SPECIAL UNITARY GROUP
# ============================================================================

test_su2_group <- function() {
  cat(rep("=", 60), "\n", sep = "")
  cat("TEST SUITE 5: SU(2) Special Unitary Group\n")
  cat(rep("=", 60), "\n", sep = "")
  
  sigma <- pauli_matrices()
  
  # Test 5.1: Pauli matrix properties
  assert_test(isTRUE(all.equal(sigma$x %*% sigma$x, diag(2))),
              "Pauli matrix: σ_x² = I")
  assert_test(isTRUE(all.equal(sigma$y %*% sigma$y, diag(2))),
              "Pauli matrix: σ_y² = I")
  assert_test(isTRUE(all.equal(sigma$z %*% sigma$z, diag(2))),
              "Pauli matrix: σ_z² = I")
  
  # Test 5.2: Pauli matrices are traceless
  assert_test(abs(sum(diag(sigma$x))) < 1e-14 && 
                abs(sum(diag(sigma$y))) < 1e-14 && 
                abs(sum(diag(sigma$z))) < 1e-14,
              "Pauli matrices are traceless")
  
  # Test 5.3: Pauli commutation relations
  comm_xy <- commutator(sigma$x, sigma$y)
  expected_xy <- 2i * sigma$z
  assert_test(max(abs(comm_xy - expected_xy)) < 1e-14,
              "Pauli commutator: [σ_x, σ_y] = 2i·σ_z")
  
  # Test 5.4: su(2) element is anti-Hermitian and traceless
  X_su2 <- 1i * sigma$x / 2
  assert_test(max(abs(X_su2 + Conj(t(X_su2)))) < 1e-14,
              "su(2) element is anti-Hermitian")
  assert_test(abs(sum(diag(X_su2))) < 1e-14,
              "su(2) element is traceless")
  
  # Test 5.5: exp(su(2)) ∈ SU(2) - unitarity
  U <- expm(X_su2)
  assert_test(max(abs(Conj(t(U)) %*% U - diag(2))) < 1e-10,
              "exp(su(2)) is unitary: U†·U = I")
  
  # Test 5.6: exp(su(2)) ∈ SU(2) - determinant = 1
  assert_test(abs(det(U) - 1) < 1e-10,
              "exp(su(2)) has det = 1",
              sprintf("det(U) = %.10f", abs(det(U))))
  
  # Test 5.7: Double cover property exp(2π·X) = -I
  exp_2pi <- expm(2*pi*X_su2)
  assert_test(max(abs(exp_2pi + diag(2))) < 1e-10,
              "SU(2) double cover: exp(2π·X) = -I")
  
  # Test 5.8: Full rotation exp(4π·X) = I
  exp_4pi <- expm(4*pi*X_su2)
  assert_test(max(abs(exp_4pi - diag(2))) < 1e-10,
              "SU(2): exp(4π·X) = I")
  
  cat("\n")
}

# ============================================================================
# TEST SUITE 6: QUANTUM ROTATION OPERATORS
# ============================================================================

test_quantum_rotations <- function() {
  cat(rep("=", 60), "\n", sep = "")
  cat("TEST SUITE 6: Quantum Rotation Operators\n")
  cat(rep("=", 60), "\n", sep = "")
  
  # Test 6.1: Rotation operator is unitary
  theta <- pi/4
  n <- c(1, 0, 0)
  U <- rotation_operator(theta, n)
  assert_test(max(abs(Conj(t(U)) %*% U - diag(2))) < 1e-10,
              "Rotation operator is unitary")
  
  # Test 6.2: Rotation operator has det = 1
  assert_test(abs(det(U) - 1) < 1e-10,
              "Rotation operator has det = 1")
  
  # Test 6.3: Rotation composition
  U1 <- rotation_operator(pi/3, n)
  U2 <- rotation_operator(pi/6, n)
  U_combined <- rotation_operator(pi/2, n)
  assert_test(max(abs(U1 %*% U2 - U_combined)) < 1e-10,
              "Rotation composition: U(θ1)·U(θ2) = U(θ1+θ2)")
  
  # Test 6.4: Identity rotation
  U_zero <- rotation_operator(0, n)
  assert_test(max(abs(U_zero - diag(2))) < 1e-10,
              "Zero rotation gives identity: U(0) = I")
  
  # Test 6.5: Full rotation around axis
  U_2pi <- rotation_operator(2*pi, n)
  assert_test(max(abs(U_2pi + diag(2))) < 1e-10,
              "2π rotation: U(2π) = -I (spinor property)")
  
  # Test 6.6: 4π rotation gives identity
  U_4pi <- rotation_operator(4*pi, n)
  assert_test(max(abs(U_4pi - diag(2))) < 1e-10,
              "4π rotation: U(4π) = I")
  
  # Test 6.7: Rotation around different axes
  n_x <- c(1, 0, 0)
  n_y <- c(0, 1, 0)
  n_z <- c(0, 0, 1)
  
  U_x <- rotation_operator(pi, n_x)
  U_y <- rotation_operator(pi, n_y)
  U_z <- rotation_operator(pi, n_z)
  
  # All should be unitary
  assert_test(max(abs(Conj(t(U_x)) %*% U_x - diag(2))) < 1e-10 &&
                max(abs(Conj(t(U_y)) %*% U_y - diag(2))) < 1e-10 &&
                max(abs(Conj(t(U_z)) %*% U_z - diag(2))) < 1e-10,
              "Rotations around x, y, z axes are all unitary")
  
  cat("\n")
}

# ============================================================================
# TEST SUITE 7: LIE BRACKET (COMMUTATOR) PROPERTIES
# ============================================================================

test_lie_bracket <- function() {
  cat(rep("=", 60), "\n", sep = "")
  cat("TEST SUITE 7: Lie Bracket Properties\n")
  cat(rep("=", 60), "\n", sep = "")
  
  set.seed(56789)
  X <- matrix(rnorm(9), 3, 3)
  Y <- matrix(rnorm(9), 3, 3)
  Z <- matrix(rnorm(9), 3, 3)
  
  # Test 7.1: Antisymmetry
  assert_test(max(abs(commutator(X, Y) + commutator(Y, X))) < 1e-13,
              "Lie bracket antisymmetry: [X,Y] = -[Y,X]")
  
  # Test 7.2: Jacobi identity
  term1 <- commutator(X, commutator(Y, Z))
  term2 <- commutator(Y, commutator(Z, X))
  term3 <- commutator(Z, commutator(X, Y))
  jacobi_sum <- term1 + term2 + term3
  assert_test(max(abs(jacobi_sum)) < 1e-12,
              "Jacobi identity: [X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0",
              sprintf("||Sum||: %.2e", max(abs(jacobi_sum))))
  
  # Test 7.3: Bilinearity (scalar multiplication)
  a <- 2.5
  b <- -1.3
  lhs <- commutator(a*X, b*Y)
  rhs <- a*b*commutator(X, Y)
  assert_test(max(abs(lhs - rhs)) < 1e-13,
              "Lie bracket bilinearity: [aX,bY] = ab[X,Y]")
  
  # Test 7.4: [X,X] = 0
  assert_test(max(abs(commutator(X, X))) < 1e-13,
              "Self-commutator is zero: [X,X] = 0")
  
  # Test 7.5: Leibniz rule (derivation property)
  # [X, YZ] = [X,Y]Z + Y[X,Z]
  YZ <- Y %*% Z
  lhs_leibniz <- commutator(X, YZ)
  rhs_leibniz <- commutator(X, Y) %*% Z + Y %*% commutator(X, Z)
  assert_test(max(abs(lhs_leibniz - rhs_leibniz)) < 1e-12,
              "Leibniz rule: [X,YZ] = [X,Y]Z + Y[X,Z]")
  
  cat("\n")
}

# ============================================================================
# TEST SUITE 8: INTEGRATION TESTS
# ============================================================================

test_integration <- function() {
  cat(rep("=", 60), "\n", sep = "")
  cat("TEST SUITE 8: Integration Tests\n")
  cat(rep("=", 60), "\n", sep = "")
  
  # Test 8.1: SO(3) rotation of spin
  set.seed(67890)
  axis <- c(0, 0, 1)  # z-axis
  theta <- pi/2
  
  # Create SO(3) rotation
  X_so3 <- skew_symmetric(axis)
  R <- expm(theta * X_so3)
  
  # Rotate a vector
  v <- c(1, 0, 0)  # x-direction
  v_rotated <- as.vector(R %*% v)
  expected <- c(0, 1, 0)  # Should point in y-direction
  
  assert_test(max(abs(v_rotated - expected)) < 1e-10,
              "SO(3): π/2 rotation around z-axis maps x→y",
              sprintf("Result: [%.4f, %.4f, %.4f]", 
                      v_rotated[1], v_rotated[2], v_rotated[3]))
  
  # Test 8.2: Consistency between SU(2) and SO(3)
  # SU(2) is double cover of SO(3)
  sigma <- pauli_matrices()
  X_su2 <- 1i * sigma$z / 2
  U <- rotation_operator(theta, c(0, 0, 1))
  
  # Both should give same physical rotation (up to phase)
  assert_test(max(abs(Conj(t(U)) %*% U - diag(2))) < 1e-10,
              "SU(2) rotation is unitary")
  
  # Test 8.3: BCH and exponential consistency
  X_small <- 0.01 * matrix(rnorm(4), 2, 2)
  Y_small <- 0.01 * matrix(rnorm(4), 2, 2)
  
  # Method 1: Direct product
  direct <- expm(X_small) %*% expm(Y_small)
  
  # Method 2: BCH approximation
  Z_bch <- bch_expansion(X_small, Y_small, order = 3)
  approx <- expm(Z_bch)
  
  error_bch <- max(abs(direct - approx))
  assert_test(error_bch < 1e-6,
              "BCH approximation accurate for small matrices",
              sprintf("Error: %.2e", error_bch))
  
  cat("\n")
}

# ============================================================================
# MAIN TEST EXECUTION
# ============================================================================

run_all_tests <- function() {
  cat("\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("STARTING COMPREHENSIVE TEST SUITE\n")
  cat("Matrix Exponential and Lie Groups Project\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("\n")
  
  # Run all test suites
  test_matrix_exponential_algorithms()
  test_theorem_properties()
  test_bch_formula()
  test_so3_group()
  test_su2_group()
  test_quantum_rotations()
  test_lie_bracket()
  test_integration()
  
  # Print final summary
  print_test_summary()
  
  # Return results
  invisible(test_counter)
}

# ============================================================================
# EXECUTE ALL TESTS
# ============================================================================

cat("\n")
cat("="*60, "\n")
cat("TEST RUNNER SCRIPT LOADED\n")
cat("="*60, "\n")
cat("\nExecuting all tests...\n")

# Run the complete test suite
test_results <- run_all_tests()

# Save results to file
results_file <- "test_results.txt"
sink(results_file)
cat("MATRIX EXPONENTIAL TEST RESULTS\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("\nTotal Tests: %d\n", test_counter$total))
cat(sprintf("Passed: %d\n", test_counter$passed))
cat(sprintf("Failed: %d\n", test_counter$failed))
cat(sprintf("Success Rate: %.1f%%\n", 100 * test_counter$passed / test_counter$total))
sink()

cat(sprintf("\n✓ Test results saved to: %s\n", results_file))

