library("bookdown")
create_gitbook(".")
file.create(".nojekyll")
rmarkdown::render_site(encoding = 'UTF-8')

# Load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)

two_pi <- 2 * pi

# -----------------------------
# Right panel: dense diagonal lines
# -----------------------------
n_lines   <- 40
n_points  <- 200
theta_vals <- seq(0, two_pi, length.out = n_points)

lines_right <- lapply(seq_len(n_lines), function(i) {
  c0 <- runif(1, 0, two_pi)
  df_base <- data.frame(theta = theta_vals, phi = theta_vals + c0)
  df_wrap <- bind_rows(
    transform(df_base, phi = phi - two_pi),
    df_base,
    transform(df_base, phi = phi + two_pi)
  )
  df_wrap[df_wrap$phi >= 0 & df_wrap$phi <= two_pi, ]
}) %>% bind_rows(.id = "line_id")

p_right <- ggplot(lines_right, aes(x = theta, y = phi, group = line_id)) +
  geom_line(linewidth = 0.3, color = "black") +
  scale_x_continuous(limits = c(0, two_pi), breaks = c(0, pi, two_pi),
                     labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(limits = c(0, two_pi), breaks = c(0, pi, two_pi),
                     labels = c("0", expression(pi), expression(2*pi))) +
  coord_fixed() +
  labs(x = expression(theta), y = expression(phi)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

# -----------------------------
# Left panel: grid + parallelogram + arrow
# -----------------------------
grid_v <- data.frame(x = c(0, pi, two_pi))
grid_h <- data.frame(y = c(0, pi, two_pi))

offsets <- c(-pi, 0, pi)
lines_left <- lapply(offsets, function(c0) {
  df_base <- data.frame(theta = theta_vals, phi = theta_vals + c0)
  df_wrap <- bind_rows(
    transform(df_base, phi = phi - two_pi),
    df_base,
    transform(df_base, phi = phi + two_pi)
  )
  df_wrap[df_wrap$phi >= 0 & df_wrap$phi <= two_pi, ]
}) %>% bind_rows(.id = "line_id")

parallelogram <- data.frame(
  theta = c(0, pi, pi, 0, 0),
  phi   = c(0, pi, 2*pi, pi, 0)
)

arrow_df <- data.frame(
  x1 = pi * 0.2,
  y1 = pi * 0.3,
  x2 = pi * 0.8,
  y2 = pi * 0.9
)

p_left <- ggplot() +
  geom_vline(data = grid_v, aes(xintercept = x),
             linetype = "dashed", color = "grey60") +
  geom_hline(data = grid_h, aes(yintercept = y),
             linetype = "dashed", color = "grey60") +
  geom_line(data = lines_left,
            aes(x = theta, y = phi, group = line_id),
            linewidth = 0.6, color = "black") +
  geom_path(data = parallelogram,
            aes(x = theta, y = phi),
            linewidth = 1.0, color = "black") +
  geom_segment(data = arrow_df,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(0.2, "cm")),
               linewidth = 0.8) +
  scale_x_continuous(limits = c(0, two_pi), breaks = c(0, pi, two_pi),
                     labels = c("0", expression(pi), expression(2*pi))) +
  scale_y_continuous(limits = c(0, two_pi), breaks = c(0, pi, two_pi),
                     labels = c("0", expression(pi), expression(2*pi))) +
  coord_fixed() +
  labs(x = expression(theta), y = expression(phi)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

# -----------------------------
# Combine panels
# -----------------------------
p_left + p_right

