# Load necessary libraries
if (!requireNamespace("ExcursionSets", quietly = TRUE)) {
  install.packages("ExcursionSets")
}
library(ExcursionSets)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

if (!requireNamespace("reshape2", quietly = TRUE)) {
  install.packages("reshape2")
}
library(reshape2)

if (!requireNamespace("excursions", quietly = TRUE)) {
  install.packages("excursions")
}
library(excursions)

# ------------------------------
# Configuration Parameters
# ------------------------------
HALF_WIDTH <- 0.5  # Set to 2 or 0.5 based on desired configuration
SIDE <- TRUE       # TRUE to choose point on boundary, FALSE otherwise
SIDE_LENGTH <- 51  # Set to 121 for the results presented in the paper (takes longer)

# Generate file name components based on parameters
fname1 <- if (HALF_WIDTH == 2) "2" else "05"
fname2 <- if (SIDE) "_side" else ""
x_val <- if (SIDE) HALF_WIDTH else 0


generator_path = paste(
  "gaussian_generator_matern_-", fname1, "to", fname1, fname2, ".RData", sep = ""
)

cond_site = data.frame("x" = x_val,"y" = 0)
if(file.exists(generator_path)){
  print(paste(
    "Loading generator:", generator_path
  ))
  load(generator_path)
} else {
  x = y = seq(-HALF_WIDTH, HALF_WIDTH, length.out = SIDE_LENGTH)
  print(paste(
    "Building the random field generator:", generator_path,
    "from", min(x), "to", max(x), "with", length(x), "entries per side."
  ))
  gen = random.field.generator(expand.grid(x, y),
                               cov.model("Matern", nu = 2.5),
                               cond_sites = cond_site)
  print("Done! Saving...")
  save(gen, file = generator_path)
}

proportion.range <- function(ep, proportion = 1) {
  left <- 1
  right <- length(ep$y)
  result <- NULL
  
  while (left <= right) {
    mid <- floor((left + right) / 2)
    
    if (ep$y[mid] < proportion) {
      result <- mid
      right <- mid - 1  # Keep searching to the left for the first occurrence
    } else {
      left <- mid + 1  # Search to the right
    }
  }
  if (is.null(result)) return(Inf)
  return(ep$x[result])
}

p = 0.99
u = qnorm(p)
N = 500
var_array = c()
epys = c()
uers = c()
lers = c()
areas = c()
perims = c()
data_path = paste(
  "statistics_-", fname1, "to", fname1, "_u99", fname2, ".RData", sep = ""
)
if (file.exists(data_path)){
  print(paste(
    "Loading data:", data_path
  ))
  load(data_path)
  rf = generateGaussian(gen, u, exceedance = TRUE) # need this for the deltas...
  ep = extent.profile(rf, u, x0 = cond_site) # need this too
} else {
  for (simulation in 1:N){
    print(paste(
      simulation, "/", N, sep = ""
    ))
    rf = generateGaussian(gen, u, exceedance = TRUE)
    var_array = cbind(var_array, rf$values)
    # plotRandomField(rf, u)
    
    ep = extent.profile(rf, u, x0 = cond_site)
    epys = cbind(epys, ep$y)
    # plot(ep)
    
    uer = upperExtremalRange(rf, u, x0 = cond_site, finite = TRUE)
    uers = c(uers, uer)
    ler = proportion.range(ep, 1)
    lers = c(lers, ler)
    
    areas = c(areas, excursionArea(rf, u))
    perims = c(perims, excursionPerimeter(rf, u))
  }
  save(lers, uers, epys, var_array, perims, areas, file = data_path)
}

# ---------------

data_path = paste(
  "bootstraps_-", fname1, "to", fname1, "_u99", fname2, ".RData", sep = ""
)
if(file.exists(data_path)){
  print(paste(
    "Loading data:", data_path
  ))
  load(data_path)
} else {
  
  delta = rf$deltas[1]
  
  N_bootstraps = 200
  df = data.frame()
  for (bs in 1:N_bootstraps){
    print(paste(
      bs, "/", N_bootstraps, sep = ""
    ))
    bs_indices = sample(N, replace = TRUE)
    
    theta_1 = median(lers[bs_indices])
    
    theta_2 = median(uers[bs_indices])
    
    ep$y = apply(epys[,bs_indices], 1, mean)
    theta_3 = proportion.range(ep, 0.5)
    # plot(ep)
    
    excurs = excursions::excursions.mc(var_array[,bs_indices], u = u, alpha = 0.5, type = ">")
    confidence_set = random.field(excurs$F, gen$grid)
    theta_4 = sqrt(excursionArea(confidence_set, 0.5) / pi)
    # plot(confidence_set, 0.5)
    
    # small_bs_indices = sample(N_small, replace = TRUE)
    theta_5 = 8/pi*sum(areas[bs_indices])/sum(perims[bs_indices])
    
    new_row = c(theta_1, theta_2, theta_3, theta_4, theta_5)
    df = rbind(df, new_row)
  }
  names(df) = c("theta_1", "theta_2", "theta_3", "theta_4", "theta_5")
  bottom = log(delta)
  save(N_bootstraps, df, bottom, file = data_path)
}

# ------------------------------------------------------------------------------
# Data visualization
# ------------------------------------------------------------------------------
# First, a boxplot of the statistics in the non-discretization area, and their standard deviations.

library(ggplot2)

x0 = c(x_val,0)

# load(generator_path)
rf = generateGaussian(gen, 1, exceedance = TRUE)
# plotRandomField(rf, 1)
ep = extent.profile(rf, 1, x0 = x0)
load(data_path)

delta = exp(bottom) # 1/range on a grid of spacing 1
furthest = if(SIDE) sqrt(5) * 60 else sqrt(2) * 60
closest = 1
df = df / delta
theta_labels <- c(expression(theta[1]), expression(theta[2]), expression(theta[3]),
                  expression(theta[4]), expression(theta[5]))

df_long <- reshape2::melt(log2(df), variable.name = "Statistic", value.name = "Value")
steps = ep$x / delta
steps = steps[-1]
break_lines = seq(0, 7, by = 1)
p = ggplot(df_long, aes(x = Statistic, y = Value)) +
  geom_hline(yintercept = log2(steps), color = "lightgrey", linetype = "solid", alpha = 0.5) + # Horizontal lines
  geom_hline(yintercept = break_lines, color = "darkgrey", linetype = "dashed", alpha = 1) + # Horizontal lines
  geom_boxplot(fill = "skyblue", color = "darkblue") +                             # Boxplot
  theme_minimal() +
  scale_y_continuous(limits = c(1.5, 7), breaks = break_lines,
                     labels = 2^break_lines) +
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels = theta_labels) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)

print("Standard deviations")
print(round(apply(log(as.matrix(df)), 2, sd), 3))
print("medians")
print(round(apply(as.matrix(df), 2, median), 1))
figure_name = paste(
  fname1, "by", fname1, fname2, ".pdf", sep = ""
)
# ggsave(figure_name, height = 4.5, width = 4.5)
# print(paste(
#   figure_name, "saved!"
# ))



