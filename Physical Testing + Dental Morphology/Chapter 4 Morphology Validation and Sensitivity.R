# =============================================================================
# Validation and sensitivity analyses for the dental morphology workflow
# =============================================================================
# Purpose
# -------
# Documents the validation, repeatability and sensitivity analyses
# used to develop the final workflow in dental_morphology_workflow.R.
#
# This script is NOT required to analyse a new specimen. It is provided to make
# the methodological development in the thesis reproducible. The example paths
# below should point to the validation datasets used for those analyses (for
# example D1/D1 if these files can be shared under the relevant permissions).
#
# Repository layout assumed:
#   R/dental_morphology_workflow.R
#   R/validation_and_sensitivity_analysis.R
#   Data/validation/...
#   Results/validation/
#
# NOTE: several Blender-validation steps require cross-section PLY files exported
# manually from Blender. Update those filenames/paths below to match the files
# included with the repository.

# =============================================================================
# USER INPUTS: thesis validation datasets
# =============================================================================
# Primary landmark/mesh dataset used for validation and sensitivity analyses
landmark_file <- "Data/validation/D1/D1.csv"
repeat_landmark_file <- "Data/validation/D1/D1_repeat.csv"
upper_mesh_file <- "Data/validation/D1/D1_Upper.stl"
lower_mesh_file <- "Data/validation/D1/D1_Lower.stl"

# Blender-exported validation cross-sections
validation_dir <- "Data/validation/Blender"
roc_validation_md_file <- file.path(validation_dir, "MD_slice_M1.ply")
roc_validation_bl_file <- file.path(validation_dir, "BL_slice_M1.ply")
csa_validation_file <- file.path(validation_dir, "Lower_I2_L_CSA.ply")

output_dir <- "Results/validation"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Check files required for the base workflow and repeated landmark analyses.
required_validation_files <- c(
  landmark_file,
  repeat_landmark_file,
  upper_mesh_file,
  lower_mesh_file,
  roc_validation_md_file,
  roc_validation_bl_file,
  csa_validation_file
)
missing_validation_files <- required_validation_files[!file.exists(required_validation_files)]
if (length(missing_validation_files) > 0) {
  stop(
    "The following validation files were not found:\n",
    paste0("  - ", missing_validation_files, collapse = "\n"),
    "\n\nEdit the USER INPUTS section before running this script."
  )
}

# Load the reusable workflow. Because the file paths above already exist in the
# environment, the workflow uses these validation files rather than its defaults.
workflow_script <- file.path("R", "dental_morphology_workflow.R")
if (!file.exists(workflow_script)) {
  # Allows the two scripts to be run from the same directory during testing.
  workflow_script <- "dental_morphology_workflow.R"
}
source(workflow_script)

# The validation sections below use additional packages.
validation_packages <- c("irr")
missing_packages <- validation_packages[
  !vapply(validation_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the following validation package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}
suppressPackageStartupMessages(library(irr))

 


# ROC and angle Blender Validation
# Validate MD and BL Blender sections together


# ---------------------------
# Export validation table for Blender plane setup
# ---------------------------

roc_blender_validation <- apex_data %>%
  rowwise() %>%
  mutate(
    MD_unit_x_3D = MD_unit_x,
    MD_unit_y_3D = MD_unit_y,
    MD_unit_z_3D = 0,
    
    BL_unit_x_3D = -MD_unit_y,
    BL_unit_y_3D =  MD_unit_x,
    BL_unit_z_3D = 0,
    
    normal_MD_x = MD_unit_y,
    normal_MD_y = -MD_unit_x,
    normal_MD_z = 0,
    
    normal_BL_x = MD_unit_x,
    normal_BL_y = MD_unit_y,
    normal_BL_z = 0,
    
    ROI_radius = 0.15 * mean(c(MD_length, BL_length), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(
    Specimen, Model, replicate, side, tooth, cusp,
    apex_x, apex_y, apex_z,
    MD_unit_x_3D, MD_unit_y_3D, MD_unit_z_3D,
    BL_unit_x_3D, BL_unit_y_3D, BL_unit_z_3D,
    normal_MD_x, normal_MD_y, normal_MD_z,
    normal_BL_x, normal_BL_y, normal_BL_z,
    ROI_radius,
    MD_length, BL_length
  )
# ---------------------------
# Load Blender-exported MD and BL section meshes
# ---------------------------

blender_mesh_list <- list(
  MD = vcgImport(roc_validation_md_file),
  BL = vcgImport(roc_validation_bl_file)
)

blender_points_list <- lapply(blender_mesh_list, function(mesh) {
  pts <- t(mesh$vb[1:3, ])
  colnames(pts) <- c("x", "y", "z")
  pts
})



# ---------------------------
# Helper: calculate subtended angle using outer points
# ---------------------------

calc_angle_from_profile_outer_points <- function(profile,
                                                 angle_cutoff = 0.60,
                                                 n_outer_points = 5) {
  
  if (nrow(profile) < (2 * n_outer_points)) {
    return(list(
      angle = NA_real_,
      line_df = data.frame(),
      left_fit_points = data.frame(),
      right_fit_points = data.frame()
    ))
  }
  
  x_limit <- angle_cutoff * max(abs(profile$position), na.rm = TRUE)
  
  angle_profile <- profile %>%
    filter(abs(position) <= x_limit)
  
  left_side <- angle_profile %>%
    filter(position < 0) %>%
    arrange(position)
  
  right_side <- angle_profile %>%
    filter(position > 0) %>%
    arrange(position)
  
  if (nrow(left_side) < n_outer_points || nrow(right_side) < n_outer_points) {
    return(list(
      angle = NA_real_,
      line_df = data.frame(),
      left_fit_points = data.frame(),
      right_fit_points = data.frame()
    ))
  }
  
  left_fit_points <- left_side %>%
    slice_head(n = n_outer_points)
  
  right_fit_points <- right_side %>%
    slice_tail(n = n_outer_points)
  
  left_fit <- lm(height ~ position, data = left_fit_points)
  right_fit <- lm(height ~ position, data = right_fit_points)
  
  m_left <- coef(left_fit)[["position"]]
  c_left <- coef(left_fit)[["(Intercept)"]]
  
  m_right <- coef(right_fit)[["position"]]
  c_right <- coef(right_fit)[["(Intercept)"]]
  
  if (!is.finite(m_left) || !is.finite(m_right) ||
      !is.finite(c_left) || !is.finite(c_right) ||
      m_left == m_right) {
    return(list(
      angle = NA_real_,
      line_df = data.frame(),
      left_fit_points = left_fit_points,
      right_fit_points = right_fit_points
    ))
  }
  
  x_intersect <- (c_right - c_left) / (m_left - m_right)
  y_intersect <- m_left * x_intersect + c_left
  
  x_left_end <- min(left_fit_points$position, na.rm = TRUE)
  y_left_end <- m_left * x_left_end + c_left
  
  x_right_end <- max(right_fit_points$position, na.rm = TRUE)
  y_right_end <- m_right * x_right_end + c_right
  
  v_left <- c(x_left_end - x_intersect, y_left_end - y_intersect)
  v_right <- c(x_right_end - x_intersect, y_right_end - y_intersect)
  
  cos_theta <- sum(v_left * v_right) /
    (sqrt(sum(v_left^2)) * sqrt(sum(v_right^2)))
  
  cos_theta <- max(min(cos_theta, 1), -1)
  
  theta_deg <- acos(cos_theta) * 180 / pi
  
  line_df <- data.frame(
    x = c(x_intersect, x_left_end, x_intersect, x_right_end),
    y = c(y_intersect, y_left_end, y_intersect, y_right_end),
    side_line = c("left", "left", "right", "right")
  )
  
  list(
    angle = theta_deg,
    line_df = line_df,
    left_fit_points = left_fit_points,
    right_fit_points = right_fit_points
  )
}







# ---------------------------
# Repeatability comparison- Tooth dimensions:
# ---------------------------

# ---------------------------
# Load both landmark datasets
# ---------------------------

data_D1 <- read.csv(landmark_file, header = TRUE)
data_D1_2 <- read.csv(repeat_landmark_file, header = TRUE)

# ---------------------------
# Calculate measurements for both datasets
# ---------------------------

results_D1 <- calc_vectors(data_D1) %>%
  mutate(dataset = "D1")

results_D1_2 <- calc_vectors(data_D1_2) %>%
  mutate(dataset = "D1.2")

# ---------------------------
# Measurements to compare
# ---------------------------

measurements <- c(
  "MD_length",
  "BL_length",
  "B_cusp_span",
  "L_cusp_span",
  "mesial_BL_cusp_distance",
  "distal_BL_cusp_distance",
  "diagonal_B1_Ldistal",
  "diagonal_Bdistal_L1"
)

# ---------------------------
# Join matching teeth between datasets
# ---------------------------

repeatability_wide <- bind_rows(results_D1, results_D1_2) %>%
  select(
    dataset,
    Specimen, Model, replicate, side, tooth,
    all_of(measurements)
  ) %>%
  pivot_longer(
    cols = all_of(measurements),
    names_to = "measurement",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = dataset,
    values_from = value
  )
# ---------------------------
# Calculate absolute and percentage differences
# ---------------------------

repeatability_differences <- repeatability_wide %>%
  mutate(
    absolute_difference = abs(D1 - `D1.2`),
    percent_difference = 100 * absolute_difference / ((D1 + `D1.2`) / 2)
  )
# ---------------------------
# ICC function
# ---------------------------

calc_icc <- function(df) {
  
  df_icc <- df %>%
    select(D1, `D1.2`) %>%
    filter(is.finite(D1), is.finite(`D1.2`))
  
  if (nrow(df_icc) < 3) {
    return(NA_real_)
  }
  
  icc_result <- irr::icc(
    df_icc,
    model = "twoway",
    type = "agreement",
    unit = "single"
  )
  
  icc_result$value
}

# ---------------------------
# Summary table by measurement
# ---------------------------

repeatability_summary <- repeatability_differences %>%
  group_by(measurement) %>%
  summarise(
    n_pairs = sum(is.finite(D1) & is.finite(`D1.2`)),
    mean_D1 = mean(D1, na.rm = TRUE),
    mean_D1_2 = mean(`D1.2`, na.rm = TRUE),
    mean_absolute_difference = mean(absolute_difference, na.rm = TRUE),
    mean_percent_difference = mean(percent_difference, na.rm = TRUE),
    max_percent_difference = max(percent_difference, na.rm = TRUE),
    ICC = calc_icc(cur_data()),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
print(repeatability_summary)



## ROC Validation



validate_blender_roc_angle <- function(validation_axis,
                                       distance_cutoff = 0.05,
                                       roc_cutoff = 0.60,
                                       angle_cutoff = 0.60,
                                       n_outer_points = 5) {
  
  blender_points <- blender_points_list[[validation_axis]]
  
  apex <- c(one_cusp$apex_x, one_cusp$apex_y, one_cusp$apex_z)
  
  MD_unit <- c(one_cusp$MD_unit_x, one_cusp$MD_unit_y, 0)
  MD_unit <- MD_unit / sqrt(sum(MD_unit^2))
  
  BL_unit <- c(-one_cusp$MD_unit_y, one_cusp$MD_unit_x, 0)
  BL_unit <- BL_unit / sqrt(sum(BL_unit^2))
  
  axis_unit <- if (validation_axis == "MD") MD_unit else BL_unit
  R_plot <- if (validation_axis == "MD") MD_plot else BL_plot
  
  V_blender <- sweep(blender_points, 2, apex)
  
  blender_profile <- data.frame(
    position = as.numeric(V_blender %*% axis_unit),
    height = V_blender[, 3]
  )
  
  if (one_cusp$Model == "Lower") {
    blender_profile$height <- -blender_profile$height
  }
  
  x_limit_roc <- roc_cutoff * max(abs(R_plot$profile$position), na.rm = TRUE)
  
  R_roc_profile <- R_plot$profile %>%
    filter(abs(position) <= x_limit_roc)
  
  blender_profile_filtered <- blender_profile %>%
    mutate(
      expected_height = predict(
        R_plot$fit,
        newdata = data.frame(position = position)
      ),
      distance_to_R_curve = abs(height - expected_height)
    ) %>%
    filter(
      position >= min(R_plot$profile$position, na.rm = TRUE),
      position <= max(R_plot$profile$position, na.rm = TRUE),
      abs(position) <= x_limit_roc,
      distance_to_R_curve <= distance_cutoff
    )
  
  blender_fit <- lm(
    height ~ poly(position, 3, raw = TRUE),
    data = blender_profile_filtered
  )
  
  blender_RoC <- calc_RoC_from_fit(blender_fit)
  blender_r2 <- summary(blender_fit)$r.squared
  
  r_RoC <- roc_results %>%
    filter(
      Specimen == one_cusp$Specimen,
      Model == one_cusp$Model,
      tooth == one_cusp$tooth,
      side == one_cusp$side,
      cusp == one_cusp$cusp
    ) %>%
    pull(if (validation_axis == "MD") MD_RoC else BL_RoC)
  
  R_angle_result <- calc_angle_from_profile_outer_points(
    R_plot$profile,
    angle_cutoff = angle_cutoff,
    n_outer_points = n_outer_points
  )
  
  Blender_angle_result <- calc_angle_from_profile_outer_points(
    blender_profile_filtered,
    angle_cutoff = angle_cutoff,
    n_outer_points = n_outer_points
  )
  
  R_angle <- R_angle_result$angle
  Blender_angle <- Blender_angle_result$angle
  
  R_angle_lines <- R_angle_result$line_df
  Blender_angle_lines <- Blender_angle_result$line_df
  
  if (nrow(R_angle_lines) > 0) R_angle_lines$Source <- "R"
  if (nrow(Blender_angle_lines) > 0) Blender_angle_lines$Source <- "Blender"
  
  angle_lines <- bind_rows(R_angle_lines, Blender_angle_lines)
  
  x_seq <- seq(
    min(R_roc_profile$position, na.rm = TRUE),
    max(R_roc_profile$position, na.rm = TRUE),
    length.out = 100
  )
  
  fit_curve <- data.frame(
    position = x_seq,
    height = predict(R_plot$fit, newdata = data.frame(position = x_seq))
  )
  
  ROI_radius <- 0.15 * mean(c(one_cusp$MD_length, one_cusp$BL_length), na.rm = TRUE)
  
  plot_subtitle <- paste0(
    "R: RoC = ", round(r_RoC, 3),
    "; angle = ", round(R_angle, 1), "°",
    "; R² = ", round(summary(R_plot$fit)$r.squared, 3),
    "; n = ", nrow(R_roc_profile),
    "\nBlender: RoC = ", round(blender_RoC, 3),
    "; angle = ", round(Blender_angle, 1), "°",
    "; R² = ", round(blender_r2, 3),
    "; n = ", nrow(blender_profile_filtered),
    "; ROI = ", round(ROI_radius, 3)
  )
  
  overlay_plot <- ggplot() +
    geom_point(
      data = R_plot$profile,
      aes(x = position, y = height),
      colour = "grey60",
      alpha = 0.3
    ) +
    geom_point(
      data = R_roc_profile,
      aes(x = position, y = height),
      colour = "black"
    ) +
    geom_line(
      data = fit_curve,
      aes(x = position, y = height),
      colour = "black"
    ) +
    geom_point(
      data = blender_profile_filtered,
      aes(x = position, y = height),
      colour = "red",
      size = 2
    ) +
    geom_point(
      data = data.frame(position = 0, height = 0),
      aes(x = position, y = height),
      shape = 4,
      size = 3,
      stroke = 1
    ) +
    geom_line(
      data = angle_lines,
      aes(x = x, y = y, group = interaction(Source, side_line), linetype = Source),
      linewidth = 1
    ) +
    theme_classic() +
    labs(
      title = paste("R and Blender RoC/angle validation:", validation_axis, "section"),
      subtitle = plot_subtitle,
      x = paste("Position along", validation_axis, "axis"),
      y = "Height relative to apex",
      linetype = "Angle lines"
    )
  
  comparison <- data.frame(
    Specimen = one_cusp$Specimen,
    Model = one_cusp$Model,
    tooth = one_cusp$tooth,
    side = one_cusp$side,
    cusp = one_cusp$cusp,
    Axis = validation_axis,
    R_RoC = r_RoC,
    Blender_RoC = blender_RoC,
    RoC_Difference = blender_RoC - r_RoC,
    RoC_Percent_Difference = 100 * (blender_RoC - r_RoC) / r_RoC,
    R_angle = R_angle,
    Blender_angle = Blender_angle,
    Angle_Difference = Blender_angle - R_angle,
    Angle_Percent_Difference = 100 * (Blender_angle - R_angle) / R_angle,
    R_R2 = summary(R_plot$fit)$r.squared,
    Blender_R2 = blender_r2,
    R_n_points = nrow(R_roc_profile),
    Blender_n_points = nrow(blender_profile_filtered)
  )
  
  list(
    comparison = comparison,
    plot = overlay_plot,
    R_angle_result = R_angle_result,
    Blender_angle_result = Blender_angle_result
  )
}

# ---------------------------
# Run for MD and BL
# ---------------------------

MD_validation <- validate_blender_roc_angle("MD")
BL_validation <- validate_blender_roc_angle("BL")

MD_validation$plot
BL_validation$plot

angle_roc_comparison_both <- bind_rows(
  MD_validation$comparison,
  BL_validation$comparison
) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
print(angle_roc_comparison_both)





# ---------------------------
# Incisor CSA Blender validation
# Example: Upper i1 L
# ---------------------------

# ---------------------------
# Choose incisor to validate
# ---------------------------

validation_row <- df_results %>%
  filter(Model == "Lower", tooth == "I2", side == "L") %>%
  slice(1)

# ---------------------------
# Load Blender-exported cross-section
# Change filename if needed
# ---------------------------

blender_csa_mesh <- vcgImport(csa_validation_file)

blender_csa_points <- t(blender_csa_mesh$vb[1:3, ])
colnames(blender_csa_points) <- c("x", "y", "z")

# ---------------------------
# Re-extract the R CSA section for the same tooth
# ---------------------------

r_extracted <- extract_incisor_intersection_2d(
  apex = c(validation_row$x_apex, validation_row$y_apex, validation_row$z_apex),
  gumline = c(validation_row$x_gumline, validation_row$y_gumline, validation_row$z_gumline),
  mesial = c(validation_row$x_mesial, validation_row$y_mesial, validation_row$z_mesial),
  distal = c(validation_row$x_distal, validation_row$y_distal, validation_row$z_distal),
  buccal = c(validation_row$x_buccal, validation_row$y_buccal, validation_row$z_buccal),
  lingual = c(validation_row$x_lingual, validation_row$y_lingual, validation_row$z_lingual),
  MD_unit = c(validation_row$MD_unit_x, validation_row$MD_unit_y),
  BL_unit = c(validation_row$BL_x, validation_row$BL_y),
  MD_length = validation_row$MD_length,
  BL_length = validation_row$BL_length,
  mesh_points = mesh_points_list[[as.character(validation_row$Model)]],
  mesh_faces = mesh_faces_list[[as.character(validation_row$Model)]],
  Model = validation_row$Model,
  section_scale = 0.80,
  bbox_buffer = 0.10
)

r_section_2d <- r_extracted$section_2d

# ---------------------------
# Project Blender section into same 2D MD/BL space
# ---------------------------

MD_unit <- c(validation_row$MD_unit_x, validation_row$MD_unit_y)
MD_unit <- MD_unit / sqrt(sum(MD_unit^2))

BL_unit <- c(validation_row$BL_x, validation_row$BL_y)
BL_unit <- BL_unit / sqrt(sum(BL_unit^2))

tooth_centre <- colMeans(rbind(
  c(validation_row$x_mesial, validation_row$y_mesial),
  c(validation_row$x_distal, validation_row$y_distal),
  c(validation_row$x_buccal, validation_row$y_buccal),
  c(validation_row$x_lingual, validation_row$y_lingual)
))

blender_xy <- blender_csa_points[, c("x", "y"), drop = FALSE]

blender_local <- sweep(blender_xy, 2, tooth_centre)

blender_section_2d <- data.frame(
  MD = as.numeric(as.matrix(blender_local) %*% MD_unit),
  BL = as.numeric(as.matrix(blender_local) %*% BL_unit)
)

# Optional: order Blender points around centre for polygon plotting
cx_b <- mean(blender_section_2d$MD, na.rm = TRUE)
cy_b <- mean(blender_section_2d$BL, na.rm = TRUE)

blender_section_2d <- blender_section_2d %>%
  mutate(angle = atan2(BL - cy_b, MD - cx_b)) %>%
  arrange(angle)

# ---------------------------
# Calculate areas
# ---------------------------

R_CSA <- if (nrow(r_section_2d) >= 3) {
  polygon_area(r_section_2d$MD, r_section_2d$BL)
} else {
  NA_real_
}

Blender_CSA <- if (nrow(blender_section_2d) >= 3) {
  polygon_area(blender_section_2d$MD, blender_section_2d$BL)
} else {
  NA_real_
}

CSA_difference <- Blender_CSA - R_CSA
CSA_percent_difference <- 100 * CSA_difference / R_CSA

csa_validation_summary <- data.frame(
  Specimen = validation_row$Specimen,
  Model = validation_row$Model,
  tooth = validation_row$tooth,
  side = validation_row$side,
  R_CSA = R_CSA,
  Blender_CSA = Blender_CSA,
  Difference = CSA_difference,
  Percent_Difference = CSA_percent_difference,
  R_n_points = nrow(r_section_2d),
  Blender_n_points = nrow(blender_section_2d),
  R_section_z = r_extracted$section_z,
  R_reason = r_extracted$reason
) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
print(csa_validation_summary)

# ---------------------------
# Overlay plot
# ---------------------------

plot_subtitle <- paste0(
  "R CSA = ", round(R_CSA, 3),
  "; n = ", nrow(r_section_2d),
  "; section z = ", round(r_extracted$section_z, 3),
  "\nBlender CSA = ", round(Blender_CSA, 3),
  "; n = ", nrow(blender_section_2d),
  "; difference = ", round(CSA_percent_difference, 2), "%"
)

ggplot() +
  geom_polygon(
    data = r_section_2d,
    aes(x = MD, y = BL),
    fill = NA,
    colour = "black",
    linewidth = 1
  ) +
  geom_point(
    data = r_section_2d,
    aes(x = MD, y = BL),
    colour = "black",
    size = 1.8
  ) +
  geom_polygon(
    data = blender_section_2d,
    aes(x = MD, y = BL),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  geom_point(
    data = blender_section_2d,
    aes(x = MD, y = BL),
    colour = "red",
    size = 1.8,
    alpha = 0.7
  ) +
  coord_equal() +
  theme_classic() +
  labs(
    title = paste(validation_row$Model, validation_row$tooth, validation_row$side,
                  "CSA validation: R vs Blender"),
    subtitle = plot_subtitle,
    x = "Mesiodistal axis",
    y = "Buccolingual axis"
  )



# ---------------------------
# Sensitivity analysis ROc and Angle


# ---------------------------
# Cutoff sensitivity
# Keep cap_scale constant
# ---------------------------

cutoff_values <- c(0.40, 0.50, 0.60, 0.70)

cutoff_sensitivity_results <- purrr::map_dfr(cutoff_values, function(cutoff) {
  
  temp_results <- apex_data %>%
    rowwise() %>%
    mutate(
      roc = list(
        calc_cusp_roc(
          apex = c(apex_x, apex_y, apex_z),
          MD_unit_x = MD_unit_x,
          MD_unit_y = MD_unit_y,
          MD_length = MD_length,
          BL_length = BL_length,
          mesh_points = mesh_points_list[[as.character(Model)]],
          mesh_faces = mesh_faces_list[[as.character(Model)]],
          Model = Model,
          cap_scale = 0.15,
          roc_cutoff = cutoff,
          angle_cutoff = cutoff,
          n_outer_points = 5
        )
      )
    ) %>%
    unnest(roc) %>%
    ungroup()
  
  tibble(
    test_type = "cutoff",
    value = cutoff,
    MD_RoC_missing = sum(is.na(temp_results$MD_RoC)),
    BL_RoC_missing = sum(is.na(temp_results$BL_RoC)),
    MD_angle_missing = sum(is.na(temp_results$MD_angle)),
    BL_angle_missing = sum(is.na(temp_results$BL_angle)),
    mean_MD_RoC = mean(temp_results$MD_RoC, na.rm = TRUE),
    mean_BL_RoC = mean(temp_results$BL_RoC, na.rm = TRUE),
    mean_MD_angle = mean(temp_results$MD_angle, na.rm = TRUE),
    mean_BL_angle = mean(temp_results$BL_angle, na.rm = TRUE)
  )
})
# ---------------------------
# Cap scale sensitivity
# Keep cutoffs constant
# ---------------------------

cap_scale_values <- c(0.10, 0.15, 0.20, 0.25)

cap_scale_sensitivity_results <- purrr::map_dfr(cap_scale_values, function(cap) {
  
  temp_results <- apex_data %>%
    rowwise() %>%
    mutate(
      roc = list(
        calc_cusp_roc(
          apex = c(apex_x, apex_y, apex_z),
          MD_unit_x = MD_unit_x,
          MD_unit_y = MD_unit_y,
          MD_length = MD_length,
          BL_length = BL_length,
          mesh_points = mesh_points_list[[as.character(Model)]],
          mesh_faces = mesh_faces_list[[as.character(Model)]],
          Model = Model,
          cap_scale = cap,
          roc_cutoff = 0.60,
          angle_cutoff = 0.60,
          n_outer_points = 5
        )
      )
    ) %>%
    unnest(roc) %>%
    ungroup()
  
  tibble(
    test_type = "cap_scale",
    value = cap,
    MD_RoC_missing = sum(is.na(temp_results$MD_RoC)),
    BL_RoC_missing = sum(is.na(temp_results$BL_RoC)),
    MD_angle_missing = sum(is.na(temp_results$MD_angle)),
    BL_angle_missing = sum(is.na(temp_results$BL_angle)),
    mean_MD_RoC = mean(temp_results$MD_RoC, na.rm = TRUE),
    mean_BL_RoC = mean(temp_results$BL_RoC, na.rm = TRUE),
    mean_MD_angle = mean(temp_results$MD_angle, na.rm = TRUE),
    mean_BL_angle = mean(temp_results$BL_angle, na.rm = TRUE)
  )
})
# ---------------------------
# Number of points used on cusp sides sensitivity
# Keep cut-offs and cap scale constant
# ---------------------------

n_outer_values <- c(3, 4, 5, 6)

outer_point_sensitivity_results <- purrr::map_dfr(n_outer_values, function(n_outer) {
  
  temp_results <- apex_data %>%
    rowwise() %>%
    mutate(
      roc = list(
        calc_cusp_roc(
          apex = c(apex_x, apex_y, apex_z),
          MD_unit_x = MD_unit_x,
          MD_unit_y = MD_unit_y,
          MD_length = MD_length,
          BL_length = BL_length,
          mesh_points = mesh_points_list[[as.character(Model)]],
          mesh_faces = mesh_faces_list[[as.character(Model)]],
          Model = Model,
          cap_scale = 0.15,
          roc_cutoff = 0.60,
          angle_cutoff = 0.60,
          n_outer_points = n_outer
        )
      )
    ) %>%
    unnest(roc) %>%
    ungroup()
  
  tibble(
    test_type = "n_outer_points",
    value = n_outer,
    MD_angle_missing = sum(is.na(temp_results$MD_angle)),
    BL_angle_missing = sum(is.na(temp_results$BL_angle)),
    mean_MD_angle = mean(temp_results$MD_angle, na.rm = TRUE),
    mean_BL_angle = mean(temp_results$BL_angle, na.rm = TRUE)
  )
})
# ---------------------------
# Combined sensitivity table
# ---------------------------

sensitivity_results <- bind_rows(
  cutoff_sensitivity_results,
  cap_scale_sensitivity_results,
  outer_point_sensitivity_results
)
# ---------------------------
# Landmark dataset sensitivity ROc and angle
# ---------------------------

# Load second landmark dataset
data_repeat <- read.csv(repeat_landmark_file, header = TRUE)

# ---------------------------
# Create apex dataset for repeat landmarks
# ---------------------------

apex_data_repeat <- data_repeat %>%
  mutate(landmark_clean = tolower(trimws(landmark))) %>%
  filter(landmark_clean %in% cusp_landmarks) %>%
  rename(
    cusp = landmark_clean,
    apex_x = x,
    apex_y = y,
    apex_z = z
  ) %>%
  left_join(
    df_results %>%
      select(
        Specimen, Model, replicate, side, tooth,
        MD_unit_x, MD_unit_y,
        MD_length, BL_length
      ),
    by = c("Specimen", "Model", "replicate", "side", "tooth")
  )

# ---------------------------
# Calculate ROC and angles for repeat landmarks
# ---------------------------

roc_results_repeat <- apex_data_repeat %>%
  rowwise() %>%
  mutate(
    roc = list(
      calc_cusp_roc(
        apex = c(apex_x, apex_y, apex_z),
        MD_unit_x = MD_unit_x,
        MD_unit_y = MD_unit_y,
        MD_length = MD_length,
        BL_length = BL_length,
        mesh_points = mesh_points_list[[as.character(Model)]],
        mesh_faces = mesh_faces_list[[as.character(Model)]],
        Model = Model,
        cap_scale = 0.15,
        roc_cutoff = 0.60,
        angle_cutoff = 0.60,
        n_outer_points = 5
      )
    )
  ) %>%
  unnest(roc) %>%
  ungroup()

# ---------------------------
# Compare datasets
# ---------------------------

landmark_sensitivity <- bind_rows(
  
  roc_results %>%
    mutate(dataset = "D1"),
  
  roc_results_repeat %>%
    mutate(dataset = "D1.2")
  
) %>%
  group_by(dataset) %>%
  summarise(
    
    MD_RoC_missing = sum(is.na(MD_RoC)),
    BL_RoC_missing = sum(is.na(BL_RoC)),
    
    MD_angle_missing = sum(is.na(MD_angle)),
    BL_angle_missing = sum(is.na(BL_angle)),
    
    mean_MD_RoC = mean(MD_RoC, na.rm = TRUE),
    mean_BL_RoC = mean(BL_RoC, na.rm = TRUE),
    
    mean_MD_angle = mean(MD_angle, na.rm = TRUE),
    mean_BL_angle = mean(BL_angle, na.rm = TRUE)
    
  )
print(landmark_sensitivity)

# ---------------------------
# CSA repeatability: D1 vs D1.2
# ---------------------------

# ---------------------------
# Load both landmark datasets
# ---------------------------

data_D1 <- read.csv(landmark_file, header = TRUE)
data_D1_2 <- read.csv(repeat_landmark_file, header = TRUE)

# ---------------------------
# Calculate vectors for both datasets
# ---------------------------

df_results_D1 <- calc_vectors(data_D1) %>%
  mutate(dataset = "D1")

df_results_D1_2 <- calc_vectors(data_D1_2) %>%
  mutate(dataset = "D1.2")

# ---------------------------
# Helper function to calculate CSA for a dataset
# ---------------------------

calc_csa_dataset <- function(df_results_input) {
  
  df_results_input %>%
    filter(tooth %in% incisor_teeth) %>%
    rowwise() %>%
    mutate(
      csa = list(
        calc_incisor_csa_intersection(
          apex = c(x_apex, y_apex, z_apex),
          gumline = c(x_gumline, y_gumline, z_gumline),
          mesial = c(x_mesial, y_mesial, z_mesial),
          distal = c(x_distal, y_distal, z_distal),
          buccal = c(x_buccal, y_buccal, z_buccal),
          lingual = c(x_lingual, y_lingual, z_lingual),
          MD_unit = c(MD_unit_x, MD_unit_y),
          BL_unit = c(BL_x, BL_y),
          MD_length = MD_length,
          BL_length = BL_length,
          mesh_points = mesh_points_list[[as.character(Model)]],
          mesh_faces = mesh_faces_list[[as.character(Model)]],
          Model = Model,
          section_scale = 0.80,
          bbox_buffer = 0.10
        )
      )
    ) %>%
    unnest(csa) %>%
    ungroup()
}

# ---------------------------
# Calculate CSA for both landmark datasets
# ---------------------------

csa_D1 <- calc_csa_dataset(df_results_D1)
csa_D1_2 <- calc_csa_dataset(df_results_D1_2)

# ---------------------------
# Join matching teeth between datasets
# ---------------------------

csa_repeatability_wide <- bind_rows(csa_D1, csa_D1_2) %>%
  select(
    dataset,
    Specimen, Model, replicate, side, tooth,
    CSA, n_points, section_z, reason
  ) %>%
  pivot_wider(
    names_from = dataset,
    values_from = c(CSA, n_points, section_z, reason)
  )
# ---------------------------
# Calculate CSA differences
# ---------------------------

csa_repeatability_differences <- csa_repeatability_wide %>%
  mutate(
    absolute_difference = abs(CSA_D1 - `CSA_D1.2`),
    percent_difference = 100 * absolute_difference / ((CSA_D1 + `CSA_D1.2`) / 2)
  )
# ---------------------------
# ICC for CSA
# ---------------------------

csa_icc_data <- csa_repeatability_differences %>%
  select(CSA_D1, `CSA_D1.2`) %>%
  filter(is.finite(CSA_D1), is.finite(`CSA_D1.2`))

csa_icc <- if (nrow(csa_icc_data) >= 3) {
  irr::icc(
    csa_icc_data,
    model = "twoway",
    type = "agreement",
    unit = "single"
  )$value
} else {
  NA_real_
}

# ---------------------------
# CSA repeatability summary
# ---------------------------

csa_repeatability_summary <- csa_repeatability_differences %>%
  summarise(
    n_pairs = sum(is.finite(CSA_D1) & is.finite(`CSA_D1.2`)),
    mean_D1_CSA = mean(CSA_D1, na.rm = TRUE),
    mean_D1_2_CSA = mean(`CSA_D1.2`, na.rm = TRUE),
    mean_absolute_difference = mean(absolute_difference, na.rm = TRUE),
    mean_percent_difference = mean(percent_difference, na.rm = TRUE),
    max_percent_difference = max(percent_difference, na.rm = TRUE),
    ICC = csa_icc
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
print(csa_repeatability_summary)

# ============================================================
# Reproducibility information
# ============================================================
# Uncomment to record the R and package versions used for a specific analysis run.
# writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))


# =============================================================================
# Save validation and sensitivity outputs
# =============================================================================
write.csv(
  repeatability_summary,
  file.path(output_dir, "tooth_dimension_repeatability_summary.csv"),
  row.names = FALSE
)
write.csv(
  angle_roc_comparison_both,
  file.path(output_dir, "RoC_angle_Blender_validation.csv"),
  row.names = FALSE
)
write.csv(
  csa_validation_summary,
  file.path(output_dir, "incisor_CSA_Blender_validation.csv"),
  row.names = FALSE
)
write.csv(
  sensitivity_results,
  file.path(output_dir, "RoC_angle_parameter_sensitivity.csv"),
  row.names = FALSE
)
write.csv(
  landmark_sensitivity,
  file.path(output_dir, "RoC_angle_landmark_sensitivity.csv"),
  row.names = FALSE
)
write.csv(
  csa_repeatability_summary,
  file.path(output_dir, "incisor_CSA_repeatability_summary.csv"),
  row.names = FALSE
)

# =============================================================================
# Reproducibility record
# =============================================================================
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo_validation.txt")
)

message("Validation and sensitivity analyses complete. Outputs directory: ", output_dir)
