### Out-lever method sensitivity analysis
#
# This script calculates out-lever lengths using the following approaches:
# 1. Direct (DI)
# 2. Alveolar plane (AV)
# 3. Frankfort Horizontal plane (FH)
# 4. Occlusal plane (OC)
# 5. The frankfurt/ alveolar average
#
# Landmark numbering matches the out-lever landmark table reported
# in the Methods chapter.
#
# IMPORTANT:
# The final out-lever distances are calculated in a 2D anatomical plane
# defined by the anteroposterior (AP) and superoinferior (SI) axes.
# Set AP_axis and SI_axis below to match the orientation of your specimens.
#
# For the specimens used in this thesis:
#   AP_axis = "y"
#   SI_axis = "x"
#   ML_axis = "z"
#
# If your specimens are oriented differently, change ONLY these assignments.

## Axis orientation ---------------------------------------------------------

AP_axis <- "y"   # anteroposterior axis
SI_axis <- "x"   # superoinferior axis
ML_axis <- "z"   # mediolateral axis; excluded from 2D out-lever calculation

valid_axes <- c("x", "y", "z")

if (!all(c(AP_axis, SI_axis, ML_axis) %in% valid_axes)) {
  stop("AP_axis, SI_axis and ML_axis must each be one of: 'x', 'y', 'z'.")
}

if (length(unique(c(AP_axis, SI_axis, ML_axis))) != 3) {
  stop("AP_axis, SI_axis and ML_axis must refer to three different coordinate axes.")
}


## Landmarks ---------------------------------------------------------------

# 1, 2  Right and left external auditory meatus (EAM)
#       Superior-most point of the external auditory meatus.
#       Used for: FH
#
# 3     Left orbit
#       Inferior-most point of the lower border of the orbit.
#       Used for: FH
#
# 4, 5  Centre of right and left condyles
#       Superior-central point of the condylar head representing the jaw joint.
#       Used for: All approaches
#
# 6, 7  Right and left alveolar margin of posterior molars
#       Posterior-most point of the alveolar margin associated with the
#       posterior-most tooth or developing tooth crypt.
#       Used for: AV
#
# 8     Alveolar incisor
#       Alveolar margin at the midline between the mandibular central incisors.
#       Used for: AV
#
# 9, 10 Centre of right and left molars
#       Superior-central point of the posterior-most fully erupted tooth.
#       Used for: OC
#
# 11    Anterior bite point
#       Centre of the medial incisors or the alveolar incisor margin.
#       Used for: OC and DI


## Read landmark data -------------------------------------------------------

table <- read.delim("Data/4Y.1000089.txt", header = FALSE, sep = "")
colnames(table) <- c("x", "y", "z")

if (nrow(table) < 11) {
  stop("Input landmark file must contain at least 11 landmark rows.")
}

head(table)


## Helper functions ---------------------------------------------------------

# 3D Euclidean distance between two landmark coordinates
distance_3D <- function(point1, point2) {
  sqrt(sum((point1 - point2)^2))
}

# 2D anatomical distance using the user-defined AP and SI axes
distance_AP_SI <- function(point1, point2, AP_axis, SI_axis) {
  sqrt(
    (point1[[AP_axis]] - point2[[AP_axis]])^2 +
      (point1[[SI_axis]] - point2[[SI_axis]])^2
  )
}

# Project a point onto a plane defined by:
# - a point lying on the plane (plane_point)
# - a normal vector to the plane (normal)
project_to_plane <- function(point, plane_point, normal) {
  parameter <- sum(normal * (plane_point - point)) / sum(normal^2)
  projected <- point + parameter * normal
  names(projected) <- c("x", "y", "z")
  projected
}


## Direct (DI) out-levers ---------------------------------------------------
# Direct 3D distances from the LEFT TMJ (landmark 5) to:
# - anterior bite point (landmark 11)
# - left posterior-most fully erupted tooth (landmark 10)

Out_Lever_DI_I <- distance_3D(table[5, ], table[11, ])
Out_Lever_DI_P <- distance_3D(table[5, ], table[10, ])

print(Out_Lever_DI_I)
print(Out_Lever_DI_P)


## Alveolar plane (AV) out-levers ------------------------------------------
# AV plane landmarks:
# 6 = Right posterior alveolar margin
# 7 = Left posterior alveolar margin
# 8 = Alveolar incisor
# 4,5 = Right and left jaw joints

# Centre of jaw rotation
JawAverage <- (table[4, ] + table[5, ]) / 2

# Define AV plane
Vector1_AV <- table[8, ] - table[6, ]
Vector2_AV <- table[8, ] - table[7, ]

Normal_AV <- c(
  (Vector1_AV[2] * Vector2_AV[3]) - (Vector1_AV[3] * Vector2_AV[2]),
  (Vector1_AV[3] * Vector2_AV[1]) - (Vector1_AV[1] * Vector2_AV[3]),
  (Vector1_AV[1] * Vector2_AV[2]) - (Vector1_AV[2] * Vector2_AV[1])
)

# Project jaw centre, anterior bite point and posterior tooth onto AV plane
Projected_J_AV <- project_to_plane(
  as.numeric(JawAverage),
  as.numeric(table[8, ]),
  Normal_AV
)

Projected_I_AV <- project_to_plane(
  as.numeric(table[11, ]),
  as.numeric(table[8, ]),
  Normal_AV
)

Projected_P_AV <- project_to_plane(
  as.numeric(table[10, ]),
  as.numeric(table[8, ]),
  Normal_AV
)

# Calculate AV out-levers in the AP-SI anatomical plane
Out_Lever_AV_I <- distance_AP_SI(
  as.list(Projected_J_AV),
  as.list(Projected_I_AV),
  AP_axis,
  SI_axis
)

Out_Lever_AV_P <- distance_AP_SI(
  as.list(Projected_J_AV),
  as.list(Projected_P_AV),
  AP_axis,
  SI_axis
)

print(Out_Lever_AV_I)
print(Out_Lever_AV_P)


## Frankfort Horizontal (FH) plane out-levers ------------------------------
# FH plane landmarks:
# 1 = Right EAM
# 2 = Left EAM
# 3 = Left orbit
# 4,5 = Right and left jaw joints

# Centre of jaw rotation
JawAverage <- (table[4, ] + table[5, ]) / 2

# Define FH plane
Vector1_FH <- table[3, ] - table[1, ]
Vector2_FH <- table[3, ] - table[2, ]

Normal_FH <- c(
  (Vector1_FH[2] * Vector2_FH[3]) - (Vector1_FH[3] * Vector2_FH[2]),
  (Vector1_FH[3] * Vector2_FH[1]) - (Vector1_FH[1] * Vector2_FH[3]),
  (Vector1_FH[1] * Vector2_FH[2]) - (Vector1_FH[2] * Vector2_FH[1])
)

# Project jaw centre, anterior bite point and posterior tooth onto FH plane
Projected_J_FH <- project_to_plane(
  as.numeric(JawAverage),
  as.numeric(table[3, ]),
  Normal_FH
)

Projected_I_FH <- project_to_plane(
  as.numeric(table[11, ]),
  as.numeric(table[3, ]),
  Normal_FH
)

Projected_P_FH <- project_to_plane(
  as.numeric(table[10, ]),
  as.numeric(table[3, ]),
  Normal_FH
)

# Calculate FH out-levers in the AP-SI anatomical plane
Out_Lever_FH_I <- distance_AP_SI(
  as.list(Projected_J_FH),
  as.list(Projected_I_FH),
  AP_axis,
  SI_axis
)

Out_Lever_FH_P <- distance_AP_SI(
  as.list(Projected_J_FH),
  as.list(Projected_P_FH),
  AP_axis,
  SI_axis
)

print(Out_Lever_FH_I)
print(Out_Lever_FH_P)


## Occlusal plane (OC) out-levers ------------------------------------------
# OC plane landmarks:
# 9  = Right posterior-most fully erupted tooth
# 10 = Left posterior-most fully erupted tooth
# 11 = Anterior bite point
# 4,5 = Right and left jaw joints

# Centre of jaw rotation
JawAverage <- (table[4, ] + table[5, ]) / 2

# Define OC plane
Vector1_OC <- table[11, ] - table[9, ]
Vector2_OC <- table[11, ] - table[10, ]

Normal_OC <- c(
  (Vector1_OC[2] * Vector2_OC[3]) - (Vector1_OC[3] * Vector2_OC[2]),
  (Vector1_OC[3] * Vector2_OC[1]) - (Vector1_OC[1] * Vector2_OC[3]),
  (Vector1_OC[1] * Vector2_OC[2]) - (Vector1_OC[2] * Vector2_OC[1])
)

# Project centre of jaw rotation onto OC plane
Projected_J_OC <- project_to_plane(
  as.numeric(JawAverage),
  as.numeric(table[11, ]),
  Normal_OC
)

# The anterior and posterior bite points already lie on the OC plane
Anterior_OC <- as.list(setNames(as.numeric(table[11, ]), c("x", "y", "z")))
Posterior_OC <- as.list(setNames(as.numeric(table[10, ]), c("x", "y", "z")))

# Calculate OC out-levers in the AP-SI anatomical plane
Out_Lever_OC_I <- distance_AP_SI(
  as.list(Projected_J_OC),
  Anterior_OC,
  AP_axis,
  SI_axis
)

Out_Lever_OC_P <- distance_AP_SI(
  as.list(Projected_J_OC),
  Posterior_OC,
  AP_axis,
  SI_axis
)

print(Out_Lever_OC_I)
print(Out_Lever_OC_P)


## Average FH and AV anterior out-lever --------------------------------------
# Calculate the mean anterior bite-point out-lever from the
# Frankfort Horizontal (FH) and Alveolar (AV) approaches.

Out_Lever_FA_I <- mean(c(Out_Lever_FH_I, Out_Lever_AV_I))

print(Out_Lever_FA_I)


## Output results -----------------------------------------------------------

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

library(openxlsx)

results_dir <- "Results"
file_name <- "Out_lever_method_sensitivity_output.xlsx"
file_path <- file.path(results_dir, file_name)

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

results <- data.frame(
  DI_Incisor_OUT = Out_Lever_DI_I,
  DI_Posterior_Molar_OUT = Out_Lever_DI_P,
  AV_Incisor_OUT = Out_Lever_AV_I,
  AV_Posterior_Molar_OUT = Out_Lever_AV_P,
  FH_Incisor_OUT = Out_Lever_FH_I,
  FH_Posterior_Molar_OUT = Out_Lever_FH_P,
  FA_Average_Incisor_OUT = Out_Lever_FA_I,
  OC_Incisor_OUT = Out_Lever_OC_I,
  OC_Posterior_Molar_OUT = Out_Lever_OC_P
)

write.xlsx(
  results,
  file = file_path,
  overwrite = TRUE
)

print(results)