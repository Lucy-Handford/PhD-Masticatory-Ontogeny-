### Chapter 3 biomechanical lever calculations
#
# This script calculates:
# - Temporalis, masseter and medial pterygoid in-levers
# - Occlusal-plane (OC) out-levers for:
#     * anterior bite point
#     * posterior-most molar
#     * first molar (M1)
# - Alveolar-plane (AV) anterior out-lever
# - Frankfort Horizontal (FH) anterior out-lever
# - Mean FH/AV anterior out-lever
#
# Landmark numbering matches the landmark table reported in Chapter 3.
#
# NOTE ON SPECIMEN ORIENTATION:
# The specimens used in this thesis were oriented such that the x-coordinate
# represented the component excluded from the 2D lever calculation.
# Therefore, the in-lever and final 2D out-lever distances are calculated
# from the y and z coordinates, consistent with the original thesis analysis.


## Landmarks ----------------------------------------------------------------

# 1, 2   L and R external auditory meatus
#        Superior-most point of the external auditory meatus.
#
# 3      Left orbit
#        Inferior-most point of the lower border of the orbit.
#
# 4, 5   L and R anterior temporalis origin
#        Central point on the orbital septum.
#
# 6, 7   L and R posterior temporalis origin
#        Anterior-most point of attachment on the glenoid fossae.
#
# 8, 9   L and R superficial masseter origin
#        Most anterior and outer corner of the zygomatic process.
#
# 10, 11 L and R deep masseter origin
#        Posterior and inferior point on the zygomatic arch, located at the
#        base of the anterior slope of the articular eminence.
#
# 12, 13 L and R superficial medial pterygoid origin
#        Anterior-most point of attachment on the pyramidal process
#        of the palatine bone.
#
# 14, 15 L and R deep medial pterygoid origin
#        Most posterior point of attachment on the pterygoid fossa
#        and medial surface of the pterygoid plate.
#
# 16, 17 L and R centre of condyles
#        Centre of condyles in the FH plane.
#
# 18, 19 L and R anterior temporalis insertion
#        Anterior and inferior-most point on the anterior border of the ramus
#        forming the lateral boundary of the retromolar triangle.
#
# 20, 21 L and R posterior temporalis insertion
#        Most inferior point of the mandibular notch.
#
# 22, 23 L and R superficial masseter insertion
#        Anterior and inferior point on the lateral surface of the mandibular corpus.
#
# 24, 25 L and R deep masseter insertion
#        Posterior-most point on the border and upper two-thirds of the ramus.
#
# 26, 27 L and R superficial medial pterygoid insertion
#        Anterior and inferior-most point on the border of the ramus.
#
# 28, 29 L and R deep medial pterygoid insertion
#        Anterior and inferior-most point on the border of the ramus.
#
# 30, 31 L and R alveolar molars
#        Posterior-most point of the posterior-most tooth or bony crypt
#        on the mandible.
#
# 32     Centre of alveolar incisors
#        Alveolar incisor margin at the centre of mandibular medial incisors.
#
# 33     Centre of medial incisors
#        Centre of the medial incisors, or alveolar incisor margin where teeth
#        are absent.
#
# 34, 35 L and R most posterior molars
#        Centre of the most posterior fully erupted teeth used to define
#        the occlusal plane.
#
# 36     First molar (M1)
#        Centre of the first molar on the left or right side of the mandible.


## Packages and input --------------------------------------------------------

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
library(openxlsx)

table <- read.delim("Data/Specimen.txt", header = FALSE, sep = "")
colnames(table) <- c("x", "y", "z")

if (nrow(table) < 36) {
  stop("Input landmark file must contain at least 36 landmark rows.")
}

head(table)


## IN-LEVERS ----------------------------------------------------------------

# Generate centre points between anterior and posterior origins and insertions
# on the LEFT side.

Centre_LTemp_O <- (table[4, ] + table[6, ]) / 2
Centre_LTemp_I <- (table[18, ] + table[20, ]) / 2

Centre_LMass_O <- (table[8, ] + table[10, ]) / 2
Centre_LMass_I <- (table[22, ] + table[24, ]) / 2

Centre_LMedP_O <- (table[12, ] + table[14, ]) / 2
Centre_LMedP_I <- (table[26, ] + table[28, ]) / 2


# Generate centre points between anterior and posterior origins and insertions
# on the RIGHT side.

Centre_RTemp_O <- (table[5, ] + table[7, ]) / 2
Centre_RTemp_I <- (table[19, ] + table[21, ]) / 2

Centre_RMass_O <- (table[9, ] + table[11, ]) / 2
Centre_RMass_I <- (table[23, ] + table[25, ]) / 2

Centre_RMedP_O <- (table[13, ] + table[15, ]) / 2
Centre_RMedP_I <- (table[27, ] + table[29, ]) / 2


# Generate an average TMJ coordinate between left and right condyles.
# The x-coordinate is excluded from the 2D lever calculation for the
# orientation used in this thesis.

Y_TMJ <- (table[16, 2] + table[17, 2]) / 2
Z_TMJ <- (table[16, 3] + table[17, 3]) / 2


### Left temporalis in-lever -------------------------------------------------

LTemp_a <- sqrt(
  (Y_TMJ - Centre_LTemp_O$y)^2 +
    (Z_TMJ - Centre_LTemp_O$z)^2
)

LTemp_b <- sqrt(
  (Y_TMJ - Centre_LTemp_I$y)^2 +
    (Z_TMJ - Centre_LTemp_I$z)^2
)

LTemp_c <- sqrt(
  (Centre_LTemp_O$y - Centre_LTemp_I$y)^2 +
    (Centre_LTemp_O$z - Centre_LTemp_I$z)^2
)

LTemp_CosB <- (
  (LTemp_a^2 + LTemp_c^2 - LTemp_b^2) /
    (2 * LTemp_a * LTemp_c)
)

LTemp_B <- acos(LTemp_CosB)

Left_Temp_In_lever <- sin(LTemp_B) * LTemp_a


### Right temporalis in-lever ------------------------------------------------

RTemp_a <- sqrt(
  (Y_TMJ - Centre_RTemp_O$y)^2 +
    (Z_TMJ - Centre_RTemp_O$z)^2
)

RTemp_b <- sqrt(
  (Y_TMJ - Centre_RTemp_I$y)^2 +
    (Z_TMJ - Centre_RTemp_I$z)^2
)

RTemp_c <- sqrt(
  (Centre_RTemp_O$y - Centre_RTemp_I$y)^2 +
    (Centre_RTemp_O$z - Centre_RTemp_I$z)^2
)

RTemp_CosB <- (
  (RTemp_a^2 + RTemp_c^2 - RTemp_b^2) /
    (2 * RTemp_a * RTemp_c)
)

RTemp_B <- acos(RTemp_CosB)

Right_Temp_In_lever <- sin(RTemp_B) * RTemp_a

AV_Temp_In_lever <- (
  Left_Temp_In_lever + Right_Temp_In_lever
) / 2


### Left masseter in-lever ---------------------------------------------------

LMass_a <- sqrt(
  (Y_TMJ - Centre_LMass_O$y)^2 +
    (Z_TMJ - Centre_LMass_O$z)^2
)

LMass_b <- sqrt(
  (Y_TMJ - Centre_LMass_I$y)^2 +
    (Z_TMJ - Centre_LMass_I$z)^2
)

LMass_c <- sqrt(
  (Centre_LMass_O$y - Centre_LMass_I$y)^2 +
    (Centre_LMass_O$z - Centre_LMass_I$z)^2
)

LMass_CosB <- (
  (LMass_a^2 + LMass_c^2 - LMass_b^2) /
    (2 * LMass_a * LMass_c)
)

LMass_B <- acos(LMass_CosB)

Left_Mass_In_lever <- sin(LMass_B) * LMass_a


### Right masseter in-lever --------------------------------------------------

RMass_a <- sqrt(
  (Y_TMJ - Centre_RMass_O$y)^2 +
    (Z_TMJ - Centre_RMass_O$z)^2
)

RMass_b <- sqrt(
  (Y_TMJ - Centre_RMass_I$y)^2 +
    (Z_TMJ - Centre_RMass_I$z)^2
)

RMass_c <- sqrt(
  (Centre_RMass_O$y - Centre_RMass_I$y)^2 +
    (Centre_RMass_O$z - Centre_RMass_I$z)^2
)

RMass_CosB <- (
  (RMass_a^2 + RMass_c^2 - RMass_b^2) /
    (2 * RMass_a * RMass_c)
)

RMass_B <- acos(RMass_CosB)

Right_Mass_In_lever <- sin(RMass_B) * RMass_a

AV_Mass_In_lever <- (
  Left_Mass_In_lever + Right_Mass_In_lever
) / 2


### Left medial pterygoid in-lever ------------------------------------------

LMedP_a <- sqrt(
  (Y_TMJ - Centre_LMedP_O$y)^2 +
    (Z_TMJ - Centre_LMedP_O$z)^2
)

LMedP_b <- sqrt(
  (Y_TMJ - Centre_LMedP_I$y)^2 +
    (Z_TMJ - Centre_LMedP_I$z)^2
)

LMedP_c <- sqrt(
  (Centre_LMedP_O$y - Centre_LMedP_I$y)^2 +
    (Centre_LMedP_O$z - Centre_LMedP_I$z)^2
)

LMedP_CosB <- (
  (LMedP_a^2 + LMedP_c^2 - LMedP_b^2) /
    (2 * LMedP_a * LMedP_c)
)

LMedP_B <- acos(LMedP_CosB)

Left_MedP_In_lever <- sin(LMedP_B) * LMedP_a


### Right medial pterygoid in-lever -----------------------------------------

RMedP_a <- sqrt(
  (Y_TMJ - Centre_RMedP_O$y)^2 +
    (Z_TMJ - Centre_RMedP_O$z)^2
)

RMedP_b <- sqrt(
  (Y_TMJ - Centre_RMedP_I$y)^2 +
    (Z_TMJ - Centre_RMedP_I$z)^2
)

RMedP_c <- sqrt(
  (Centre_RMedP_O$y - Centre_RMedP_I$y)^2 +
    (Centre_RMedP_O$z - Centre_RMedP_I$z)^2
)

RMedP_CosB <- (
  (RMedP_a^2 + RMedP_c^2 - RMedP_b^2) /
    (2 * RMedP_a * RMedP_c)
)

RMedP_B <- acos(RMedP_CosB)

Right_MedP_In_lever <- sin(RMedP_B) * RMedP_a

AV_MedP_In_lever <- (
  Left_MedP_In_lever + Right_MedP_In_lever
) / 2


## OUT-LEVERS ---------------------------------------------------------------

## Alveolar plane (AV): ANTERIOR bite point only ----------------------------

# AV plane:
# 30,31 = left and right posterior alveolar molar landmarks
# 32    = alveolar incisor
# 16,17 = left and right condyles
# 33    = anterior bite point

JawAverage <- (table[16, ] + table[17, ]) / 2

Vector1_AV <- table[32, ] - table[30, ]
Vector2_AV <- table[32, ] - table[31, ]

crossproductX_AV <- (Vector1_AV[2] * Vector2_AV[3]) -
  (Vector1_AV[3] * Vector2_AV[2])

crossproductY_AV <- (Vector1_AV[3] * Vector2_AV[1]) -
  (Vector1_AV[1] * Vector2_AV[3])

crossproductZ_AV <- (Vector1_AV[1] * Vector2_AV[2]) -
  (Vector1_AV[2] * Vector2_AV[1])


# Project centre of jaw rotation onto AV plane

J_AV <- (
  crossproductX_AV * (table[32, 1] - JawAverage[1]) +
    crossproductY_AV * (table[32, 2] - JawAverage[2]) +
    crossproductZ_AV * (table[32, 3] - JawAverage[3])
) / (
  crossproductX_AV^2 +
    crossproductY_AV^2 +
    crossproductZ_AV^2
)

ProjectPointX_J_AV <- JawAverage[1] + J_AV * crossproductX_AV
ProjectPointY_J_AV <- JawAverage[2] + J_AV * crossproductY_AV
ProjectPointZ_J_AV <- JawAverage[3] + J_AV * crossproductZ_AV

ProjectedPoint3D_J_AV <- data.frame(
  x = ProjectPointX_J_AV,
  y = ProjectPointY_J_AV,
  z = ProjectPointZ_J_AV
)


# Project anterior bite point onto AV plane

I_AV <- (
  crossproductX_AV * (table[32, 1] - table[33, 1]) +
    crossproductY_AV * (table[32, 2] - table[33, 2]) +
    crossproductZ_AV * (table[32, 3] - table[33, 3])
) / (
  crossproductX_AV^2 +
    crossproductY_AV^2 +
    crossproductZ_AV^2
)

ProjectPointX_I_AV <- table[33, 1] + I_AV * crossproductX_AV
ProjectPointY_I_AV <- table[33, 2] + I_AV * crossproductY_AV
ProjectPointZ_I_AV <- table[33, 3] + I_AV * crossproductZ_AV

ProjectedPoint3D_I_AV <- data.frame(
  x = ProjectPointX_I_AV,
  y = ProjectPointY_I_AV,
  z = ProjectPointZ_I_AV
)


# Calculate anterior AV out-lever in the y-z plane

Out_Lever_AV_I <- sqrt(
  (ProjectedPoint3D_J_AV$y - ProjectedPoint3D_I_AV$y)^2 +
    (ProjectedPoint3D_J_AV$z - ProjectedPoint3D_I_AV$z)^2
)

print(Out_Lever_AV_I)


## Frankfort Horizontal (FH): ANTERIOR bite point only ----------------------

# FH plane:
# 1,2   = left and right EAM
# 3     = left orbit
# 16,17 = left and right condyles
# 33    = anterior bite point

JawAverage <- (table[16, ] + table[17, ]) / 2

Vector1_FH <- table[3, ] - table[1, ]
Vector2_FH <- table[3, ] - table[2, ]

crossproductX_FH <- (Vector1_FH[2] * Vector2_FH[3]) -
  (Vector1_FH[3] * Vector2_FH[2])

crossproductY_FH <- (Vector1_FH[3] * Vector2_FH[1]) -
  (Vector1_FH[1] * Vector2_FH[3])

crossproductZ_FH <- (Vector1_FH[1] * Vector2_FH[2]) -
  (Vector1_FH[2] * Vector2_FH[1])


# Project centre of jaw rotation onto FH plane

J_FH <- (
  crossproductX_FH * (table[3, 1] - JawAverage[1]) +
    crossproductY_FH * (table[3, 2] - JawAverage[2]) +
    crossproductZ_FH * (table[3, 3] - JawAverage[3])
) / (
  crossproductX_FH^2 +
    crossproductY_FH^2 +
    crossproductZ_FH^2
)

ProjectPointX_J_FH <- JawAverage[1] + J_FH * crossproductX_FH
ProjectPointY_J_FH <- JawAverage[2] + J_FH * crossproductY_FH
ProjectPointZ_J_FH <- JawAverage[3] + J_FH * crossproductZ_FH

ProjectedPoint3D_J_FH <- data.frame(
  x = ProjectPointX_J_FH,
  y = ProjectPointY_J_FH,
  z = ProjectPointZ_J_FH
)


# Project anterior bite point onto FH plane

I_FH <- (
  crossproductX_FH * (table[3, 1] - table[33, 1]) +
    crossproductY_FH * (table[3, 2] - table[33, 2]) +
    crossproductZ_FH * (table[3, 3] - table[33, 3])
) / (
  crossproductX_FH^2 +
    crossproductY_FH^2 +
    crossproductZ_FH^2
)

ProjectPointX_I_FH <- table[33, 1] + I_FH * crossproductX_FH
ProjectPointY_I_FH <- table[33, 2] + I_FH * crossproductY_FH
ProjectPointZ_I_FH <- table[33, 3] + I_FH * crossproductZ_FH

ProjectedPoint3D_I_FH <- data.frame(
  x = ProjectPointX_I_FH,
  y = ProjectPointY_I_FH,
  z = ProjectPointZ_I_FH
)


# Calculate anterior FH out-lever in the y-z plane

Out_Lever_FH_I <- sqrt(
  (ProjectedPoint3D_J_FH$y - ProjectedPoint3D_I_FH$y)^2 +
    (ProjectedPoint3D_J_FH$z - ProjectedPoint3D_I_FH$z)^2
)

print(Out_Lever_FH_I)


## Mean FH/AV anterior out-lever --------------------------------------------

Out_Lever_FA_I <- mean(
  c(Out_Lever_FH_I, Out_Lever_AV_I)
)

print(Out_Lever_FA_I)


## Occlusal plane (OC): ALL bite points ------------------------------------

# OC plane:
# 34,35 = left and right posterior-most fully erupted molars
# 33    = anterior bite point
# 16,17 = left and right condyles
# 36    = first molar (M1)
#
# The OC plane is defined using landmarks 33, 34 and 35.
# The anterior bite point and posterior-most molar lie on this plane.
# The M1 landmark is projected to the OC plane before its out-lever
# distance is calculated.

JawAverage <- (table[16, ] + table[17, ]) / 2

Vector1_OC <- table[33, ] - table[34, ]
Vector2_OC <- table[33, ] - table[35, ]

crossproductX_OC <- (Vector1_OC[2] * Vector2_OC[3]) -
  (Vector1_OC[3] * Vector2_OC[2])

crossproductY_OC <- (Vector1_OC[3] * Vector2_OC[1]) -
  (Vector1_OC[1] * Vector2_OC[3])

crossproductZ_OC <- (Vector1_OC[1] * Vector2_OC[2]) -
  (Vector1_OC[2] * Vector2_OC[1])


# Project centre of jaw rotation onto OC plane

J_OC <- (
  crossproductX_OC * (table[33, 1] - JawAverage[1]) +
    crossproductY_OC * (table[33, 2] - JawAverage[2]) +
    crossproductZ_OC * (table[33, 3] - JawAverage[3])
) / (
  crossproductX_OC^2 +
    crossproductY_OC^2 +
    crossproductZ_OC^2
)

ProjectPointX_J_OC <- JawAverage[1] + J_OC * crossproductX_OC
ProjectPointY_J_OC <- JawAverage[2] + J_OC * crossproductY_OC
ProjectPointZ_J_OC <- JawAverage[3] + J_OC * crossproductZ_OC

ProjectedPoint3D_J_OC <- data.frame(
  x = ProjectPointX_J_OC,
  y = ProjectPointY_J_OC,
  z = ProjectPointZ_J_OC
)


# Project M1 (landmark 36) onto OC plane

M1_OC <- (
  crossproductX_OC * (table[33, 1] - table[36, 1]) +
    crossproductY_OC * (table[33, 2] - table[36, 2]) +
    crossproductZ_OC * (table[33, 3] - table[36, 3])
) / (
  crossproductX_OC^2 +
    crossproductY_OC^2 +
    crossproductZ_OC^2
)

ProjectPointX_M1_OC <- table[36, 1] + M1_OC * crossproductX_OC
ProjectPointY_M1_OC <- table[36, 2] + M1_OC * crossproductY_OC
ProjectPointZ_M1_OC <- table[36, 3] + M1_OC * crossproductZ_OC

ProjectedPoint3D_M1_OC <- data.frame(
  x = ProjectPointX_M1_OC,
  y = ProjectPointY_M1_OC,
  z = ProjectPointZ_M1_OC
)


# Bite-point coordinates on the OC plane

Incisor_y <- table[33, 2]
Incisor_z <- table[33, 3]

Posterior_y <- table[34, 2]
Posterior_z <- table[34, 3]


# Calculate OC anterior out-lever

Out_Lever_OC_I <- sqrt(
  (ProjectedPoint3D_J_OC$y - Incisor_y)^2 +
    (ProjectedPoint3D_J_OC$z - Incisor_z)^2
)


# Calculate OC posterior-most molar out-lever

Out_Lever_OC_P <- sqrt(
  (ProjectedPoint3D_J_OC$y - Posterior_y)^2 +
    (ProjectedPoint3D_J_OC$z - Posterior_z)^2
)


# Calculate OC M1 out-lever using projected M1 coordinate

Out_Lever_OC_M1 <- sqrt(
  (ProjectedPoint3D_J_OC$y - ProjectedPoint3D_M1_OC$y)^2 +
    (ProjectedPoint3D_J_OC$z - ProjectedPoint3D_M1_OC$z)^2
)

print(Out_Lever_OC_I)
print(Out_Lever_OC_P)
print(Out_Lever_OC_M1)


## Output data --------------------------------------------------------------

results_dir <- "Results"
file_name <- "Chapter_3_lever_calculations_output.xlsx"
file_path <- file.path(results_dir, file_name)

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

results <- data.frame(
  Left_Temporalis_IN = Left_Temp_In_lever,
  Right_Temporalis_IN = Right_Temp_In_lever,
  Average_Temporalis_IN = AV_Temp_In_lever,
  Left_Masseter_IN = Left_Mass_In_lever,
  Right_Masseter_IN = Right_Mass_In_lever,
  Average_Masseter_IN = AV_Mass_In_lever,
  Left_Medial_Pterygoid_IN = Left_MedP_In_lever,
  Right_Medial_Pterygoid_IN = Right_MedP_In_lever,
  Average_Medial_Pterygoid_IN = AV_MedP_In_lever,
  AV_Anterior_OUT = Out_Lever_AV_I,
  FH_Anterior_OUT = Out_Lever_FH_I,
  FH_AV_Average_Anterior_OUT = Out_Lever_FA_I,
  OC_Anterior_OUT = Out_Lever_OC_I,
  OC_Posterior_Molar_OUT = Out_Lever_OC_P,
  OC_M1_OUT = Out_Lever_OC_M1
)

write.xlsx(
  results,
  file = file_path,
  overwrite = TRUE
)

print(results)