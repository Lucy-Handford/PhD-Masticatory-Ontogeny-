# Muscle CSA plane sensitivity analysis
#
# Purpose:
# This script evaluates alternative anatomical landmarks for defining
# cross-sectional planes used to quantify masticatory muscle CSA.
#
# It:
# 1. Defines the Frankfort Horizontal (FH) plane.
# 2. Defines the zygomatic (ZG) plane.
# 3. Calculates perpendicular distances between these planes and
#    candidate muscle/anatomical landmarks.
#
# The analysis was used to assess the suitability of landmark-based
# plane definitions across the ontogenetic sample.

##Landmarks
#1	R EAM (Right External Auditory Meatus) 
#2	L EAM (Left External Auditory Meatus) 
#3	L ORB (Left Orbit)
#4	Most anterior point of temporalis origin attachment 
#5	Most medial part of the infratemporal crest
#6	R Most anterior and inferior part of the zygomatic bone
#7	L Most anterior and inferior part of the zygomatic bone
#8	L most posterior and inferior part of the zygomatic bone
#9 Angle of the mandible (most posterior and inferior point of masseter and medial pterygoid muscle attachment)
#10	Base of the Lingula


table <- read.delim('Data/Adult.txt', header=FALSE, sep="")
colnames(table) <- c("x", "y", "z")
head(table)


##Distance above FH plane temporalis----
#Read Table where 1 = R EAM, 2 = L EAM, 3 = R Orb, 4 = 	Most anterior point of temporalis origin attachment, 5 = Most medial part of the infra temporal crest 


# Ensure table is numeric
if (!is.matrix(table) && !is.numeric(table)) {
  table <- as.matrix(table)  # Convert to matrix
  table <- apply(table, 2, as.numeric)  # Ensure all columns are numeric
}

# Get Vectors to define the FH plane using landmarks
Vector1 <- table[3, ] - table[1, ]
Vector2 <- table[3, ] - table[2, ]

# Compute cross product (normal vector of the plane)
N_FH <- c(
  (Vector1[2] * Vector2[3]) - (Vector1[3] * Vector2[2]),  # X component
  (Vector1[3] * Vector2[1]) - (Vector1[1] * Vector2[3]),  # Y component
  (Vector1[1] * Vector2[2]) - (Vector1[2] * Vector2[1])   # Z component
)

# Point on the FH plane
P0 <- table[1, ]

# Function to calculate the perpendicular distance from point L to the plane
perpendicular_distance <- function(L, P0, N) {
  numerator <- abs(sum((L - P0) * N))
  denominator <- sqrt(sum(N * N))
  distance <- numerator / denominator
  return(distance)
}

# Calculate the perpendicular distance from point 4 to the FH plane
L4 <- table[4, ]
distance_L4 <- perpendicular_distance(L4, P0, N_FH)
print(distance_L4)

# Calculate the perpendicular distance from point 5 to the FH plane
L5 <- table[5, ]
distance_L5 <- perpendicular_distance(L5, P0, N_FH)
print(distance_L5)





## Get Vectors to define the ZG plane using landmarks
Vector3 <- table[7, ] - table[6, ]
Vector4 <- table[7, ] - table[8, ]

# Compute cross product (normal vector of the plane)
N_FH <- c(
  (Vector3[2] * Vector4[3]) - (Vector3[3] * Vector4[2]),  # X component
  (Vector3[3] * Vector4[1]) - (Vector3[1] * Vector4[3]),  # Y component
  (Vector3[1] * Vector4[2]) - (Vector3[2] * Vector4[1])   # Z component
)

# Point on the ZG plane 
P0 <- table[7, ]

# Function to calculate the perpendicular distance from point L to the plane
perpendicular_distance <- function(L, P0, N) {
  numerator <- abs(sum((L - P0) * N))
  denominator <- sqrt(sum(N * N))
  distance <- numerator / denominator
  return(distance)
}

# Calculate the perpendicular distance from point 4 to the FH plane
L9 <- table[9, ]
distance_L9 <- perpendicular_distance(L9, P0, N_FH)
print(distance_L9)

# Calculate the perpendicular distance from point 5 to the FH plane
L10 <- table[10, ]
distance_L10 <- perpendicular_distance(L10, P0, N_FH)
print(distance_L10)



## Make graphs

library(ggplot2)

Data <- read.delim('Data/Temporalis2.txt')

ggplot(Data, aes(x = Distance)) +
  geom_point(aes(y = Adult_CSA, color = "Adult CSA"), size = 1.2) + 
  geom_point(aes(y = Juvenile_CSA, color = "Juvenile CSA"), size = 1.2) +
  labs(title = "Variation in Temporalis CSA",
       x = "Distance above FH plane (mm)",
       y = "CSA (mm²)") +
  theme_minimal() +
  scale_color_manual(name = "Legend", values = c("Adult CSA" = "blue", "Juvenile CSA" = "red")) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))  # Center the title



library(ggplot2)

Data <- read.delim('Data/Masseter.txt')

ggplot(Data, aes(x = Distance)) +
  geom_point(aes(y = Adult_CSA, color = "Adult CSA"), size = 1.2) + 
  geom_point(aes(y = Juvenile_CSA, color = "Juvenile CSA"), size = 1.2) +
  labs(title = "Variation in Masseter CSA",
       x = "Distance below ZG plane (mm)",
       y = "CSA (mm²)") +
  theme_minimal() +
  scale_color_manual(name = "Legend", values = c("Adult CSA" = "blue", "Juvenile CSA" = "red")) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))  # Center the title



library(ggplot2)

Data <- read.delim('Data/MP.txt')

ggplot(Data, aes(x = Distance)) +
  geom_point(aes(y = Adult_CSA, color = "Adult CSA"), size = 1.2) + 
  geom_point(aes(y = Juvenile_CSA, color = "Juvenile CSA"), size = 1.2) +
  labs(title = "Variation in Medial Pterygoid CSA",
       x = "Distance below ZG plane (mm)",
       y = "CSA (mm²)") +
  theme_minimal() +
  scale_color_manual(name = "Legend", values = c("Adult CSA" = "blue", "Juvenile CSA" = "red")) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))  # Center the title





