library(sampling)
library(randomForest)
library(nnet)
library(data.table)


# The data can be downloaded via http://amelia.uni-trier.de/?page_id=121
AMELIA_files <- paste0(path, list.files(path))

for(i in 1:6){
  load(AMELIA_files[i])
}

# Households bigger than size 6 are subsumed under one category
HHS <- ifelse(HHS >= 6, 6, HHS)
# Create categories for age; this is because one of the later used strategies only works
# with catgeorical data
AGE_cat <- rep("77+x", length(AGE))
AGE_cat[AGE < 18] <- "0-17"
AGE_cat[AGE > 17 & AGE <= 28] <- "18-28"
AGE_cat[AGE > 28 & AGE <= 39] <- "29-39"
AGE_cat[AGE > 40 & AGE <= 50] <- "40-50"
AGE_cat[AGE > 50 & AGE <= 61] <- "51-61"
AGE_cat[AGE > 61 & AGE <= 72] <- "62-72"
AGE_cat[AGE > 72 & AGE <= 77] <- "73-77"

# "AMELIA" with PIDs in decreasing order of age
AMELIA <- data.table(HID, HHS, AGE, AGE_cat, MST, SEX, INC) 

AMELIA[order(HID, AGE)] # Sort according to age within a household
# Create person ID within a Household
AMELIA[, PID := order(AGE, decreasing = TRUE), HID] 
AMELIA[, AGE_cat:=as.factor(AGE_cat)]
AMELIA[, MST:=as.factor(MST)]
AMELIA[, SEX:=as.factor(SEX)]
AMELIA[, HHS:=as.factor(HHS)]

# Remove unused variables
rm(AGE, MST, SEX, HHS, HID, PID, EDU, INC)

# Create a unique person ID for every single individual
AMELIA$PIDX <- rep(1:nrow(AMELIA))

# AMELIA_HH_wide will be our true household dataset, from which we can subset different 
# different groups of households. It is based on the individual data and their unqiue identifier
AMELIA_HH_wide <- dcast(AMELIA, HID + HHS ~ PID, value.var = c("AGE_cat", "AGE", "MST", "SEX", "INC", "PIDX"))

# Draw Single-Stage-Cluster Sample (This one sample serves as basis for every single analysis; no
# variance estimation via Monte-Carlo simulation was carried out)

SI_SAMPLE <- srswor(37812, nrow(AMELIA_HH_wide)) # 1-percent sample, see for example german microcensus
AMELIA_wide_sample <- AMELIA_HH_wide[SI_SAMPLE>0, ]
head(AMELIA_wide_sample) # AMELIA_wide_sample is now the sample to refer to in the code

# At this point it must be noted that we can now proceed differently for the modelling. On the one hand, 
# we can regard the members in a household as a specific combination of members and characteristics for
# whose occurrence a certain probability is estimated via a model (combination approach). 
# However, it is also possible to consider the characteristics of 
# one household member as a dependent variable, which is explained by the characteristics of the 
# other members (modeling approach). Both approaches will be exmined below.

###################
###################
# COMBINATION APPROCH STARTS HERE
################################

# Overall, the method used can be explained quite simply: the rows coded
# as "true combinations" based on a Single stage cluster sample of households (PSUs) and their
# members (SSUS) are set against as many random combinations as possible, based on the individuals
# from the same sample which are randomly combined with other individuals into new random households.
# They are coded with "0".
# The subsequent application of the model works by randomly creating households and let
# the model evaluate their quality. 

# We can improve model quality by setting the true values against a larger number of 
# possible random combinations. 

# The following function takes a sample of individuals and returns a set of random combinations
# of individuals in households. Like explained above, this function allows to build the model 
# data 

ind_to_random_hid <- function(x, household_size, drop, synthetic_hhs) {
  random_1 <- x
  if(household_size == 2) {
  random_1$HID_rand <- rep(1:(nrow(random_1)/household_size), times= household_size, 10^8, 
                           replace=FALSE)[1:nrow(random_1)]
  random_1 <- random_1[order(HID_rand, AGE)] 
  random_1 <- random_1[, PID := order(AGE, decreasing = TRUE), HID_rand]
  random_1 <- dcast(random_1, HID_rand + synthetic_hhs ~ PID, 
                    value.var = c("AGE_cat", "AGE", "MST", "SEX", "INC", "PIDX"))
  random_1$AGE_mean <- (random_1$AGE_1 + random_1$AGE_2) / 2
  random_1$AGE_diff <- sqrt(((random_1$AGE_mean - random_1$AGE_1)^2 + (random_1$AGE_mean - random_1$AGE_2)^2) /2)
  random_1$SEX_diff <- interaction(random_1$SEX_1, random_1$SEX_2, drop = drop)
  random_1$INC <- random_1$INC_1 + random_1$INC_2
  f1 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2))
  f2 <- pmax(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2))
  random_1$MST_group <- droplevels(interaction(f1, f2, sep=""))
  }
  else if(household_size == 3) {
    random_1$HID_rand <- rep(1:(nrow(random_1)/household_size), times= household_size, 10^8, 
                             replace=FALSE)[1:nrow(random_1)]
    random_1 <- random_1[order(HID_rand, AGE)] 
    random_1 <- random_1[, PID := order(AGE, decreasing = TRUE), HID_rand]
    random_1 <- dcast(random_1, HID_rand + synthetic_hhs ~ PID, 
                      value.var = c("AGE_cat", "AGE", "MST", "SEX", "INC", "PIDX"))
    random_1$AGE_mean <- (random_1$AGE_1 + random_1$AGE_2 + random_1$AGE_3) / 3
    random_1$AGE_diff <- sqrt(((random_1$AGE_mean - random_1$AGE_1)^2 + (random_1$AGE_mean - random_1$AGE_2)^2 +
      (random_1$AGE_mean - random_1$AGE_3)^2) / 3) # Standardabweichung von Hand berechnet! Und da sage noch einer, Grundlagen seien "trocken"...
    random_1$SEX_diff <- interaction(random_1$SEX_1, random_1$SEX_2, random_1$SEX_3, drop = drop)
    random_1$INC <- random_1$INC_1 + random_1$INC_2 + random_1$INC_3
    f1 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2))
    f2 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3))
    f3 <- pmax(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3))
    random_1$MST_group <- droplevels(interaction(f1, f2, f3, sep = ""))
  }
  else if(household_size == 4) {
    random_1$HID_rand <- rep(1:(nrow(random_1)/household_size), times= household_size, 10^8, 
                             replace=FALSE)[1:nrow(random_1)]
    random_1 <- random_1[order(HID_rand, AGE)] 
    random_1 <- random_1[, PID := order(AGE, decreasing = TRUE), HID_rand]
    random_1 <- dcast(random_1, HID_rand + synthetic_hhs ~ PID, 
                      value.var = c("AGE_cat", "AGE", "MST", "SEX", "INC", "PIDX"))
    random_1$AGE_mean <- (random_1$AGE_1 + random_1$AGE_2 + random_1$AGE_3 + random_1$AGE_4) / 4
    random_1$AGE_diff <- sqrt(((random_1$AGE_mean - random_1$AGE_1)^2 + (random_1$AGE_mean - random_1$AGE_2)^2 +
                                 (random_1$AGE_mean - random_1$AGE_3)^2 + (random_1$AGE_mean - random_1$AGE_4)^2) / 4) # Standardabweichung von Hand berechnet! Und da sage noch einer, Grundlagen seien "trocken"...
    random_1$SEX_diff <- interaction(random_1$SEX_1, random_1$SEX_2, random_1$SEX_3, random_1$SEX_4, 
                                     drop = drop)
    random_1$INC <- random_1$INC_1 + random_1$INC_2 + random_1$INC_3 + random_1$INC_4
    f1 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2))
    f2 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3))
    f3 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3), 
               as.numeric(random_1$MST_4))
    f4 <- pmax(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3),
               as.numeric(random_1$MST_4))
    random_1$MST_group <- droplevels(interaction(f1, f2, f3, f4, sep = ""))
  }
    else if(household_size == 5) {
      random_1$HID_rand <- rep(1:(nrow(random_1)/household_size), times= household_size, 10^8, 
                               replace=FALSE)[1:nrow(random_1)]
      random_1 <- random_1[order(HID_rand, AGE)] 
      random_1 <- random_1[, PID := order(AGE, decreasing = TRUE), HID_rand]
      random_1 <- dcast(random_1, HID_rand + synthetic_hhs ~ PID, 
                        value.var = c("AGE_cat", "AGE", "MST", "SEX", "INC", "PIDX"))
      random_1$AGE_mean <- (random_1$AGE_1 + random_1$AGE_2 + random_1$AGE_3 + random_1$AGE_4 +
                              random_1$AGE_5) / 5
      random_1$AGE_diff <- sqrt(((random_1$AGE_mean - random_1$AGE_1)^2 + (random_1$AGE_mean - random_1$AGE_2)^2 +
                                   (random_1$AGE_mean - random_1$AGE_3)^2 + (random_1$AGE_mean - random_1$AGE_4)^2 +
                                  (random_1$AGE_mean - random_1$AGE_5)^2) / 5) # Standardabweichung von Hand berechnet! Und da sage noch einer, Grundlagen seien "trocken"...
      random_1$SEX_diff <- interaction(random_1$SEX_1, random_1$SEX_2, random_1$SEX_3, random_1$SEX_4,
                                       random_1$SEX_5, drop = drop)
      random_1$INC <- random_1$INC_1 + random_1$INC_2 + random_1$INC_3 + random_1$INC_4 + 
        random_1$INC_5
      f1 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2))
      f2 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3))
      f3 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3), 
                 as.numeric(random_1$MST_4))
      f4 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3),
                 as.numeric(random_1$MST_4), as.numeric(random_1$MST_5))
      f5 <- pmax(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3),
                 as.numeric(random_1$MST_4), as.numeric(random_1$MST_5))
      random_1$MST_group <- droplevels(interaction(f1, f2, f3, f4, f5, sep = ""))
    } 
  else if(household_size== 6) {
        random_1$HID_rand <- rep(1:(nrow(random_1)/household_size), times= household_size, 10^8, 
                                 replace=FALSE)[1:nrow(random_1)]
        random_1 <- random_1[order(HID_rand, AGE)] 
        random_1 <- random_1[, PID := order(AGE, decreasing = TRUE), HID_rand]
        random_1 <- dcast(random_1, HID_rand + synthetic_hhs ~ PID, 
                          value.var = c("AGE_cat", "AGE", "MST", "SEX", "INC", "PIDX"))
        random_1$AGE_mean <- (random_1$AGE_1 + random_1$AGE_2 + random_1$AGE_3 + random_1$AGE_4 +
                                random_1$AGE_5 + random_1$AGE_6) / 6
        random_1$AGE_diff <- sqrt(((random_1$AGE_mean - random_1$AGE_1)^2 + (random_1$AGE_mean - random_1$AGE_2)^2 +
                                     (random_1$AGE_mean - random_1$AGE_3)^2 + (random_1$AGE_mean - random_1$AGE_4)^2 +
                                     (random_1$AGE_mean - random_1$AGE_5)^2 +
                                     (random_1$AGE_mean - random_1$AGE_6)^2) / 6) # Standardabweichung von Hand berechnet! Und da sage noch einer, Grundlagen seien "trocken"...
        random_1$INC <- random_1$INC_1 + random_1$INC_2 + random_1$INC_3 + random_1$INC_4 + 
          random_1$INC_5 + random_1$INC_6
        f1 <- pmin(as.numeric(random_1$SEX_1), as.numeric(random_1$SEX_2))
        f2 <- pmin(as.numeric(random_1$SEX_1), as.numeric(random_1$SEX_2), as.numeric(random_1$SEX_3))
        f3 <- pmin(as.numeric(random_1$SEX_1), as.numeric(random_1$SEX_2), as.numeric(random_1$SEX_3),
                   as.numeric(random_1$SEX_4))
        f4 <- pmin(as.numeric(random_1$SEX_1), as.numeric(random_1$SEX_2), as.numeric(random_1$SEX_3),
                   as.numeric(random_1$SEX_4), as.numeric(random_1$SEX_5))
        f5 <- pmin(as.numeric(random_1$SEX_1), as.numeric(random_1$SEX_2), as.numeric(random_1$SEX_3),
                   as.numeric(random_1$SEX_4), as.numeric(random_1$SEX_5), as.numeric(random_1$SEX_6))
        f6 <- pmax(as.numeric(random_1$SEX_1), as.numeric(random_1$SEX_2), as.numeric(random_1$SEX_3),
                   as.numeric(random_1$SEX_4), as.numeric(random_1$SEX_5), as.numeric(random_1$SEX_6))
        random_1$SEX_diff <- droplevels(interaction(f1, f2, f3, f4, f5, f6, sep = ""))
        f1 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2))
        f2 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3))
        f3 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3), 
                   as.numeric(random_1$MST_4))
        f4 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3),
                   as.numeric(random_1$MST_4), as.numeric(random_1$MST_5))
        f5 <- pmin(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3),
                   as.numeric(random_1$MST_4), as.numeric(random_1$MST_5), as.numeric(random_1$MST_6))
        f6 <- pmax(as.numeric(random_1$MST_1), as.numeric(random_1$MST_2), as.numeric(random_1$MST_3),
                   as.numeric(random_1$MST_4), as.numeric(random_1$MST_5), as.numeric(random_1$MST_6))
        random_1$MST_group <- droplevels(interaction(f1, f2, f3, f4, f5, f6, sep = ""))
  }
  hid_data <- random_1
  return(hid_data)
}

### Let´s subset the true household distribution for every household size

hhs_2_true <- AMELIA_HH_wide[,c("HID", "HHS", "AGE_cat_1", "AGE_cat_2", "AGE_1", "AGE_2","SEX_1", 
                                "SEX_2", "MST_1", "MST_2", "INC_1", "INC_2", "PIDX_1", "PIDX_2")][AMELIA_HH_wide$HHS == 2]

hhs_2_true$AGE_mean <- (hhs_2_true$AGE_1 + hhs_2_true$AGE_2) / 2
hhs_2_true$AGE_diff <- sqrt(((hhs_2_true$AGE_mean - hhs_2_true$AGE_1)^2 + (hhs_2_true$AGE_mean - hhs_2_true$AGE_2)^2) /2)

hhs_2_true$INC <- hhs_2_true$INC_1 + hhs_2_true$INC_2 

# Using combinations instead of perumutations. Therefore, we order the factor categories
f1 <- pmin(as.numeric(hhs_2_true$MST_1), as.numeric(hhs_2_true$MST_2))
f2 <- pmax(as.numeric(hhs_2_true$MST_1), as.numeric(hhs_2_true$MST_2))
hhs_2_true$MST_group <- droplevels(interaction(f1, f2, sep=""))

f1 <- pmin(as.numeric(hhs_2_true$AGE_cat_1), as.numeric(hhs_2_true$AGE_cat_2))
f2 <- pmax(as.numeric(hhs_2_true$AGE_cat_1), as.numeric(hhs_2_true$AGE_cat_2))
hhs_2_true$AGE_comb <- droplevels(interaction(f1, f2, sep=""))

f1 <- pmin(as.numeric(hhs_2_true$SEX_1), as.numeric(hhs_2_true$SEX_2))
f2 <- pmax(as.numeric(hhs_2_true$SEX_1), as.numeric(hhs_2_true$SEX_2))
hhs_2_true$SEX_diff <- droplevels(interaction(f1, f2, sep=""))


hhs_3_true <- AMELIA_HH_wide[,c("HID", "HHS", "AGE_cat_1", "AGE_cat_2", "AGE_cat_3", "AGE_1", "AGE_2",
                                "AGE_3", "SEX_1", "SEX_2", "SEX_3", "MST_1", "MST_2", "MST_3", 
                                "INC_1", "INC_2", "INC_3", "PIDX_1", 
                                "PIDX_2", "PIDX_3")][AMELIA_HH_wide$HHS == 3]

hhs_3_true$AGE_mean <- (hhs_3_true$AGE_1 + hhs_3_true$AGE_2 + hhs_3_true$AGE_3) / 3
hhs_3_true$AGE_diff <- sqrt(((hhs_3_true$AGE_mean - hhs_3_true$AGE_1)^2 + (hhs_3_true$AGE_mean - hhs_3_true$AGE_2)^2 +
                             (hhs_3_true$AGE_mean - hhs_3_true$AGE_3)^2) / 3) # Standardabweichung von Hand berechnet! Und da sage noch einer, Grundlagen seien "trocken"...
hhs_3_true$INC <- hhs_3_true$INC_1 + hhs_3_true$INC_2 + hhs_3_true$INC_3

f1 <- pmin(as.numeric(hhs_3_true$MST_1), as.numeric(hhs_3_true$MST_2))
f2 <- pmin(as.numeric(hhs_3_true$MST_1), as.numeric(hhs_3_true$MST_2), as.numeric(hhs_3_true$MST_3))
f3 <- pmax(as.numeric(hhs_3_true$MST_1), as.numeric(hhs_3_true$MST_2), as.numeric(hhs_3_true$MST_3))

hhs_3_true$MST_group <- droplevels(interaction(f1, f2, f3, sep = ""))

f1 <- pmin(as.numeric(hhs_3_true$AGE_cat_1), as.numeric(hhs_3_true$AGE_cat_2))
f2 <- pmin(as.numeric(hhs_3_true$AGE_cat_1), as.numeric(hhs_3_true$AGE_cat_2), as.numeric(hhs_3_true$AGE_cat_3))
f3 <- pmax(as.numeric(hhs_3_true$AGE_cat_1), as.numeric(hhs_3_true$AGE_cat_2), as.numeric(hhs_3_true$AGE_cat_3))

hhs_3_true$AGE_comb <- droplevels(interaction(f1, f2, f3, sep = ""))

f1 <- pmin(as.numeric(hhs_3_true$SEX_1), as.numeric(hhs_3_true$SEX_2))
f2 <- pmin(as.numeric(hhs_3_true$SEX_1), as.numeric(hhs_3_true$SEX_2), as.numeric(hhs_3_true$SEX_3))
f3 <- pmax(as.numeric(hhs_3_true$SEX_1), as.numeric(hhs_3_true$SEX_2), as.numeric(hhs_3_true$SEX_3))
hhs_3_true$SEX_diff <- droplevels(interaction(f1, f2, f3, sep = ""))

# Same for households of size 4

hhs_4_true <- AMELIA_HH_wide[,c("HID", "HHS", "AGE_cat_1", "AGE_cat_2", "AGE_cat_3", "AGE_cat_4",
                                "AGE_1", "AGE_2","AGE_3","AGE_4", "SEX_1", "SEX_2", "SEX_3", "SEX_4", 
                                "MST_1", "MST_2", "MST_3", "MST_4", "INC_1", "INC_2", "INC_3", "INC_4",
                                "PIDX_1", "PIDX_2", "PIDX_3", "PIDX_4")][AMELIA_HH_wide$HHS == 4]

hhs_4_true$AGE_mean <- (hhs_4_true$AGE_1 + hhs_4_true$AGE_2 + hhs_4_true$AGE_3 + hhs_4_true$AGE_4) / 4
hhs_4_true$AGE_diff <- sqrt(((hhs_4_true$AGE_mean - hhs_4_true$AGE_1)^2 + (hhs_4_true$AGE_mean - hhs_4_true$AGE_2)^2 +
                               (hhs_4_true$AGE_mean - hhs_4_true$AGE_3)^2 +
                               (hhs_4_true$AGE_mean - hhs_4_true$AGE_4)^2)) / 4 # Standardabweichung von Hand berechnet! Und da sage noch einer, Grundlagen seien "trocken"...
hhs_4_true$INC <- hhs_4_true$INC_1 + hhs_4_true$INC_2 + hhs_4_true$INC_3 + hhs_4_true$INC_4

f1 <- pmin(as.numeric(hhs_4_true$MST_1), as.numeric(hhs_4_true$MST_2))
f2 <- pmin(as.numeric(hhs_4_true$MST_1), as.numeric(hhs_4_true$MST_2), as.numeric(hhs_4_true$MST_3))
f3 <- pmin(as.numeric(hhs_4_true$MST_1), as.numeric(hhs_4_true$MST_2), as.numeric(hhs_4_true$MST_3),
           as.numeric(hhs_4_true$MST_4))
f4 <- pmax(as.numeric(hhs_4_true$MST_1), as.numeric(hhs_4_true$MST_2), as.numeric(hhs_4_true$MST_3),
           as.numeric(hhs_4_true$MST_4))

hhs_4_true$MST_group <- droplevels(interaction(f1, f2, f3, f4, sep = ""))

f1 <- pmin(as.numeric(hhs_4_true$AGE_cat_1), as.numeric(hhs_4_true$AGE_cat_2))
f2 <- pmin(as.numeric(hhs_4_true$AGE_cat_1), as.numeric(hhs_4_true$AGE_cat_2), as.numeric(hhs_4_true$AGE_cat_3))
f3 <- pmin(as.numeric(hhs_4_true$AGE_cat_1), as.numeric(hhs_4_true$AGE_cat_2), as.numeric(hhs_4_true$AGE_cat_3),
           as.numeric(hhs_4_true$AGE_cat_4))
f4 <- pmax(as.numeric(hhs_4_true$AGE_cat_1), as.numeric(hhs_4_true$AGE_cat_2), as.numeric(hhs_4_true$AGE_cat_3),
           as.numeric(hhs_4_true$AGE_cat_4))

hhs_4_true$AGE_comb <- droplevels(interaction(f1, f2, f3, f4, sep = ""))

f1 <- pmin(as.numeric(hhs_4_true$SEX_1), as.numeric(hhs_4_true$SEX_2))
f2 <- pmin(as.numeric(hhs_4_true$SEX_1), as.numeric(hhs_4_true$SEX_2), as.numeric(hhs_4_true$SEX_3))
f3 <- pmin(as.numeric(hhs_4_true$SEX_1), as.numeric(hhs_4_true$SEX_2), as.numeric(hhs_4_true$SEX_3),
           as.numeric(hhs_4_true$SEX_4))
f4 <- pmax(as.numeric(hhs_4_true$SEX_1), as.numeric(hhs_4_true$SEX_2), as.numeric(hhs_4_true$SEX_3),
           as.numeric(hhs_4_true$SEX_4))

hhs_4_true$SEX_diff <- droplevels(interaction(f1, f2, f3, f4, sep = ""))

# Same for households of size 5

hhs_5_true <- AMELIA_HH_wide[,c("HID", "HHS", "AGE_cat_1", "AGE_cat_2", "AGE_cat_3", "AGE_cat_4",
                                "AGE_cat_5",  "AGE_1", "AGE_2","AGE_3","AGE_4", "AGE_5", "SEX_1", "SEX_2", 
                                "SEX_3","SEX_4", "SEX_5",
                                "MST_1", "MST_2", "MST_3", "MST_4", "MST_5", "INC_1", "INC_2",
                                "INC_3", "INC_4", "INC_5", "PIDX_1", 
                                "PIDX_2", "PIDX_3", "PIDX_4", "PIDX_5")][AMELIA_HH_wide$HHS == 5]

hhs_5_true$AGE_mean <- (hhs_5_true$AGE_1 + hhs_5_true$AGE_2 + hhs_5_true$AGE_3 + hhs_5_true$AGE_4 +
                          hhs_5_true$AGE_5) / 5
hhs_5_true$AGE_diff <- sqrt(((hhs_5_true$AGE_mean - hhs_5_true$AGE_1)^2 + (hhs_5_true$AGE_mean - hhs_5_true$AGE_2)^2 +
                               (hhs_5_true$AGE_mean - hhs_5_true$AGE_3)^2 +
                               (hhs_5_true$AGE_mean - hhs_5_true$AGE_4)^2 + 
                               (hhs_5_true$AGE_mean - hhs_5_true$AGE_5)^2)) / 5 
hhs_5_true$INC <- hhs_5_true$INC_1 + hhs_5_true$INC_2 + hhs_5_true$INC_3 + hhs_5_true$INC_4 + hhs_5_true$INC_5

f1 <- pmin(as.numeric(hhs_5_true$MST_1), as.numeric(hhs_5_true$MST_2))
f2 <- pmin(as.numeric(hhs_5_true$MST_1), as.numeric(hhs_5_true$MST_2), as.numeric(hhs_5_true$MST_3))
f3 <- pmin(as.numeric(hhs_5_true$MST_1), as.numeric(hhs_5_true$MST_2), as.numeric(hhs_5_true$MST_3),
           as.numeric(hhs_5_true$MST_4))
f4 <- pmin(as.numeric(hhs_5_true$MST_1), as.numeric(hhs_5_true$MST_2), as.numeric(hhs_5_true$MST_3),
           as.numeric(hhs_5_true$MST_4), as.numeric(hhs_5_true$MST_5))
f5 <- pmax(as.numeric(hhs_5_true$MST_1), as.numeric(hhs_5_true$MST_2), as.numeric(hhs_5_true$MST_3),
                 as.numeric(hhs_5_true$MST_4), as.numeric(hhs_5_true$MST_5))

hhs_5_true$MST_group <- droplevels(interaction(f1, f2, f3, f4, f5, sep = ""))

f1 <- pmin(as.numeric(hhs_5_true$AGE_cat_1), as.numeric(hhs_5_true$AGE_cat_2))
f2 <- pmin(as.numeric(hhs_5_true$AGE_cat_1), as.numeric(hhs_5_true$AGE_cat_2), as.numeric(hhs_5_true$AGE_cat_3))
f3 <- pmin(as.numeric(hhs_5_true$AGE_cat_1), as.numeric(hhs_5_true$AGE_cat_2), as.numeric(hhs_5_true$AGE_cat_3),
           as.numeric(hhs_5_true$AGE_cat_4))
f4 <- pmin(as.numeric(hhs_5_true$AGE_cat_1), as.numeric(hhs_5_true$AGE_cat_2), as.numeric(hhs_5_true$AGE_cat_3),
           as.numeric(hhs_5_true$AGE_cat_4), as.numeric(hhs_5_true$AGE_cat_5))
f5 <- pmax(as.numeric(hhs_5_true$AGE_cat_1), as.numeric(hhs_5_true$AGE_cat_2), as.numeric(hhs_5_true$AGE_cat_3),
           as.numeric(hhs_5_true$AGE_cat_4), as.numeric(hhs_5_true$AGE_cat_5))

hhs_5_true$AGE_comb <- droplevels(interaction(f1, f2, f3, f4, f5, sep = ""))

f1 <- pmin(as.numeric(hhs_5_true$SEX_1), as.numeric(hhs_5_true$SEX_2))
f2 <- pmin(as.numeric(hhs_5_true$SEX_1), as.numeric(hhs_5_true$SEX_2), as.numeric(hhs_5_true$SEX_3))
f3 <- pmin(as.numeric(hhs_5_true$SEX_1), as.numeric(hhs_5_true$SEX_2), as.numeric(hhs_5_true$SEX_3),
           as.numeric(hhs_5_true$SEX_4))
f4 <- pmin(as.numeric(hhs_5_true$SEX_1), as.numeric(hhs_5_true$SEX_2), as.numeric(hhs_5_true$SEX_3),
           as.numeric(hhs_5_true$SEX_4), as.numeric(hhs_5_true$SEX_5))
f5 <- pmax(as.numeric(hhs_5_true$SEX_1), as.numeric(hhs_5_true$SEX_2), as.numeric(hhs_5_true$SEX_3),
           as.numeric(hhs_5_true$SEX_4), as.numeric(hhs_5_true$SEX_5))
hhs_5_true$SEX_diff <- droplevels(interaction(f1, f2, f3, f4, f5, sep = ""))

# Same for household size 6
hhs_6_true <- AMELIA_HH_wide[,c("HID", "HHS", "AGE_cat_1", "AGE_cat_2", "AGE_cat_3", "AGE_cat_4",
                                "AGE_cat_5", "AGE_cat_6",  "AGE_1", "AGE_2","AGE_3","AGE_4", "AGE_5", 
                                "AGE_6", "SEX_1", "SEX_2", 
                                "SEX_3","SEX_4", "SEX_5", "SEX_6",
                                "MST_1", "MST_2", "MST_3", "MST_4", "MST_5", "MST_6", 
                                "INC_1", "INC_2", "INC_3", "INC_4", "INC_5", "INC_6",
                                "PIDX_1", 
                                "PIDX_2", "PIDX_3", "PIDX_4", "PIDX_5", "PIDX_6")][AMELIA_HH_wide$HHS == 6]


hhs_6_true$AGE_mean <- (hhs_6_true$AGE_1 + hhs_6_true$AGE_2 + hhs_6_true$AGE_3 + hhs_6_true$AGE_4 +
                          hhs_6_true$AGE_5 + hhs_6_true$AGE_6) / 6
hhs_6_true$AGE_diff <- sqrt(((hhs_6_true$AGE_mean - hhs_6_true$AGE_1)^2 + 
                               (hhs_6_true$AGE_mean - hhs_6_true$AGE_2)^2 +
                               (hhs_6_true$AGE_mean - hhs_6_true$AGE_3)^2 +
                               (hhs_6_true$AGE_mean - hhs_6_true$AGE_4)^2 + 
                               (hhs_6_true$AGE_mean - hhs_6_true$AGE_5)^2 +
                               (hhs_6_true$AGE_mean - hhs_6_true$AGE_6)^2) / 6)
hhs_6_true$INC <- hhs_6_true$INC_1 + hhs_6_true$INC_2 + hhs_6_true$INC_3 + hhs_6_true$INC_4 + hhs_6_true$INC_5 +
  hhs_6_true$INC_6
f1 <- pmin(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2))
f2 <- pmin(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2), as.numeric(hhs_6_true$SEX_3))
f3 <- pmin(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2), as.numeric(hhs_6_true$SEX_3),
           as.numeric(hhs_6_true$SEX_4))
f4 <- pmin(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2), as.numeric(hhs_6_true$SEX_3),
           as.numeric(hhs_6_true$SEX_4), as.numeric(hhs_6_true$SEX_5))
f5 <- pmin(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2), as.numeric(hhs_6_true$SEX_3),
           as.numeric(hhs_6_true$SEX_4), as.numeric(hhs_6_true$SEX_5), as.numeric(hhs_6_true$SEX_6))
f6 <- pmax(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2), as.numeric(hhs_6_true$SEX_3),
           as.numeric(hhs_6_true$SEX_4), as.numeric(hhs_6_true$SEX_5), as.numeric(hhs_6_true$SEX_6))

f1 <- pmin(as.numeric(hhs_6_true$MST_1), as.numeric(hhs_6_true$MST_2))
f2 <- pmin(as.numeric(hhs_6_true$MST_1), as.numeric(hhs_6_true$MST_2), as.numeric(hhs_6_true$MST_3))
f3 <- pmin(as.numeric(hhs_6_true$MST_1), as.numeric(hhs_6_true$MST_2), as.numeric(hhs_6_true$MST_3),
           as.numeric(hhs_6_true$MST_4))
f4 <- pmin(as.numeric(hhs_6_true$MST_1), as.numeric(hhs_6_true$MST_2), as.numeric(hhs_6_true$MST_3),
           as.numeric(hhs_6_true$MST_4), as.numeric(hhs_6_true$MST_5))
f5 <- pmin(as.numeric(hhs_6_true$MST_1), as.numeric(hhs_6_true$MST_2), as.numeric(hhs_6_true$MST_3),
           as.numeric(hhs_6_true$MST_4), as.numeric(hhs_6_true$MST_5), as.numeric(hhs_6_true$MST_6))
f6 <- pmax(as.numeric(hhs_6_true$MST_1), as.numeric(hhs_6_true$MST_2), as.numeric(hhs_6_true$MST_3),
           as.numeric(hhs_6_true$MST_4), as.numeric(hhs_6_true$MST_5), as.numeric(hhs_6_true$MST_6))

hhs_6_true$MST_group <- droplevels(interaction(f1, f2, f3, f4, f5, f6, sep = ""))

f1 <- pmin(as.numeric(hhs_6_true$AGE_cat_1), as.numeric(hhs_6_true$AGE_cat_2))
f2 <- pmin(as.numeric(hhs_6_true$AGE_cat_1), as.numeric(hhs_6_true$AGE_cat_2), as.numeric(hhs_6_true$AGE_cat_3))
f3 <- pmin(as.numeric(hhs_6_true$AGE_cat_1), as.numeric(hhs_6_true$AGE_cat_2), as.numeric(hhs_6_true$AGE_cat_3),
           as.numeric(hhs_6_true$AGE_cat_4))
f4 <- pmin(as.numeric(hhs_6_true$AGE_cat_1), as.numeric(hhs_6_true$AGE_cat_2), as.numeric(hhs_6_true$AGE_cat_3),
           as.numeric(hhs_6_true$AGE_cat_4), as.numeric(hhs_6_true$AGE_cat_5))
f5 <- pmin(as.numeric(hhs_6_true$AGE_cat_1), as.numeric(hhs_6_true$AGE_cat_2), as.numeric(hhs_6_true$AGE_cat_3),
           as.numeric(hhs_6_true$AGE_cat_4), as.numeric(hhs_6_true$AGE_cat_5), as.numeric(hhs_6_true$AGE_cat_6))
f6 <- pmax(as.numeric(hhs_6_true$AGE_cat_1), as.numeric(hhs_6_true$AGE_cat_2), as.numeric(hhs_6_true$AGE_cat_3),
           as.numeric(hhs_6_true$AGE_cat_4), as.numeric(hhs_6_true$AGE_cat_5), as.numeric(hhs_6_true$AGE_cat_6))

hhs_6_true$AGE_comb <- droplevels(interaction(f1, f2, f3, f4, f5, f6, sep = ""))

f1 <- pmin(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2))
f2 <- pmin(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2), as.numeric(hhs_6_true$SEX_3))
f3 <- pmin(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2), as.numeric(hhs_6_true$SEX_3),
           as.numeric(hhs_6_true$SEX_4))
f4 <- pmin(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2), as.numeric(hhs_6_true$SEX_3),
           as.numeric(hhs_6_true$SEX_4), as.numeric(hhs_6_true$SEX_5))
f5 <- pmin(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2), as.numeric(hhs_6_true$SEX_3),
           as.numeric(hhs_6_true$SEX_4), as.numeric(hhs_6_true$SEX_5), as.numeric(hhs_6_true$SEX_6))
f6 <- pmax(as.numeric(hhs_6_true$SEX_1), as.numeric(hhs_6_true$SEX_2), as.numeric(hhs_6_true$SEX_3),
           as.numeric(hhs_6_true$SEX_4), as.numeric(hhs_6_true$SEX_5), as.numeric(hhs_6_true$SEX_6))

hhs_6_true$SEX_diff <- droplevels(interaction(f1, f2, f3, f4, f5, f6, sep = ""))

# The following function draws a random value for the variable of an individual, based on
# its distribution function observed in the sample (see Chapter 4, 'Combination approach')
create_synth <- function(hhs, model, x, cluster_sample) {
  repeat{
    # Household size 2
    if(hhs == 2) {
      probs_AGE_1 <- density(cluster_sample$AGE_1)
      probs_AGE_2 <- density(cluster_sample$AGE_2)
      probs_SEX_1 <- as.vector(prop.table(table((cluster_sample$SEX_1))))
      probs_SEX_2 <- as.vector(prop.table(table((cluster_sample$SEX_2))))
      probs_MST_1 <- as.vector(prop.table(table((cluster_sample$MST_1))))
      probs_MST_2 <- as.vector(prop.table(table((cluster_sample$MST_2))))
      
      AGE_1 <- sample(probs_AGE_1$x, 1, prob = probs_AGE_1$y)
      AGE_2 <- sample(probs_AGE_2$x, 1, prob = probs_AGE_2$y)
      # SEX_1 1/2 SEX_2 1/2
      SEX_1 <- sample(1:2, 1, prob = probs_SEX_1)
      SEX_2 <- sample(1:2, 1, prob = probs_SEX_2)
      # MST_1 1-5
      MST_1 <- sample(1:5, 1, prob = probs_MST_1)
      MST_2 <- sample(1:5, 1, prob = probs_MST_2)
      ###
      choice <- as.data.frame(cbind(AGE_1, AGE_2, SEX_1, SEX_2, MST_1, MST_2))
      choice$SEX_1 <- factor(choice$SEX_1, levels = levels(hhs_2_true$SEX_1))
      choice$SEX_2 <- factor(choice$SEX_2, levels = levels(hhs_2_true$SEX_2))
      choice$MST_1 <- factor(choice$MST_1, levels = levels(hhs_2_true$MST_1))
      choice$MST_2 <- factor(choice$MST_2, levels = levels(hhs_2_true$MST_2))
    }
    if(hhs == 3) { 
      probs_AGE_1 <- density(cluster_sample$AGE_1)
      probs_AGE_2 <- density(cluster_sample$AGE_2)
      probs_AGE_3 <- density(cluster_sample$AGE_3)

      probs_SEX_1 <- as.vector(prop.table(table((cluster_sample$SEX_1))))
      probs_SEX_2 <- as.vector(prop.table(table((cluster_sample$SEX_2))))
      probs_SEX_3 <- as.vector(prop.table(table((cluster_sample$SEX_3))))
      
      probs_MST_1 <- as.vector(prop.table(table((cluster_sample$MST_1))))
      probs_MST_2 <- as.vector(prop.table(table((cluster_sample$MST_2))))
      probs_MST_3 <- as.vector(prop.table(table((cluster_sample$MST_3))))

      AGE_1 <- sample(probs_AGE_1$x, 1, prob = probs_AGE_1$y)
      AGE_2 <- sample(probs_AGE_2$x, 1, prob = probs_AGE_2$y)
      AGE_3 <- sample(probs_AGE_3$x, 1, prob = probs_AGE_3$y)
      # SEX_1 1/2 SEX_2 1/2
      SEX_1 <- sample(1:2, 1, prob = probs_SEX_1)
      SEX_2 <- sample(1:2, 1, prob = probs_SEX_2)
      SEX_3 <- sample(1:2, 1, prob = probs_SEX_3)
      # MST_1 1-5
      MST_1 <- sample(1:5, 1, prob = probs_MST_1)
      MST_2 <- sample(1:5, 1, prob = probs_MST_2)
      MST_3 <- sample(1:5, 1, prob = probs_MST_3)
      ###
      choice <- as.data.frame(cbind(AGE_1, AGE_2, AGE_3, SEX_1, SEX_2, SEX_3, MST_1, MST_2, MST_3))
      choice$SEX_1 <- factor(choice$SEX_1, levels = levels(hhs_3_true$SEX_1))
      choice$SEX_2 <- factor(choice$SEX_2, levels = levels(hhs_3_true$SEX_2))
      choice$SEX_3 <- factor(choice$SEX_3, levels = levels(hhs_3_true$SEX_3))
      choice$MST_1 <- factor(choice$MST_1, levels = levels(hhs_3_true$MST_1))
      choice$MST_2 <- factor(choice$MST_2, levels = levels(hhs_3_true$MST_2))
      choice$MST_3 <- factor(choice$MST_3, levels = levels(hhs_3_true$MST_3))
    }
    if(hhs == 4) {    
      probs_AGE_1 <- density(cluster_sample$AGE_1)
      probs_AGE_2 <- density(cluster_sample$AGE_2)
      probs_AGE_3 <- density(cluster_sample$AGE_3)
      probs_AGE_4 <- density(cluster_sample$AGE_4)

      probs_SEX_1 <- as.vector(prop.table(table((cluster_sample$SEX_1))))
      probs_SEX_2 <- as.vector(prop.table(table((cluster_sample$SEX_2))))
      probs_SEX_3 <- as.vector(prop.table(table((cluster_sample$SEX_3))))
      probs_SEX_4 <- as.vector(prop.table(table((cluster_sample$SEX_4))))

      probs_MST_1 <- as.vector(prop.table(table((cluster_sample$MST_1))))
      probs_MST_2 <- as.vector(prop.table(table((cluster_sample$MST_2))))
      probs_MST_3 <- as.vector(prop.table(table((cluster_sample$MST_3))))
      probs_MST_4 <- as.vector(prop.table(table((cluster_sample$MST_4))))

      AGE_1 <- sample(probs_AGE_1$x, 1, prob = probs_AGE_1$y)
      AGE_2 <- sample(probs_AGE_2$x, 1, prob = probs_AGE_2$y)
      AGE_3 <- sample(probs_AGE_3$x, 1, prob = probs_AGE_3$y)
      AGE_4 <- sample(probs_AGE_4$x, 1, prob = probs_AGE_4$y)
      # SEX_1 1/2 SEX_2 1/2
      SEX_1 <- sample(1:2, 1, prob = probs_SEX_1)
      SEX_2 <- sample(1:2, 1, prob = probs_SEX_2)
      SEX_3 <- sample(1:2, 1, prob = probs_SEX_3)
      SEX_4 <- sample(1:2, 1, prob = probs_SEX_4)
      # MST_1 1-5
      MST_1 <- sample(1:5, 1, prob = probs_MST_1)
      MST_2 <- sample(1:5, 1, prob = probs_MST_2)
      MST_3 <- sample(1:5, 1, prob = probs_MST_3)
      MST_4 <- sample(1:5, 1, prob = probs_MST_4)
      ###
      choice <- as.data.frame(cbind(AGE_1, AGE_2, AGE_3, AGE_4, SEX_1, SEX_2, SEX_3, SEX_4, 
                                    MST_1, MST_2, MST_3, MST_4))
      choice$SEX_1 <- factor(choice$SEX_1, levels = levels(hhs_4_true$SEX_1))
      choice$SEX_2 <- factor(choice$SEX_2, levels = levels(hhs_4_true$SEX_2))
      choice$SEX_3 <- factor(choice$SEX_3, levels = levels(hhs_4_true$SEX_3))
      choice$SEX_4 <- factor(choice$SEX_4, levels = levels(hhs_4_true$SEX_4))
      choice$MST_1 <- factor(choice$MST_1, levels = levels(hhs_4_true$MST_1))
      choice$MST_2 <- factor(choice$MST_2, levels = levels(hhs_4_true$MST_2))
      choice$MST_3 <- factor(choice$MST_3, levels = levels(hhs_4_true$MST_3))
      choice$MST_4 <- factor(choice$MST_4, levels = levels(hhs_4_true$MST_4))
    }
    if(hhs == 5) {   
      probs_AGE_1 <- density(cluster_sample$AGE_1)
      probs_AGE_2 <- density(cluster_sample$AGE_2)
      probs_AGE_3 <- density(cluster_sample$AGE_3)
      probs_AGE_4 <- density(cluster_sample$AGE_4)
      probs_AGE_5 <- density(cluster_sample$AGE_5)

      probs_SEX_1 <- as.vector(prop.table(table((cluster_sample$SEX_1))))
      probs_SEX_2 <- as.vector(prop.table(table((cluster_sample$SEX_2))))
      probs_SEX_3 <- as.vector(prop.table(table((cluster_sample$SEX_3))))
      probs_SEX_4 <- as.vector(prop.table(table((cluster_sample$SEX_4))))
      probs_SEX_5 <- as.vector(prop.table(table((cluster_sample$SEX_5))))
      
      probs_MST_1 <- as.vector(prop.table(table((cluster_sample$MST_1))))
      probs_MST_2 <- as.vector(prop.table(table((cluster_sample$MST_2))))
      probs_MST_3 <- as.vector(prop.table(table((cluster_sample$MST_3))))
      probs_MST_4 <- as.vector(prop.table(table((cluster_sample$MST_4))))
      probs_MST_5 <- as.vector(prop.table(table((cluster_sample$MST_5))))

      AGE_1 <- sample(probs_AGE_1$x, 1, prob = probs_AGE_1$y)
      AGE_2 <- sample(probs_AGE_2$x, 1, prob = probs_AGE_2$y)
      AGE_3 <- sample(probs_AGE_3$x, 1, prob = probs_AGE_3$y)
      AGE_4 <- sample(probs_AGE_4$x, 1, prob = probs_AGE_4$y)
      AGE_5 <- sample(probs_AGE_5$x, 1, prob = probs_AGE_5$y)
      # SEX_1 1/2 SEX_2 1/2
      SEX_1 <- sample(1:2, 1, prob = probs_SEX_1)
      SEX_2 <- sample(1:2, 1, prob = probs_SEX_2)
      SEX_3 <- sample(1:2, 1, prob = probs_SEX_3)
      SEX_4 <- sample(1:2, 1, prob = probs_SEX_4)
      SEX_5 <- sample(1:2, 1, prob = probs_SEX_5)
      # MST_1 1-5
      MST_1 <- sample(1:5, 1, prob = probs_MST_1)
      MST_2 <- sample(1:5, 1, prob = probs_MST_2)
      MST_3 <- sample(1:5, 1, prob = probs_MST_3)
      MST_4 <- sample(1:5, 1, prob = probs_MST_4)
      MST_5 <- sample(1:5, 1, prob = probs_MST_5)
      ###
      choice <- as.data.frame(cbind(AGE_1, AGE_2, AGE_3, AGE_4, AGE_5, SEX_1, SEX_2, SEX_3, SEX_4, SEX_5, 
                                    MST_1, MST_2, MST_3, MST_4, MST_5))
      choice$SEX_1 <- factor(choice$SEX_1, levels = levels(hhs_5_true$SEX_1))
      choice$SEX_2 <- factor(choice$SEX_2, levels = levels(hhs_5_true$SEX_2))
      choice$SEX_3 <- factor(choice$SEX_3, levels = levels(hhs_5_true$SEX_3))
      choice$SEX_4 <- factor(choice$SEX_4, levels = levels(hhs_5_true$SEX_4))
      choice$SEX_5 <- factor(choice$SEX_5, levels = levels(hhs_5_true$SEX_5))
      choice$MST_1 <- factor(choice$MST_1, levels = levels(hhs_5_true$MST_1))
      choice$MST_2 <- factor(choice$MST_2, levels = levels(hhs_5_true$MST_2))
      choice$MST_3 <- factor(choice$MST_3, levels = levels(hhs_5_true$MST_3))
      choice$MST_4 <- factor(choice$MST_4, levels = levels(hhs_5_true$MST_4))
      choice$MST_5 <- factor(choice$MST_5, levels = levels(hhs_5_true$MST_5))
    }      
    if(hhs == 6) {   
      probs_AGE_1 <- density(cluster_sample$AGE_1)
      probs_AGE_2 <- density(cluster_sample$AGE_2)
      probs_AGE_3 <- density(cluster_sample$AGE_3)
      probs_AGE_4 <- density(cluster_sample$AGE_4)
      probs_AGE_5 <- density(cluster_sample$AGE_5)
      probs_AGE_6 <- density(cluster_sample$AGE_6)
      
      probs_SEX_1 <- as.vector(prop.table(table((cluster_sample$SEX_1))))
      probs_SEX_2 <- as.vector(prop.table(table((cluster_sample$SEX_2))))
      probs_SEX_3 <- as.vector(prop.table(table((cluster_sample$SEX_3))))
      probs_SEX_4 <- as.vector(prop.table(table((cluster_sample$SEX_4))))
      probs_SEX_5 <- as.vector(prop.table(table((cluster_sample$SEX_5))))
      probs_SEX_6 <- as.vector(prop.table(table((cluster_sample$SEX_6))))
      
      probs_MST_1 <- as.vector(prop.table(table((cluster_sample$MST_1))))
      probs_MST_2 <- as.vector(prop.table(table((cluster_sample$MST_2))))
      probs_MST_3 <- as.vector(prop.table(table((cluster_sample$MST_3))))
      probs_MST_4 <- as.vector(prop.table(table((cluster_sample$MST_4))))
      probs_MST_5 <- as.vector(prop.table(table((cluster_sample$MST_5))))
      probs_MST_6 <- as.vector(prop.table(table((cluster_sample$MST_6))))
      AGE_1 <- sample(probs_AGE_1$x, 1, prob = probs_AGE_1$y)
      AGE_2 <- sample(probs_AGE_2$x, 1, prob = probs_AGE_2$y)
      AGE_3 <- sample(probs_AGE_3$x, 1, prob = probs_AGE_3$y)
      AGE_4 <- sample(probs_AGE_4$x, 1, prob = probs_AGE_4$y)
      AGE_5 <- sample(probs_AGE_5$x, 1, prob = probs_AGE_5$y)
      AGE_6 <- sample(probs_AGE_6$x, 1, prob = probs_AGE_6$y)
      # SEX_1 1/2 SEX_2 1/2
      SEX_1 <- sample(1:2, 1, prob = probs_SEX_1)
      SEX_2 <- sample(1:2, 1, prob = probs_SEX_2)
      SEX_3 <- sample(1:2, 1, prob = probs_SEX_3)
      SEX_4 <- sample(1:2, 1, prob = probs_SEX_4)
      SEX_5 <- sample(1:2, 1, prob = probs_SEX_5)
      SEX_6 <- sample(1:2, 1, prob = probs_SEX_6)
      # MST_1 1-5
      MST_1 <- sample(1:5, 1, prob = probs_MST_1)
      MST_2 <- sample(1:5, 1, prob = probs_MST_2)
      MST_3 <- sample(1:5, 1, prob = probs_MST_3)
      MST_4 <- sample(1:5, 1, prob = probs_MST_4)
      MST_5 <- sample(1:5, 1, prob = probs_MST_5)
      MST_6 <- sample(1:5, 1, prob = probs_MST_6)
      ###
      choice <- as.data.frame(cbind(AGE_1, AGE_2, AGE_3, AGE_4, AGE_5, AGE_6, SEX_1, SEX_2, SEX_3,
                                    SEX_4, SEX_5, SEX_6,
                                    MST_1, MST_2, MST_3, MST_4, MST_5, MST_6))
      choice$SEX_1 <- factor(choice$SEX_1, levels = levels(hhs_6_true$SEX_1))
      choice$SEX_2 <- factor(choice$SEX_2, levels = levels(hhs_6_true$SEX_2))
      choice$SEX_3 <- factor(choice$SEX_3, levels = levels(hhs_6_true$SEX_3))
      choice$SEX_4 <- factor(choice$SEX_4, levels = levels(hhs_6_true$SEX_4))
      choice$SEX_5 <- factor(choice$SEX_5, levels = levels(hhs_6_true$SEX_5))
      choice$SEX_6 <- factor(choice$SEX_6, levels = levels(hhs_6_true$SEX_6))
      choice$MST_1 <- factor(choice$MST_1, levels = levels(hhs_6_true$MST_1))
      choice$MST_2 <- factor(choice$MST_2, levels = levels(hhs_6_true$MST_2))
      choice$MST_3 <- factor(choice$MST_3, levels = levels(hhs_6_true$MST_3))
      choice$MST_4 <- factor(choice$MST_4, levels = levels(hhs_6_true$MST_4))
      choice$MST_5 <- factor(choice$MST_5, levels = levels(hhs_6_true$MST_5))
      choice$MST_6 <- factor(choice$MST_6, levels = levels(hhs_6_true$MST_6))
    }
    
    if(model == "NN") {      # Neural net   
      pr_synthetic_hhs <- predict(x, newdata = choice, type = "raw")
    }
    if(model == "RF") {     # Random Forest   
      pr_synthetic_hhs <- predict(x, newdata = choice, type = "prob")
    }
    if(model == "Logit") {  # Logit Model      
      pr_synthetic_hhs <- predict(x, newdata = choice, type = "response")
    }  
    r <- runif(1) # Acceptance-rejecting sampling
    if(r < pr_synthetic_hhs) {
      synth_nn <- choice
      synth_nn$prob <- pr_synthetic_hhs
      break
    }} 
  return(synth_nn)
}

# Create a function to simulate the creation of synthetic households with the combination
# approach. The function allows to vary the sample size of the indivudals, the number of 
# random combinations set against the
# true combinations, the size of the synthetic population (e.g the number of households in a certain class)
#
Monte_Carlo_Simulation <- function(n_model_comb, household_size,
                                   True_hhs, AMELIA_HHS, hhs_true, NP) {
  out <- list()  
  
   ### Our Single-Stage-Cluster-Sample
  # Convert wide format back to individual format; this enables us to recombine the individuals
  Sample_x <-  AMELIA_wide_sample[which(AMELIA_wide_sample$HHS == household_size),]
  full_synthetic_simulation <- AMELIA[HID %in% Sample_x$HID]
  
  datalist = list()
  for(i in 1:n_model_comb) { # Number of random combinations should bebetween 30 and 100
                             # Sample = 1000 # Model Comb = 100 // == 1000*100 = 1000000 "0 coded
                             # random combinations 1000 "true" 1 codes households. Size 101000
    dat <- ind_to_random_hid(full_synthetic_simulation, household_size = household_size,
                             synthetic_hhs = full_synthetic_simulation$HHS,
                             drop = FALSE)
    datalist[[i]] <- dat
    print(i)
  }
  
  ra <- do.call(rbind, datalist)
  ra$TRUE_ <- 0 # Code random households with "0
  Sample_x$TRUE_ <- 1 # Code "true" (obserseved) households from sample with "1"
  model_data_cx <- rbind(Sample_x, ra, fill = TRUE) # Bind them together
  model_data_cx$TRUE_ <- as.factor(model_data_cx$TRUE_)
  model_data_cx$SEX_diff <- as.factor(model_data_cx$SEX_diff)
  # Apply models based on household size
  if(household_size == 2) {
    hhs_nn_cx <- nnet(TRUE_ ~ AGE_1 + AGE_2 + SEX_1 + SEX_2 + MST_1 + MST_2, data = model_data_cx, size = 4,
                      maxit = 500)
    hhs_rf_cx <- randomForest(TRUE_ ~ AGE_1 + AGE_2 + SEX_1 + SEX_2 + MST_1 + MST_2, data = model_data_cx, type = "prob")
    hhs_logit_cx <- glm(TRUE_ ~ AGE_1 + AGE_2 + SEX_1 + SEX_2 + MST_1 + MST_2, data = model_data_cx, family = "binomial")
  }
  # As we can see, bigger household size leads to more members, leads to more covariates
  if(household_size == 3) {
    hhs_nn_cx <- nnet(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + SEX_1 + SEX_2 + SEX_3 + MST_1 + MST_2 + MST_3, data = model_data_cx, 
                      size = 4, maxit = 500)
    hhs_rf_cx <- randomForest(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + SEX_1 + SEX_2 + SEX_3 + MST_1 + MST_2 + MST_3, data = model_data_cx, type = "prob")
    hhs_logit_cx <- glm(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + SEX_1 + SEX_2 + SEX_3 + MST_1 + MST_2 + MST_3, data = model_data_cx, family = "binomial")
  }
  
  if(household_size == 4) {
    hhs_nn_cx <- nnet(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + MST_1 + MST_2 + MST_3 + MST_4, data = model_data_cx, 
                      size = 4, maxit = 500)
    hhs_rf_cx <- randomForest(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + MST_1 + MST_2 + MST_3 + MST_4, data = model_data_cx, type = "prob")
    hhs_logit_cx <- glm(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + MST_1 + MST_2 + MST_3 + MST_4, data = model_data_cx, family = "binomial")
  }
  
  if(household_size == 5) {
    hhs_nn_cx <- nnet(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5, data = model_data_cx, size = 4,
                      maxit = 500)
    hhs_rf_cx <- randomForest(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5, data = model_data_cx, type = "prob")
    hhs_logit_cx <- glm(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5, data = model_data_cx, family = "binomial")
  }
  
  if(household_size == 6) {
    hhs_nn_cx <- nnet(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + AGE_6 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + SEX_6 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5 + MST_6, data = model_data_cx, size = 4,
                      maxit = 500)
    hhs_rf_cx <- randomForest(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + AGE_6 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + SEX_6 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5 + MST_6, data = model_data_cx, type = "prob")
    hhs_logit_cx <- glm(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + AGE_6 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + SEX_6 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5 + MST_6, data = model_data_cx, family = "binomial")
  }
  
  result <- list() # Store every accepted household in a list
  for(i in 1:NP) {
    resultX <- create_synth(hhs = household_size, model = "NN", x = hhs_nn_cx, 
                            cluster_sample = Sample_x)
    result[[length(result)+1]] <- resultX
    print(i)
  }
  synthetic_hhs_nn <- result # Save synthetic population created by NN
  result <- list()
  for(i in 1:NP) {
    resultX <- create_synth(hhs = household_size, model = "RF", x = hhs_rf_cx, 
                            cluster_sample = Sample_x)
    result[[length(result)+1]] <- resultX
    print(i)
  }
  synthetic_hhs_rf <- result
  result <- list()
  for(i in 1:NP) {
    resultX <- create_synth(hhs = household_size, model = "Logit", 
                            x = hhs_logit_cx, 
                            cluster_sample = Sample_x)
    result[[length(result)+1]] <- resultX
    print(i)
  }  
  synthetic_hhs_logit <- result # Save synthetic population created by RF
  out <- list(synthetic_hhs_logit, synthetic_hhs_nn, synthetic_hhs_rf)
  return(out)
}
# Create a function to subsequently create age categories
AGE_cat_creator <- function(x) {
x$AGE_cat_1 <- "83+x"
x$AGE_cat_1[which(x$AGE_1 < 18)] <- "0-17"
x$AGE_cat_1[which(x$AGE_1 > 17 & x$AGE_1 <= 28)] <- "18-28"
x$AGE_cat_1[which(x$AGE_1 > 28 & x$AGE_1 <= 39)] <- "29-39"
x$AGE_cat_1[which(x$AGE_1 > 40 & x$AGE_1 <= 50)] <- "40-50"
x$AGE_cat_1[which(x$AGE_1 > 50 & x$AGE_1 <= 61)] <- "51-61"
x$AGE_cat_1[which(x$AGE_1 > 61 & x$AGE_1 <= 72)] <- "62-72"
x$AGE_cat_1[which(x$AGE_1 > 72 & x$AGE_1 <= 83)] <- "73-83"

x$AGE_cat_2 <- "83+x"
x$AGE_cat_2[which(x$AGE_2 < 18)] <- "0-17"
x$AGE_cat_2[which(x$AGE_2 > 17 & x$AGE_2 <= 28)] <- "18-28"
x$AGE_cat_2[which(x$AGE_2 > 28 & x$AGE_2 <= 39)] <- "29-39"
x$AGE_cat_2[which(x$AGE_2 > 40 & x$AGE_2 <= 50)] <- "40-50"
x$AGE_cat_2[which(x$AGE_2 > 50 & x$AGE_2 <= 61)] <- "51-61"
x$AGE_cat_2[which(x$AGE_2 > 61 & x$AGE_2 <= 72)] <- "62-72"
x$AGE_cat_2[which(x$AGE_2 > 72 & x$AGE_2 <= 83)] <- "73-83"

x$AGE_cat_3 <- "83+x"
x$AGE_cat_3[which(x$AGE_3 < 18)] <- "0-17"
x$AGE_cat_3[which(x$AGE_3 > 17 & x$AGE_3 <= 28)] <- "18-28"
x$AGE_cat_3[which(x$AGE_3 > 28 & x$AGE_3 <= 39)] <- "29-39"
x$AGE_cat_3[which(x$AGE_3 > 40 & x$AGE_3 <= 50)] <- "40-50"
x$AGE_cat_3[which(x$AGE_3 > 50 & x$AGE_3 <= 61)] <- "51-61"
x$AGE_cat_3[which(x$AGE_3 > 61 & x$AGE_3 <= 72)] <- "62-72"
x$AGE_cat_3[which(x$AGE_3 > 72 & x$AGE_3 <= 83)] <- "73-83"

x$AGE_cat_4 <- "83+x"
x$AGE_cat_4[which(x$AGE_4 < 18)] <- "0-17"
x$AGE_cat_4[which(x$AGE_4 > 17 & x$AGE_4 <= 28)] <- "18-28"
x$AGE_cat_4[which(x$AGE_4 > 28 & x$AGE_4 <= 39)] <- "29-39"
x$AGE_cat_4[which(x$AGE_4 > 40 & x$AGE_4 <= 50)] <- "40-50"
x$AGE_cat_4[which(x$AGE_4 > 50 & x$AGE_4 <= 61)] <- "51-61"
x$AGE_cat_4[which(x$AGE_4 > 61 & x$AGE_4 <= 72)] <- "62-72"
x$AGE_cat_4[which(x$AGE_4 > 72 & x$AGE_4 <= 83)] <- "73-83"

x$AGE_cat_5 <- "83+x"
x$AGE_cat_5[which(x$AGE_5 < 18)] <- "0-17"
x$AGE_cat_5[which(x$AGE_5 > 17 & x$AGE_5 <= 28)] <- "18-28"
x$AGE_cat_5[which(x$AGE_5 > 28 & x$AGE_5 <= 39)] <- "29-39"
x$AGE_cat_5[which(x$AGE_5 > 40 & x$AGE_5 <= 50)] <- "40-50"
x$AGE_cat_5[which(x$AGE_5 > 50 & x$AGE_5 <= 61)] <- "51-61"
x$AGE_cat_5[which(x$AGE_5 > 61 & x$AGE_5 <= 72)] <- "62-72"
x$AGE_cat_5[which(x$AGE_5 > 72 & x$AGE_5 <= 83)] <- "73-83"

x$AGE_cat_6 <- "83+x"
x$AGE_cat_6[which(x$AGE_6 < 18)] <- "0-17"
x$AGE_cat_6[which(x$AGE_6 > 17 & x$AGE_6 <= 28)] <- "18-28"
x$AGE_cat_6[which(x$AGE_6 > 28 & x$AGE_6 <= 39)] <- "29-39"
x$AGE_cat_6[which(x$AGE_6 > 40 & x$AGE_6 <= 50)] <- "40-50"
x$AGE_cat_6[which(x$AGE_6 > 50 & x$AGE_6 <= 61)] <- "51-61"
x$AGE_cat_6[which(x$AGE_6 > 61 & x$AGE_6 <= 72)] <- "62-72"
x$AGE_cat_6[which(x$AGE_6 > 72 & x$AGE_6 <= 83)] <- "73-83"
return(x)
}
# Create a function to obtain variables which give the "distance" between household
# members, like mean of age, difference of age or combinations of categorical variables
# (Two men living in a Household == Sex_Comb = "11") # See Figure 4 in the thesis
calc_difference <- function(x, hhs_size) {
  AMELIA_wide_start_c <- x 
  if(hhs_size == 2) { 
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2))
    f2 <- pmax(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2))
    AMELIA_wide_start_c$MST_group <- droplevels(interaction(f1, f2, sep=""))
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2))
    f2 <- pmax(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2))
    AMELIA_wide_start_c$AGE_comb <- droplevels(interaction(f1, f2, sep=""))
    AMELIA_wide_start_c$AGE_mean <- (AMELIA_wide_start_c$AGE_1 + AMELIA_wide_start_c$AGE_2) / 2
    AMELIA_wide_start_c$AGE_diff <- sqrt(((AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_1)^2 + (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_2)^2) /2)
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2))
    f2 <- pmax(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2))
    AMELIA_wide_start_c$SEX_diff <- droplevels(interaction(f1, f2, sep=""))
  } else if(hhs_size == 3) {
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), as.numeric(AMELIA_wide_start_c$MST_3))
    f3 <- pmax(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), as.numeric(AMELIA_wide_start_c$MST_3))
    
    AMELIA_wide_start_c$MST_group <- droplevels(interaction(f1, f2, f3, sep = ""))
    
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), as.numeric(AMELIA_wide_start_c$AGE_cat_3))
    f3 <- pmax(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), as.numeric(AMELIA_wide_start_c$AGE_cat_3))
    
    AMELIA_wide_start_c$AGE_comb <- droplevels(interaction(f1, f2, f3, sep = ""))
    AMELIA_wide_start_c$AGE_mean <- (AMELIA_wide_start_c$AGE_1 + AMELIA_wide_start_c$AGE_2 +
                                       AMELIA_wide_start_c$AGE_3) / 3
    AMELIA_wide_start_c$AGE_diff <- sqrt(((AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_1)^2 + 
                                            (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_2)^2 +
                                           (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_3)^2)/3)
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), as.numeric(AMELIA_wide_start_c$SEX_3))
    f3 <- pmax(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), as.numeric(AMELIA_wide_start_c$SEX_3))
    AMELIA_wide_start_c$SEX_diff <- droplevels(interaction(f1, f2, f3, sep = ""))
  } else if(hhs_size == 4) {
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), as.numeric(AMELIA_wide_start_c$MST_3))
    f3 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), as.numeric(AMELIA_wide_start_c$MST_3),
               as.numeric(AMELIA_wide_start_c$MST_4))
    f4 <- pmax(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), as.numeric(AMELIA_wide_start_c$MST_3),
               as.numeric(AMELIA_wide_start_c$MST_4))
    
    AMELIA_wide_start_c$MST_group <- droplevels(interaction(f1, f2, f3, f4, sep = ""))
    
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), as.numeric(AMELIA_wide_start_c$AGE_cat_3))
    f3 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), as.numeric(AMELIA_wide_start_c$AGE_cat_3),
               as.numeric(AMELIA_wide_start_c$AGE_cat_4))
    f4 <- pmax(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), as.numeric(AMELIA_wide_start_c$AGE_cat_3),
               as.numeric(AMELIA_wide_start_c$AGE_cat_4))
    
    AMELIA_wide_start_c$AGE_comb <- droplevels(interaction(f1, f2, f3, f4, sep = ""))
    AMELIA_wide_start_c$AGE_mean <- (AMELIA_wide_start_c$AGE_1 + AMELIA_wide_start_c$AGE_2 +
                                       AMELIA_wide_start_c$AGE_3 + AMELIA_wide_start_c$AGE_4) / 4
    AMELIA_wide_start_c$AGE_diff <- sqrt(((AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_1)^2 + 
                                            (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_2)^2 +
                                           (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_3)^2 +
                                          (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_4)^2) /4)
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), as.numeric(AMELIA_wide_start_c$SEX_3))
    f3 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), as.numeric(AMELIA_wide_start_c$SEX_3),
               as.numeric(AMELIA_wide_start_c$SEX_4))
    f4 <- pmax(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), as.numeric(AMELIA_wide_start_c$SEX_3),
               as.numeric(AMELIA_wide_start_c$SEX_4))
    AMELIA_wide_start_c$SEX_diff <- droplevels(interaction(f1, f2, f3, f4, sep = ""))
  } else if(hhs_size == 5) {
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), as.numeric(AMELIA_wide_start_c$MST_3))
    f3 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), as.numeric(AMELIA_wide_start_c$MST_3),
               as.numeric(AMELIA_wide_start_c$MST_4))
    f4 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), as.numeric(AMELIA_wide_start_c$MST_3),
               as.numeric(AMELIA_wide_start_c$MST_4), as.numeric(AMELIA_wide_start_c$MST_5))
    f5 <- pmax(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), as.numeric(AMELIA_wide_start_c$MST_3),
               as.numeric(AMELIA_wide_start_c$MST_4), as.numeric(AMELIA_wide_start_c$MST_5))
    
    AMELIA_wide_start_c$MST_group <- droplevels(interaction(f1, f2, f3, f4, f5, sep = ""))
    
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), as.numeric(AMELIA_wide_start_c$AGE_cat_3))
    f3 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), as.numeric(AMELIA_wide_start_c$AGE_cat_3),
               as.numeric(AMELIA_wide_start_c$AGE_cat_4))
    f4 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), as.numeric(AMELIA_wide_start_c$AGE_cat_3),
               as.numeric(AMELIA_wide_start_c$AGE_cat_4), as.numeric(AMELIA_wide_start_c$AGE_cat_5))
    f5 <- pmax(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), as.numeric(AMELIA_wide_start_c$AGE_cat_3),
               as.numeric(AMELIA_wide_start_c$AGE_cat_4), as.numeric(AMELIA_wide_start_c$AGE_cat_5))
    
    AMELIA_wide_start_c$AGE_comb <- droplevels(interaction(f1, f2, f3, f4, f5, sep = ""))
    AMELIA_wide_start_c$AGE_mean <- (AMELIA_wide_start_c$AGE_1 + AMELIA_wide_start_c$AGE_2 +
                                       AMELIA_wide_start_c$AGE_3 + AMELIA_wide_start_c$AGE_4 +
                                       AMELIA_wide_start_c$AGE_5) / 5
    AMELIA_wide_start_c$AGE_diff <- sqrt(((AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_1)^2 + 
                                            (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_2)^2 +
                                            (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_3)^2 +
                                            (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_4)^2 +
                                            (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_5)^2) /5)
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), as.numeric(AMELIA_wide_start_c$SEX_3))
    f3 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), as.numeric(AMELIA_wide_start_c$SEX_3),
               as.numeric(AMELIA_wide_start_c$SEX_4))
    f4 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), as.numeric(AMELIA_wide_start_c$SEX_3),
               as.numeric(AMELIA_wide_start_c$SEX_4), as.numeric(AMELIA_wide_start_c$SEX_5))
    f5 <- pmax(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), as.numeric(AMELIA_wide_start_c$SEX_3),
               as.numeric(AMELIA_wide_start_c$SEX_4), as.numeric(AMELIA_wide_start_c$SEX_5))
    AMELIA_wide_start_c$SEX_diff <- droplevels(interaction(f1, f2, f3, f4, f5, sep = ""))
  } else if(hhs_size == 6) {
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), 
               as.numeric(AMELIA_wide_start_c$MST_3))
    f3 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), 
               as.numeric(AMELIA_wide_start_c$MST_3),
               as.numeric(AMELIA_wide_start_c$MST_4))
    f4 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2),
               as.numeric(AMELIA_wide_start_c$MST_3),
               as.numeric(AMELIA_wide_start_c$MST_4), as.numeric(AMELIA_wide_start_c$MST_5))
    f5 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), 
               as.numeric(AMELIA_wide_start_c$MST_3),
               as.numeric(AMELIA_wide_start_c$MST_4), as.numeric(AMELIA_wide_start_c$MST_5),
               as.numeric(AMELIA_wide_start_c$MST_6))
    f6 <- pmax(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), 
               as.numeric(AMELIA_wide_start_c$MST_3),
               as.numeric(AMELIA_wide_start_c$MST_4), as.numeric(AMELIA_wide_start_c$MST_5), 
               as.numeric(AMELIA_wide_start_c$MST_6))
    
    AMELIA_wide_start_c$MST_group <- droplevels(interaction(f1, f2, f3, f4, f5, f6, sep = ""))
    
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), 
               as.numeric(AMELIA_wide_start_c$AGE_cat_3))
    f3 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), 
               as.numeric(AMELIA_wide_start_c$AGE_cat_3),
               as.numeric(AMELIA_wide_start_c$AGE_cat_4))
    f4 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2),
               as.numeric(AMELIA_wide_start_c$AGE_cat_3),
               as.numeric(AMELIA_wide_start_c$AGE_cat_4), as.numeric(AMELIA_wide_start_c$AGE_cat_5))
    f5 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), 
               as.numeric(AMELIA_wide_start_c$AGE_cat_3),
               as.numeric(AMELIA_wide_start_c$AGE_cat_4), as.numeric(AMELIA_wide_start_c$AGE_cat_5),
               as.numeric(AMELIA_wide_start_c$AGE_cat_6))
    f6 <- pmax(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), 
               as.numeric(AMELIA_wide_start_c$AGE_cat_3),
               as.numeric(AMELIA_wide_start_c$AGE_cat_4), as.numeric(AMELIA_wide_start_c$AGE_cat_5), 
               as.numeric(AMELIA_wide_start_c$AGE_cat_6))
    
    AMELIA_wide_start_c$AGE_comb <- droplevels(interaction(f1, f2, f3, f4, f5, f6, sep = ""))
    AMELIA_wide_start_c$AGE_mean <- (AMELIA_wide_start_c$AGE_1 + AMELIA_wide_start_c$AGE_2 +
                                       AMELIA_wide_start_c$AGE_3 + AMELIA_wide_start_c$AGE_4 +
                                       AMELIA_wide_start_c$AGE_5 + AMELIA_wide_start_c$AGE_6) / 6
    AMELIA_wide_start_c$AGE_diff <- sqrt(((AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_1)^2 + 
                                            (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_2)^2 +
                                            (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_3)^2 +
                                            (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_4)^2 +
                                            (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_5)^2 +
                                            (AMELIA_wide_start_c$AGE_mean - AMELIA_wide_start_c$AGE_6)^2) /5)
    
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), 
               as.numeric(AMELIA_wide_start_c$SEX_3))
    f3 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), 
               as.numeric(AMELIA_wide_start_c$SEX_3),
               as.numeric(AMELIA_wide_start_c$SEX_4))
    f4 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2),
               as.numeric(AMELIA_wide_start_c$SEX_3),
               as.numeric(AMELIA_wide_start_c$SEX_4), as.numeric(AMELIA_wide_start_c$SEX_5))
    f5 <- pmin(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), 
               as.numeric(AMELIA_wide_start_c$SEX_3),
               as.numeric(AMELIA_wide_start_c$SEX_4), as.numeric(AMELIA_wide_start_c$SEX_5),
               as.numeric(AMELIA_wide_start_c$SEX_6))
    f6 <- pmax(as.numeric(AMELIA_wide_start_c$SEX_1), as.numeric(AMELIA_wide_start_c$SEX_2), 
               as.numeric(AMELIA_wide_start_c$SEX_3),
               as.numeric(AMELIA_wide_start_c$SEX_4), as.numeric(AMELIA_wide_start_c$SEX_5), 
               as.numeric(AMELIA_wide_start_c$SEX_6))
    AMELIA_wide_start_c$SEX_diff <- droplevels(interaction(f1, f2, f3, f4, f5, f6, sep = ""))
  }
  return(AMELIA_wide_start_c)
}
# Create a function to calculate the deviation between the synthetic population and the true
# population for 1-Dimensional  aggregated data
calc_bias <- function(est, true) {
  r <- list()
  est$AGE_comb <- factor(est$AGE_comb, levels = levels(true$AGE_comb))
  est$MST_group <- factor(est$MST_group, levels = levels(true$MST_group))
  est$SEX_diff <- factor(est$SEX_diff, levels = levels(true$SEX_diff))
  r$deviation_AGE <- sum(abs(prop.table(table(est$AGE_comb)) - prop.table(table(true$AGE_comb)))) / length(levels(true$AGE_comb))
  r$deviation_MST <- sum(abs(prop.table(table(est$MST_group)) - prop.table(table(true$MST_group)))) / length(levels(true$MST_group))
  r$deviation_SEX <- sum(abs(prop.table(table(est$SEX_diff)) - prop.table(table(true$SEX_diff)))) / length(levels(true$SEX_diff))
  return(r)
}

# For the full population we need 
# 1178764 households of size 2
# 703453 households of size 3
# 678141 households of size 4
# 249307 households of size 5
# 110497 households of size 6

# Create synthetic households of size 2
# 300.000 Size of Model data = 60 * 5000
test_2 <- Monte_Carlo_Simulation(n_model_comb = 30, NP = nrow(hhs_2_true), hhs_true = hhs_2_true,
                               household_size = 2)

# and so on...

# Convert to readeable format
test_logit_2 <- data.frame(matrix(unlist(test_2[[1]]), nrow(hhs_2_true), byrow=T))
test_nn_2 <- data.frame(matrix(unlist(test_2[[2]]), nrow(hhs_2_true), byrow=T))
test_rf_2 <- data.frame(matrix(unlist(test_2[[3]]), nrow(hhs_2_true), byrow=T))

# Create age cat
test_logit_2 <- AGE_cat_creator(test_logit_2)
test_nn_2 <- AGE_cat_creator(test_nn_2)
test_rf_2 <- AGE_cat_creator(test_rf_2)

# Calc difference in aggreate data for HHS2

test_logit_2 <- calc_difference(test_logit_2, hhs_size = 2)
test_nn_2 <- calc_difference(test_nn_2, hhs_size = 2)
test_rf_2 <- calc_difference(test_rf_2, hhs_size = 2)

# and so on..

# Calc the deviation in aggregate data as given in Figure 4
calc_bias(AMELIA_wide_start_logit_2, hhs_2_true) # 
# and so on for the other households

# Load result given in .csv file (you find it in seafile or on USB)

aggregate_simulation <- read.csv("Simulation_results_thesis.csv", header = TRUE, sep = ";")
cols <- c("Strategy", "Method", "Household.size", "Aggregate.Variable")
# Cleaning the data
aggregate_simulation[cols] <- lapply(aggregate_simulation[cols], factor)
aggregate_simulation$Deviation..in... <- as.numeric(sub(",", ".", aggregate_simulation$Deviation..in..., fixed = TRUE))

# Creating the plots for figure 4
SEX_diff_plot <- ggplot(aggregate_simulation[which(aggregate_simulation$Aggregate.Variable == "SEX_diff"),]
                        ,aes(col = Strategy,
                             shape = Method,
                              y= Deviation..in..., 
                              x=Household.size)) + 
  geom_point() +
  geom_line() +
  theme(legend.position="bottom") +
  xlab("Household size") +
  ylab("Deviation (%)") +
  theme_bw()
SEX_diff_plot

MST_group_plot <- ggplot(aggregate_simulation[which(aggregate_simulation$Aggregate.Variable == 
                                                      "MST_group"),]
                        ,aes(col = Strategy,
                             shape = Method,
                             y= Deviation..in..., 
                             x=Household.size)) + 
  geom_point() +
  geom_line() +
  theme(legend.position="bottom") +
  xlab("Household size") +
  ylab("Deviation (%)") +
  theme_bw()
MST_group_plot

AGE_comb_plot <- ggplot(aggregate_simulation[which(aggregate_simulation$Aggregate.Variable == 
                                                      "AGE_Comb"),]
                         ,aes(col = Strategy,
                              shape = Method,
                              y= Deviation..in..., 
                              x=Household.size)) + 
  geom_point() +
  geom_line() +
  theme(legend.position="bottom") +
  xlab("Household size") +
  ylab("Deviation (%)") +
  theme_bw()
AGE_comb_plot

###########################
# HERE BEGINS THE APPROACH USING MODELS FOR THE CHARACTERISTICS OF A HOUSEHOLD MEMBER BASED BY 
# OTHER MEMBERS (Model approach)

# Subset the household size classes from the sample. Why? Because we need a model for every
# household size

AMELIA_PID2_wide_sample <- AMELIA_wide_sample[which(AMELIA_wide_sample$HHS == 2),]
AMELIA_PID3_wide_sample <- AMELIA_wide_sample[which(AMELIA_wide_sample$HHS == 3),]
AMELIA_PID4_wide_sample <- AMELIA_wide_sample[which(AMELIA_wide_sample$HHS == 4),]
AMELIA_PID5_wide_sample <- AMELIA_wide_sample[which(AMELIA_wide_sample$HHS == 5),]
AMELIA_PID6_wide_sample <- AMELIA_wide_sample[which(AMELIA_wide_sample$HHS == 6),]


# To be able to apply our model later on, we need a starting point. Since the characteristics of
# the individuals are estimated like a pyramid (starting with the oldest one) it has to be determined
# who is the oldest person in a household. Therefore, we will start with a table which gives
# the true values of the oldest person. In reality, this would have to be estimated or created
# via resampling.

# Create starting points
AMELIA_wide_start_nnet_2 <- hhs_2_true[, c("AGE_1", "AGE_cat_1", "MST_1", "SEX_1")]
AMELIA_wide_start_rf_2 <- hhs_2_true[, c("AGE_1", "AGE_cat_1", "MST_1", "SEX_1")]
AMELIA_wide_start_logit_2 <- hhs_2_true[, c("AGE_1", "AGE_cat_1", "MST_1", "SEX_1")]

AMELIA_wide_start_nnet_3 <- hhs_3_true[, c("AGE_1","AGE_cat_1", "MST_1", "SEX_1")]
AMELIA_wide_start_rf_3 <- hhs_3_true[, c("AGE_1","AGE_cat_1", "MST_1", "SEX_1")]
AMELIA_wide_start_logit_3 <- hhs_3_true[, c("AGE_1", "AGE_cat_1", "MST_1", "SEX_1")]

AMELIA_wide_start_nnet_4 <- hhs_4_true[, c("AGE_1","AGE_cat_1", "MST_1", "SEX_1")]
AMELIA_wide_start_rf_4 <- hhs_4_true[, c("AGE_1","AGE_cat_1", "MST_1", "SEX_1")]
AMELIA_wide_start_logit_4 <- hhs_4_true[, c("AGE_1", "AGE_cat_1", "MST_1", "SEX_1")]

AMELIA_wide_start_nnet_5 <-  hhs_5_true[, c("AGE_1","AGE_cat_1", "MST_1", "SEX_1")]
AMELIA_wide_start_rf_5 <- hhs_5_true[, c("AGE_1","AGE_cat_1", "MST_1", "SEX_1")]
AMELIA_wide_start_logit_5 <- hhs_5_true[, c("AGE_1", "AGE_cat_1", "MST_1", "SEX_1")]

AMELIA_wide_start_nnet_6 <-  hhs_6_true[, c("AGE_1","AGE_cat_1", "MST_1", "SEX_1")]
AMELIA_wide_start_rf_6 <-  hhs_6_true[, c("AGE_1","AGE_cat_1", "MST_1", "SEX_1")]
AMELIA_wide_start_logit_6 <- hhs_6_true[, c("AGE_1", "AGE_cat_1", "MST_1", "SEX_1")]
#######################
# 2 # RUN THE MODEL APPROACH (takes some time)

nn_AGE_cat_2 <- nnet(AGE_cat_2 ~ AGE_1 + AGE_cat_1 + MST_1 + SEX_1, size = 4, 
                     data = AMELIA_PID2_wide_sample, maxit = 500)
p_nn_AGE_cat_2 <- predict(nn_AGE_cat_2, newdata = AMELIA_wide_start_nnet_2, type = "raw")
AGE_cat_2 <- apply(p_nn_AGE_cat_2, 1, function(x)sample(colnames(p_nn_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_nnet_2[, AGE_cat_2 := as.factor(AGE_cat_2)] 

rf_AGE_cat_2 <- randomForest(AGE_cat_2 ~ AGE_cat_1 + MST_1 + SEX_1, 
                             data = AMELIA_PID2_wide_sample)
p_rf_AGE_cat_2 <- predict(rf_AGE_cat_2, newdata = AMELIA_wide_start_rf_2, type = "prob")
AGE_cat_2 <- apply(p_rf_AGE_cat_2, 1, function(x)sample(colnames(p_rf_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_rf_2[, AGE_cat_2 := as.factor(AGE_cat_2)] 

logit_AGE_cat_2 <- multinom(AGE_cat_2 ~ AGE_cat_1 + MST_1 + SEX_1, 
                             data = AMELIA_PID2_wide_sample)
p_logit_AGE_cat_2 <- predict(logit_AGE_cat_2, newdata = AMELIA_wide_start_logit_2, type = "prob")
AGE_cat_2 <- apply(p_logit_AGE_cat_2, 1, function(x)sample(colnames(p_logit_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_logit_2[, AGE_cat_2 := as.factor(AGE_cat_2)] 

nn_SEX_2 <- nnet(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, size = 4, 
                 maxit = 500, data = AMELIA_PID2_wide_sample)
p_nn_SEX_2 <- predict(nn_SEX_2, newdata = AMELIA_wide_start_nnet_2, type = "raw")
p_nn_SEX_2_2 <- 1 - p_nn_SEX_2
p_nn_SEX_2 <- cbind(p_nn_SEX_2, p_nn_SEX_2_2)
colnames(p_nn_SEX_2) <- c(2, 1)
SEX_2 <- apply(p_nn_SEX_2, 1, function(x)sample(colnames(p_nn_SEX_2), size = 1, prob = x))
AMELIA_wide_start_nnet_2[, SEX_2 := as.factor(SEX_2)] 

rf_SEX_2 <- randomForest(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, data = AMELIA_PID2_wide_sample)
p_rf_SEX_2 <- predict(rf_SEX_2, newdata = AMELIA_wide_start_rf_2, type = "prob")
SEX_2 <- apply(p_rf_SEX_2, 1, function(x)sample(colnames(p_rf_SEX_2), size = 1, prob = x))
AMELIA_wide_start_rf_2[, SEX_2 := as.factor(SEX_2)] 

logit_SEX_2 <- multinom(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1,
                        data = AMELIA_PID2_wide_sample)
p_logit_SEX_2 <- predict(logit_SEX_2, newdata = AMELIA_wide_start_logit_2, type = "prob")
p_logit_SEX_2_2 <- 1 - p_logit_SEX_2
p_logit_SEX_2 <- cbind(p_logit_SEX_2, p_logit_SEX_2_2)
colnames(p_logit_SEX_2) <- c(2, 1)
SEX_2 <- apply(p_logit_SEX_2, 1, function(x)sample(colnames(p_logit_SEX_2), size = 1, prob = x))
AMELIA_wide_start_logit_2[, SEX_2 := as.factor(SEX_2)] 

nn_MST_2 <- nnet(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, size = 4,
                 maxit = 500, data = AMELIA_PID2_wide_sample)
p_nn_MST_2 <- predict(nn_MST_2, newdata = AMELIA_wide_start_nnet_2, type = "raw")
MST_2 <- apply(p_nn_MST_2, 1, function(x)sample(colnames(p_nn_MST_2), size = 1, prob = x))
AMELIA_wide_start_nnet_2[, MST_2 := as.factor(MST_2)] 

rf_MST_2 <- randomForest(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, data = AMELIA_PID2_wide_sample)
p_rf_MST_2 <- predict(rf_MST_2, newdata = AMELIA_wide_start_rf_2, type = "prob")
MST_2_rf <- apply(p_rf_MST_2, 1, function(x)sample(colnames(p_rf_MST_2), size = 1, prob = x))
AMELIA_wide_start_rf_2[, MST_2 := as.factor(MST_2_rf)] 

logit_MST_2 <- multinom(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, data = AMELIA_PID2_wide_sample)
p_logit_MST_2 <- predict(logit_MST_2, newdata = AMELIA_wide_start_logit_2, type = "prob")
MST_2 <- apply(p_logit_MST_2, 1, function(x)sample(colnames(p_logit_MST_2), size = 1, prob = x))
AMELIA_wide_start_logit_2[, MST_2 := as.factor(MST_2)] 
############### 3
##########
nn_AGE_cat_2 <- nnet(AGE_cat_2 ~ AGE_1 + AGE_cat_1 + MST_1 + SEX_1, size = 4, 
                     data = AMELIA_PID3_wide_sample, maxit = 500)
p_nn_AGE_cat_2 <- predict(nn_AGE_cat_2, newdata = AMELIA_wide_start_nnet_3, type = "raw")
AGE_cat_2 <- apply(p_nn_AGE_cat_2, 1, function(x)sample(colnames(p_nn_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_nnet_3[, AGE_cat_2 := as.factor(AGE_cat_2)] 

rf_AGE_cat_2 <- randomForest(AGE_cat_2 ~ AGE_cat_1 + MST_1 + SEX_1, 
                             data = AMELIA_PID3_wide_sample)
p_rf_AGE_cat_2 <- predict(rf_AGE_cat_2, newdata = AMELIA_wide_start_rf_3, type = "prob")
AGE_cat_2 <- apply(p_rf_AGE_cat_2, 1, function(x)sample(colnames(p_rf_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_rf_3[, AGE_cat_2 := as.factor(AGE_cat_2)] 

logit_AGE_cat_2 <- multinom(AGE_cat_2 ~ AGE_cat_1 + MST_1 + SEX_1, 
                            data = AMELIA_PID3_wide_sample)
p_logit_AGE_cat_2 <- predict(logit_AGE_cat_2, newdata = AMELIA_wide_start_logit_3, type = "prob")
AGE_cat_2 <- apply(p_logit_AGE_cat_2, 1, function(x)sample(colnames(p_logit_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_logit_3[, AGE_cat_2 := as.factor(AGE_cat_2)] 

nn_SEX_2 <- nnet(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, size = 4,
                 data = AMELIA_PID3_wide_sample, maxit = 500)
p_nn_SEX_2 <- predict(nn_SEX_2, newdata = AMELIA_wide_start_nnet_3, type = "raw")
p_nn_SEX_2_2 <- 1 - p_nn_SEX_2
p_nn_SEX_2 <- cbind(p_nn_SEX_2, p_nn_SEX_2_2)
colnames(p_nn_SEX_2) <- c(2, 1)
SEX_2 <- apply(p_nn_SEX_2, 1, function(x)sample(colnames(p_nn_SEX_2), size = 1, prob = x))
AMELIA_wide_start_nnet_3[, SEX_2 := as.factor(SEX_2)] 

rf_SEX_2 <- randomForest(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, data = AMELIA_PID3_wide_sample)
p_rf_SEX_2 <- predict(rf_SEX_2, newdata = AMELIA_wide_start_rf_3, type = "prob")
SEX_2 <- apply(p_rf_SEX_2, 1, function(x)sample(colnames(p_rf_SEX_2), size = 1, prob = x))
AMELIA_wide_start_rf_3[, SEX_2 := as.factor(SEX_2)] 

logit_SEX_2 <- multinom(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1,
                        data = AMELIA_PID3_wide_sample)
p_logit_SEX_2 <- predict(logit_SEX_2, newdata = AMELIA_wide_start_logit_3, type = "prob")
p_logit_SEX_2_2 <- 1 - p_logit_SEX_2
p_logit_SEX_2 <- cbind(p_logit_SEX_2, p_logit_SEX_2_2)
colnames(p_logit_SEX_2) <- c(2, 1)
SEX_2 <- apply(p_logit_SEX_2, 1, function(x)sample(colnames(p_logit_SEX_2), size = 1, prob = x))
AMELIA_wide_start_logit_3[, SEX_2 := as.factor(SEX_2)] 

nn_MST_2 <- nnet(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, 
                 size = 4, data = AMELIA_PID3_wide_sample, maxit = 500)
p_nn_MST_2 <- predict(nn_MST_2, newdata = AMELIA_wide_start_nnet_3, type = "raw")
MST_2 <- apply(p_nn_MST_2, 1, function(x)sample(colnames(p_nn_MST_2), size = 1, prob = x))
AMELIA_wide_start_nnet_3[, MST_2 := as.factor(MST_2)] 

rf_MST_2 <- randomForest(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, data = AMELIA_PID3_wide_sample)
p_rf_MST_2 <- predict(rf_MST_2, newdata = AMELIA_wide_start_rf_3, type = "prob")
MST_2 <- apply(p_rf_MST_2, 1, function(x)sample(colnames(p_rf_MST_2), size = 1, prob = x))
AMELIA_wide_start_rf_3[, MST_2 := as.factor(MST_2)] 

logit_MST_2 <- multinom(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, data = AMELIA_PID3_wide_sample)
p_logit_MST_2 <- predict(logit_MST_2, newdata = AMELIA_wide_start_logit_3, type = "prob")
MST_2 <- apply(p_logit_MST_2, 1, function(x)sample(colnames(p_logit_MST_2), size = 1, prob = x))
AMELIA_wide_start_logit_3[, MST_2 := as.factor(MST_2)] 

nn_AGE_cat_3 <- nnet(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                     size = 4, data = AMELIA_PID3_wide_sample, maxit = 500)
p_nn_AGE_cat_3 <- predict(nn_AGE_cat_3, newdata = AMELIA_wide_start_nnet_3, type = "raw")
AGE_cat_3 <- apply(p_nn_AGE_cat_3, 1, function(x)sample(colnames(p_nn_AGE_cat_3), size = 1, prob = x))
AMELIA_wide_start_nnet_3[, AGE_cat_3 := as.factor(AGE_cat_3)] 

rf_AGE_cat_3 <- randomForest(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                     data = AMELIA_PID3_wide_sample)
p_rf_AGE_cat_3 <- predict(rf_AGE_cat_3, newdata = AMELIA_wide_start_rf_3, type = "prob")
AGE_cat_3 <- apply(p_rf_AGE_cat_3, 1, function(x)sample(colnames(p_rf_AGE_cat_3), size = 1, prob = x))
AMELIA_wide_start_rf_3[, AGE_cat_3 := as.factor(AGE_cat_3)] 

logit_AGE_cat_3 <- randomForest(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                             data = AMELIA_PID3_wide_sample)
p_logit_AGE_cat_3 <- predict(logit_AGE_cat_3, newdata = AMELIA_wide_start_logit_3, type = "prob")
AGE_cat_3 <- apply(p_logit_AGE_cat_3, 1, function(x)sample(colnames(p_logit_AGE_cat_3), size = 1, prob = x))
AMELIA_wide_start_logit_3[, AGE_cat_3 := as.factor(AGE_cat_3)] 

nn_SEX_3 <- nnet(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                              size = 4, data = AMELIA_PID3_wide_sample, maxit = 500)
p_nn_SEX_3 <- predict(nn_SEX_3, newdata = AMELIA_wide_start_nnet_3, type = "raw")
p_nn_SEX_3_2 <- 1 - p_nn_SEX_3
p_nn_SEX_3 <- cbind(p_nn_SEX_3, p_nn_SEX_3_2)
colnames(p_nn_SEX_3) <- c(2, 1)
SEX_3 <- apply(p_nn_SEX_3, 1, function(x)sample(colnames(p_nn_SEX_3), size = 1, prob = x))
AMELIA_wide_start_nnet_3[, SEX_3 := as.factor(SEX_3)]

rf_SEX_3 <- randomForest(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
              data = AMELIA_PID3_wide_sample)
p_rf_SEX_3 <- predict(rf_SEX_3, newdata = AMELIA_wide_start_rf_3, type = "prob")
SEX_3 <- apply(p_rf_SEX_3, 1, function(x)sample(colnames(p_rf_SEX_3), size = 1, prob = x))
AMELIA_wide_start_rf_3[, SEX_3 := as.factor(SEX_3)] 

logit_SEX_3 <- multinom(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                 data = AMELIA_PID3_wide_sample)
p_logit_SEX_3 <- predict(logit_SEX_3, newdata = AMELIA_wide_start_nnet_3, type = "prob")
p_logit_SEX_3_2 <- 1 - p_logit_SEX_3
p_logit_SEX_3 <- cbind(p_logit_SEX_3, p_logit_SEX_3_2)
colnames(p_logit_SEX_3) <- c(2, 1)
SEX_3 <- apply(p_logit_SEX_3, 1, function(x)sample(colnames(p_logit_SEX_3), size = 1, prob = x))
AMELIA_wide_start_logit_3[, SEX_3 := as.factor(SEX_3)]

nn_MST_3 <- nnet(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                 size = 4, data = AMELIA_PID3_wide_sample, maxit = 500)
p_nn_MST_3 <- predict(nn_MST_3, newdata = AMELIA_wide_start_nnet_3, type = "raw")
MST_3 <- apply(p_nn_MST_3, 1, function(x)sample(colnames(p_nn_MST_3), size = 1, prob = x))
AMELIA_wide_start_nnet_3[, MST_3 := as.factor(MST_3)]

rf_MST_3 <- randomForest(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                         data = AMELIA_PID3_wide_sample)
p_rf_MST_3 <- predict(rf_MST_3, newdata = AMELIA_wide_start_rf_3, type = "prob")
MST_3 <- apply(p_rf_MST_3, 1, function(x)sample(colnames(p_rf_MST_3), size = 1, prob = x))
AMELIA_wide_start_rf_3[, MST_3 := as.factor(MST_3)] 

logit_MST_3 <- multinom(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                         data = AMELIA_PID3_wide_sample)
p_logit_MST_3 <- predict(logit_MST_3, newdata = AMELIA_wide_start_logit_3, type = "prob")
MST_3 <- apply(p_logit_MST_3, 1, function(x)sample(colnames(p_logit_MST_3), size = 1, prob = x))
AMELIA_wide_start_logit_3[, MST_3 := as.factor(MST_3)] 

################### 4
###############

nn_AGE_cat_2 <- nnet(AGE_cat_2 ~ AGE_1 + AGE_cat_1 + MST_1 + SEX_1, size = 4, 
                     data = AMELIA_PID4_wide_sample, maxit = 500)
p_nn_AGE_cat_2 <- predict(nn_AGE_cat_2, newdata = AMELIA_wide_start_nnet_4, type = "raw")
AGE_cat_2 <- apply(p_nn_AGE_cat_2, 1, function(x)sample(colnames(p_nn_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_nnet_4[, AGE_cat_2 := as.factor(AGE_cat_2)] 

rf_AGE_cat_2 <- randomForest(AGE_cat_2 ~ AGE_cat_1 + MST_1 + SEX_1, 
                             data = AMELIA_PID4_wide_sample)
p_rf_AGE_cat_2 <- predict(rf_AGE_cat_2, newdata = AMELIA_wide_start_rf_4, type = "prob")
AGE_cat_2_rf <- apply(p_rf_AGE_cat_2, 1, function(x)sample(colnames(p_rf_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_rf_4[, AGE_cat_2 := as.factor(AGE_cat_2_rf)] 

logit_AGE_cat_2 <- multinom(AGE_cat_2 ~ AGE_cat_1 + MST_1 + SEX_1, 
                             data = AMELIA_PID4_wide_sample)
p_logit_AGE_cat_2 <- predict(logit_AGE_cat_2, newdata = AMELIA_wide_start_logit_4, type = "prob")
AGE_cat_2_logit <- apply(p_logit_AGE_cat_2, 1, function(x)sample(colnames(p_logit_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_logit_4[, AGE_cat_2 := as.factor(AGE_cat_2_rf)] 

nn_SEX_2 <- nnet(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, size = 4, 
                 data = AMELIA_PID4_wide_sample, maxit = 500)
p_nn_SEX_2 <- predict(nn_SEX_2, newdata = AMELIA_wide_start_nnet_4, type = "raw")
p_nn_SEX_2_2 <- 1 - p_nn_SEX_2
p_nn_SEX_2 <- cbind(p_nn_SEX_2, p_nn_SEX_2_2)
colnames(p_nn_SEX_2) <- c(2, 1)
SEX_2 <- apply(p_nn_SEX_2, 1, function(x)sample(colnames(p_nn_SEX_2), size = 1, prob = x))
AMELIA_wide_start_nnet_4[, SEX_2 := as.factor(SEX_2)] 

rf_SEX_2 <- randomForest(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, data = AMELIA_PID4_wide_sample)
p_rf_SEX_2 <- predict(rf_SEX_2, newdata = AMELIA_wide_start_rf_4, type = "prob")
SEX_2 <- apply(p_rf_SEX_2, 1, function(x)sample(colnames(p_rf_SEX_2), size = 1, prob = x))
AMELIA_wide_start_rf_4[, SEX_2 := as.factor(SEX_2)] 

logit_SEX_2 <- nnet(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, size = 4, 
                 data = AMELIA_PID4_wide_sample, maxit = 500)
p_logit_SEX_2 <- predict(logit_SEX_2, newdata = AMELIA_wide_start_logit_4, type = "raw")
p_logit_SEX_2_2 <- 1 - p_logit_SEX_2
p_logit_SEX_2 <- cbind(p_logit_SEX_2, p_logit_SEX_2_2)
colnames(p_logit_SEX_2) <- c(2, 1)
SEX_2 <- apply(p_logit_SEX_2, 1, function(x)sample(colnames(p_logit_SEX_2), size = 1, prob = x))
AMELIA_wide_start_logit_4[, SEX_2 := as.factor(SEX_2)] 

nn_MST_2 <- nnet(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, size = 4, 
                 data = AMELIA_PID4_wide_sample, maxit = 500)
p_nn_MST_2 <- predict(nn_MST_2, newdata = AMELIA_wide_start_nnet_4, type = "raw")
MST_2 <- apply(p_nn_MST_2, 1, function(x)sample(colnames(p_nn_MST_2), size = 1, prob = x))
AMELIA_wide_start_nnet_4[, MST_2 := as.factor(MST_2)] 

rf_MST_2 <- randomForest(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, data = AMELIA_PID4_wide_sample)
p_rf_MST_2 <- predict(rf_MST_2, newdata = AMELIA_wide_start_rf_4, type = "prob")
MST_2 <- apply(p_rf_MST_2, 1, function(x)sample(colnames(p_rf_MST_2), size = 1, prob = x))
AMELIA_wide_start_rf_4[, MST_2 := as.factor(MST_2)] 

logit_MST_2 <- multinom(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, 
                        data = AMELIA_PID4_wide_sample)
p_logit_MST_2 <- predict(logit_MST_2, newdata = AMELIA_wide_start_logit_4, 
                         type = "prob")
MST_2 <- apply(p_logit_MST_2, 1, function(x)sample(colnames(p_logit_MST_2), size = 1, prob = x))
AMELIA_wide_start_logit_4[, MST_2 := as.factor(MST_2)] 

nn_AGE_cat_3 <- nnet(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                     size = 4, data = AMELIA_PID4_wide_sample, maxit = 500)
p_nn_AGE_cat_3 <- predict(nn_AGE_cat_3, newdata = AMELIA_wide_start_nnet_4, type = "raw")
AGE_cat_3 <- apply(p_nn_AGE_cat_3, 1, function(x)sample(colnames(p_nn_AGE_cat_3), size = 1, prob = x))
AMELIA_wide_start_nnet_4[, AGE_cat_3 := as.factor(AGE_cat_3)] 

rf_AGE_cat_3 <- randomForest(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                             data = AMELIA_PID4_wide_sample)
p_rf_AGE_cat_3 <- predict(rf_AGE_cat_3, newdata = AMELIA_wide_start_rf_4, type = "prob")
AGE_cat_3 <- apply(p_rf_AGE_cat_3, 1, function(x)sample(colnames(p_rf_AGE_cat_3), size = 1, prob = x))
AMELIA_wide_start_rf_4[, AGE_cat_3 := as.factor(AGE_cat_3)] 

logit_AGE_cat_3 <- multinom(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                             data = AMELIA_PID4_wide_sample)
p_logit_AGE_cat_3 <- predict(logit_AGE_cat_3, newdata = AMELIA_wide_start_logit_4, type = "prob")
AGE_cat_3 <- apply(p_logit_AGE_cat_3, 1, function(x)sample(colnames(p_logit_AGE_cat_3), size = 1, prob = x))
AMELIA_wide_start_logit_4[, AGE_cat_3 := as.factor(AGE_cat_3)]

nn_SEX_3 <- nnet(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                 size = 4, data = AMELIA_PID4_wide_sample, maxit = 500)
p_nn_SEX_3 <- predict(nn_SEX_3, newdata = AMELIA_wide_start_nnet_4, type = "raw")
p_nn_SEX_3_2 <- 1 - p_nn_SEX_3
p_nn_SEX_3 <- cbind(p_nn_SEX_3, p_nn_SEX_3_2)
colnames(p_nn_SEX_3) <- c(2, 1)
SEX_3 <- apply(p_nn_SEX_3, 1, function(x)sample(colnames(p_nn_SEX_3), size = 1, prob = x))
AMELIA_wide_start_nnet_4[, SEX_3 := as.factor(SEX_3)]

rf_SEX_3 <- randomForest(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                         data = AMELIA_PID4_wide_sample)
p_rf_SEX_3 <- predict(rf_SEX_3, newdata = AMELIA_wide_start_rf_4, type = "prob")
SEX_3 <- apply(p_rf_SEX_3, 1, function(x)sample(colnames(p_rf_SEX_3), size = 1, prob = x))
AMELIA_wide_start_rf_4[, SEX_3 := as.factor(SEX_3)] 

logit_SEX_3 <- nnet(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                 size = 4, data = AMELIA_PID4_wide_sample)
p_logit_SEX_3 <- predict(logit_SEX_3, newdata = AMELIA_wide_start_logit_4, type = "raw")
p_logit_SEX_3_2 <- 1 - p_logit_SEX_3
p_logit_SEX_3 <- cbind(p_logit_SEX_3, p_logit_SEX_3_2)
colnames(p_logit_SEX_3) <- c(2, 1)
SEX_3 <- apply(p_logit_SEX_3, 1, function(x)sample(colnames(p_logit_SEX_3), size = 1, prob = x))
AMELIA_wide_start_logit_4[, SEX_3 := as.factor(SEX_3)]

nn_MST_3 <- nnet(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                 size = 4, data = AMELIA_PID4_wide_sample, maxit = 500)
p_nn_MST_3 <- predict(nn_MST_3, newdata = AMELIA_wide_start_nnet_4, type = "raw")
MST_3 <- apply(p_nn_MST_3, 1, function(x)sample(colnames(p_nn_MST_3), size = 1, prob = x))
AMELIA_wide_start_nnet_4[, MST_3 := as.factor(MST_3)]

rf_MST_3 <- randomForest(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                         data = AMELIA_PID4_wide_sample)
p_rf_MST_3 <- predict(rf_MST_3, newdata = AMELIA_wide_start_rf_4, type = "prob")
MST_3 <- apply(p_rf_MST_3, 1, function(x)sample(colnames(p_rf_MST_3), size = 1, prob = x))
AMELIA_wide_start_rf_4[, MST_3 := as.factor(MST_3)] 

logit_MST_3 <- multinom(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                         data = AMELIA_PID4_wide_sample)
p_logit_MST_3 <- predict(logit_MST_3, newdata = AMELIA_wide_start_logit_4, type = "prob")
MST_3 <- apply(p_logit_MST_3, 1, function(x)sample(colnames(p_logit_MST_3), size = 1, prob = x))
AMELIA_wide_start_logit_4[, MST_3 := as.factor(MST_3)] 

nn_AGE_cat_4 <- nnet(AGE_cat_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + MST_3 +
                       SEX_1 + SEX_2 + SEX_3, 
                       size = 4, data = AMELIA_PID4_wide_sample, maxit = 500)
p_nn_AGE_cat_4 <- predict(nn_AGE_cat_4, newdata = AMELIA_wide_start_nnet_4, type = "raw")
AGE_cat_4 <- apply(p_nn_AGE_cat_4, 1, function(x)sample(colnames(p_nn_AGE_cat_4), size = 1, prob = x))
AMELIA_wide_start_nnet_4[, AGE_cat_4 := as.factor(AGE_cat_4)]

rf_AGE_cat_4 <- randomForest(AGE_cat_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + MST_3 +
                       SEX_1 + SEX_2 +  SEX_3, 
                       data = AMELIA_PID4_wide_sample)
p_rf_AGE_cat_4 <- predict(rf_AGE_cat_4, newdata = AMELIA_wide_start_rf_4, type = "prob")
AGE_cat_4 <- apply(p_rf_AGE_cat_4, 1, function(x)sample(colnames(p_rf_AGE_cat_4), size = 1, prob = x))
AMELIA_wide_start_rf_4[, AGE_cat_4 := as.factor(AGE_cat_4)] 

logit_AGE_cat_4 <- multinom(AGE_cat_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + MST_3 +
                               SEX_1 + SEX_2 +  SEX_3, 
                             data = AMELIA_PID4_wide_sample)
p_logit_AGE_cat_4 <- predict(logit_AGE_cat_4, newdata = AMELIA_wide_start_logit_4, type = "prob")
AGE_cat_4 <- apply(p_logit_AGE_cat_4, 1, function(x)sample(colnames(p_logit_AGE_cat_4), size = 1, prob = x))
AMELIA_wide_start_logit_4[, AGE_cat_4 := as.factor(AGE_cat_4)] 

nn_SEX_4 <- nnet(SEX_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                       SEX_1 + SEX_2 +  SEX_3, 
                     size = 4, data = AMELIA_PID4_wide_sample, maxit = 500)
p_nn_SEX_4 <- predict(nn_SEX_4, newdata = AMELIA_wide_start_nnet_4, type = "raw")
p_nn_SEX_4_2 <- 1 - p_nn_SEX_4
p_nn_SEX_4 <- cbind(p_nn_SEX_4, p_nn_SEX_4_2)
colnames(p_nn_SEX_4) <- c(2, 1)
SEX_4 <- apply(p_nn_SEX_4, 1, function(x)sample(colnames(p_nn_SEX_4), size = 1, prob = x))
AMELIA_wide_start_nnet_4[, SEX_4 := as.factor(SEX_4)]

rf_SEX_4 <- randomForest(SEX_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                   SEX_1 + SEX_2 +  SEX_3, 
                   data = AMELIA_PID4_wide_sample)
p_rf_SEX_4 <- predict(rf_SEX_4, newdata = AMELIA_wide_start_rf_4, type = "prob")
SEX_4 <- apply(p_rf_SEX_4, 1, function(x)sample(colnames(p_rf_SEX_4), size = 1, prob = x))
AMELIA_wide_start_rf_4[, SEX_4 := as.factor(SEX_4)] 

logit_SEX_4 <- multinom(SEX_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                   SEX_1 + SEX_2 +  SEX_3, 
                   data = AMELIA_PID4_wide_sample)
p_logit_SEX_4 <- predict(logit_SEX_4, newdata = AMELIA_wide_start_logit_4, type = "prob")
p_logit_SEX_4_2 <- 1 - p_logit_SEX_4
p_logit_SEX_4 <- cbind(p_logit_SEX_4, p_logit_SEX_4_2)
colnames(p_logit_SEX_4) <- c(2, 1)
SEX_4 <- apply(p_logit_SEX_4, 1, function(x)sample(colnames(p_logit_SEX_4), size = 1, prob = x))
AMELIA_wide_start_logit_4[, SEX_4 := as.factor(SEX_4)]

nn_MST_4 <- nnet(MST_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                   SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                 size = 4, data = AMELIA_PID4_wide_sample, maxit = 500)
p_nn_MST_4 <- predict(nn_MST_4, newdata = AMELIA_wide_start_nnet_4, type = "raw")
MST_4 <- apply(p_nn_MST_4, 1, function(x)sample(colnames(p_nn_MST_4), size = 1, prob = x))
AMELIA_wide_start_nnet_4[, MST_4 := as.factor(MST_4)]

rf_MST_4 <- randomForest(MST_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                           SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                   data = AMELIA_PID4_wide_sample)
p_rf_MST_4 <- predict(rf_MST_4, newdata = AMELIA_wide_start_rf_4, type = "prob")
MST_4 <- apply(p_rf_MST_4, 1, function(x)sample(colnames(p_rf_MST_4), size = 1, prob = x))
AMELIA_wide_start_rf_4[, MST_4 := as.factor(MST_4)] 

logit_MST_4 <- multinom(MST_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                           SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                         data = AMELIA_PID4_wide_sample)
p_logit_MST_4 <- predict(rf_MST_4, newdata = AMELIA_wide_start_logit_4, type = "prob")
MST_4 <- apply(p_logit_MST_4, 1, function(x)sample(colnames(p_logit_MST_4), size = 1, prob = x))
AMELIA_wide_start_logit_4[, MST_4 := as.factor(MST_4)] 

########################
########################5
nn_AGE_cat_2 <- nnet(AGE_cat_2 ~ AGE_1 + AGE_cat_1 + MST_1 + SEX_1, size = 4, 
                     data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_AGE_cat_2 <- predict(nn_AGE_cat_2, newdata = AMELIA_wide_start_nnet_5, type = "raw")
AGE_cat_2 <- apply(p_nn_AGE_cat_2, 1, function(x)sample(colnames(p_nn_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, AGE_cat_2 := as.factor(AGE_cat_2)] 

rf_AGE_cat_2 <- randomForest(AGE_cat_2 ~ AGE_cat_1 + MST_1 + SEX_1, 
                             data = AMELIA_PID5_wide_sample)
p_rf_AGE_cat_2 <- predict(rf_AGE_cat_2, newdata = AMELIA_wide_start_rf_5, type = "prob")
AGE_cat_2 <- apply(p_rf_AGE_cat_2, 1, function(x)sample(colnames(p_rf_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_rf_5[, AGE_cat_2 := as.factor(AGE_cat_2)] 

logit_AGE_cat_2 <- multinom(AGE_cat_2 ~ AGE_cat_1 + MST_1 + SEX_1, 
                             data = AMELIA_PID5_wide_sample)
p_logit_AGE_cat_2 <- predict(logit_AGE_cat_2, newdata = AMELIA_wide_start_logit_5, type = "prob")
AGE_cat_2 <- apply(p_logit_AGE_cat_2, 1, function(x)sample(colnames(p_logit_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_logit_5[, AGE_cat_2 := as.factor(AGE_cat_2)] 

nn_SEX_2 <- nnet(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, size = 4, 
                 data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_SEX_2 <- predict(nn_SEX_2, newdata = AMELIA_wide_start_nnet_5, type = "raw")
p_nn_SEX_2_2 <- 1 - p_nn_SEX_2
p_nn_SEX_2 <- cbind(p_nn_SEX_2, p_nn_SEX_2_2)
colnames(p_nn_SEX_2) <- c(2, 1)
SEX_2 <- apply(p_nn_SEX_2, 1, function(x)sample(colnames(p_nn_SEX_2), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, SEX_2 := as.factor(SEX_2)] 

rf_SEX_2 <- randomForest(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, 
                         data = AMELIA_PID5_wide_sample)
p_rf_SEX_2 <- predict(rf_SEX_2, newdata = AMELIA_wide_start_rf_5, type = "prob")
SEX_2 <- apply(p_rf_SEX_2, 1, function(x)sample(colnames(p_rf_SEX_2), size = 1, prob = x))
AMELIA_wide_start_rf_5[, SEX_2 := as.factor(SEX_2)] 

logit_SEX_2 <- multinom(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, 
                 data = AMELIA_PID5_wide_sample)
p_logit_SEX_2 <- predict(logit_SEX_2, newdata = AMELIA_wide_start_logit_5, 
                         type = "prob")
p_logit_SEX_2_2 <- 1 - p_logit_SEX_2
p_logit_SEX_2 <- cbind(p_logit_SEX_2, p_logit_SEX_2_2)
colnames(p_logit_SEX_2) <- c(2, 1)
SEX_2 <- apply(p_logit_SEX_2, 1, function(x)sample(colnames(p_logit_SEX_2), size = 1, prob = x))
AMELIA_wide_start_logit_5[, SEX_2 := as.factor(SEX_2)] 

nn_MST_2 <- nnet(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, size = 4, 
                 data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_MST_2 <- predict(nn_MST_2, newdata = AMELIA_wide_start_nnet_5, type = "raw")
MST_2 <- apply(p_nn_MST_2, 1, function(x)sample(colnames(p_nn_MST_2), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, MST_2 := as.factor(MST_2)] 

rf_MST_2 <- randomForest(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, data = AMELIA_PID5_wide_sample)
p_rf_MST_2 <- predict(rf_MST_2, newdata = AMELIA_wide_start_rf_5, type = "prob")
MST_2_rf <- apply(p_rf_MST_2, 1, function(x)sample(colnames(p_rf_MST_2), size = 1, prob = x))
AMELIA_wide_start_rf_5[, MST_2 := as.factor(MST_2_rf)] 

logit_MST_2 <- multinom(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, data = AMELIA_PID5_wide_sample)
p_logit_MST_2 <- predict(logit_MST_2, newdata = AMELIA_wide_start_logit_5, type = "prob")
MST_2_logit <- apply(p_logit_MST_2, 1, function(x)sample(colnames(p_logit_MST_2), size = 1, prob = x))
AMELIA_wide_start_logit_5[, MST_2 := as.factor(MST_2_logit)] 

nn_AGE_cat_3 <- nnet(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                     size = 4, data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_AGE_cat_3 <- predict(nn_AGE_cat_3, newdata = AMELIA_wide_start_nnet_5, type = "raw")
AGE_cat_3 <- apply(p_nn_AGE_cat_3, 1, function(x)sample(colnames(p_nn_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, AGE_cat_3 := as.factor(AGE_cat_3)] 

rf_AGE_cat_3 <- randomForest(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                             data = AMELIA_PID5_wide_sample)
p_rf_AGE_cat_3 <- predict(rf_AGE_cat_3, newdata = AMELIA_wide_start_rf_5, type = "prob")
AGE_cat_3 <- apply(p_rf_AGE_cat_3, 1, function(x)sample(colnames(p_rf_AGE_cat_3), size = 1, prob = x))
AMELIA_wide_start_rf_5[, AGE_cat_3 := as.factor(AGE_cat_3)] 

logit_AGE_cat_3 <- multinom(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                             data = AMELIA_PID5_wide_sample)
p_logit_AGE_cat_3 <- predict(logit_AGE_cat_3, newdata = AMELIA_wide_start_logit_5, type = "prob")
AGE_cat_3 <- apply(p_logit_AGE_cat_3, 1, function(x)sample(colnames(p_logit_AGE_cat_3), size = 1, prob = x))
AMELIA_wide_start_logit_5[, AGE_cat_3 := as.factor(AGE_cat_3)] 

nn_SEX_3 <- nnet(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                 size = 4, data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_SEX_3 <- predict(nn_SEX_3, newdata = AMELIA_wide_start_nnet_5, type = "raw")
p_nn_SEX_3_2 <- 1 - p_nn_SEX_3
p_nn_SEX_3 <- cbind(p_nn_SEX_3, p_nn_SEX_3_2)
colnames(p_nn_SEX_3) <- c(2, 1)
SEX_3 <- apply(p_nn_SEX_3, 1, function(x)sample(colnames(p_nn_SEX_3), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, SEX_3 := as.factor(SEX_3)]

rf_SEX_3 <- randomForest(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                         data = AMELIA_PID5_wide_sample)
p_rf_SEX_3 <- predict(rf_SEX_3, newdata = AMELIA_wide_start_rf_5, type = "prob")
SEX_3 <- apply(p_rf_SEX_3, 1, function(x)sample(colnames(p_rf_SEX_3), size = 1, prob = x))
AMELIA_wide_start_rf_5[, SEX_3 := as.factor(SEX_3)] 

logit_SEX_3 <- multinom(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                 data = AMELIA_PID5_wide_sample)
p_logit_SEX_3 <- predict(logit_SEX_3, newdata = AMELIA_wide_start_logit_5, 
                         type = "prob")
p_logit_SEX_3_2 <- 1 - p_logit_SEX_3
p_logit_SEX_3 <- cbind(p_logit_SEX_3, p_logit_SEX_3_2)
colnames(p_logit_SEX_3) <- c(2, 1)
SEX_3 <- apply(p_logit_SEX_3, 1, function(x)sample(colnames(p_logit_SEX_3), size = 1, prob = x))
AMELIA_wide_start_logit_5[, SEX_3 := as.factor(SEX_3)]

nn_MST_3 <- nnet(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                 size = 4, data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_MST_3 <- predict(nn_MST_3, newdata = AMELIA_wide_start_nnet_5, type = "raw")
MST_3 <- apply(p_nn_MST_3, 1, function(x)sample(colnames(p_nn_MST_3), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, MST_3 := as.factor(MST_3)]

rf_MST_3 <- randomForest(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                         data = AMELIA_PID5_wide_sample)
p_rf_MST_3 <- predict(rf_MST_3, newdata = AMELIA_wide_start_rf_5, type = "prob")
MST_3_rf <- apply(p_rf_MST_3, 1, function(x)sample(colnames(p_rf_MST_3), size = 1, prob = x))
AMELIA_wide_start_rf_5[, MST_3 := as.factor(MST_3_rf)] 

logit_MST_3 <- multinom(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                         data = AMELIA_PID5_wide_sample)
p_logit_MST_3 <- predict(logit_MST_3, newdata = AMELIA_wide_start_logit_5, type = "prob")
MST_3 <- apply(p_logit_MST_3, 1, function(x)sample(colnames(p_logit_MST_3), size = 1, prob = x))
AMELIA_wide_start_logit_5[, MST_3 := as.factor(MST_3)] 

nn_AGE_cat_4 <- nnet(AGE_cat_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + MST_3 +
                       SEX_1 + SEX_2 + SEX_3, 
                     size = 4, data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_AGE_cat_4 <- predict(nn_AGE_cat_4, newdata = AMELIA_wide_start_nnet_5, type = "raw")
AGE_cat_4 <- apply(p_nn_AGE_cat_4, 1, function(x)sample(colnames(p_nn_AGE_cat_4), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, AGE_cat_4 := as.factor(AGE_cat_4)]

rf_AGE_cat_4 <- randomForest(AGE_cat_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + MST_3 +
                               SEX_1 + SEX_2 +  SEX_3, 
                             data = AMELIA_PID5_wide_sample)
p_rf_AGE_cat_4 <- predict(rf_AGE_cat_4, newdata = AMELIA_wide_start_rf_5, type = "prob")
AGE_cat_4 <- apply(p_rf_AGE_cat_4, 1, function(x)sample(colnames(p_rf_AGE_cat_4), size = 1, prob = x))
AMELIA_wide_start_rf_5[, AGE_cat_4 := as.factor(AGE_cat_4)] 

logit_AGE_cat_4 <- multinom(AGE_cat_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + MST_3 +
                               SEX_1 + SEX_2 +  SEX_3, 
                             data = AMELIA_PID5_wide_sample)
p_logit_AGE_cat_4 <- predict(logit_AGE_cat_4, newdata = AMELIA_wide_start_logit_5, type = "prob")
AGE_cat_4 <- apply(p_logit_AGE_cat_4, 1, function(x)sample(colnames(p_logit_AGE_cat_4), size = 1, prob = x))
AMELIA_wide_start_logit_5[, AGE_cat_4 := as.factor(AGE_cat_4)] 

nn_SEX_4 <- nnet(SEX_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                   SEX_1 + SEX_2 +  SEX_3, 
                 size = 4, data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_SEX_4 <- predict(nn_SEX_4, newdata = AMELIA_wide_start_nnet_5, type = "raw")
p_nn_SEX_4_2 <- 1 - p_nn_SEX_4
p_nn_SEX_4 <- cbind(p_nn_SEX_4, p_nn_SEX_4_2)
colnames(p_nn_SEX_4) <- c(2, 1)
SEX_4 <- apply(p_nn_SEX_4, 1, function(x)sample(colnames(p_nn_SEX_4), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, SEX_4 := as.factor(SEX_4)]

rf_SEX_4 <- randomForest(SEX_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                           SEX_1 + SEX_2 +  SEX_3, 
                         data = AMELIA_PID5_wide_sample)
p_rf_SEX_4 <- predict(rf_SEX_4, newdata = AMELIA_wide_start_rf_5, type = "prob")
SEX_4 <- apply(p_rf_SEX_4, 1, function(x)sample(colnames(p_rf_SEX_4), size = 1, prob = x))
AMELIA_wide_start_rf_5[, SEX_4 := as.factor(SEX_4)] 

logit_SEX_4 <- multinom(SEX_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                   SEX_1 + SEX_2 +  SEX_3, 
                   data = AMELIA_PID5_wide_sample)
p_logit_SEX_4 <- predict(logit_SEX_4, newdata = AMELIA_wide_start_logit_5, type = "prob")
p_logit_SEX_4_2 <- 1 - p_logit_SEX_4
p_logit_SEX_4 <- cbind(p_logit_SEX_4, p_logit_SEX_4_2)
colnames(p_logit_SEX_4) <- c(2, 1)
SEX_4 <- apply(p_logit_SEX_4, 1, function(x)sample(colnames(p_logit_SEX_4), size = 1, prob = x))
AMELIA_wide_start_logit_5[, SEX_4 := as.factor(SEX_4)]

nn_MST_4 <- nnet(MST_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                   SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                 size = 4, data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_MST_4 <- predict(nn_MST_4, newdata = AMELIA_wide_start_nnet_5, type = "raw")
MST_4 <- apply(p_nn_MST_4, 1, function(x)sample(colnames(p_nn_MST_4), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, MST_4 := as.factor(MST_4)]

rf_MST_4 <- randomForest(MST_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                           SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                         data = AMELIA_PID5_wide_sample)
p_rf_MST_4 <- predict(rf_MST_4, newdata = AMELIA_wide_start_rf_5, type = "prob")
MST_4_rf <- apply(p_rf_MST_4, 1, function(x)sample(colnames(p_rf_MST_4), size = 1, prob = x))
AMELIA_wide_start_rf_5[, MST_4 := as.factor(MST_4_rf)] 

logit_MST_4 <- multinom(MST_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                           SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                         data = AMELIA_PID5_wide_sample)
p_logit_MST_4 <- predict(logit_MST_4, newdata = AMELIA_wide_start_logit_5, type = "prob")
MST_4_logit <- apply(p_logit_MST_4, 1, function(x)sample(colnames(p_logit_MST_4), size = 1, prob = x))
AMELIA_wide_start_logit_5[, MST_4 := as.factor(MST_4_logit)] 

nn_AGE_cat_5 <- nnet(AGE_cat_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                       MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                     size = 4, data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_AGE_cat_5 <- predict(nn_AGE_cat_5, newdata = AMELIA_wide_start_nnet_5, type = "raw")
AGE_cat_5 <- apply(p_nn_AGE_cat_5, 1, function(x)sample(colnames(p_nn_AGE_cat_5), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, AGE_cat_5 := as.factor(AGE_cat_5)]

rf_AGE_cat_5 <- randomForest(AGE_cat_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                       MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                       data = AMELIA_PID5_wide_sample)
p_rf_AGE_cat_5 <- predict(rf_AGE_cat_5, newdata = AMELIA_wide_start_rf_5, type = "prob")
AGE_cat_5 <- apply(p_rf_AGE_cat_5, 1, function(x)sample(colnames(p_rf_AGE_cat_5), size = 1, prob = x))
AMELIA_wide_start_rf_5[, AGE_cat_5 := as.factor(AGE_cat_5)] 

logit_AGE_cat_5 <- multinom(AGE_cat_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                               MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                             data = AMELIA_PID5_wide_sample)
p_logit_AGE_cat_5 <- predict(logit_AGE_cat_5, newdata = AMELIA_wide_start_logit_5, type = "prob")
AGE_cat_5 <- apply(p_logit_AGE_cat_5, 1, function(x)sample(colnames(p_logit_AGE_cat_5), size = 1, prob = x))
AMELIA_wide_start_logit_5[, AGE_cat_5 := as.factor(AGE_cat_5)] 

nn_SEX_5 <- nnet(SEX_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                       MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                     size = 4, data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_SEX_5 <- predict(nn_SEX_5, newdata = AMELIA_wide_start_nnet_5, type = "raw")
p_nn_SEX_5_2 <- 1 - p_nn_SEX_5
p_nn_SEX_5 <- cbind(p_nn_SEX_5, p_nn_SEX_5_2)
colnames(p_nn_SEX_5) <- c(2, 1)
SEX_5 <- apply(p_nn_SEX_5, 1, function(x)sample(colnames(p_nn_SEX_5), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, SEX_5 := as.factor(SEX_5)]

rf_SEX_5 <- randomForest(SEX_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                           MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                   data = AMELIA_PID5_wide_sample)
p_rf_SEX_5 <- predict(rf_SEX_5, newdata = AMELIA_wide_start_rf_5, type = "prob")
SEX_5 <- apply(p_rf_SEX_5, 1, function(x)sample(colnames(p_rf_SEX_5), size = 1, prob = x))
AMELIA_wide_start_rf_5[, SEX_5 := as.factor(SEX_5)] 

logit_SEX_5 <- multinom(SEX_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                   MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                   data = AMELIA_PID5_wide_sample)
p_logit_SEX_5 <- predict(logit_SEX_5, newdata = AMELIA_wide_start_logit_5, type = "prob")
p_logit_SEX_5_2 <- 1 - p_logit_SEX_5
p_logit_SEX_5 <- cbind(p_logit_SEX_5, p_logit_SEX_5_2)
colnames(p_logit_SEX_5) <- c(2, 1)
SEX_5 <- apply(p_logit_SEX_5, 1, function(x)sample(colnames(p_logit_SEX_5), 
size = 1, prob = x))
AMELIA_wide_start_logit_5[, SEX_5 := as.factor(SEX_5)]

nn_MST_5 <- nnet(MST_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                   MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                 size = 4, data = AMELIA_PID5_wide_sample, maxit = 500)
p_nn_MST_5 <- predict(nn_MST_5, newdata = AMELIA_wide_start_nnet_5, type = "raw")
MST_5 <- apply(p_nn_MST_5, 1, function(x)sample(colnames(p_nn_MST_5), size = 1, prob = x))
AMELIA_wide_start_nnet_5[, MST_5 := as.factor(MST_5)]

rf_MST_5 <- randomForest(MST_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                           MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                   data = AMELIA_PID5_wide_sample)
p_rf_MST_5 <- predict(rf_MST_5, newdata = AMELIA_wide_start_rf_5, type = "prob")
MST_5_rf <- apply(p_rf_MST_5, 1, function(x)sample(colnames(p_rf_MST_5), size = 1, prob = x))
AMELIA_wide_start_rf_5[, MST_5 := as.factor(MST_5_rf)] 

logit_MST_5 <- multinom(MST_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                           MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                         data = AMELIA_PID5_wide_sample)
p_logit_MST_5 <- predict(rf_MST_5, newdata = AMELIA_wide_start_logit_5, type = "prob")
MST_5 <- apply(p_logit_MST_5, 1, function(x)sample(colnames(p_logit_MST_5), size = 1, prob = x))
AMELIA_wide_start_logit_5[, MST_5 := as.factor(MST_5)] 

##############################
##############################
#############################6
nn_AGE_cat_2 <- nnet(AGE_cat_2 ~ AGE_1 + AGE_cat_1 + MST_1 + SEX_1, size = 4, 
                     data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_AGE_cat_2 <- predict(nn_AGE_cat_2, newdata = AMELIA_wide_start_nnet_6, type = "raw")
AGE_cat_2 <- apply(p_nn_AGE_cat_2, 1, function(x)sample(colnames(p_nn_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, AGE_cat_2 := as.factor(AGE_cat_2)] 

rf_AGE_cat_2 <- randomForest(AGE_cat_2 ~ AGE_cat_1 + MST_1 + SEX_1, 
                             data = AMELIA_PID6_wide_sample)
p_rf_AGE_cat_2 <- predict(rf_AGE_cat_2, newdata = AMELIA_wide_start_rf_6, type = "prob")
AGE_cat_2 <- apply(p_rf_AGE_cat_2, 1, function(x)sample(colnames(p_rf_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_rf_6[, AGE_cat_2 := as.factor(AGE_cat_2)] 

logit_AGE_cat_2 <- multinom(AGE_cat_2 ~ AGE_cat_1 + MST_1 + SEX_1, 
                             data = AMELIA_PID6_wide_sample)
p_logit_AGE_cat_2 <- predict(logit_AGE_cat_2, newdata = AMELIA_wide_start_logit_6, type = "prob")
AGE_cat_2 <- apply(p_logit_AGE_cat_2, 1, function(x)sample(colnames(p_logit_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_logit_6[, AGE_cat_2 := as.factor(AGE_cat_2)] 

nn_SEX_2 <- nnet(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, size = 4, 
                 data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_SEX_2 <- predict(nn_SEX_2, newdata = AMELIA_wide_start_nnet_6, type = "raw")
p_nn_SEX_2_2 <- 1 - p_nn_SEX_2
p_nn_SEX_2 <- cbind(p_nn_SEX_2, p_nn_SEX_2_2)
colnames(p_nn_SEX_2) <- c(2, 1)
SEX_2 <- apply(p_nn_SEX_2, 1, function(x)sample(colnames(p_nn_SEX_2), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, SEX_2 := as.factor(SEX_2)] 

rf_SEX_2 <- randomForest(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, data = AMELIA_PID6_wide_sample)
p_rf_SEX_2 <- predict(rf_SEX_2, newdata = AMELIA_wide_start_rf_6, type = "prob")
SEX_2 <- apply(p_rf_SEX_2, 1, function(x)sample(colnames(p_rf_SEX_2), size = 1, prob = x))
AMELIA_wide_start_rf_6[, SEX_2 := as.factor(SEX_2)] 

logit_SEX_2 <- multinom(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1, 
                 data = AMELIA_PID6_wide_sample)
p_logit_SEX_2 <- predict(logit_SEX_2, newdata = AMELIA_wide_start_logit_6, 
                         type = "prob")
p_logit_SEX_2_2 <- 1 - p_logit_SEX_2
p_logit_SEX_2 <- cbind(p_logit_SEX_2, p_logit_SEX_2_2)
colnames(p_logit_SEX_2) <- c(2, 1)
SEX_2 <- apply(p_logit_SEX_2, 1, function(x)sample(colnames(p_logit_SEX_2), size = 1, prob = x))
AMELIA_wide_start_logit_6[, SEX_2 := as.factor(SEX_2)] 

nn_MST_2 <- nnet(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, size = 4,
                 data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_MST_2 <- predict(nn_MST_2, newdata = AMELIA_wide_start_nnet_6, type = "raw")
MST_2 <- apply(p_nn_MST_2, 1, function(x)sample(colnames(p_nn_MST_2), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, MST_2 := as.factor(MST_2)] 

rf_MST_2 <- randomForest(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, 
                         data = AMELIA_PID6_wide_sample)
p_rf_MST_2 <- predict(rf_MST_2, newdata = AMELIA_wide_start_rf_6, type = "prob")
MST_2_rf <- apply(p_rf_MST_2, 1, function(x)sample(colnames(p_rf_MST_2), size = 1, prob = x))
AMELIA_wide_start_rf_6[, MST_2 := as.factor(MST_2_rf)] 

logit_MST_2 <- multinom(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2, 
                         data = AMELIA_PID6_wide_sample)
p_logit_MST_2 <- predict(logit_MST_2, newdata = AMELIA_wide_start_logit_6, type = "prob")
MST_2 <- apply(p_logit_MST_2, 1, function(x)sample(colnames(p_logit_MST_2), size = 1, prob = x))
AMELIA_wide_start_logit_6[, MST_2 := as.factor(MST_2)] 

nn_AGE_cat_3 <- nnet(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                     size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_AGE_cat_3 <- predict(nn_AGE_cat_3, newdata = AMELIA_wide_start_nnet_6, type = "raw")
AGE_cat_3 <- apply(p_nn_AGE_cat_3, 1, function(x)sample(colnames(p_nn_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, AGE_cat_3 := as.factor(AGE_cat_3)] 

rf_AGE_cat_3 <- randomForest(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                             data = AMELIA_PID6_wide_sample)
p_rf_AGE_cat_3 <- predict(rf_AGE_cat_3, newdata = AMELIA_wide_start_rf_6, type = "prob")
AGE_cat_3 <- apply(p_rf_AGE_cat_3, 1, function(x)sample(colnames(p_rf_AGE_cat_3), size = 1, prob = x))
AMELIA_wide_start_rf_6[, AGE_cat_3 := as.factor(AGE_cat_3)] 

logit_AGE_cat_3 <- multinom(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                             data = AMELIA_PID6_wide_sample)
p_logit_AGE_cat_3 <- predict(logit_AGE_cat_3, newdata = AMELIA_wide_start_logit_6, type = "prob")
AGE_cat_3 <- apply(p_logit_AGE_cat_3, 1, function(x)sample(colnames(p_logit_AGE_cat_3), size = 1, prob = x))
AMELIA_wide_start_logit_6[, AGE_cat_3 := as.factor(AGE_cat_3)] 

nn_SEX_3 <- nnet(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                 size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_SEX_3 <- predict(nn_SEX_3, newdata = AMELIA_wide_start_nnet_6, type = "raw")
p_nn_SEX_3_2 <- 1 - p_nn_SEX_3
p_nn_SEX_3 <- cbind(p_nn_SEX_3, p_nn_SEX_3_2)
colnames(p_nn_SEX_3) <- c(2, 1)
SEX_3 <- apply(p_nn_SEX_3, 1, function(x)sample(colnames(p_nn_SEX_3), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, SEX_3 := as.factor(SEX_3)]

rf_SEX_3 <- randomForest(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                         data = AMELIA_PID6_wide_sample)
p_rf_SEX_3 <- predict(rf_SEX_3, newdata = AMELIA_wide_start_rf_6, type = "prob")
SEX_3 <- apply(p_rf_SEX_3, 1, function(x)sample(colnames(p_rf_SEX_3), size = 1, prob = x))
AMELIA_wide_start_rf_6[, SEX_3 := as.factor(SEX_3)] 

logit_SEX_3 <- multinom(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2, 
                 data = AMELIA_PID6_wide_sample)
p_logit_SEX_3 <- predict(logit_SEX_3, newdata = AMELIA_wide_start_logit_6, type = "prob")
p_logit_SEX_3_2 <- 1 - p_logit_SEX_3
p_logit_SEX_3 <- cbind(p_logit_SEX_3, p_logit_SEX_3_2)
colnames(p_logit_SEX_3) <- c(2, 1)
SEX_3 <- apply(p_logit_SEX_3, 1, function(x)sample(colnames(p_logit_SEX_3), size = 1, prob = x))
AMELIA_wide_start_logit_6[, SEX_3 := as.factor(SEX_3)]

nn_MST_3 <- nnet(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                 size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_MST_3 <- predict(nn_MST_3, newdata = AMELIA_wide_start_nnet_6, type = "raw")
MST_3 <- apply(p_nn_MST_3, 1, function(x)sample(colnames(p_nn_MST_3), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, MST_3 := as.factor(MST_3)]

rf_MST_3 <- randomForest(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                         data = AMELIA_PID6_wide_sample)
p_rf_MST_3 <- predict(rf_MST_3, newdata = AMELIA_wide_start_rf_6, type = "prob")
MST_3_rf <- apply(p_rf_MST_3, 1, function(x)sample(colnames(p_rf_MST_3), size = 1, prob = x))
AMELIA_wide_start_rf_6[, MST_3 := as.factor(MST_3_rf)] 

logit_MST_3 <- multinom(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3, 
                         data = AMELIA_PID6_wide_sample)
p_logit_MST_3 <- predict(logit_MST_3, newdata = AMELIA_wide_start_logit_6, type = "prob")
MST_3 <- apply(p_logit_MST_3, 1, function(x)sample(colnames(p_logit_MST_3), size = 1, prob = x))
AMELIA_wide_start_logit_6[, MST_3 := as.factor(MST_3)] 

nn_AGE_cat_4 <- nnet(AGE_cat_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + MST_3 +
                       SEX_1 + SEX_2 + SEX_3, 
                     size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_AGE_cat_4 <- predict(nn_AGE_cat_4, newdata = AMELIA_wide_start_nnet_6, type = "raw")
AGE_cat_4 <- apply(p_nn_AGE_cat_4, 1, function(x)sample(colnames(p_nn_AGE_cat_4), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, AGE_cat_4 := as.factor(AGE_cat_4)]

rf_AGE_cat_4 <- randomForest(AGE_cat_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + MST_3 +
                               SEX_1 + SEX_2 +  SEX_3, 
                             data = AMELIA_PID6_wide_sample)
p_rf_AGE_cat_4 <- predict(rf_AGE_cat_4, newdata = AMELIA_wide_start_rf_6, type = "prob")
AGE_cat_4 <- apply(p_rf_AGE_cat_4, 1, function(x)sample(colnames(p_rf_AGE_cat_4), size = 1, prob = x))
AMELIA_wide_start_rf_6[, AGE_cat_4 := as.factor(AGE_cat_4)] 

logit_AGE_cat_4 <- multinom(AGE_cat_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + MST_3 +
                               SEX_1 + SEX_2 +  SEX_3, 
                             data = AMELIA_PID6_wide_sample)
p_logit_AGE_cat_4 <- predict(logit_AGE_cat_4, newdata = AMELIA_wide_start_logit_6, type = "prob")
AGE_cat_4 <- apply(p_logit_AGE_cat_4, 1, function(x)sample(colnames(p_logit_AGE_cat_4), size = 1, prob = x))
AMELIA_wide_start_logit_6[, AGE_cat_4 := as.factor(AGE_cat_4)] 

nn_SEX_4 <- nnet(SEX_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                   SEX_1 + SEX_2 +  SEX_3, 
                 size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_SEX_4 <- predict(nn_SEX_4, newdata = AMELIA_wide_start_nnet_6, type = "raw")
p_nn_SEX_4_2 <- 1 - p_nn_SEX_4
p_nn_SEX_4 <- cbind(p_nn_SEX_4, p_nn_SEX_4_2)
colnames(p_nn_SEX_4) <- c(2, 1)
SEX_4 <- apply(p_nn_SEX_4, 1, function(x)sample(colnames(p_nn_SEX_4), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, SEX_4 := as.factor(SEX_4)]

rf_SEX_4 <- randomForest(SEX_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                           SEX_1 + SEX_2 +  SEX_3, 
                         data = AMELIA_PID6_wide_sample)
p_rf_SEX_4 <- predict(rf_SEX_4, newdata = AMELIA_wide_start_rf_6, type = "prob")
SEX_4 <- apply(p_rf_SEX_4, 1, function(x)sample(colnames(p_rf_SEX_4), size = 1, prob = x))
AMELIA_wide_start_rf_6[, SEX_4 := as.factor(SEX_4)] 

logit_SEX_4 <- multinom(SEX_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                   SEX_1 + SEX_2 +  SEX_3, 
                   data = AMELIA_PID6_wide_sample)
p_logit_SEX_4 <- predict(logit_SEX_4, newdata = AMELIA_wide_start_logit_6, type = "prob")
p_logit_SEX_4_2 <- 1 - p_logit_SEX_4
p_logit_SEX_4 <- cbind(p_logit_SEX_4, p_logit_SEX_4_2)
colnames(p_logit_SEX_4) <- c(2, 1)
SEX_4 <- apply(p_logit_SEX_4, 1, function(x)sample(colnames(p_logit_SEX_4), size = 1, prob = x))
AMELIA_wide_start_logit_6[, SEX_4 := as.factor(SEX_4)]

nn_MST_4 <- nnet(MST_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                   SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                 size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_MST_4 <- predict(nn_MST_4, newdata = AMELIA_wide_start_nnet_6, type = "raw")
MST_4 <- apply(p_nn_MST_4, 1, function(x)sample(colnames(p_nn_MST_4), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, MST_4 := as.factor(MST_4)]

rf_MST_4 <- randomForest(MST_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                           SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                         data = AMELIA_PID6_wide_sample)
p_rf_MST_4 <- predict(rf_MST_4, newdata = AMELIA_wide_start_rf_6, type = "prob")
MST_4_rf <- apply(p_rf_MST_4, 1, function(x)sample(colnames(p_rf_MST_4), size = 1, prob = x))
AMELIA_wide_start_rf_6[, MST_4 := as.factor(MST_4_rf)] 

logit_MST_4 <- multinom(MST_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                           SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                         data = AMELIA_PID6_wide_sample)
p_logit_MST_4 <- predict(logit_MST_4, newdata = AMELIA_wide_start_logit_6, type = "prob")
MST_4 <- apply(p_logit_MST_4, 1, function(x)sample(colnames(p_logit_MST_4), size = 1, prob = x))
AMELIA_wide_start_logit_6[, MST_4 := as.factor(MST_4)] 

nn_AGE_cat_5 <- nnet(AGE_cat_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                       MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                     size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_AGE_cat_5 <- predict(nn_AGE_cat_5, newdata = AMELIA_wide_start_nnet_6, type = "raw")
AGE_cat_5 <- apply(p_nn_AGE_cat_5, 1, function(x)sample(colnames(p_nn_AGE_cat_5), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, AGE_cat_5 := as.factor(AGE_cat_5)]

rf_AGE_cat_5 <- randomForest(AGE_cat_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                               MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                             data = AMELIA_PID6_wide_sample)
p_rf_AGE_cat_5 <- predict(rf_AGE_cat_5, newdata = AMELIA_wide_start_rf_6, type = "prob")
AGE_cat_5 <- apply(p_rf_AGE_cat_5, 1, function(x)sample(colnames(p_rf_AGE_cat_5), size = 1, prob = x))
AMELIA_wide_start_rf_6[, AGE_cat_5 := as.factor(AGE_cat_5)] 

logit_AGE_cat_5 <- multinom(AGE_cat_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                               MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                             data = AMELIA_PID6_wide_sample)
p_logit_AGE_cat_5 <- predict(logit_AGE_cat_5, newdata = AMELIA_wide_start_logit_6, type = "prob")
AGE_cat_5 <- apply(p_logit_AGE_cat_5, 1, function(x)sample(colnames(p_logit_AGE_cat_5), size = 1, prob = x))
AMELIA_wide_start_logit_6[, AGE_cat_5 := as.factor(AGE_cat_5)] 

nn_SEX_5 <- nnet(SEX_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                   MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                 size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_SEX_5 <- predict(nn_SEX_5, newdata = AMELIA_wide_start_nnet_6, type = "raw")
p_nn_SEX_5_2 <- 1 - p_nn_SEX_5
p_nn_SEX_5 <- cbind(p_nn_SEX_5, p_nn_SEX_5_2)
colnames(p_nn_SEX_5) <- c(2, 1)
SEX_5 <- apply(p_nn_SEX_5, 1, function(x)sample(colnames(p_nn_SEX_5), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, SEX_5 := as.factor(SEX_5)]

rf_SEX_5 <- randomForest(SEX_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 +
                           MST_1 + MST_2 + MST_3 +
                           MST_4 + SEX_1 + SEX_2 + SEX_3 + SEX_4, 
                         data = AMELIA_PID6_wide_sample)
p_rf_SEX_5 <- predict(rf_SEX_5, newdata = AMELIA_wide_start_rf_6, type = "prob")
SEX_5 <- apply(p_rf_SEX_5, 1, function(x)sample(colnames(p_rf_SEX_5), size = 1, prob = x))
AMELIA_wide_start_rf_6[, SEX_5 := as.factor(SEX_5)] 

logit_SEX_5 <- multinom(SEX_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                   MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4, 
                   data = AMELIA_PID6_wide_sample)
p_logit_SEX_5 <- predict(logit_SEX_5, newdata = AMELIA_wide_start_logit_6, type = "prob")
p_logit_SEX_5_2 <- 1 - p_logit_SEX_5
p_logit_SEX_5 <- cbind(p_logit_SEX_5, p_logit_SEX_5_2)
colnames(p_logit_SEX_5) <- c(2, 1)
SEX_5 <- apply(p_logit_SEX_5, 1, function(x)sample(colnames(p_logit_SEX_5), size = 1, prob = x))
AMELIA_wide_start_logit_6[, SEX_5 := as.factor(SEX_5)]

nn_MST_5 <- nnet(MST_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                   MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                 size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_MST_5 <- predict(nn_MST_5, newdata = AMELIA_wide_start_nnet_6, type = "raw")
MST_5 <- apply(p_nn_MST_5, 1, function(x)sample(colnames(p_nn_MST_5), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, MST_5 := as.factor(MST_5)]

rf_MST_5 <- randomForest(MST_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                           MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                         data = AMELIA_PID6_wide_sample)
p_rf_MST_5 <- predict(rf_MST_5, newdata = AMELIA_wide_start_rf_6, type = "prob")
MST_5_rf <- apply(p_rf_MST_5, 1, function(x)sample(colnames(p_rf_MST_5), size = 1, prob = x))
AMELIA_wide_start_rf_6[, MST_5 := as.factor(MST_5_rf)] 

logit_MST_5 <- multinom(MST_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                           MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                         data = AMELIA_PID6_wide_sample)
p_logit_MST_5 <- predict(logit_MST_5, newdata = AMELIA_wide_start_logit_6, type = "prob")
MST_5 <- apply(p_logit_MST_5, 1, function(x)sample(colnames(p_logit_MST_5), size = 1, prob = x))
AMELIA_wide_start_logit_6[, MST_5 := as.factor(MST_5)] 

nn_AGE_cat_6 <- nnet(AGE_cat_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 +
                       MST_1 + MST_2 + MST_3 +
                       MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                     size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_AGE_cat_6 <- predict(nn_AGE_cat_6, newdata = AMELIA_wide_start_nnet_6, type = "raw")
AGE_cat_6 <- apply(p_nn_AGE_cat_6, 1, function(x)sample(colnames(p_nn_AGE_cat_6), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, AGE_cat_6 := as.factor(AGE_cat_6)]

rf_AGE_cat_6 <- randomForest(AGE_cat_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 +
                               MST_1 + MST_2 + MST_3 +
                               MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                             data = AMELIA_PID6_wide_sample)
p_rf_AGE_cat_6 <- predict(rf_AGE_cat_6, newdata = AMELIA_wide_start_rf_6, type = "prob")
AGE_cat_6 <- apply(p_rf_AGE_cat_6, 1, function(x)sample(colnames(p_rf_AGE_cat_6), size = 1, prob = x))
AMELIA_wide_start_rf_6[, AGE_cat_6 := as.factor(AGE_cat_6)] 

logit_AGE_cat_6 <- multinom(AGE_cat_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 +
                               MST_1 + MST_2 + MST_3 +
                               MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                             data = AMELIA_PID6_wide_sample)
p_logit_AGE_cat_6 <- predict(logit_AGE_cat_6, newdata = AMELIA_wide_start_logit_6, type = "prob")
AGE_cat_6 <- apply(p_logit_AGE_cat_6, 1, function(x)sample(colnames(p_logit_AGE_cat_6), size = 1, prob = x))
AMELIA_wide_start_logit_6[, AGE_cat_6 := as.factor(AGE_cat_6)] 

nn_SEX_6 <- nnet(SEX_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + AGE_cat_6 +
                       MST_1 + MST_2 + MST_3 +
                       MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                     size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_SEX_6 <- predict(nn_SEX_6, newdata = AMELIA_wide_start_nnet_6, type = "raw")
p_nn_SEX_6_2 <- 1 - p_nn_SEX_6
p_nn_SEX_6 <- cbind(p_nn_SEX_6, p_nn_SEX_6_2)
colnames(p_nn_SEX_6) <- c(2, 1)
SEX_6 <- apply(p_nn_SEX_6, 1, function(x)sample(colnames(p_nn_SEX_6), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, SEX_6 := as.factor(SEX_6)]

rf_SEX_6 <- randomForest(SEX_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 +
                           AGE_cat_6 + MST_1 + MST_2 + MST_3 +
                           MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                         data = AMELIA_PID6_wide_sample)
p_rf_SEX_6 <- predict(rf_SEX_6, newdata = AMELIA_wide_start_rf_6, type = "prob")
SEX_6 <- apply(p_rf_SEX_6, 1, function(x)sample(colnames(p_rf_SEX_6), size = 1, prob = x))
AMELIA_wide_start_rf_6[, SEX_6 := as.factor(SEX_6)] 

logit_SEX_6 <- multinom(SEX_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + AGE_cat_6 +
                   MST_1 + MST_2 + MST_3 +
                   MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5, 
                   data = AMELIA_PID6_wide_sample)
p_logit_SEX_6 <- predict(logit_SEX_6, newdata = AMELIA_wide_start_logit_6, type = "prob")
p_logit_SEX_6_2 <- 1 - p_logit_SEX_6
p_logit_SEX_6 <- cbind(p_logit_SEX_6, p_logit_SEX_6_2)
colnames(p_logit_SEX_6) <- c(2, 1)
SEX_6 <- apply(p_logit_SEX_6, 1, function(x)sample(colnames(p_logit_SEX_6), size = 1, prob = x))
AMELIA_wide_start_logit_6[, SEX_6 := as.factor(SEX_6)]

nn_MST_6 <- nnet(MST_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + AGE_cat_6 +
                   MST_1 + MST_2 + MST_3 +
                   MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5 + SEX_6, 
                 size = 4, data = AMELIA_PID6_wide_sample, maxit = 500)
p_nn_MST_6 <- predict(nn_MST_6, newdata = AMELIA_wide_start_nnet_6, type = "raw")
MST_6 <- apply(p_nn_MST_6, 1, function(x)sample(colnames(p_nn_MST_6), size = 1, prob = x))
AMELIA_wide_start_nnet_6[, MST_6 := as.factor(MST_6)]

rf_MST_6 <- randomForest(MST_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + AGE_cat_6 +
                           MST_1 + MST_2 + MST_3 +
                           MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5 + SEX_6, 
                         data = AMELIA_PID6_wide_sample)
p_rf_MST_6 <- predict(rf_MST_6, newdata = AMELIA_wide_start_rf_6, type = "prob")
MST_6_rf <- apply(p_rf_MST_6, 1, function(x)sample(colnames(p_rf_MST_6), size = 1, prob = x))
AMELIA_wide_start_rf_6[, MST_6 := as.factor(MST_6_rf)] 

logit_MST_6 <- multinom(MST_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + AGE_cat_6 +
                           MST_1 + MST_2 + MST_3 +
                           MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5 + SEX_6, 
                         data = AMELIA_PID6_wide_sample)
p_logit_MST_6 <- predict(logit_MST_6, newdata = AMELIA_wide_start_logit_6, type = "prob")
MST_6 <- apply(p_logit_MST_6, 1, function(x)sample(colnames(p_logit_MST_6), size = 1, prob = x))
AMELIA_wide_start_logit_6[, MST_6 := as.factor(MST_6)] 
#########################################
#########################################
# Generate categories to create a graph like the one in Kolb (2013). Therefore,
# let´s create a counter for certain relative frequencies which are based on connected
# attributes over several variables

colb_creater_2 <- function(x) {
# 1: All older then 18, mixed gender, both married
x$def_1 <- 0
x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1 
                                  & x$MST_group == "22" & x$SEX_diff == "12"] <- 1
# 2: All older then 18, same gender, both married
x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1 
        & x$MST_group == "22" & x$SEX_diff != "12"] <- 2
# 3: All older then 18, mixed gender, both never married
x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1 
        & x$MST_group == "11" & x$SEX_diff == "12"] <- 3
# 4: One person below 18, mixed gender, both never married
x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) < 2 
        & x$MST_group == "11" & x$SEX_diff == "12"] <- 4
# 5: One person below 18, same gender, both never married
x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) < 2 
        & x$MST_group == "11" & x$SEX_diff != "12"] <- 5
# 6: One person below 18, mixed gender, one divorced
x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) < 2 
        & x$MST_group == "15" & x$SEX_diff == "12"] <- 6
# 7: One person below 18, same gender , one divorced
x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) < 2 
        & x$MST_group == "15" & x$SEX_diff != "12"] <- 7
# 8: One person below 18, same gender, one seperated
x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) < 2 
        & x$MST_group == "14" & x$SEX_diff != "12"] <- 8
return(x)
}

colb_creater_3 <- function(x) {
  # 1 All older then 18, mixed gender, one married, rest never married
  x$def_1 <- 0
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1 &
            as.numeric(x$AGE_cat_3) > 1 & x$SEX_diff != "111" & x$SEX_diff != "222" &
            x$MST_group == "112"] <- 1
  #2 One person below 18 years, mixed gender, one married, rest never married
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1 &
            as.numeric(x$AGE_cat_3) < 2 & x$SEX_diff != "111" & x$SEX_diff != "222" &
            x$MST_group == "112"] <- 2
  #3 All older then 18, mixed gender, two married, rest never married
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1 &
            as.numeric(x$AGE_cat_3) > 1 & x$SEX_diff != "111" & x$SEX_diff != "222" &
            x$MST_group == "212"] <- 3
  #4 One person below 18 years, mixed gender, two married, rest never married
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1 &
            as.numeric(x$AGE_cat_3) < 2 & x$SEX_diff != "111" & x$SEX_diff != "222" &
            x$MST_group == "212"] <- 4
  #5 Two persons below 18 years, mixed gender, one married, rest never married
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) < 2 &
            as.numeric(x$AGE_cat_3) < 2 & x$SEX_diff != "111" & x$SEX_diff != "222" &
            x$MST_group == "112"] <- 5
  #6 All older then 18, mixed gender, all married
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1 &
            as.numeric(x$AGE_cat_3) > 1 & x$SEX_diff != "111" & x$SEX_diff != "222" &
            x$MST_group == "222"] <- 6
  #7 Only men, all older then 18
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1 &
            as.numeric(x$AGE_cat_3) > 1 & x$SEX_diff == "111"] <- 7
  # 8 Only woman, all older then 18
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1 &
            as.numeric(x$AGE_cat_3) > 1 & x$SEX_diff == "222"] <- 8
  return(x)
}

colb_creater_4 <- function(x) {
  # Two men, two woman, all at least 18 years old
  x$def_1 <- 0
  x$def_1[as.numeric(x$AGE_cat_1) > 1  & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
            x$SEX_diff == "2112"] <- 1
  # Three man, one woman, all at least 18 years old
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
            x$SEX_diff == "1112"] <- 2
  # Only men or only woman, all at least 18 years old
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
            x$SEX_diff == "1111"|as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_1) > 1 & 
            as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 
          & x$SEX_diff == "2222"] <- 3
  # Three adults one person below 17 years
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) < 2 ] <- 4
  # One adult and three persons below 17 years
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) < 2
          & as.numeric(x$AGE_cat_3) < 2 & as.numeric(x$AGE_cat_4) < 2 ] <- 5
  # Two adults and two persons below 17 years
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) < 2 & as.numeric(x$AGE_cat_4) < 2 ] <- 6
  # All persons below 18 years
  x$def_1[as.numeric(x$AGE_cat_1) < 2 & as.numeric(x$AGE_cat_2) > 2
          & as.numeric(x$AGE_cat_3) < 2 & as.numeric(x$AGE_cat_4) < 2 ] <- 7
  return(x)
}

colb_creater_5 <- function(x) {
  #1 Only men, all older then 17
  x$def_1 <- 0
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 & as.numeric(x$AGE_cat_5)  > 1 &
            x$SEX_diff == "11111"] <- 1
  #2 4 Men, 1 Woman, all older then 17
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 & as.numeric(x$AGE_cat_5)  > 1 &
            x$SEX_diff == "11112"] <- 2
  #3 3 Men, 2 Woman, all older then 17
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 & as.numeric(x$AGE_cat_5)  > 1 &
            x$SEX_diff == "21112"] <- 3
  #4 2 Men, two Woman, all older then 17
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 & as.numeric(x$AGE_cat_5)  > 1 &
            x$SEX_diff == "22112"] <- 4
  #5 Only woman, all older then 17
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 & as.numeric(x$AGE_cat_5)  > 1 &
            x$SEX_diff == "22222"] <- 5
  #6 All older then 17, expect for one person
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
            as.numeric(x$AGE_cat_5) < 2] <- 6
  #7 All older then 17, exept for two persons
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) < 2 &
            as.numeric(x$AGE_cat_5) < 2] <- 7
  #8 Two persons older then 17, three below 17
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) < 2 & as.numeric(x$AGE_cat_4) < 2 &
            as.numeric(x$AGE_cat_5) < 2] <- 8
  return(x)
}

colb_creater_6 <- function(x) {
  x$def_1 <- 0
  #1 Five men, one woman, all at least 18
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
            & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
              as.numeric(x$AGE_cat_5) > 1 & as.numeric(x$AGE_cat_6) & x$SEX_diff == "111112"] <- 1
  #2 Four men, two man, all at least 18
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
            as.numeric(x$AGE_cat_5) > 1 & as.numeric(x$AGE_cat_6) & x$SEX_diff == "211112"] <- 2
  #3 Three man, three woman, all at least 18
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
            as.numeric(x$AGE_cat_5) > 1 & as.numeric(x$AGE_cat_6) & x$SEX_diff == "221112"] <- 3
  #4 Five man, one woman, one person below 18 years
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
  & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
    as.numeric(x$AGE_cat_5) > 1 & as.numeric(x$AGE_cat_6) < 2 & x$SEX_diff == "111112"] <- 4
  #5 Five man, one woman, two persons below 18 years
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
            as.numeric(x$AGE_cat_5) < 2 & as.numeric(x$AGE_cat_6) < 2 & x$SEX_diff == "111112"] <- 5
  #6 Five man, one woman, three persons below 18 years
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) < 2 &
            as.numeric(x$AGE_cat_5) < 2 & as.numeric(x$AGE_cat_6) < 2 & x$SEX_diff == "111112"] <- 6
  #7 Four man, two woman, three persons below 18 years
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) < 2 &
            as.numeric(x$AGE_cat_5) < 2 & as.numeric(x$AGE_cat_6) < 2 & x$SEX_diff == "211112"] <- 7
  #8 Four man, two woman, three persons below 18 years
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) < 2 & as.numeric(x$AGE_cat_4) < 2 &
            as.numeric(x$AGE_cat_5) < 2 & as.numeric(x$AGE_cat_6) < 2 & x$SEX_diff == "211112"] <- 8
return(x)
}

# Apply colb creator (needed for Figure 11)

hhs_2_true <- colb_creater_2(hhs_2_true)
synthetic_hhs2_logit_colb <- colb_creater_2(test_logit_2)
synthetic_hhs2_nn_colb <- colb_creater_2(test_nn_2)
synthetic_hhs2_rf_colb <- colb_creater_2(test_rf_2)
synthetic_hhs2_logit_mod_colb <- colb_creater_2(AMELIA_wide_start_logit_2)
synthetic_hhs2_nn_mod_colb <- colb_creater_2(AMELIA_wide_start_nnet_2)
synthetic_hhs2_rf_mod_colb <- colb_creater_2(AMELIA_wide_start_rf_2)


hhs_3_true <- colb_creater_3(hhs_3_true)
synthetic_hhs3_logit_colb <- colb_creater_3(test_logit_3)
synthetic_hhs3_nn_colb <- colb_creater_3(test_nn_3)
synthetic_hhs3_rf_colb <- colb_creater_3(test_rf_3)
synthetic_hhs3_logit_mod_colb <- colb_creater_3(AMELIA_wide_start_logit_3)
synthetic_hhs3_nn_mod_colb <- colb_creater_3(AMELIA_wide_start_nnet_3)
synthetic_hhs3_rf_mod_colb <- colb_creater_3(AMELIA_wide_start_rf_3)

hhs_4_true <- colb_creater_4(hhs_4_true)
synthetic_hhs4_logit_colb <- colb_creater_4(test_logit_4)
synthetic_hhs4_nn_colb <- colb_creater_4(test_nn_4)
synthetic_hhs4_rf_colb <- colb_creater_4(test_rf_4)
synthetic_hhs4_logit_mod_colb <- colb_creater_4(AMELIA_wide_start_logit_4)
synthetic_hhs4_nn_mod_colb <- colb_creater_4(AMELIA_wide_start_nnet_4)
synthetic_hhs4_rf_mod_colb <- colb_creater_4(AMELIA_wide_start_rf_4)

hhs_5_true <- colb_creater_5(hhs_5_true)
synthetic_hhs5_logit_colb <- colb_creater_5(test_logit_5)
synthetic_hhs5_nn_colb <- colb_creater_5(test_nn_5)
synthetic_hhs5_rf_colb <- colb_creater_5(test_rf_5)
synthetic_hhs5_logit_mod_colb <- colb_creater_5(AMELIA_wide_start_logit_5)
synthetic_hhs5_nn_mod_colb <- colb_creater_5(AMELIA_wide_start_nnet_5)
synthetic_hhs5_rf_mod_colb <- colb_creater_5(AMELIA_wide_start_rf_5)

hhs_6_true <- colb_creater_6(hhs_6_true)
synthetic_hhs6_logit_colb <- colb_creater_6(test_logit_6)
synthetic_hhs6_nn_colb <- colb_creater_6(test_nn_6)
synthetic_hhs6_rf_colb <- colb_creater_6(test_rf_6)
synthetic_hhs6_logit_mod_colb <- colb_creater_6(AMELIA_wide_start_logit_6)
synthetic_hhs6_nn_mod_colb <- colb_creater_6(AMELIA_wide_start_nnet_6)
synthetic_hhs6_rf_mod_colb <- colb_creater_6(AMELIA_wide_start_rf_6)

hhs2_t <- as.vector(prop.table(table(hhs_2_true$def_1)))
hhs2_c_log <- as.vector(prop.table(table(synthetic_hhs2_logit_colb$def_1)))
hhs2_c_nn <- as.vector(prop.table(table(synthetic_hhs2_nn_colb$def_1)))
hhs2_c_rf <- as.vector(prop.table(table(synthetic_hhs2_rf_colb$def_1)))
hhs2_m_log <- as.vector(prop.table(table(synthetic_hhs2_logit_mod_colb$def_1)))
hhs2_m_nn <- as.vector(prop.table(table(synthetic_hhs2_nn_mod_colb$def_1)))
hhs2_m_rf <- as.vector(prop.table(table(synthetic_hhs2_rf_mod_colb$def_1)))

colb_2 <- cbind(hhs2_t, hhs2_c_log, hhs2_c_nn, hhs2_c_rf, hhs2_m_log, 
                hhs2_m_nn, hhs2_m_rf)

colb_2 <- melt(colb_2)

###
# and so on...
# Give names for multivariate household type

#
colb_2$Var1 <- factor(colb_2$X1, labels = c("Other",
                                              "All older then 18 years, mixed gender, both married",
                                              "All older then 18 years, same gender, both married",
                                              "All older then 18 years, mixed gender, both never married",
                                              "One person below 18 years, mixed gender, both never married",
                                              "One person below 18 years, same gender, both never married",
                                              "One person below 18 years, mixed gender, one divorced",
                                              "One person below 18 years, same gender , one divorced",
                                              "One person below 18 years, same gender, one seperated"))

colb_3$Var1 <- factor(colb_3$X1, labels = c("Other",
                                              "All at least 18 years, mixed gender, one married, rest never married",
                                              "One person below 18 years, mixed gender, one married, rest never married",
                                              "All at least 18 years, mixed gender, two married, rest never married",
                                              "One person below 18 years, mixed gender, two married, rest never married",
                                              "Two persons below 18 years, mixed gender, one married, rest never married",
                                              "All at least 18 years, mixed gender, all married",
                                              "Only men, all older then 18 years",
                                              "Only woman, all older then 18 years"))

colb_4$Var1 <- factor(colb_4$X1, labels = c("Other",
                                              "Two men, two woman, all at least 18 years old",
                                              "Three man, one woman, all at least 18 years old",
                                              "Only men or only woman, all at least 18 years old",
                                              "Three adults one person below 17 years",
                                              "One adult and three persons below 17 years",
                                              "Two adults and two persons below 17 years",
                                              "All persons below 18 years"))

colb_5$Var1 <- factor(colb_5$X1, labels = c("Other",
                                              "Only men, all at least 18 years",
                                              "Four Men, one Woman, all at least then 18 years",
                                              "Three Men, two Woman, all at least 18 years",
                                              "Two Men, two Woman, all at least 18 years",
                                              "Only woman, all at least 18 years",
                                              "All older then 17 years, except for one person",
                                              "All older then 17 years, exept for two persons",
                                              "Two persons older then 17 years, three below 17"))
colb_6$Var1 <- factor(colb_6$X1, labels = c("Other",
                                              "Five men, one woman, all at least 18 years",
                                              "Four men, two man, all at least 18 years",
                                              "Three man, three woman, all at least 18 years",
                                              "Five man, one woman, one person below 18 years",
                                              "Five man, one woman, two persons below 18 years",
                                              "Five man, one woman, three persons below 18 years",
                                              "Four man, two woman, three persons below 18 years",
                                              "Four man, two woman, three persons below 18 years"))

# Create plots for figure 11
kolb_plot_2 <-  ggplot(colb_2, aes(y = value, x = X2,
                                   fill = Var1)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  theme(legend.position="bottom") +
  scale_fill_brewer(palette="Spectral") +
  xlab("Strategy and Method") +
  ylab("Frequency") +
  labs(fill="Household Type") +
  theme_bw()
kolb_plot_2

kolb_plot_3 <-  ggplot(colb_3, aes(y = value, x = X2,
                                   fill = Var1)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  theme(legend.position="bottom") +
  scale_fill_brewer(palette="Spectral") +
  xlab("Strategy and Method") +
  ylab("Frequency") +
  labs(fill="Household Type") +
  theme_bw()
kolb_plot_3

kolb_plot_4 <-  ggplot(colb_4, aes(y = value, x = X2,
                                   fill = Var1)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  theme(legend.position="bottom") +
  scale_fill_brewer(palette="Spectral") +
  xlab("Strategy and Method") +
  ylab("Frequency") +
  labs(fill="Household Type") +
  theme_bw()
kolb_plot_4

kolb_plot_5 <-  ggplot(colb_5, aes(y = value, x = X2,
                                   fill = Var1)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  theme(legend.position="bottom") +
  scale_fill_brewer(palette="Spectral") +
  xlab("Strategy and Method") +
  ylab("Frequency") +
  labs(fill="Household Type") +
  theme_bw()
kolb_plot_5

kolb_plot_6 <-  ggplot(colb_6, aes(y = value, x = X2,
                                   fill = Var1)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  theme(legend.position="bottom") +
  scale_fill_brewer(palette="Spectral") +
  xlab("Strategy and Method") +
  ylab("Frequency") +
  labs(fill="Household Type") +
  theme_bw()
kolb_plot_6

######################
#####################
intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
}

intersect_all(a,b,c)
# Number of unique Households in every synthetic population; compared with resampling uniquness
# See Table 5
#2 
unique(hhs_2_true[,c("AGE_cat_1", "AGE_cat_2", "SEX_1", "SEX_2",
                     "MST_1", "MST_2")]) 
nrow(unique(AMELIA_PID2_wide_sample[,c("AGE_cat_1", "AGE_cat_2", "SEX_1", "SEX_2",
                                       "MST_1", "MST_2")]))
nrow(unique(synthetic_hhs2_nn_colb[,c("AGE_cat_1", "AGE_cat_2", "SEX_1", "SEX_2",
                                      "MST_1", "MST_2")])) 
nrow(unique(synthetic_hhs2_rf_colb[,c("AGE_cat_1", "AGE_cat_2", "SEX_1", "SEX_2",
                                      "MST_1", "MST_2")])) 

###
# and so on..

# How often did a unique household appear? See Table 6
library(dplyr) # dplyr package needed for 'pipe operator' %>%
synthetic_hhs6_rf_colb %>% 
  select(AGE_cat_1, AGE_cat_2, AGE_cat_3, AGE_cat_4,
         AGE_cat_5, AGE_cat_6,
         SEX_1, SEX_2, SEX_3, SEX_4, SEX_5, SEX_6,
         MST_1, MST_2, MST_3, MST_4, MST_5, MST_6) %>%
  group_by(AGE_cat_1, AGE_cat_2, AGE_cat_3, AGE_cat_4,
           AGE_cat_5, AGE_cat_6,
           SEX_1, SEX_2, SEX_3, SEX_4, SEX_5, SEX_6,
           MST_1, MST_2, MST_3, MST_4, MST_5, MST_6) %>% 
  mutate(count = n()) %>%
  arrange(desc(count)) %>%
  unique() %>%
  View()

## 
# and so on...
  
# Now we create the mosaic plots (Figure 6-10)             

par(mfrow = c(1, 7))
####### Give pairwise comparison of Age for hhs2
#TRUE
title("True")
plot(hhs_2_true$AGE_cat_1, hhs_2_true$AGE_cat_2, xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))

# Sample
title("Sample")
plot(AMELIA_PID2_wide_sample$AGE_cat_1, AMELIA_PID2_wide_sample$AGE_cat_2, xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
#C-Log
title("C-LOG")
plot(synthetic_hhs2_logit_colb$AGE_cat_1, synthetic_hhs2_logit_colb$AGE_cat_2,
     xlab = "Age1", ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
#C-NN
title("C-NN")
plot(synthetic_hhs2_nn_colb$AGE_cat_1, synthetic_hhs2_nn_colb$AGE_cat_2,
     xlab = "Age1", ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
#C-RF
title("C-RF")
plot(synthetic_hhs2_rf_colb$AGE_cat_1, synthetic_hhs2_nn_colb$AGE_cat_2,
     xlab = "Age1", ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
#M-LOG
title("M-LOG")
plot(synthetic_hhs2_logit_mod_colb$AGE_cat_1, synthetic_hhs2_logit_mod_colb$AGE_cat_2,
     xlab = "Age1", ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
#M-NN  
title("M-NN")
plot(synthetic_hhs2_nn_mod_colb$AGE_cat_1, synthetic_hhs2_nn_mod_colb$AGE_cat_2,
     xlab = "Age1", ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
#M-RF 
title("M-RF")
plot(synthetic_hhs2_rf_mod_colb$AGE_cat_1, synthetic_hhs2_rf_mod_colb$AGE_cat_2,
     xlab = "Age1", ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))

par(mfrow = c(3, 7))
####### Give pairwise comparison of Age for hhs3
#TRUE
title("True")
plot(hhs_3_true$AGE_cat_1, hhs_3_true$AGE_cat_2, xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_3_true$AGE_cat_1, hhs_3_true$AGE_cat_3, xlab = "Age1",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_3_true$AGE_cat_2, hhs_3_true$AGE_cat_3, xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#C-LOG
title("C-LOG")
plot(synthetic_hhs3_logit_colb$AGE_cat_1, synthetic_hhs3_logit_colb$AGE_cat_2, 
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_colb$AGE_cat_1, synthetic_hhs3_logit_colb$AGE_cat_3,
     xlab = "Age1",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_colb$AGE_cat_2, synthetic_hhs3_logit_colb$AGE_cat_3, 
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#C-NN
title("C-NN")
plot(synthetic_hhs3_nn_colb$AGE_cat_1, synthetic_hhs3_nn_colb$AGE_cat_2, 
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_colb$AGE_cat_1, synthetic_hhs3_nn_colb$AGE_cat_3,
     xlab = "Age1",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_colb$AGE_cat_2, synthetic_hhs3_nn_colb$AGE_cat_3, 
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#C-RF
title("C-RF")
plot(synthetic_hhs3_rf_colb$AGE_cat_1, synthetic_hhs3_nn_colb$AGE_cat_2, 
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_colb$AGE_cat_1, synthetic_hhs3_nn_colb$AGE_cat_3, 
     xlab = "Age1",
     ylab = "Age3",col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_colb$AGE_cat_2, synthetic_hhs3_nn_colb$AGE_cat_3, 
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#M-LOG
title("M-LOG")
plot(synthetic_hhs3_logit_mod_colb$AGE_cat_1, synthetic_hhs3_logit_mod_colb$AGE_cat_2, 
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_mod_colb$AGE_cat_1, synthetic_hhs3_logit_mod_colb$AGE_cat_3, 
     xlab = "Age1",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_mod_colb$AGE_cat_2, synthetic_hhs3_logit_mod_colb$AGE_cat_3,
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#M-NN  
title("M-NN")
plot(synthetic_hhs3_nn_mod_colb$AGE_cat_1, synthetic_hhs3_nn_mod_colb$AGE_cat_2, 
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_mod_colb$AGE_cat_1, synthetic_hhs3_nn_mod_colb$AGE_cat_3, 
     xlab = "Age1",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_mod_colb$AGE_cat_2, synthetic_hhs3_nn_mod_colb$AGE_cat_3,
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#M-RF 
title("M-RF")
plot(synthetic_hhs3_rf_mod_colb$AGE_cat_1, synthetic_hhs3_rf_mod_colb$AGE_cat_2, 
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_mod_colb$AGE_cat_1, synthetic_hhs3_rf_mod_colb$AGE_cat_3,
     xlab = "Age1",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_mod_colb$AGE_cat_2, synthetic_hhs3_rf_mod_colb$AGE_cat_3, 
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#
# MST
####### Give pairwise comparison of Age for hhs2
#TRUE
title("True")
plot(hhs_2_true$SEX_1, hhs_2_true$SEX_2, xlab = "Sex1",
     ylab = "Sex2", col = brewer.pal(name = "Spectral", n = 10))

#C-Log
title("C-LOG")
plot(synthetic_hhs2_logit_colb$SEX_1, 
     synthetic_hhs2_logit_colb$SEX_2,
     xlab = "Sex1", ylab = "Sex2", col = brewer.pal(name = "Spectral", n = 10))
#C-NN
title("C-NN")
plot(synthetic_hhs2_nn_colb$SEX_1, 
     synthetic_hhs2_nn_colb$SEX_2,
     xlab = "Sex1", ylab = "Sex2", col = brewer.pal(name = "Spectral", n = 10))
#C-RF
title("C-RF")
plot(synthetic_hhs2_rf_colb$SEX_1,
     synthetic_hhs2_rf_colb$SEX_2,
     xlab = "Sex1", ylab = "Sex2", col = brewer.pal(name = "Spectral", n = 10))
#M-LOG
title("M-LOG")
plot(synthetic_hhs2_logit_mod_colb$SEX_1, 
     synthetic_hhs2_logit_mod_colb$SEX_2,
     xlab = "Sex1", ylab = "Sex2", col = brewer.pal(name = "Spectral", n = 10))
#M-NN  
title("M-NN")
plot(synthetic_hhs2_nn_mod_colb$SEX_1, 
     synthetic_hhs2_nn_mod_colb$SEX_2,
     xlab = "Sex1", ylab = "Sex2", col = brewer.pal(name = "Spectral", n = 10))
#M-RF 
title("M-RF")
plot(synthetic_hhs2_rf_mod_colb$SEX_1, 
     synthetic_hhs2_rf_mod_colb$SEX_2,
     xlab = "Sex1", ylab = "Sex2", col = brewer.pal(name = "Spectral", n = 10))

par(mfrow = c(3, 7))
####### Give pairwise comparison of Age for hhs3
#TRUE
title("True")
plot(hhs_3_true$MST_1, hhs_3_true$MST_2, xlab = "MST1",
     ylab = "MST2", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_3_true$MST_1, hhs_3_true$MST_3, xlab = "MST1",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_3_true$MST_2, hhs_3_true$MST_3, xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#C-LOG
title("C-LOG")
plot(synthetic_hhs3_logit_colb$MST_1, synthetic_hhs3_logit_colb$MST_2, 
     xlab = "MST1",
     ylab = "MST2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_colb$MST_1, synthetic_hhs3_logit_colb$MST_3,
     xlab = "MST1",
     ylab = "A3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_colb$AGE_cat_2, synthetic_hhs3_logit_colb$AGE_cat_3, 
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#C-NN
title("C-NN")
plot(synthetic_hhs3_nn_colb$AGE_cat_1, synthetic_hhs3_nn_colb$AGE_cat_2, 
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_colb$AGE_cat_1, synthetic_hhs3_nn_colb$AGE_cat_3,
     xlab = "Age1",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_colb$AGE_cat_2, synthetic_hhs3_nn_colb$AGE_cat_3, 
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#C-RF
title("C-RF")
plot(synthetic_hhs3_rf_colb$AGE_cat_1, synthetic_hhs3_nn_colb$AGE_cat_2, 
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_colb$AGE_cat_1, synthetic_hhs3_nn_colb$AGE_cat_3, 
     xlab = "Age1",
     ylab = "Age3",col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_colb$AGE_cat_2, synthetic_hhs3_nn_colb$AGE_cat_3, 
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#M-LOG
title("M-LOG")
plot(synthetic_hhs3_logit_mod_colb$AGE_cat_1, synthetic_hhs3_logit_mod_colb$AGE_cat_2, 
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_mod_colb$AGE_cat_1, synthetic_hhs3_logit_mod_colb$AGE_cat_3, 
     xlab = "Age1",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_mod_colb$AGE_cat_2, synthetic_hhs3_logit_mod_colb$AGE_cat_3,
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#M-NN  
title("M-NN")
plot(synthetic_hhs3_nn_mod_colb$AGE_cat_1, synthetic_hhs3_nn_mod_colb$AGE_cat_2, 
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_mod_colb$AGE_cat_1, synthetic_hhs3_nn_mod_colb$AGE_cat_3, 
     xlab = "Age1",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_mod_colb$AGE_cat_2, synthetic_hhs3_nn_mod_colb$AGE_cat_3,
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
#M-RF 
title("M-RF")
plot(synthetic_hhs3_rf_mod_colb$AGE_cat_1, synthetic_hhs3_rf_mod_colb$AGE_cat_2, 
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_mod_colb$AGE_cat_1, synthetic_hhs3_rf_mod_colb$AGE_cat_3,
     xlab = "Age1",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_mod_colb$AGE_cat_2, synthetic_hhs3_rf_mod_colb$AGE_cat_3, 
     xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))

####### Give pairwise comparison of Age for hhs4
#TRUE
title("True")
plot(hhs_4_true$AGE_cat_1, hhs_4_true$AGE_cat_2,
     xlab = "Age1",
     ylab = "Age2", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_4_true$AGE_cat_1, hhs_4_true$AGE_cat_3,xlab = "Age1",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_4_true$AGE_cat_1, hhs_4_true$AGE_cat_4, xlab = "Age1",
     ylab = "Age4", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_4_true$AGE_cat_2, hhs_4_true$AGE_cat_3, xlab = "Age2",
     ylab = "Age3", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_4_true$AGE_cat_2, hhs_4_true$AGE_cat_4, xlab = "Age2",
     ylab = "Age4", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_4_true$AGE_cat_3, hhs_4_true$AGE_cat_4, xlab = "Age3",
     ylab = "Age4", col = brewer.pal(name = "Spectral", n = 10))
#M-LOG 
title("M-LOG")
plot(synthetic_hhs4_logit_mod_colb$AGE_cat_1, 
     synthetic_hhs4_logit_mod_colb$AGE_cat_2,
     xlab = "Age1", ylab = "Age2", 
     col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs4_logit_mod_colb$AGE_cat_1, 
     synthetic_hhs4_logit_mod_colb$AGE_cat_3,
     xlab = "Age1", ylab = "Age3", 
     col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs4_logit_mod_colb$AGE_cat_1, 
     synthetic_hhs4_logit_mod_colb$AGE_cat_4,
     xlab = "Age1", ylab = "Age4", 
     col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs4_logit_mod_colb$AGE_cat_2, 
     synthetic_hhs4_logit_mod_colb$AGE_cat_3,
     xlab = "Age2", ylab = "Age3", 
     col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs4_logit_mod_colb$AGE_cat_2, 
     synthetic_hhs4_logit_mod_colb$AGE_cat_4,
     xlab = "Age2", ylab = "Age4", 
     col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs4_logit_mod_colb$AGE_cat_3, 
     synthetic_hhs4_logit_mod_colb$AGE_cat_4,
     xlab = "Age3", ylab = "Age4", 
     col = brewer.pal(name = "Spectral", n = 10))

#M-NN 
title("M-NN")
plot(AMELIA_wide_start_nnet_4$AGE_cat_1, 
     AMELIA_wide_start_nnet_4$AGE_cat_2,
     xlab = "Age1", ylab = "Age2", 
     col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs4_nn_mod_colb$AGE_cat_1, 
     synthetic_hhs4_nn_mod_colb$AGE_cat_3,
     xlab = "Age1", ylab = "Age3", 
     col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs4_nn_mod_colb$AGE_cat_1, 
     synthetic_hhs4_nn_mod_colb$AGE_cat_4,
     xlab = "Age1", ylab = "Age4", 
     col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs4_nn_mod_colb$AGE_cat_2, 
     synthetic_hhs4_nn_mod_colb$AGE_cat_3,
     xlab = "Age2", ylab = "Age3", 
     col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs4_nn_mod_colb$AGE_cat_2, 
     synthetic_hhs4_nn_mod_colb$AGE_cat_4,
     xlab = "Age2", ylab = "Age4", 
     col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs4_nn_mod_colb$AGE_cat_3, 
     synthetic_hhs4_nn_mod_colb$AGE_cat_4,
     xlab = "Age3", ylab = "Age4", 
     col = brewer.pal(name = "Spectral", n = 10))

par(mfrow = c(3, 7))
##########
####### Give pairwise comparison of MST for hhs3
#TRUE
title("True")
plot(hhs_3_true$MST_1, hhs_3_true$MST_2, xlab = "MST1",
     ylab = "MST2", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_3_true$MST_1, hhs_3_true$MST_3, xlab = "MST1",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_3_true$MST_2, hhs_3_true$MST_3, xlab = "MST2",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
#C_LOG
title("C-LOG")
plot(synthetic_hhs3_logit_colb$MST_1, synthetic_hhs3_logit_colb$MST_2, xlab = "MST1",
     ylab = "MST2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_colb$MST_1, synthetic_hhs3_logit_colb$MST_3, xlab = "MST1",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_colb$MST_2, synthetic_hhs3_logit_colb$MST_3, xlab = "MST2",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
#C-NN
title("C-NN")
plot(synthetic_hhs3_nn_colb$MST_1, synthetic_hhs3_nn_colb$MST_2, 
     xlab = "MST1",
     ylab = "MST2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_colb$MST_1, synthetic_hhs3_nn_colb$MST_3,
     xlab = "MST1",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_colb$MST_2, synthetic_hhs3_nn_colb$MST_3, 
     xlab = "MST2",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
#C-RF
title("C-RF")
plot(synthetic_hhs3_rf_colb$MST_1, synthetic_hhs3_rf_colb$MST_2, 
     xlab = "MST1",
     ylab = "MST2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_colb$MST_1, synthetic_hhs3_rf_colb$MST_3,
     xlab = "MST1",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_colb$MST_2, synthetic_hhs3_rf_colb$MST_3, 
     xlab = "MST2",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
#M-LOG
title("M-LOG")
plot(synthetic_hhs3_logit_mod_colb$MST_1, synthetic_hhs3_logit_mod_colb$MST_2, 
     xlab = "MST1",
     ylab = "MST2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_mod_colb$MST_1, synthetic_hhs3_logit_mod_colb$MST_3, 
     xlab = "MST1",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_logit_mod_colb$MST_2, synthetic_hhs3_logit_mod_colb$MST_3,
     xlab = "MST2",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
#M-NN  
title("M-NN")
plot(synthetic_hhs3_nn_mod_colb$MST_1, synthetic_hhs3_nn_mod_colb$MST_2, 
     xlab = "MST1",
     ylab = "MST2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_mod_colb$MST_1, synthetic_hhs3_nn_mod_colb$MST_3, 
     xlab = "MST1",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_nn_mod_colb$MST_2, synthetic_hhs3_nn_mod_colb$MST_3,
     xlab = "MST2",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
#M-RF 
title("M-RF")
plot(synthetic_hhs3_rf_mod_colb$MST_1, synthetic_hhs3_rf_mod_colb$MST_2, 
     xlab = "MST1",
     ylab = "MST2", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_mod_colb$MST_1, synthetic_hhs3_rf_mod_colb$MST_3, 
     xlab = "MST1",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs3_rf_mod_colb$MST_2, synthetic_hhs3_rf_mod_colb$MST_3,
     xlab = "MST2",
     ylab = "MST3", col = brewer.pal(name = "Spectral", n = 10))
par(mfrow = c(4,7))
######
# MST_hhs6_last
title("True")
plot(hhs_6_true$MST_3, hhs_6_true$MST_4, xlab = "MST3",
     ylab = "MST4", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_6_true$MST_3, hhs_6_true$MST_5, xlab = "MST3",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_6_true$MST_3, hhs_6_true$MST_6, xlab = "MST3",
     ylab = "MST6", col = brewer.pal(name = "Spectral", n = 10))
plot(hhs_6_true$MST_4, hhs_6_true$MST_5, xlab = "MST4",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))

#C_LOG
title("C-LOG")
plot(synthetic_hhs6_logit_colb$MST_3, synthetic_hhs6_logit_colb$MST_4, xlab = "MST3",
     ylab = "MST4", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_logit_colb$MST_3, synthetic_hhs6_logit_colb$MST_5, xlab = "MST3",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_logit_colb$MST_3, synthetic_hhs6_logit_colb$MST_6, xlab = "MST3",
     ylab = "MST6", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_logit_colb$MST_4, synthetic_hhs6_logit_colb$MST_5, xlab = "MST4",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
#C-NN
title("C-NN")
plot(synthetic_hhs6_nn_colb$MST_3, synthetic_hhs6_nn_colb$MST_4, xlab = "MST3",
     ylab = "MST4", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_nn_colb$MST_3, synthetic_hhs6_nn_colb$MST_5, xlab = "MST3",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_nn_colb$MST_3, synthetic_hhs6_nn_colb$MST_6, xlab = "MST3",
     ylab = "MST6", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_nn_colb$MST_4, synthetic_hhs6_nn_colb$MST_5, xlab = "MST4",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
#C-RF
title("C-RF")
plot(synthetic_hhs6_rf_colb$MST_3, synthetic_hhs6_rf_colb$MST_4, xlab = "MST3",
     ylab = "MST4", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_rf_colb$MST_3, synthetic_hhs6_rf_colb$MST_5, xlab = "MST3",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_rf_colb$MST_3, synthetic_hhs6_rf_colb$MST_6, xlab = "MST3",
     ylab = "MST6", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_rf_colb$MST_4, synthetic_hhs6_rf_colb$MST_5, xlab = "MST4",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
#M-LOG
title("M-LOG")
plot(synthetic_hhs6_logit_mod_colb$MST_3, synthetic_hhs6_logit_mod_colb$MST_4, xlab = "MST3",
     ylab = "MST4", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_logit_mod_colb$MST_3, synthetic_hhs6_logit_mod_colb$MST_5, xlab = "MST3",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_logit_mod_colb$MST_3, synthetic_hhs6_logit_mod_colb$MST_6, xlab = "MST3",
     ylab = "MST6", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_logit_mod_colb$MST_4, synthetic_hhs6_logit_mod_colb$MST_5, xlab = "MST4",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
#M-NN  
title("M-NN")
plot(synthetic_hhs6_nn_mod_colb$MST_3, synthetic_hhs6_nn_mod_colb$MST_4, xlab = "MST3",
     ylab = "MST4", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_nn_mod_colb$MST_3, synthetic_hhs6_nn_mod_colb$MST_5, xlab = "MST3",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_nn_mod_colb$MST_3, synthetic_hhs6_nn_mod_colb$MST_6, xlab = "MST3",
     ylab = "MST6", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_nn_mod_colb$MST_4, synthetic_hhs6_nn_mod_colb$MST_5, xlab = "MST4",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
title("M-RF")
plot(synthetic_hhs6_rf_mod_colb$MST_3, synthetic_hhs6_rf_mod_colb$MST_4, xlab = "MST3",
     ylab = "MST4", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_rf_mod_colb$MST_3, synthetic_hhs6_rf_mod_colb$MST_5, xlab = "MST3",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_rf_mod_colb$MST_3, synthetic_hhs6_rf_mod_colb$MST_6, xlab = "MST3",
     ylab = "MST6", col = brewer.pal(name = "Spectral", n = 10))
plot(synthetic_hhs6_rf_mod_colb$MST_4, synthetic_hhs6_rf_mod_colb$MST_5, xlab = "MST4",
     ylab = "MST5", col = brewer.pal(name = "Spectral", n = 10))
#############
# 
par(mfrow = c(1,1))

# Create figure 5; therefore, extract AGE_diff
x <- data.frame(True = hhs_6_true$AGE_diff,
                C_LOG = synthetic_hhs6_logit_colb$AGE_diff,
                C_NN = synthetic_hhs6_nn_colb$AGE_diff, 
                C_RF = synthetic_hhs6_rf_colb$AGE_diff)
library(ggplot2)
library(reshape2)
data<- melt(x)
AGE_diff_2_compare <- ggplot(data,aes(x=value, fill=variable)) + 
  geom_density(alpha = 0.25) +
  theme_bw() +
  xlab("Individual age deviation from age mean in household") +
  ylab("Density") +
  labs(fill="Dataset") +
  ggtitle("HHS2")

AGE_diff_3_compare <- ggplot(data,aes(x=value, fill=variable)) + 
  geom_density(alpha = 0.25) +
  theme_bw() +
  xlab("Individual age deviation from age mean in household") +
  ylab("Density") +
  labs(fill="Dataset") +
  ggtitle("HHS3")

AGE_diff_4_compare <- ggplot(data,aes(x=value, fill=variable)) + 
  geom_density(alpha = 0.25) +
  theme_bw() +
  xlab("Individual age deviation from age mean in household") +
  ylab("Density") + 
  labs(fill="Dataset") +
  ggtitle("HHS4")

AGE_diff_5_compare <- ggplot(data,aes(x=value, fill=variable)) + 
  geom_density(alpha = 0.25) +
  theme_bw() +
  xlab("Individual age deviation from age mean in household") +
  ylab("Density") +
  labs(fill="Dataset") +
  ggtitle("HHS5")

AGE_diff_6_compare <- ggplot(data,aes(x=value, fill=variable)) + 
  geom_density(alpha = 0.25) +
  theme_bw() +
  xlab("Individual age deviation from age mean in household") +
  ylab("Density") + 
  labs(fill="Dataset") +
  ggtitle("HHS6")

# This is figure 5
grid.arrange(AGE_diff_2_compare, AGE_diff_3_compare, AGE_diff_4_compare,
             AGE_diff_5_compare, AGE_diff_6_compare)

##########
# Find structural zeros, according to pre-defined rules
# 1 = FR (First rule as given in chapter 5.3)
# 2 = SR (Second rule as given in chapter 5.3)

sz_hhs2 <- function(x) {
  x$sc <- 0
  x$sc[as.numeric(x$AGE_cat_1) < 2 & as.numeric(x$AGE_cat_2) < 2] <- 1
  x$sc[as.numeric(x$AGE_cat_1) < 2 & x$MST_1 != 1|
         as.numeric(x$AGE_cat_2) < 2 & x$MST_2 != 1] <- 2
  return(x)
}

sz_hhs3 <- function(x) {
  x$sc <- 0
  x$sc[as.numeric(x$AGE_cat_1) < 2 & as.numeric(x$AGE_cat_2) < 2 
       & as.numeric(x$AGE_cat_3) < 2] <- 1
  x$sc[as.numeric(x$AGE_cat_1) < 2 & x$MST_1 != 1|
         as.numeric(x$AGE_cat_2) < 2 & x$MST_2 != 1|
         as.numeric(x$AGE_cat_3) < 2 & x$MST_3 != 1] <- 2
  return(x)
}

sz_hhs4 <- function(x) {
  x$sc <- 0
  x$sc[as.numeric(x$AGE_cat_1) < 2 & as.numeric(x$AGE_cat_2) < 2 
       & as.numeric(x$AGE_cat_3) & as.numeric(x$AGE_cat_4) < 2] <- 1
  x$sc[as.numeric(x$AGE_cat_1) < 2 & x$MST_1 != 1|
         as.numeric(x$AGE_cat_2) < 2 & x$MST_2 != 1|
         as.numeric(x$AGE_cat_3) < 2 & x$MST_3 != 1|
         as.numeric(x$AGE_cat_4) < 2 & x$MST_4 != 1] <- 2
  return(x)
}

sz_hhs5 <- function(x) {
  x$sc <- 0
  x$sc[as.numeric(x$AGE_cat_1) < 2 & as.numeric(x$AGE_cat_2) < 2 
       & as.numeric(x$AGE_cat_3) & as.numeric(x$AGE_cat_4) < 2
       & as.numeric(x$AGE_cat_5) < 2] <- 1
  x$sc[as.numeric(x$AGE_cat_1) < 2 & x$MST_1 != 1|
         as.numeric(x$AGE_cat_2) < 2 & x$MST_2 != 1|
         as.numeric(x$AGE_cat_3) < 2 & x$MST_3 != 1|
         as.numeric(x$AGE_cat_4) < 2 & x$MST_4 != 1|
         as.numeric(x$AGE_cat_5) < 2 & x$MST_5 != 1] <- 2
  return(x)
}

sz_hhs6 <- function(x) {
  x$sc <- 0
  x$sc[as.numeric(x$AGE_cat_1) < 2 & as.numeric(x$AGE_cat_2) < 2 
       & as.numeric(x$AGE_cat_3) & as.numeric(x$AGE_cat_4) < 2
       & as.numeric(x$AGE_cat_5) & as.numeric(x$AGE_cat_6) < 2] <- 1
  x$sc[as.numeric(x$AGE_cat_1) < 2 & x$MST_1 != 1|
         as.numeric(x$AGE_cat_2) < 2 & x$MST_2 != 1|
         as.numeric(x$AGE_cat_3) < 2 & x$MST_3 != 1|
         as.numeric(x$AGE_cat_4) < 2 & x$MST_4 != 1|
         as.numeric(x$AGE_cat_5) < 2 & x$MST_5 != 1|
         as.numeric(x$AGE_cat_6) < 2 & x$MST_6 != 1] <- 2
  return(x)
}

synthetic_hhs2_logit_colb <- sz_hhs2(synthetic_hhs2_logit_colb)
synthetic_hhs2_nn_colb <- sz_hhs2(synthetic_hhs2_nn_colb)
synthetic_hhs2_rf_colb <- sz_hhs2(synthetic_hhs2_rf_colb)

# and so on...

#################
############
#Pairwise correlation between household members; age

cor(synthetic_hhs4_rf_mod_colb[,c("AGE_1", "AGE_2", "AGE_3", "AGE_4")])

# Cramers V matrix

library(vcd)

# Initialize empty matrix to store coefficients
empty_m <- matrix(ncol = length(synthetic_hhs2_rf_mod_colb[,c("SEX_1", "SEX_2", "MST_1", "MST_2")]),
                  nrow = length(synthetic_hhs2_rf_mod_colb[,c("SEX_1", "SEX_2", "MST_1", "MST_2")]),
                  dimnames = list(names(synthetic_hhs2_rf_mod_colb[,c("SEX_1", "SEX_2", "MST_1", "MST_2")]), 
                                  names(synthetic_hhs2_rf_mod_colb[,c("SEX_1", "SEX_2", "MST_1", "MST_2")])))
# Function that accepts matrix for coefficients and data and returns a correlation matrix
calculate_cramer <- function(m, df) {
  for (r in seq(nrow(m))){
    for (c in seq(ncol(m))){
      m[[r, c]] <- assocstats(table(df[[r]], df[[c]]))$cramer
    }
  }
  return(m)
}