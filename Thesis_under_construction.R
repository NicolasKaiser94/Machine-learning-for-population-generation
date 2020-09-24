library(sampling)
library(DescTools)
library(plyr)
library(forcats)
library(randomForest)
library(nnet)
library(data.table)
library(moments)
library(rpart)

AMELIA_files <- paste0(path, list.files(path))

for(i in 1:6){
  load(AMELIA_files[i])
}

# Create a data frame for personal level information
HHS <- ifelse(HHS >= 6, 6, HHS)
# Alter klassieren 
AGE_cat <- rep("83+x", length(AGE))
AGE_cat[AGE < 18] <- "0-17"
AGE_cat[AGE > 17 & AGE <= 28] <- "18-28"
AGE_cat[AGE > 28 & AGE <= 39] <- "29-39"
AGE_cat[AGE > 40 & AGE <= 50] <- "40-50"
AGE_cat[AGE > 50 & AGE <= 61] <- "51-61"
AGE_cat[AGE > 61 & AGE <= 72] <- "62-72"
AGE_cat[AGE > 72 & AGE <= 83] <- "73-83"

AMELIA <- data.table(HID, HHS, AGE, AGE_cat, MST, SEX, INC) # "AMELIA" with PIDs in decreasing order of age

AMELIA[order(HID, AGE)] # Innerhalb eines HH nach Alter sortieren
AMELIA[, PID := order(AGE, decreasing = TRUE), HID] # Innerhalb eines HH Personen-ID generieren, älteste Person ist 1
AMELIA[, AGE_cat:=as.factor(AGE_cat)]
AMELIA[, MST:=as.factor(MST)]
AMELIA[, SEX:=as.factor(SEX)]
AMELIA[, HHS:=as.factor(HHS)]

rm(AGE, MST, SEX, HHS, HID, PID, EDU, INC)

# We create personal IDs which indicate the order of age within a household. It´s possible
# to not have this information in reality. Nevertheless, it is necessary to apply
# one of the modelling aproaches we want to test. In this approach, we model the characterics of,
# for instance, the second oldest person in a two person household by the characteristics of the oldest
# person. We do this, for every single characteristic. 
# For the moment, let´s assume that we HAVE this information. 

AMELIA$PIDX <- rep(1:nrow(AMELIA))

# Create a vector for the relative proportion of a household class. Might be useful for
# more sophisticated weighting approaches
AMELIA$inc_prop_hhs <- 0
AMELIA$inc_prop_hhs[AMELIA$HHS == 1] <- 0.2277337
AMELIA$inc_prop_hhs[AMELIA$HHS == 2] <- 0.311736
AMELIA$inc_prop_hhs[AMELIA$HHS == 3] <- 0.1860352
AMELIA$inc_prop_hhs[AMELIA$HHS == 4] <- 0.1793412
AMELIA$inc_prop_hhs[AMELIA$HHS == 5] <- 0.06593175
AMELIA$inc_prop_hhs[AMELIA$HHS == 6] <- 0.02922205

AMELIA_HH_wide <- dcast(AMELIA, HID + HHS ~ PID, value.var = c("AGE_cat", "AGE", "MST", "SEX", "INC", "PIDX"))

# AMELIA_HH_wide will be our true household dataset, from which we can subset different groups.

# At this point it must be noted that we can now proceed differently for the modelling. On the one hand, 
# we can regard the members in a household as a specific combination of members and characteristics for
# whose occurrence a certain probability is estimated via a model. For this purpose. However, 
# it is also possible to consider the characteristics of 
# one household member as a dependent variable, which is explained by the characteristics of the other 
# members. Both approaches will be exmined below.

###################
###################
# COMBINATION APPROCH STARTS HERE
################################

# Overall, the method used can be explained quite simply: the rows coded
# as "true combinations" (based on a sample of true combinations) of household members are set against
# as many random combinations as possible. 
# The subsequent application of the model works according to the same principle: random combinations are 
# created from a sample of individuals until our model recognizes them as plausible.
# We can improve model quality by setting the true values against a larger number of possible random 
# combinations. 

# The following function takes a sample of individuals and returns a set of random combinations
# of individuals. Like explained above, this function hels to build the model data and to create
# synthetic households

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
hhs_2_true$SEX_diff <- interaction(hhs_2_true$SEX_1, hhs_2_true$SEX_2, drop = FALSE)
hhs_2_true$INC <- hhs_2_true$INC_1 + hhs_2_true$INC_2 

# Using combinations instead of perumutations. Therefore, we order the factor categories
f1 <- pmin(as.numeric(hhs_2_true$MST_1), as.numeric(hhs_2_true$MST_2))
f2 <- pmax(as.numeric(hhs_2_true$MST_1), as.numeric(hhs_2_true$MST_2))
hhs_2_true$MST_group <- droplevels(interaction(f1, f2, sep=""))

f1 <- pmin(as.numeric(hhs_2_true$AGE_cat_1), as.numeric(hhs_2_true$AGE_cat_2))
f2 <- pmax(as.numeric(hhs_2_true$AGE_cat_1), as.numeric(hhs_2_true$AGE_cat_2))
hhs_2_true$AGE_comb <- droplevels(interaction(f1, f2, sep=""))


hhs_3_true <- AMELIA_HH_wide[,c("HID", "HHS", "AGE_cat_1", "AGE_cat_2", "AGE_cat_3", "AGE_1", "AGE_2",
                                "AGE_3", "SEX_1", "SEX_2", "SEX_3", "MST_1", "MST_2", "MST_3", 
                                "INC_1", "INC_2", "INC_3", "PIDX_1", 
                                "PIDX_2", "PIDX_3")][AMELIA_HH_wide$HHS == 3]

hhs_3_true$AGE_mean <- (hhs_3_true$AGE_1 + hhs_3_true$AGE_2 + hhs_3_true$AGE_3) / 3
hhs_3_true$AGE_diff <- sqrt(((hhs_3_true$AGE_mean - hhs_3_true$AGE_1)^2 + (hhs_3_true$AGE_mean - hhs_3_true$AGE_2)^2 +
                             (hhs_3_true$AGE_mean - hhs_3_true$AGE_3)^2) / 3) # Standardabweichung von Hand berechnet! Und da sage noch einer, Grundlagen seien "trocken"...
hhs_3_true$SEX_diff <- interaction(hhs_3_true$SEX_1, hhs_3_true$SEX_2, hhs_3_true$SEX_3)
hhs_3_true$INC <- hhs_3_true$INC_1 + hhs_3_true$INC_2 + hhs_3_true$INC_3

f1 <- pmin(as.numeric(hhs_3_true$MST_1), as.numeric(hhs_3_true$MST_2))
f2 <- pmin(as.numeric(hhs_3_true$MST_1), as.numeric(hhs_3_true$MST_2), as.numeric(hhs_3_true$MST_3))
f3 <- pmax(as.numeric(hhs_3_true$MST_1), as.numeric(hhs_3_true$MST_2), as.numeric(hhs_3_true$MST_3))

hhs_3_true$MST_group <- droplevels(interaction(f1, f2, f3, sep = ""))

f1 <- pmin(as.numeric(hhs_3_true$AGE_cat_1), as.numeric(hhs_3_true$AGE_cat_2))
f2 <- pmin(as.numeric(hhs_3_true$AGE_cat_1), as.numeric(hhs_3_true$AGE_cat_2), as.numeric(hhs_3_true$AGE_cat_3))
f3 <- pmax(as.numeric(hhs_3_true$AGE_cat_1), as.numeric(hhs_3_true$AGE_cat_2), as.numeric(hhs_3_true$AGE_cat_3))

hhs_3_true$AGE_comb <- droplevels(interaction(f1, f2, f3, sep = ""))

# Same for households of size 4

hhs_4_true <- AMELIA_HH_wide[,c("HID", "HHS", "AGE_cat_1", "AGE_cat_2", "AGE_cat_3", "AGE_cat_4",
                                "AGE_1", "AGE_2","AGE_3","AGE_4", "SEX_1", "SEX_2", "SEX_3", "SEX_4", 
                                "MST_1", "MST_2", "MST_3", "MST_4", "INC_1", "INC_2", "INC_3", "INC_4",
                                "PIDX_1", "PIDX_2", "PIDX_3", "PIDX_4")][AMELIA_HH_wide$HHS == 4]

hhs_4_true$AGE_mean <- (hhs_4_true$AGE_1 + hhs_4_true$AGE_2 + hhs_4_true$AGE_3 + hhs_4_true$AGE_4) / 4
hhs_4_true$AGE_diff <- sqrt(((hhs_4_true$AGE_mean - hhs_4_true$AGE_1)^2 + (hhs_4_true$AGE_mean - hhs_4_true$AGE_2)^2 +
                               (hhs_4_true$AGE_mean - hhs_4_true$AGE_3)^2 +
                               (hhs_4_true$AGE_mean - hhs_4_true$AGE_4)^2)) / 4 # Standardabweichung von Hand berechnet! Und da sage noch einer, Grundlagen seien "trocken"...
hhs_4_true$SEX_diff <- interaction(hhs_4_true$SEX_1, hhs_4_true$SEX_2, hhs_4_true$SEX_3,
                                   hhs_4_true$SEX_4)
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
                               (hhs_5_true$AGE_mean - hhs_5_true$AGE_5)^2)) / 5 # Standardabweichung von Hand berechnet! Und da sage noch einer, Grundlagen seien "trocken"...
hhs_5_true$SEX_diff <- interaction(hhs_5_true$SEX_1, hhs_5_true$SEX_2, hhs_5_true$SEX_3,
                                   hhs_5_true$SEX_4, hhs_5_true$SEX_5)
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
hhs_6_true$SEX_diff <- droplevels(interaction(f1, f2, f3, f4, f5, f6, sep = ""))

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


# Create subsets of individuals belonging to a certain household size for modelling
AMELIA_HHS2 <- AMELIA[HHS == 2]
AMELIA_HHS3 <- AMELIA[HHS == 3] 
AMELIA_HHS4 <- AMELIA[HHS == 4]
AMELIA_HHS5 <- AMELIA[HHS == 5]
AMELIA_HHS6 <- AMELIA[HHS == 6]

# Draw a random numer which determines the categorical output of a variable
synth_nn <- list()
create_synth <- function(hhs, model, x) {
  repeat{
    # Household size 2
    # #AGE_1 0-90 AGE_2
    if(hhs == 2) {
      AGE_1 <- sample(0:90, 1)
      AGE_2 <- sample(0:90, 1)
      # SEX_1 1/2 SEX_2 1/2
      SEX_1 <- sample(1:2, 1)
      SEX_2 <- sample(1:2, 1)
      # MST_1 1-5
      MST_1 <- sample(1:5, 1)
      MST_2 <- sample(1:5, 1)
      ###
      choice <- as.data.frame(cbind(AGE_1, AGE_2, SEX_1, SEX_2, MST_1, MST_2))
      choice$SEX_1 <- factor(choice$SEX_1, levels = levels(hhs_2_true$SEX_1))
      choice$SEX_2 <- factor(choice$SEX_2, levels = levels(hhs_2_true$SEX_2))
      choice$MST_1 <- factor(choice$MST_1, levels = levels(hhs_2_true$MST_1))
      choice$MST_2 <- factor(choice$MST_2, levels = levels(hhs_2_true$MST_2))
    }
    if(hhs == 3) {    
      AGE_1 <- sample(0:90, 1)
      AGE_2 <- sample(0:90, 1)
      AGE_3 <- sample(0:90, 1)
      # SEX_1 1/2 SEX_2 1/2
      SEX_1 <- sample(1:2, 1)
      SEX_2 <- sample(1:2, 1)
      SEX_3 <- sample(1:2, 1)
      # MST_1 1-5
      MST_1 <- sample(1:5, 1)
      MST_2 <- sample(1:5, 1)
      MST_3 <- sample(1:5, 1)
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
      AGE_1 <- sample(0:90, 1)
      AGE_2 <- sample(0:90, 1)
      AGE_3 <- sample(0:90, 1)
      AGE_4 <- sample(0:90, 1)
      # SEX_1 1/2 SEX_2 1/2
      SEX_1 <- sample(1:2, 1)
      SEX_2 <- sample(1:2, 1)
      SEX_3 <- sample(1:2, 1)
      SEX_4 <- sample(1:2, 1)
      # MST_1 1-5
      MST_1 <- sample(1:5, 1)
      MST_2 <- sample(1:5, 1)
      MST_3 <- sample(1:5, 1)
      MST_4 <- sample(1:5, 1)
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
      AGE_1 <- sample(0:90, 1)
      AGE_2 <- sample(0:90, 1)
      AGE_3 <- sample(0:90, 1)
      AGE_4 <- sample(0:90, 1)
      AGE_5 <- sample(0:90, 1)
      # SEX_1 1/2 SEX_2 1/2
      SEX_1 <- sample(1:2, 1)
      SEX_2 <- sample(1:2, 1)
      SEX_3 <- sample(1:2, 1)
      SEX_4 <- sample(1:2, 1)
      SEX_5 <- sample(1:2, 1)
      # MST_1 1-5
      MST_1 <- sample(1:5, 1)
      MST_2 <- sample(1:5, 1)
      MST_3 <- sample(1:5, 1)
      MST_4 <- sample(1:5, 1)
      MST_5 <- sample(1:5, 1)
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
      AGE_1 <- sample(0:90, 1)
      AGE_2 <- sample(0:90, 1)
      AGE_3 <- sample(0:90, 1)
      AGE_4 <- sample(0:90, 1)
      AGE_5 <- sample(0:90, 1)
      AGE_6 <- sample(0:90, 1)
      # SEX_1 1/2 SEX_2 1/2
      SEX_1 <- sample(1:2, 1)
      SEX_2 <- sample(1:2, 1)
      SEX_3 <- sample(1:2, 1)
      SEX_4 <- sample(1:2, 1)
      SEX_5 <- sample(1:2, 1)
      SEX_6 <- sample(1:2, 1)
      # MST_1 1-5
      MST_1 <- sample(1:5, 1)
      MST_2 <- sample(1:5, 1)
      MST_3 <- sample(1:5, 1)
      MST_4 <- sample(1:5, 1)
      MST_5 <- sample(1:5, 1)
      MST_6 <- sample(1:5, 1)
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
    
    if(model == "NN") {        
      pr_synthetic_hhs <- predict(x, newdata = choice, type = "raw")
    }
    if(model == "RF") {        
      pr_synthetic_hhs <- predict(x, newdata = choice, type = "prob")
    }
    r <- runif(1)
    if(r < pr_synthetic_hhs) {
      synth_nn <- choice
      break
    }} 
  return(synth_nn)
}

# Create a function to simulate the creation of synthetic households with the combination
# approach. The function allows to vary the sample size of the indivudals, the sample size
# of the true households used for the model, the number of random combinations set against the
# true combinations, the number of combinations the model has to pick the realistic ones from,
# an argument "cut_categories" to collapse factors which occur seldom into one factor.
#
Monte_Carlo_Simulation <- function(HID_sample, n_model_comb, household_size,
                                   True_hhs, AMELIA_HHS, hhs_true, NP) {
  out <- list()  
  
  SI_SAMPLE <- srswor(n = HID_sample, nrow(hhs_true)) ### Our sample.
  AMELIA_sample <- hhs_true[SI_SAMPLE>0, ]
  head(AMELIA_sample)
  
  full_synthetic_simulation <- AMELIA[HID %in% AMELIA_sample$HID]
  
  datalist = list()
  for(i in 1:n_model_comb) { # 100 Random household allocations based on sample
    dat <- ind_to_random_hid(full_synthetic_simulation, household_size = household_size,
                             synthetic_hhs = full_synthetic_simulation$HHS,
                             drop = FALSE)
    datalist[[i]] <- dat
    print(i)
  }
  
  ra <- do.call(rbind, datalist)
  ra$TRUE_ <- 0
  AMELIA_sample$TRUE_ <- 1
  model_data_cx <- rbind(AMELIA_sample, ra, fill = TRUE)
  model_data_cx$TRUE_ <- as.factor(model_data_cx$TRUE_)
  model_data_cx$SEX_diff <- as.factor(model_data_cx$SEX_diff)
  
  if(household_size == 2) {
    hhs_nn_cx <- nnet(TRUE_ ~ AGE_1 + AGE_2 + SEX_1 + SEX_2 + MST_1 + MST_2, data = model_data_cx, size = 2)
    hhs_rf_cx <- randomForest(TRUE_ ~ AGE_1 + AGE_2 + SEX_1 + SEX_2 + MST_1 + MST_2, data = model_data_cx, type = "prob")
    hhs_logit_cx <- glm(TRUE_ ~ AGE_1 + AGE_2 + SEX_1 + SEX_2 + MST_1 + MST_2, data = model_data_cx, family = "binomial")
  }
  
  if(household_size == 3) {
    hhs_nn_cx <- nnet(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + SEX_1 + SEX_2 + SEX_3 + MST_1 + MST_2 + MST_3, data = model_data_cx, size = 2)
    hhs_rf_cx <- randomForest(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + SEX_1 + SEX_2 + SEX_3 + MST_1 + MST_2 + MST_3, data = model_data_cx, type = "prob")
    hhs_logit_cx <- glm(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + SEX_1 + SEX_2 + SEX_3 + MST_1 + MST_2 + MST_3, data = model_data_cx, family = "binomial")
  }
  
  if(household_size == 4) {
    hhs_nn_cx <- nnet(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + MST_1 + MST_2 + MST_3 + MST_4, data = model_data_cx, size = 2)
    hhs_rf_cx <- randomForest(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + MST_1 + MST_2 + MST_3 + MST_4, data = model_data_cx, type = "prob")
    hhs_logit_cx <- glm(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + MST_1 + MST_2 + MST_3 + MST_4, data = model_data_cx, family = "binomial")
  }
  
  if(household_size == 5) {
    hhs_nn_cx <- nnet(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5, data = model_data_cx, size = 2)
    hhs_rf_cx <- randomForest(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5, data = model_data_cx, type = "prob")
    hhs_logit_cx <- glm(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5, data = model_data_cx, family = "binomial")
  }
  
  if(household_size == 6) {
    hhs_nn_cx <- nnet(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + AGE_6 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + SEX_6 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5 + MST_6, data = model_data_cx, size = 2)
    hhs_rf_cx <- randomForest(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + AGE_6 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + SEX_6 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5 + MST_6, data = model_data_cx, type = "prob")
    hhs_logit_cx <- glm(TRUE_ ~ AGE_1 + AGE_2 + AGE_3 + AGE_4 + AGE_5 + AGE_6 + SEX_1 + SEX_2 + SEX_3 + SEX_4 + SEX_5 + SEX_6 + MST_1 + MST_2 + MST_3 + MST_4 + MST_5 + MST_6, data = model_data_cx, family = "binomial")
  }
  
  result <- list()
  for(i in 1:NP) {
    resultX <- create_synth(hhs = household_size, model = "NN", x = hhs_nn_cx)
    result[[length(result)+1]] <- resultX
    print(i)
  }
  synthetic_hhs_nn <- result
  result <- list()
  for(i in 1:NP) {
    resultX <- create_synth(hhs = household_size, model = "RF", x = hhs_rf_cx)
    result[[length(result)+1]] <- resultX
    print(i)
  }
  synthetic_hhs_rf <- result
  out <- list(synthetic_hhs_nn, synthetic_hhs_rf)
  return(out)
}

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
    AMELIA_wide_start_c$SEX_diff <- interaction(AMELIA_wide_start_c$SEX_1, AMELIA_wide_start_c$SEX_2, drop = FALSE)
  } else if(hhs_size == 3) {
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), as.numeric(AMELIA_wide_start_c$MST_3))
    f3 <- pmax(as.numeric(AMELIA_wide_start_c$MST_1), as.numeric(AMELIA_wide_start_c$MST_2), as.numeric(AMELIA_wide_start_c$MST_3))
    
    AMELIA_wide_start_c$MST_group <- droplevels(interaction(f1, f2, f3, sep = ""))
    
    f1 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2))
    f2 <- pmin(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), as.numeric(AMELIA_wide_start_c$AGE_cat_3))
    f3 <- pmax(as.numeric(AMELIA_wide_start_c$AGE_cat_1), as.numeric(AMELIA_wide_start_c$AGE_cat_2), as.numeric(AMELIA_wide_start_c$AGE_cat_3))
    
    AMELIA_wide_start_c$AGE_comb <- droplevels(interaction(f1, f2, f3, sep = ""))
    AMELIA_wide_start_c$SEX_diff <- interaction(AMELIA_wide_start_c$SEX_1, AMELIA_wide_start_c$SEX_2,
                                                AMELIA_wide_start_c$SEX_3, drop = FALSE)
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
    AMELIA_wide_start_c$SEX_diff <- interaction(AMELIA_wide_start_c$SEX_1, AMELIA_wide_start_c$SEX_2,
                                                AMELIA_wide_start_c$SEX_3,
                                                AMELIA_wide_start_c$SEX_4, drop = FALSE)
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
    AMELIA_wide_start_c$SEX_diff <- interaction(AMELIA_wide_start_c$SEX_1, AMELIA_wide_start_c$SEX_2,
                                                AMELIA_wide_start_c$SEX_3,
                                                AMELIA_wide_start_c$SEX_4,
                                                AMELIA_wide_start_c$SEX_5, drop = FALSE)
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
    AMELIA_wide_start_c$SEX_diff <- interaction(AMELIA_wide_start_c$SEX_1,
                                                AMELIA_wide_start_c$SEX_2,
                                                AMELIA_wide_start_c$SEX_3,
                                                AMELIA_wide_start_c$SEX_4,
                                                AMELIA_wide_start_c$SEX_5, 
                                                AMELIA_wide_start_c$SEX_6, drop = FALSE)
  }
  return(AMELIA_wide_start_c)
}

calc_bias <- function(est, true) {
  r <- list()
  est$AGE_comb <- factor(est$AGE_comb, levels = levels(true$AGE_comb))
  est$MST_group <- factor(est$MST_group, levels = levels(true$MST_group))
  est$SEX_diff <- factor(est$SEX_diff, levels = levels(true$SEX_diff))
  r$deviation_AGE <- sum((prop.table(table(est$AGE_comb)) - prop.table(table(true$AGE_comb)))^2) *100
  r$deviation_MST <- sum((prop.table(table(est$MST_group)) - prop.table(table(true$MST_group)))^2) *100
  r$deviation_SEX <- sum((prop.table(table(est$SEX_diff)) - prop.table(table(true$SEX_diff)))^2) *100
  return(r)
}

# For the full population we need 
# 1178764 households of size 2
# 703453 households of size 3
# 678141 households of size 4
# 249307 households of size 5
# 110497 households of size 6

test_2 <- Monte_Carlo_Simulation(1000, n_model_comb = 50, NP = 100, hhs_true = hhs_2_true,
                               household_size = 2)

test_nn <- do.call(rbind,test_2[[1]])
test_rf <- do.call(rbind,test_2[[2]])

test_nn$HHS <- 2
test_nn$HID <- rep(1:nrow(test_nn))
test_rf$HHS <- 2
test_rf$HID <- rep(1:nrow(test_rf))

test_nn <- AGE_cat_creator(test_nn)
test_rf <- AGE_cat_creator(test_rf)

test_nn$AGE_cat_1 <- as.factor(test_nn$AGE_cat_1)
test_nn$AGE_cat_2 <- as.factor(test_nn$AGE_cat_2)
test_nn$AGE_cat_3 <- as.factor(test_nn$AGE_cat_3)
test_nn$AGE_cat_4 <- as.factor(test_nn$AGE_cat_4)
test_nn$AGE_cat_5 <- as.factor(test_nn$AGE_cat_5)
test_nn$AGE_cat_6 <- as.factor(test_nn$AGE_cat_6)

test_rf$AGE_cat_1 <- as.factor(test_rf$AGE_cat_1)
test_rf$AGE_cat_2 <- as.factor(test_rf$AGE_cat_2)
test_rf$AGE_cat_3 <- as.factor(test_rf$AGE_cat_3)
test_rf$AGE_cat_4 <- as.factor(test_rf$AGE_cat_4)
test_rf$AGE_cat_5 <- as.factor(test_rf$AGE_cat_5)
test_rf$AGE_cat_6 <- as.factor(test_rf$AGE_cat_6)

test_nn <- calc_difference(test_nn, hhs_size = 2)
test_rf <- calc_difference(test_rf, hhs_size = 2)

plot(AMELIA_wide_start_nnet_2$MST_group)
plot(hhs_2_true$MST_group)
calc_bias(test_nn, hhs_2_true)
calc_bias(test_rf, hhs_2_true)
calc_bias(AMELIA_wide_start_nnet_2, hhs_2_true)
calc_bias(AMELIA_wide_start_rf_2, hhs_2_true)
# For the random forest, it is a problem to have a predictor with more than 53 categories.
# This problem can be solved by merging factor categories with very low numbers (1, 2 etc.)

# We have 5 synthetic populations for every household size under optimal conditions (perfect household size allocation)
# and such under the realistic case of estimating the household size. Now: How do we compare the
# resulting datasets with the true population?

###########################
# HERE BEGINS THE APPROACH USING MODELS FOR THE CHARACTERISTICS OF A HOUSEHOLD MEMBER BASED BY 
# OTHER MEMBERS

AMELIA_PID1 <- AMELIA[PID == 1] # AMELIA-Datensatz mit ältester Person im HH
AMELIA_PID1_wide <- dcast(AMELIA_PID1, HID + HHS ~ PID, value.var = c("AGE", "AGE_cat", "MST", "SEX")) 
# Transformation ins wide-Format

start_table <- table(AMELIA_PID1_wide[, c("AGE_1", "AGE_cat_1", "MST_1", "SEX_1")]) # Diese Tabelle könnte als Basis dienen
AMELIA_wide_start <- AMELIA_PID1_wide[, c("HID", "AGE_1", "AGE_cat_1", "MST_1", "SEX_1")]

AMELIA_PID2 <- AMELIA[PID == 1 | PID == 2] # AMELIA-Datensatz mit ältester und zweitältester Person im HH
# Transformation ins wide-Format
AMELIA_PID2_wide <- dcast(AMELIA_PID2, HID + HHS ~ PID, value.var = c("AGE", "AGE_cat", "MST", "SEX")) 
AMELIA_PID2_wide <- AMELIA_PID2_wide[as.numeric(HHS) > 1] # Nur mindestens ZWeipersonenhaushalte

AMELIA_PID3 <- AMELIA[PID == 1 | PID == 2 | PID == 3]
# Transformation ins wide-Format
AMELIA_PID3_wide <- dcast(AMELIA_PID3, HID + HHS ~ PID, value.var = c("AGE", "AGE_cat", "MST", "SEX")) 
AMELIA_PID3_wide <- AMELIA_PID3_wide[as.numeric(HHS) > 2] 

AMELIA_PID4 <- AMELIA[PID == 1 | PID == 2 | PID == 3 | PID == 4]
# Transformation ins wide-Format
AMELIA_PID4_wide <- dcast(AMELIA_PID4, HID + HHS ~ PID, value.var = c("AGE", "AGE_cat", "MST", "SEX")) 
AMELIA_PID4_wide <- AMELIA_PID4_wide[as.numeric(HHS) > 3] 

AMELIA_PID5 <- AMELIA[PID == 1 | PID == 2 | PID == 3 | PID == 4 | PID == 5]
# Transformation ins wide-Format
AMELIA_PID5_wide <- dcast(AMELIA_PID5, HID + HHS ~ PID, value.var = c("AGE", "AGE_cat", "MST", "SEX")) 
AMELIA_PID5_wide <- AMELIA_PID5_wide[as.numeric(HHS) > 4] 

AMELIA_PID6 <- AMELIA[PID == 1 | PID == 2 | PID == 3 | PID == 4 | PID == 5 | PID == 6]
# Transformation ins wide-Format
AMELIA_PID6_wide <- dcast(AMELIA_PID6, HID + HHS ~ PID, value.var = c("AGE", "AGE_cat", "MST", "SEX")) 
AMELIA_PID6_wide <- AMELIA_PID6_wide[as.numeric(HHS) > 5] 

SI_SAMPLE <- srswor(5000, nrow(AMELIA)) 
AMELIA_sample <- AMELIA[SI_SAMPLE>0, ]
head(AMELIA_sample)

SI_SAMPLE <- srswor(5000, nrow(AMELIA_PID1_wide)) 
AMELIA_PID1_wide_sample <- AMELIA_PID1_wide[SI_SAMPLE>0, ]
head(AMELIA_PID1_wide_sample)

SI_SAMPLE <- srswor(5000, nrow(AMELIA_PID2_wide)) 
AMELIA_PID2_wide_sample <- AMELIA_PID2_wide[SI_SAMPLE>0, ]
head(AMELIA_PID2_wide_sample)

SI_SAMPLE <- srswor(5000, nrow(AMELIA_PID3_wide)) 
AMELIA_PID3_wide_sample <- AMELIA_PID3_wide[SI_SAMPLE>0, ]
head(AMELIA_PID3_wide_sample)

SI_SAMPLE <- srswor(5000, nrow(AMELIA_PID4_wide)) 
AMELIA_PID4_wide_sample <- AMELIA_PID4_wide[SI_SAMPLE>0, ]
head(AMELIA_PID4_wide_sample)

SI_SAMPLE <- srswor(5000, nrow(AMELIA_PID5_wide)) 
AMELIA_PID5_wide_sample <- AMELIA_PID5_wide[SI_SAMPLE>0, ]
head(AMELIA_PID5_wide_sample)

SI_SAMPLE <- srswor(5000, nrow(AMELIA_PID6_wide)) 
AMELIA_PID6_wide_sample <- AMELIA_PID6_wide[SI_SAMPLE>0, ]
head(AMELIA_PID6_wide_sample)

# To be able to apply our model later on, we need a starting point. Since the chatacericts of
# the individuals are estimated like a pyramid (starting with the oldest one) it has to be determined
# who is the oldest person in a household. Therefore, we will create an extra model which
# predicts who is the oldest person in a household

# Generate dichotomous variable 
AMELIA$PIDE <- ifelse(AMELIA$PID == 1, 1, 0)
AMELIA$PIDE <- as.factor(AMELIA$PIDE)

nn_PIDE <- nnet(PIDE ~ AGE + SEX + MST, size = 2, data = AMELIA_sample)
rf_PIDE <- randomForest(PIDE ~ AGE + SEX + MST, data = AMELIA_sample)
plot(rf_PIDE)

# Build model

start_table <- table(AMELIA_PID1_wide[, c("AGE_1", "AGE_cat_1", "MST_1", "SEX_1")]) 
AMELIA_wide_start <- AMELIA_PID1_wide[, c("HID", "AGE_1", "AGE_cat_1", "MST_1", "SEX_1")]

AMELIA_wide_start_nnet <- AMELIA_wide_start
AMELIA_wide_start_rpart <- AMELIA_wide_start
AMELIA_wide_start_rf <- AMELIA_wide_start

# HHS: Haushaltsgröße
nn_hhs <- nnet(HHS ~ AGE_1 + AGE_cat_1 + MST_1 + SEX_1, size = 2, data = AMELIA_PID1_wide_sample)
p_nn_hhs <- predict(nn_hhs, newdata = AMELIA_wide_start_nnet, type = "raw")
HHS <- apply(p_nn_hhs, 1, function(x)sample(colnames(p_nn_hhs), size = 1, prob = x))
AMELIA_wide_start_nnet[, HHS := as.factor(HHS)] 

rf_hhs <- randomForest(HHS ~ AGE_1 + AGE_cat_1 + MST_1 + SEX_1, data = AMELIA_PID1_wide_sample)
p_rf_hhs <- predict(rf_hhs, newdata = AMELIA_wide_start_nnet, type = "prob")
HHS <- apply(p_rf_hhs, 1, function(x)sample(colnames(p_rf_hhs), size = 1, prob = x))
AMELIA_wide_start_rf[, HHS := as.factor(HHS)] 

AMELIA_wide_start_nnet$HHS_True <- AMELIA_PID1$HHS
AMELIA_wide_start_rf$HHS_True <- AMELIA_PID1$HHS

plot(AMELIA_PID1$HHS)
plot(as.factor(AMELIA_wide_start_nnet$HHS))
plot(as.factor(AMELIA_wide_start_rf$HHS))

nn_AGE_cat_2 <- nnet(AGE_cat_2 ~ AGE_1 + AGE_cat_1 + MST_1 + SEX_1 + HHS, size = 2, 
                     data = AMELIA_PID2_wide_sample)
p_nn_AGE_cat_2 <- predict(nn_AGE_cat_2, newdata = AMELIA_wide_start_nnet, type = "raw")
AGE_cat_2 <- apply(p_nn_AGE_cat_2, 1, function(x)sample(colnames(p_nn_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_nnet[, AGE_cat_2 := as.factor(AGE_cat_2)] 

rf_AGE_cat_2 <- randomForest(AGE_cat_2 ~ AGE_cat_1 + MST_1 + SEX_1 + HHS, 
                             data = AMELIA_PID2_wide_sample)
p_rf_AGE_cat_2 <- predict(rf_AGE_cat_2, newdata = AMELIA_wide_start_rf, type = "prob")
AGE_cat_2 <- apply(p_rf_AGE_cat_2, 1, function(x)sample(colnames(p_rf_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_rf[, AGE_cat_2 := as.factor(AGE_cat_2)] 

plot(rf_AGE_cat_2)

nn_SEX_2 <- nnet(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + HHS, size = 2, data = AMELIA_PID2_wide_sample)
p_nn_SEX_2 <- predict(nn_SEX_2, newdata = AMELIA_wide_start_nnet, type = "raw")
p_nn_SEX_2_2 <- 1 - p_nn_SEX_2
p_nn_SEX_2 <- cbind(p_nn_SEX_2, p_nn_SEX_2_2)
colnames(p_nn_SEX_2) <- c(1, 2)
SEX_2 <- apply(p_nn_SEX_2, 1, function(x)sample(colnames(p_nn_SEX_2), size = 1, prob = x))
AMELIA_wide_start_nnet[, SEX_2 := as.factor(SEX_2)] 

rf_SEX_2 <- randomForest(SEX_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + HHS, data = AMELIA_PID2_wide_sample)
p_rf_SEX_2 <- predict(rf_SEX_2, newdata = AMELIA_wide_start_rf, type = "prob")
SEX_2 <- apply(p_rf_SEX_2, 1, function(x)sample(colnames(p_rf_SEX_2), size = 1, prob = x))
AMELIA_wide_start_rf[, SEX_2 := as.factor(SEX_2)] 

nn_MST_2 <- nnet(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2 + HHS, size = 2, data = AMELIA_PID2_wide_sample)
p_nn_MST_2 <- predict(nn_MST_2, newdata = AMELIA_wide_start_nnet, type = "raw")
MST_2 <- apply(p_nn_MST_2, 1, function(x)sample(colnames(p_nn_MST_2), size = 1, prob = x))
AMELIA_wide_start_nnet[, MST_2 := as.factor(MST_2)] 

rf_MST_2 <- randomForest(MST_2 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + SEX_1 + SEX_2 + HHS, data = AMELIA_PID2_wide_sample)
p_rf_MST_2 <- predict(rf_MST_2, newdata = AMELIA_wide_start_rf, type = "prob")
MST_2_rf <- apply(p_rf_MST_2, 1, function(x)sample(colnames(p_rf_MST_2), size = 1, prob = x))
AMELIA_wide_start_rf[, MST_2 := as.factor(MST_2_rf)] 
# 3

nn_AGE_cat_3 <- nnet(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2 + HHS, 
                     size = 2, data = AMELIA_PID3_wide_sample)
p_nn_AGE_cat_3 <- predict(nn_AGE_cat_3, newdata = AMELIA_wide_start_nnet, type = "raw")
AGE_cat_3 <- apply(p_nn_AGE_cat_3, 1, function(x)sample(colnames(p_nn_AGE_cat_2), size = 1, prob = x))
AMELIA_wide_start_nnet[, AGE_cat_3 := as.factor(AGE_cat_3)] 

rf_AGE_cat_3 <- randomForest(AGE_cat_3 ~ AGE_cat_1 + AGE_cat_2 + MST_1 + MST_2 + SEX_1 + SEX_2 + HHS, 
                     data = AMELIA_PID3_wide_sample)
p_rf_AGE_cat_3 <- predict(rf_AGE_cat_3, newdata = AMELIA_wide_start_rf, type = "prob")
AGE_cat_3 <- apply(p_rf_AGE_cat_3, 1, function(x)sample(colnames(p_rf_AGE_cat_3), size = 1, prob = x))
AMELIA_wide_start_rf[, AGE_cat_3 := as.factor(AGE_cat_3)] 

nn_SEX_3 <- nnet(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + HHS, 
                              size = 2, data = AMELIA_PID3_wide_sample)
p_nn_SEX_3 <- predict(nn_SEX_3, newdata = AMELIA_wide_start_nnet, type = "raw")
p_nn_SEX_3 <- predict(nn_SEX_3, newdata = AMELIA_wide_start_nnet, type = "raw")
p_nn_SEX_3_2 <- 1 - p_nn_SEX_3
p_nn_SEX_3 <- cbind(p_nn_SEX_3, p_nn_SEX_3_2)
colnames(p_nn_SEX_3) <- c(1, 2)
SEX_3 <- apply(p_nn_SEX_3, 1, function(x)sample(colnames(p_nn_SEX_3), size = 1, prob = x))
AMELIA_wide_start_nnet[, SEX_3 := as.factor(SEX_3)]

rf_SEX_3 <- randomForest(SEX_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + HHS, 
              data = AMELIA_PID3_wide_sample)
p_rf_SEX_3 <- predict(rf_SEX_3, newdata = AMELIA_wide_start_rf, type = "prob")
SEX_3 <- apply(p_rf_SEX_3, 1, function(x)sample(colnames(p_rf_SEX_3), size = 1, prob = x))
AMELIA_wide_start_rf[, SEX_3 := as.factor(SEX_3)] 


nn_MST_3 <- nnet(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3 + HHS, 
                 size = 2, data = AMELIA_PID3_wide_sample)
p_nn_MST_3 <- predict(nn_MST_3, newdata = AMELIA_wide_start_nnet, type = "raw")
MST_3 <- apply(p_nn_MST_3, 1, function(x)sample(colnames(p_nn_MST_3), size = 1, prob = x))
AMELIA_wide_start_nnet[, MST_3 := as.factor(MST_3)]

rf_MST_3 <- randomForest(MST_3 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + SEX_1 + SEX_2 + SEX_3 + HHS, 
                         data = AMELIA_PID3_wide_sample)
p_rf_MST_3 <- predict(rf_MST_3, newdata = AMELIA_wide_start_rf, type = "prob")
MST_3_rf <- apply(p_rf_MST_3, 1, function(x)sample(colnames(p_rf_MST_3), size = 1, prob = x))
AMELIA_wide_start_rf[, MST_3 := as.factor(MST_3_rf)] 

# 4

nn_AGE_cat_4 <- nnet(AGE_cat_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + MST_3 +
                       SEX_1 + SEX_2 + SEX_3 + HHS, 
                       size = 2, data = AMELIA_PID4_wide_sample)
p_nn_AGE_cat_4 <- predict(nn_AGE_cat_4, newdata = AMELIA_wide_start_nnet, type = "raw")
AGE_cat_4 <- apply(p_nn_AGE_cat_4, 1, function(x)sample(colnames(p_nn_AGE_cat_4), size = 1, prob = x))
AMELIA_wide_start_nnet[, AGE_cat_4 := as.factor(AGE_cat_4)]

rf_AGE_cat_4 <- randomForest(AGE_cat_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + MST_1 + MST_2 + MST_3 +
                       SEX_1 + SEX_2 +  SEX_3 + HHS, 
                       data = AMELIA_PID4_wide_sample)
p_rf_AGE_cat_4 <- predict(rf_AGE_cat_4, newdata = AMELIA_wide_start_rf, type = "prob")
AGE_cat_4 <- apply(p_rf_AGE_cat_4, 1, function(x)sample(colnames(p_rf_AGE_cat_4), size = 1, prob = x))
AMELIA_wide_start_rf[, AGE_cat_4 := as.factor(AGE_cat_4)] 

nn_SEX_4 <- nnet(SEX_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                       SEX_1 + SEX_2 +  SEX_3 + HHS, 
                     size = 2, data = AMELIA_PID4_wide_sample)
p_nn_SEX_4 <- predict(nn_SEX_4, newdata = AMELIA_wide_start_nnet, type = "raw")
p_nn_SEX_4_2 <- 1 - p_nn_SEX_4
p_nn_SEX_4 <- cbind(p_nn_SEX_4, p_nn_SEX_4_2)
colnames(p_nn_SEX_4) <- c(1, 2)
SEX_4 <- apply(p_nn_SEX_4, 1, function(x)sample(colnames(p_nn_SEX_4), size = 1, prob = x))
AMELIA_wide_start_nnet[, SEX_4 := as.factor(SEX_4)]

rf_SEX_4 <- randomForest(SEX_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                   SEX_1 + SEX_2 +  SEX_3 + HHS, 
                   data = AMELIA_PID4_wide_sample)
p_rf_SEX_4 <- predict(rf_SEX_4, newdata = AMELIA_wide_start_rf, type = "prob")
SEX_4 <- apply(p_rf_SEX_4, 1, function(x)sample(colnames(p_rf_SEX_4), size = 1, prob = x))
AMELIA_wide_start_rf[, SEX_4 := as.factor(SEX_4)] 

nn_MST_4 <- nnet(MST_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                   SEX_1 + SEX_2 +  SEX_3 + SEX_4 + HHS, 
                 size = 2, data = AMELIA_PID4_wide_sample)
p_nn_MST_4 <- predict(nn_MST_4, newdata = AMELIA_wide_start_nnet, type = "raw")
MST_4 <- apply(p_nn_MST_4, 1, function(x)sample(colnames(p_nn_MST_4), size = 1, prob = x))
AMELIA_wide_start_nnet[, MST_4 := as.factor(MST_4)]

rf_MST_4 <- randomForest(MST_4 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                           SEX_1 + SEX_2 +  SEX_3 + SEX_4 + HHS, 
                   data = AMELIA_PID4_wide_sample)
p_rf_MST_4 <- predict(rf_MST_4, newdata = AMELIA_wide_start_rf, type = "prob")
MST_4_rf <- apply(p_rf_MST_4, 1, function(x)sample(colnames(p_rf_MST_4), size = 1, prob = x))
AMELIA_wide_start_rf[, MST_4 := as.factor(MST_4_rf)] 

# 5

nn_AGE_cat_5 <- nnet(AGE_cat_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                       MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + HHS, 
                     size = 2, data = AMELIA_PID5_wide_sample)
p_nn_AGE_cat_5 <- predict(nn_AGE_cat_5, newdata = AMELIA_wide_start_nnet, type = "raw")
AGE_cat_5 <- apply(p_nn_AGE_cat_5, 1, function(x)sample(colnames(p_nn_AGE_cat_5), size = 1, prob = x))
AMELIA_wide_start_nnet[, AGE_cat_5 := as.factor(AGE_cat_5)]

rf_AGE_cat_5 <- randomForest(AGE_cat_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + MST_1 + MST_2 + MST_3 +
                       MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + HHS, 
                       data = AMELIA_PID5_wide_sample)
p_rf_AGE_cat_5 <- predict(rf_AGE_cat_5, newdata = AMELIA_wide_start_rf, type = "prob")
AGE_cat_5 <- apply(p_rf_AGE_cat_5, 1, function(x)sample(colnames(p_rf_AGE_cat_5), size = 1, prob = x))
AMELIA_wide_start_rf[, AGE_cat_5 := as.factor(AGE_cat_5)] 

nn_SEX_5 <- nnet(SEX_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                       MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + HHS, 
                     size = 2, data = AMELIA_PID5_wide_sample)
p_nn_SEX_5 <- predict(nn_SEX_5, newdata = AMELIA_wide_start_nnet, type = "raw")
p_nn_SEX_5_2 <- 1 - p_nn_SEX_5
p_nn_SEX_5 <- cbind(p_nn_SEX_5, p_nn_SEX_5_2)
colnames(p_nn_SEX_5) <- c(1, 2)
SEX_5 <- apply(p_nn_SEX_5, 1, function(x)sample(colnames(p_nn_SEX_5), size = 1, prob = x))
AMELIA_wide_start_nnet[, SEX_5 := as.factor(SEX_5)]

rf_SEX_5 <- randomForest(SEX_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                           MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + HHS, 
                   data = AMELIA_PID5_wide_sample)
p_rf_SEX_5 <- predict(rf_SEX_5, newdata = AMELIA_wide_start_rf, type = "prob")
SEX_5 <- apply(p_rf_SEX_5, 1, function(x)sample(colnames(p_rf_SEX_5), size = 1, prob = x))
AMELIA_wide_start_rf[, SEX_5 := as.factor(SEX_5)] 

nn_MST_5 <- nnet(MST_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                   MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5 + HHS, 
                 size = 2, data = AMELIA_PID5_wide_sample)
p_nn_MST_5 <- predict(nn_MST_5, newdata = AMELIA_wide_start_nnet, type = "raw")
MST_5 <- apply(p_nn_MST_5, 1, function(x)sample(colnames(p_nn_MST_5), size = 1, prob = x))
AMELIA_wide_start_nnet[, MST_5 := as.factor(MST_5)]

rf_MST_5 <- randomForest(MST_5 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + MST_1 + MST_2 + MST_3 +
                           MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5 + HHS, 
                   data = AMELIA_PID5_wide_sample)
p_rf_MST_5 <- predict(rf_MST_5, newdata = AMELIA_wide_start_rf, type = "prob")
MST_5_rf <- apply(p_rf_MST_5, 1, function(x)sample(colnames(p_rf_MST_5), size = 1, prob = x))
AMELIA_wide_start_rf[, MST_5 := as.factor(MST_5_rf)] 

# 6

nn_AGE_cat_6 <- nnet(AGE_cat_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 +
                       MST_1 + MST_2 + MST_3 +
                       MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5 + HHS, 
                     size = 2, data = AMELIA_PID6_wide_sample)
p_nn_AGE_cat_6 <- predict(nn_AGE_cat_6, newdata = AMELIA_wide_start_nnet, type = "raw")
AGE_cat_6 <- apply(p_nn_AGE_cat_6, 1, function(x)sample(colnames(p_nn_AGE_cat_6), size = 1, prob = x))
AMELIA_wide_start_nnet[, AGE_cat_6 := as.factor(AGE_cat_6)]

rf_AGE_cat_6 <- randomForest(AGE_cat_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 +
                               MST_1 + MST_2 + MST_3 +
                               MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5 + HHS, 
                             data = AMELIA_PID6_wide_sample)
p_rf_AGE_cat_6 <- predict(rf_AGE_cat_6, newdata = AMELIA_wide_start_rf, type = "prob")
AGE_cat_6 <- apply(p_rf_AGE_cat_6, 1, function(x)sample(colnames(p_rf_AGE_cat_6), size = 1, prob = x))
AMELIA_wide_start_rf[, AGE_cat_6 := as.factor(AGE_cat_6)] 

nn_SEX_6 <- nnet(SEX_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + AGE_cat_6 +
                       MST_1 + MST_2 + MST_3 +
                       MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5 + HHS, 
                     size = 2, data = AMELIA_PID6_wide_sample)
p_nn_SEX_6 <- predict(nn_SEX_6, newdata = AMELIA_wide_start_nnet, type = "raw")
p_nn_SEX_6_2 <- 1 - p_nn_SEX_6
p_nn_SEX_6 <- cbind(p_nn_SEX_6, p_nn_SEX_6_2)
colnames(p_nn_SEX_6) <- c(1, 2)
SEX_6 <- apply(p_nn_SEX_6, 1, function(x)sample(colnames(p_nn_SEX_6), size = 1, prob = x))
AMELIA_wide_start_nnet[, SEX_6 := as.factor(SEX_6)]

rf_SEX_6 <- randomForest(SEX_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 +
                           AGE_cat_6 + MST_1 + MST_2 + MST_3 +
                           MST_4 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5 + HHS, 
                         data = AMELIA_PID6_wide_sample)
p_rf_SEX_6 <- predict(rf_SEX_6, newdata = AMELIA_wide_start_rf, type = "prob")
SEX_6 <- apply(p_rf_SEX_6, 1, function(x)sample(colnames(p_rf_SEX_6), size = 1, prob = x))
AMELIA_wide_start_rf[, SEX_6 := as.factor(SEX_6)] 

nn_MST_6 <- nnet(MST_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + AGE_cat_6 +
                   MST_1 + MST_2 + MST_3 +
                   MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5 + SEX_6 + HHS, 
                 size = 2, data = AMELIA_PID6_wide_sample)
p_nn_MST_6 <- predict(nn_MST_6, newdata = AMELIA_wide_start_nnet, type = "raw")
MST_6 <- apply(p_nn_MST_6, 1, function(x)sample(colnames(p_nn_MST_6), size = 1, prob = x))
AMELIA_wide_start_nnet[, MST_6 := as.factor(MST_6)]

rf_MST_6 <- randomForest(MST_6 ~ AGE_cat_1 + AGE_cat_2 + AGE_cat_3 + AGE_cat_4 + AGE_cat_5 + AGE_cat_6 +
                           MST_1 + MST_2 + MST_3 +
                           MST_4 + MST_5 + SEX_1 + SEX_2 +  SEX_3 + SEX_4 + SEX_5 + SEX_6 + HHS, 
                         data = AMELIA_PID6_wide_sample)
p_rf_MST_6 <- predict(rf_MST_6, newdata = AMELIA_wide_start_rf, type = "prob")
MST_6_rf <- apply(p_rf_MST_6, 1, function(x)sample(colnames(p_rf_MST_6), size = 1, prob = x))
AMELIA_wide_start_rf[, MST_6 := as.factor(MST_6_rf)] 

AMELIA_wide_start_nnet_1 <- AMELIA_wide_start_nnet[HHS == 1]
AMELIA_wide_start_nnet_2 <- AMELIA_wide_start_nnet[HHS == 2]
AMELIA_wide_start_nnet_3 <- AMELIA_wide_start_nnet[HHS == 3]
AMELIA_wide_start_nnet_4 <- AMELIA_wide_start_nnet[HHS == 4]
AMELIA_wide_start_nnet_5 <- AMELIA_wide_start_nnet[HHS == 5]
AMELIA_wide_start_nnet_6 <- AMELIA_wide_start_nnet[HHS == 6]

AMELIA_wide_start_rf_1 <- AMELIA_wide_start_rf[HHS == 1]
AMELIA_wide_start_rf_2 <- AMELIA_wide_start_rf[HHS == 2]
AMELIA_wide_start_rf_3 <- AMELIA_wide_start_rf[HHS == 3]
AMELIA_wide_start_rf_4 <- AMELIA_wide_start_rf[HHS == 4]
AMELIA_wide_start_rf_5 <- AMELIA_wide_start_rf[HHS == 5]
AMELIA_wide_start_rf_6 <- AMELIA_wide_start_rf[HHS == 6]

AMELIA_wide_start_nnet_2 <- calc_difference(AMELIA_wide_start_nnet_2, hhs_size = 2)
AMELIA_wide_start_nnet_3 <- calc_difference(AMELIA_wide_start_nnet_3, hhs_size = 3)
AMELIA_wide_start_nnet_4 <- calc_difference(AMELIA_wide_start_nnet_4, hhs_size = 4)
AMELIA_wide_start_nnet_5 <- calc_difference(AMELIA_wide_start_nnet_5, hhs_size = 5)
AMELIA_wide_start_nnet_6 <- calc_difference(AMELIA_wide_start_nnet_6, hhs_size = 6)

AMELIA_wide_start_rf_2 <- calc_difference(AMELIA_wide_start_rf_2, hhs_size = 2)
AMELIA_wide_start_rf_3 <- calc_difference(AMELIA_wide_start_rf_3, hhs_size = 3)
AMELIA_wide_start_rf_4 <- calc_difference(AMELIA_wide_start_rf_4, hhs_size = 4)
AMELIA_wide_start_rf_5 <- calc_difference(AMELIA_wide_start_rf_5, hhs_size = 5)
AMELIA_wide_start_rf_6 <- calc_difference(AMELIA_wide_start_rf_6, hhs_size = 6)

calc_bias(AMELIA_wide_start_nnet_2, true = hhs_2_true)
calc_bias(AMELIA_wide_start_nnet_3, true = hhs_3_true)
calc_bias(AMELIA_wide_start_nnet_4, true = hhs_4_true)
calc_bias(AMELIA_wide_start_nnet_5, true = hhs_5_true)
calc_bias(AMELIA_wide_start_nnet_6, true = hhs_6_true)

calc_bias(AMELIA_wide_start_rf_2, true = hhs_2_true)
calc_bias(AMELIA_wide_start_rf_3, true = hhs_3_true)
calc_bias(AMELIA_wide_start_rf_4, true = hhs_4_true)
calc_bias(AMELIA_wide_start_rf_5, true = hhs_5_true)
calc_bias(AMELIA_wide_start_rf_6, true = hhs_6_true)

# Generate categories to create a graph like the one in Kolb (2012). Therefore,
# let´s create a counter for certain relative frequencies which are based on connected
# attributes over several variables

colb_creater <- function(x) {
# Two married adults, both below age 65
x$def_1 <- 0
x$def_1[as.numeric(x$AGE_cat_1) <= 6 & as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) <= 6
                                  & x$MST_1 == 2 & x$MST_2 == 2] <- 1


# Two married adults, at least one above 65

x$def_1[as.numeric(x$AGE_cat_1) >= 6 & as.numeric(x$AGE_cat_2) > 1
                 & x$MST_1 == 2 & x$MST_2 == 2] <- 2

# Two adults, both never married below age 65
x$def_1[as.numeric(x$AGE_cat_1) < 6 & as.numeric(x$AGE_cat_2) < 6
                 & x$MST_1 == 1 & x$MST_2 == 1] <- 3

# One adult, one person below 18 years
x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) < 2] <- 4

# Two adults, one divorced, one never married

x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1 &
                   x$MST_group == "15"] <- 5

return(x)
}

colb_creater_4 <- function(x) {
  # Two men, two woman, all at least 18 years old
  x$def_1 <- 0
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
            x$SEX_diff == "1.1.2.2"] <- 1
  
  
  # Three man, one woman, all at least 18 years old
  
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
            x$SEX_diff == "1.1.1.2"] <- 2
  
  # Only men or only woman, all at least 18 years old
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) > 1 &
            x$SEX_diff == "1.1.1.1"|x$SEX_diff == "2.2.2.2"] <- 3
  
  # Three adults one person below 17 years
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) > 1 & as.numeric(x$AGE_cat_4) < 2 ] <- 4
  
  # One adult and three persons below 17 years
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) < 2
          & as.numeric(x$AGE_cat_3) < 2 & as.numeric(x$AGE_cat_4) < 2 ] <- 5
  
  # Two adults and two persons below 17 years
  x$def_1[as.numeric(x$AGE_cat_1) > 1 & as.numeric(x$AGE_cat_2) > 1
          & as.numeric(x$AGE_cat_3) < 2 & as.numeric(x$AGE_cat_4) < 2 ] <- 6
  
  
  
  return(x)
}

test <- colb_creater_4(hhs_4_true)


hhs_4_true <- colb_creater_4(hhs_4_true)
synthetic_hhs4_nn_colb <- colb_creater_4(test_4_vt$synthetic_hhs_nn)
synthetic_hhs4_rf_colb <- colb_creater_4(test_4_vt$synthetic_hhs_rf)
synthetic_hhs4_nn_mod_colb <- colb_creater_4(AMELIA_wide_start_nnet_4)
synthetic_hhs4_rf_mod_colb <- colb_creater_4(AMELIA_wide_start_rf_4)

hhs4_t <- as.data.frame(prop.table(table(hhs_4_true$def_1)))
prop.table(table(synthetic_hhs4_nn_colb$def_1))
prop.table(table(synthetic_hhs4_rf_colb$def_1))
prop.table(table(synthetic_hhs4_nn_mod_colb$def_1))
prop.table(table(synthetic_hhs4_rf_mod_colb$def_1))

prop.table(table(test$def_1))

hhs_2_true <- colb_creater(hhs_2_true)
synthetic_hhs2_nn_colb <- colb_creater(test_2_t$synthetic_hhs_nn)
synthetic_hhs2_rf_colb <- colb_creater(test_2_t$synthetic_hhs_rf)
synthetic_hhs2_nn_mod_colb <- colb_creater(AMELIA_wide_start_nnet_2)
synthetic_hhs2_rf_mod_colb <- colb_creater(AMELIA_wide_start_rf_2)

prop.table(table(hhs_2_true$def_1))
prop.table(table(synthetic_hhs2_nn_colb$def_1))
prop.table(table(synthetic_hhs2_rf_colb$def_1))
prop.table(table(synthetic_hhs2_nn_mod_colb$def_1))
prop.table(table(synthetic_hhs2_rf_mod_colb$def_1))

# library
library(ggplot2)

hhs2_t <- data.frame(prop.table(table(hhs_2_true$def_1)))
hhs2_c_nn <- data.frame(prop.table(table(synthetic_hhs2_nn_colb$def_1)))
hhs2_c_rf <- data.frame(prop.table(table(synthetic_hhs2_rf_colb$def_1)))
hhs2_m_nn <- data.frame(prop.table(table(synthetic_hhs2_nn_mod_colb$def_1)))
hhs2_m_rf <- data.frame(prop.table(table(synthetic_hhs2_rf_mod_colb$def_1)))

hhs4_t <- data.frame(prop.table(table(hhs_4_true$def_1)))
hhs4_c_nn <- data.frame(prop.table(table(synthetic_hhs4_nn_colb$def_1)))
hhs4_c_rf <- data.frame(prop.table(table(synthetic_hhs4_rf_colb$def_1)))
hhs4_m_nn <- data.frame(prop.table(table(synthetic_hhs4_nn_mod_colb$def_1)))
hhs4_m_rf <- data.frame(prop.table(table(synthetic_hhs4_rf_mod_colb$def_1)))



data <- rbind(hhs2_t, hhs2_c_nn, hhs2_c_rf, hhs2_m_nn, hhs2_m_rf)
data_2 <- rbind(hhs4_t, hhs4_c_nn, hhs4_c_rf, hhs4_m_nn, hhs4_m_rf)
data_2$Method[29:35] <- c("hhs4_m_rf")
data_2$Method <- as.factor(data_2$Method)

data$Method[25:30] <- c("hhs2_m_rf")
data$Method <- as.factor(data$Method)

colnames(data) <- c("hhs2_t", "Freq", "hhs2_c_nn", "Freq", "hhs_c_rf", "Freq", 
                    "hhs2_m_nn", "Freq", "hhs2_m_rf", "Freq")


#
data_2$Var1 <- factor(data_2$Var1, labels = c("Other", 
                                              "Two men, two woman, all at least 18 years old",
                                              "Three man, one woman, all at least 18 years old",
                                              "Only men or only woman, all at least 18 years old",
                                          "Three adults one person below 17 years",
                                          "One adult and three persons below 17 years",
                                          "Two adults and two persons below 17 years"))
colnames(data_2) <- c("Household_Type", "Frequency", "Method")

kolb_plot_1 <-  ggplot(data, aes(fill=Household_Type, y=Frequency, x=Method)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  theme(legend.position="bottom") +
  theme_bw()

kolb_plot_2 <- ggplot(data_2, aes(fill=Household_Type, y=Frequency, x=Method)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  theme(legend.position="bottom") +
  theme_bw()

kolb_plot_1
kolb_plot_2
# Some ideas for later about fine tuning the random forest
set.seed(100)
train <- sample(nrow(model_data_c), 0.7*nrow(model_data_c), replace = FALSE)
TrainSet <- model_data_c[train,]
ValidSet <- model_data_c[-train,]
summary(TrainSet)
summary(ValidSet)

# Create a Random Forest model with default parameters
model1 <- randomForest(TRUE_ ~ AGE_diff + SEX_diff + MST_group, data = TrainSet, importance = TRUE)
model1

# Fine tuning parameters of Random Forest model
model2 <- randomForest(TRUE_ ~ AGE_diff + SEX_diff + MST_group, data = TrainSet, ntree = 1000,
                       importance = TRUE)
model2

# Predicting on train set
predTrain <- predict(model2, TrainSet, type = "class")
# Checking classification accuracy
table(predTrain, TrainSet$TRUE_)  

# Predicting on Validation set
predValid <- predict(model2, ValidSet, type = "class")
# Checking classification accuracy
mean(predValid == ValidSet$TRUE_)                    
table(predValid,ValidSet$TRUE_)

# To check important variables
importance(model2)        
varImpPlot(model2) 

