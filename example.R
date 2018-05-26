# Bootstrap Magnitude Based Predictions (BMBP)
# Copyright 2018 Mladen Jovanovic
#
# To reference plese use:
#
#      Jovanovic, M. (2018, May, 25). Bootstrap Magnitude Based Predictions. 
#          Complementary Training. Retrieved from http://www.complementarytraining.com/bmbp
#
#

# This is Example file
# Using Vertical Jump Height as a test

source("bMBP.r")

set.seed(6667)
# Generate data
# Vertical jump height

N = 25
athletes <- sprintf("Athlete_%03d", seq(1, N))
verticalJumpReal <- rnorm(n = N, mean = 40, sd = 8)

# Two devices
verticalJumpTestPre1 <- verticalJumpReal + rnorm(n = N, mean = 0, sd = 0.5)  
verticalJumpTestPre2 <- verticalJumpReal + rnorm(n = N, mean = 0, sd = 0.5)

df <- data.frame(Athlete = athletes,
                 Pre1 = verticalJumpTestPre1,
                 Pre2 = verticalJumpTestPre2)

vjRel.two <- bootReliability(df, "Pre1", "Pre2", samples = 200, plot = TRUE)


# Multiple devices
verticalJumpTestPre1 <- verticalJumpReal + rnorm(n = N, mean = 0, sd = 0.5)  
verticalJumpTestPre2 <- verticalJumpReal + rnorm(n = N, mean = 0, sd = 0.5)
verticalJumpTestPre3 <- verticalJumpReal + rnorm(n = N, mean = 3, sd = 0) 
verticalJumpTestPre4 <- verticalJumpReal + rnorm(n = N, mean = 0, sd = 2) 

df <- data.frame(Athlete = athletes,
                 Pre1 = verticalJumpTestPre1,
                 Pre2 = verticalJumpTestPre2,
                 Pre3 = verticalJumpTestPre3,
                 Pre4 = verticalJumpTestPre4)


vjRel.multiple <- bootReliabilityMultiple(df[-1], samples = 200, plot = TRUE)

# Multiple days
verticalJumpTestDay1 <- verticalJumpReal + rnorm(n = N, mean = 0, sd = 2)  
verticalJumpTestDay2 <- verticalJumpTestDay1 + rnorm(n = N, mean = 2, sd = 1.5)
verticalJumpTestDay3 <- verticalJumpTestDay2 + rnorm(n = N, mean = 1, sd = 1) 
verticalJumpTestDay4 <- verticalJumpTestDay3 + rnorm(n = N, mean = 0, sd = 0.5) 

df <- data.frame(Athlete = athletes,
                 Day1 = verticalJumpTestDay1,
                 Day2 = verticalJumpTestDay2,
                 Day3 = verticalJumpTestDay3,
                 Day4 = verticalJumpTestDay4)


vjRel.multiple <- bootReliabilityMultiple(df[-1], type = "sequential",
                                          samples = 200, plot = TRUE)

# Test on Hopkins Reliability sample data
df <- read.csv("hopkins-sample-reliability.csv", header = TRUE)

hopkins.multiple <- bootReliabilityMultiple(df, type = "sequential",
                                          samples = 200, plot = TRUE)

overall <- select(hopkins.multiple$overallSummary, Estimate, Value)
variable <- select(hopkins.multiple$variableSummary, Estimate, Variable, Value)
pairwise <- select(hopkins.multiple$pairwiseSummary, Estimate, Variable, Value)

################
# Change after interventions
verticalJumpTestPre <- rnorm(n = N, mean = 40, sd = 8)
verticalJumpTestPost <- verticalJumpTestPre + rnorm(n = N, mean = 5, sd = 6)

df <- data.frame(Athlete = athletes,
                 Pre = verticalJumpTestPre,
                 Post = verticalJumpTestPost)
vjChange <- bootChange(df, pre = "Pre", post = "Post", samples = 200, plot = TRUE, na.rm = TRUE)
