################################################################################
# Script: WormSimulation.R
# Author: Isaac Harris
#
# Purpose:
#   Simulate the dynamics of two strategies in a worm population:
#     1) A "transgenerational" strategy (strat = 1) that experiences
#        transgenerational fitness effects after starvation.
#     2) A "control" strategy (strat = 0) with no such effects.
#
#   The population is structured into cohorts (same age, size, strategy,
#   and starvation history). Starvation events trigger dispersal and founding
#   of a new population with a random subset of worms, carrying their
#   transgenerational state forward.
#
# Output:
#   - A data frame `dataframepop` describing population sizes and strategy
#     proportions over time.
#   - A ggplot object `Plot` showing the proportion of each strategy over time.
################################################################################

#############################
# 1. Load required packages #
#############################

library(ggplot2)
library(progress)

set.seed(123)

################################################################################
# 2. Parameters                                                                 #
################################################################################

# Initial cohorts: a list of lists, each representing a cohort of worms.
# Fields:
#   age        : age in days
#   size       : number of worms in the cohort
#   starvation : how many starvation events the lineage has experienced
#   strat      : strategy (1 = transgenerational, 0 = control)
cohorts <- list(
  list(age = 0, size = 1, starvation = 1, strat = 1),  # transgenerational founder
  list(age = 0, size = 1, starvation = 0, strat = 0)   # control founder
)

# Food supply
Food  <- 1.0e+8          # initial number of bacterial cells
Foods <- c(1.0e+8)       # (optional) track food over time if desired

# Food eaten per worm per day by age class:
#   indices: 1 = eggs/L1/L2, 2 = L3/L4, 3 = reproductive, 4 = post-reproductive
Eat <- c(1, 2, 4, 1)

# Transgenerational effects on reproduction for starvation generations:
#   Effects[1] : F1 (after 1 starvation) multiplier
#   Effects[2] : F2
#   Effects[3] : F3+
Effects <- c(2, 1.0, 0.5)

# Population tracking
Pop      <- c(0)  # total population over time
transPop <- 0     # transgenerational subpopulation over time
conPop   <- 0     # control subpopulation over time

# Age-specific reproduction schedule:
# indices correspond to age in days, e.g.
#   R[1] = L1/L2, R[2] = L3/L4, R[3] = Day 1 adult, etc.
R <- c(0, 0, 41, 125, 80, 18, 5) 

# Simulation horizon (days)
n_days <- 100

# Progress bar
pb <- progress_bar$new(total = n_days)

################################################################################
# 3. Simulation loop                                                            #
################################################################################

for (day in 1:n_days) {
  pb$tick()
  Sys.sleep(1 / 100)  # purely cosmetic for the progress bar
  
  # Track changes made this iteration
  new_cohorts   <- list()  # cohorts produced this day
  tpop          <- 0       # total population this day
  food_consumed <- 0       # total food consumed this day
  t_count       <- 0       # number of transgenerational worms
  c_count       <- 0       # number of control worms
  
  # Ageing and population calculations
  for (i in seq_along(cohorts)) {
    cohort <- cohorts[[i]]              # current cohort
    cohort$age <- cohort$age + 1        # all worms age by 1 day
    
    # Update total population
    tpop <- tpop + cohort$size
    
    # Split by strategy
    if (cohort$strat == 1) {
      t_count <- t_count + cohort$size
    } else {
      c_count <- c_count + cohort$size
    }
    
    # Food consumption by age class
    if (cohort$age == 1) {
      # Eggs/L1/L2
      food_consumed <- food_consumed + cohort$size * Eat[1]
    } else if (cohort$age == 2) {
      # L3/L4
      food_consumed <- food_consumed + cohort$size * Eat[2]
    } else if (cohort$age >= 3 && cohort$age <= 7) {
      # Reproductive worms
      food_consumed <- food_consumed + cohort$size * Eat[3]
    } else {
      # Post-reproductive worms
      food_consumed <- food_consumed + cohort$size * Eat[4]
    }
    
    # Reproduction (age-specific)
    if (cohort$age > 2 && cohort$age <= 7) {  # reproductive age window
      
      # Default modifier (no transgenerational effect)
      repro_modifier <- 1
      
      # Apply transgenerational effects for strategy 1 only
      if (cohort$starvation == 1 && cohort$strat == 1) {
        repro_modifier <- Effects[1]  # F1 effect
      } else if (cohort$starvation == 2 && cohort$strat == 1) {
        repro_modifier <- Effects[2]  # F2 effect
      } else if (cohort$starvation >= 3 && cohort$strat == 1) {
        repro_modifier <- Effects[3]  # F3+ effect
      }
      
      # Eggs laid per worm at this age
      rate <- R[cohort$age] * repro_modifier
      
      # New cohort of age 0; offspring inherit parental starvation history + 1
      new_cohorts <- append(
        new_cohorts,
        list(
          list(
            age        = 0,
            size       = floor(rate * cohort$size),
            starvation = cohort$starvation + 1,
            strat      = cohort$strat
          )
        )
      )
    }
    
    # Update cohort in the global list
    cohorts[[i]] <- cohort
  }
  
  # Update food supply
  Food <- Food - food_consumed
  
  # Add newly produced cohorts to population
  cohorts <- c(cohorts, new_cohorts)
  
  # Store total and strategy-specific population sizes
  Pop[day]      <- tpop
  transPop[day] <- t_count
  conPop[day]   <- c_count
  
  # Store population data in a data frame (updated each iteration)
  dataframepop <- data.frame(
    Populationsize   = Pop,
    Transgenerational = transPop,
    control          = conPop
  )
  
  # Optionally store food for debugging:
  # Foods[day] <- Food 
  
  # Starvation event: food depleted
  if (Food <= 0) {
    
    # Draw a random founding population size between 10 and 100
    n_founders <- sample(10:100, 1)
    
    # Proportion of founders that are transgenerational vs control
    nt <- n_founders * dataframepop$Transgenerational[day] /
      dataframepop$Populationsize[day]
    nc <- n_founders * dataframepop$control[day] /
      dataframepop$Populationsize[day]
    
    # New founding cohorts (age 1 because founding occurs at larval arrest)
    cohorts <- list(
      list(age = 1, size = round(nt), starvation = 1, strat = 1),
      list(age = 1, size = round(nc), starvation = 0, strat = 0)
    )
    
    # Randomise the new food supply
    Food <- rnorm(mean = 1e+8, sd = 1e+7, n = 1)
  }
}

################################################################################
# 4. Graphs                                                                     #
################################################################################

# Calculate proportions of each strategy
dataframepop$tpro <- dataframepop$Transgenerational / dataframepop$Populationsize
dataframepop$cpro <- dataframepop$control          / dataframepop$Populationsize
dataframepop$day  <- 1:n_days

# Plot proportion of each strategy over time
Plot <- ggplot(dataframepop, aes(x = day)) +
  geom_point(aes(y = tpro, colour = "Transgenerational"), size = 3) +
  geom_line(aes(y = tpro, colour = "Transgenerational"), size = 1.5) +
  geom_point(aes(y = cpro, colour = "Control"), size = 3) +
  geom_line(aes(y = cpro, colour = "Control"), size = 1.5) +
  scale_color_manual(
    name   = "Strategy",
    values = c("Transgenerational" = "#0072B2",
               "Control"           = "#D55E00")
  ) +
  xlab("Time (days)") +
  ylab("Proportion of population") +
  ylim(c(0, 1)) +
  theme_classic() +
  theme(
    axis.title   = element_text(size = 30),
    axis.text    = element_text(size = 25),
    panel.grid.major = element_line(colour = "grey"),
    panel.grid.minor = element_line(colour = "grey"),
    legend.title = element_text(size = 25),
    legend.text  = element_text(size = 25)
  )

Plot
# Optionally save:
# ggsave(filename = here::here("outputs", "strategy_proportions.png"), Plot,
#        width = 8, height = 6, dpi = 300)

################################################################################
# 5. 1000 day Simulation loop                                                            #
################################################################################
n_days = 1000

# Progress bar
pb <- progress_bar$new(total = n_days)

for (day in 1:n_days) {
  pb$tick()
  Sys.sleep(1 / 1000)  # purely cosmetic for the progress bar
  
  # Track changes made this iteration
  new_cohorts   <- list()  # cohorts produced this day
  tpop          <- 0       # total population this day
  food_consumed <- 0       # total food consumed this day
  t_count       <- 0       # number of transgenerational worms
  c_count       <- 0       # number of control worms
  
  # Ageing and population calculations
  for (i in seq_along(cohorts)) {
    cohort <- cohorts[[i]]              # current cohort
    cohort$age <- cohort$age + 1        # all worms age by 1 day
    
    # Update total population
    tpop <- tpop + cohort$size
    
    # Split by strategy
    if (cohort$strat == 1) {
      t_count <- t_count + cohort$size
    } else {
      c_count <- c_count + cohort$size
    }
    
    # Food consumption by age class
    if (cohort$age == 1) {
      # Eggs/L1/L2
      food_consumed <- food_consumed + cohort$size * Eat[1]
    } else if (cohort$age == 2) {
      # L3/L4
      food_consumed <- food_consumed + cohort$size * Eat[2]
    } else if (cohort$age >= 3 && cohort$age <= 7) {
      # Reproductive worms
      food_consumed <- food_consumed + cohort$size * Eat[3]
    } else {
      # Post-reproductive worms
      food_consumed <- food_consumed + cohort$size * Eat[4]
    }
    
    # Reproduction (age-specific)
    if (cohort$age > 2 && cohort$age <= 7) {  # reproductive age window
      
      # Default modifier (no transgenerational effect)
      repro_modifier <- 1
      
      # Apply transgenerational effects for strategy 1 only
      if (cohort$starvation == 1 && cohort$strat == 1) {
        repro_modifier <- Effects[1]  # F1 effect
      } else if (cohort$starvation == 2 && cohort$strat == 1) {
        repro_modifier <- Effects[2]  # F2 effect
      } else if (cohort$starvation >= 3 && cohort$strat == 1) {
        repro_modifier <- Effects[3]  # F3+ effect
      }
      
      # Eggs laid per worm at this age
      rate <- R[cohort$age] * repro_modifier
      
      # New cohort of age 0; offspring inherit parental starvation history + 1
      new_cohorts <- append(
        new_cohorts,
        list(
          list(
            age        = 0,
            size       = floor(rate * cohort$size),
            starvation = cohort$starvation + 1,
            strat      = cohort$strat
          )
        )
      )
    }
    
    # Update cohort in the global list
    cohorts[[i]] <- cohort
  }
  
  # Update food supply
  Food <- Food - food_consumed
  
  # Add newly produced cohorts to population
  cohorts <- c(cohorts, new_cohorts)
  
  # Store total and strategy-specific population sizes
  Pop[day]      <- tpop
  transPop[day] <- t_count
  conPop[day]   <- c_count
  
  # Store population data in a data frame (updated each iteration)
  dataframepop <- data.frame(
    Populationsize   = Pop,
    Transgenerational = transPop,
    control          = conPop
  )
  
  # Optionally store food for debugging:
  # Foods[day] <- Food 
  
  # Starvation event: food depleted
  if (Food <= 0) {
    
    # Draw a random founding population size between 10 and 100
    n_founders <- sample(10:100, 1)
    
    # Proportion of founders that are transgenerational vs control
    nt <- n_founders * dataframepop$Transgenerational[day] /
      dataframepop$Populationsize[day]
    nc <- n_founders * dataframepop$control[day] /
      dataframepop$Populationsize[day]
    
    # New founding cohorts (age 1 because founding occurs at larval arrest)
    cohorts <- list(
      list(age = 1, size = round(nt), starvation = 1, strat = 1),
      list(age = 1, size = round(nc), starvation = 0, strat = 0)
    )
    
    # Randomise the new food supply
    Food <- rnorm(mean = 1e+8, sd = 1e+7, n = 1)
  }
}

################################################################################
# 6.1000 day Graphs                                                                     #
################################################################################

# Calculate proportions of each strategy
dataframepop$tpro <- dataframepop$Transgenerational / dataframepop$Populationsize
dataframepop$cpro <- dataframepop$control          / dataframepop$Populationsize
dataframepop$day  <- 1:n_days

# Plot proportion of each strategy over time
Plot2 <- ggplot(dataframepop, aes(x = day)) +
  geom_point(aes(y = tpro, colour = "Transgenerational"), size = 3) +
  geom_line(aes(y = tpro, colour = "Transgenerational"), size = 1.5) +
  geom_point(aes(y = cpro, colour = "Control"), size = 3) +
  geom_line(aes(y = cpro, colour = "Control"), size = 1.5) +
  scale_color_manual(
    name   = "Strategy",
    values = c("Transgenerational" = "#0072B2",
               "Control"           = "#D55E00")
  ) +
  xlab("Time (days)") +
  ylab("Proportion of population") +
  ylim(c(0, 1)) +
  theme_classic() +
  theme(
    axis.title   = element_text(size = 30),
    axis.text    = element_text(size = 25),
    panel.grid.major = element_line(colour = "grey"),
    panel.grid.minor = element_line(colour = "grey"),
    legend.title = element_text(size = 25),
    legend.text  = element_text(size = 25)
  )

Plot2
# Optionally save:
# ggsave(filename = here::here("outputs", "strategy_proportions.png"), Plot,
#        width = 8, height = 6, dpi = 300)



################################################################################
#Heatmap of parameters                                                         #
################################################################################

library(parallel)

run_worm_sim <- function(f1_effect, f3_effect) {
  # Re-initialize local variables
  Effects <- c(f1_effect, 1.0, f3_effect)
  Food <- 1.0e+8
  n_days <- 1000 # Shorter duration for the grid search
  
  cohorts <- list(
    list(age = 0, size = 1, starvation = 1, strat = 1),
    list(age = 0, size = 1, starvation = 0, strat = 0)
  )
  
  for (day in 1:n_days) {

    # Track changes made this iteration
    new_cohorts   <- list()  # cohorts produced this day
    tpop          <- 0       # total population this day
    food_consumed <- 0       # total food consumed this day
    t_count       <- 0       # number of transgenerational worms
    c_count       <- 0       # number of control worms
    
    # Ageing and population calculations
    for (i in seq_along(cohorts)) {
      cohort <- cohorts[[i]]              # current cohort
      cohort$age <- cohort$age + 1        # all worms age by 1 day
      
      # Update total population
      tpop <- tpop + cohort$size
      
      # Split by strategy
      if (cohort$strat == 1) {
        t_count <- t_count + cohort$size
      } else {
        c_count <- c_count + cohort$size
      }
      
      # Food consumption by age class
      if (cohort$age == 1) {
        # Eggs/L1/L2
        food_consumed <- food_consumed + cohort$size * Eat[1]
      } else if (cohort$age == 2) {
        # L3/L4
        food_consumed <- food_consumed + cohort$size * Eat[2]
      } else if (cohort$age >= 3 && cohort$age <= 7) {
        # Reproductive worms
        food_consumed <- food_consumed + cohort$size * Eat[3]
      } else {
        # Post-reproductive worms
        food_consumed <- food_consumed + cohort$size * Eat[4]
      }
      
      # Reproduction (age-specific)
      if (cohort$age > 2 && cohort$age <= 7) {  # reproductive age window
        
        # Default modifier (no transgenerational effect)
        repro_modifier <- 1
        
        # Apply transgenerational effects for strategy 1 only
        if (cohort$starvation == 1 && cohort$strat == 1) {
          repro_modifier <- Effects[1]  # F1 effect
        } else if (cohort$starvation == 2 && cohort$strat == 1) {
          repro_modifier <- Effects[2]  # F2 effect
        } else if (cohort$starvation >= 3 && cohort$strat == 1) {
          repro_modifier <- Effects[3]  # F3+ effect
        }
        
        # Eggs laid per worm at this age
        rate <- R[cohort$age] * repro_modifier
        
        # New cohort of age 0; offspring inherit parental starvation history + 1
        new_cohorts <- append(
          new_cohorts,
          list(
            list(
              age        = 0,
              size       = floor(rate * cohort$size),
              starvation = cohort$starvation + 1,
              strat      = cohort$strat
            )
          )
        )
      }
      
      # Update cohort in the global list
      cohorts[[i]] <- cohort
    }
    
    # Update food supply
    Food <- Food - food_consumed
    
    # Add newly produced cohorts to population
    cohorts <- c(cohorts, new_cohorts)
    
    # Store total and strategy-specific population sizes
    Pop[day]      <- tpop
    transPop[day] <- t_count
    conPop[day]   <- c_count
    
    # Store population data in a data frame (updated each iteration)
    dataframepop <- data.frame(
      Populationsize   = Pop,
      Transgenerational = transPop,
      control          = conPop
    )
    
    # Optionally store food for debugging:
    # Foods[day] <- Food 
    
    # Starvation event: food depleted
    if (Food <= 0) {
      
      # Draw a random founding population size between 10 and 100
      n_founders <- sample(10:100, 1)
      
      # Proportion of founders that are transgenerational vs control
      nt <- n_founders * dataframepop$Transgenerational[day] /
        dataframepop$Populationsize[day]
      nc <- n_founders * dataframepop$control[day] /
        dataframepop$Populationsize[day]
      
      # New founding cohorts (age 1 because founding occurs at larval arrest)
      cohorts <- list(
        list(age = 1, size = round(nt), starvation = 1, strat = 1),
        list(age = 1, size = round(nc), starvation = 0, strat = 0)
      )
      
      # Randomise the new food supply
      Food <- rnorm(mean = 1e+8, sd = 1e+7, n = 1)
    }
  }
  
  # Return final proportion of transgenerational strategy
  final_t_pro <- t_count / (t_count + c_count)
  return(final_t_pro)
}



# Define ranges to test
f1_range <- seq(1, 3, by = 0.02)
f3_range <- seq(0.1, 1, by = 0.01)

# Create the grid
param_grid <- expand.grid(f1 = f1_range, f3 = f3_range)

# Run the simulation for every combination
param_grid$outcome <- mapply(run_worm_sim, param_grid$f1, param_grid$f3)

# Define original parameters
orig_f1 <- 2.0
orig_f3 <- 0.5
