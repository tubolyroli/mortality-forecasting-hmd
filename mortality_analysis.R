# Computational Statistics Group Homework
# Topic05 - Mortality Forecasting

# Members of the team:
# Bendegúz Birkmayer
# Bojta Rácz
# Kristóf Légrádi
# Roland Tuboly

##### 0. PREREQUISITES & WD #####

# Loading packages
library(dplyr)
library(tidyr)
library(MortCast)
library(forecast)
library(zoo)
library(ggplot2)
library(ggrepel)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# Set working directory
# setwd()

##### 1. DATA LOADING & PREPROCESSING #####

# Get a list of all country sub-folders (e.g., data/australia, data/sweden)
country_folders <- list.dirs(path = "data", full.names = TRUE, recursive = FALSE)

# Initialize lists to store dataframes for each country
list_deaths <- list()
list_exposures <- list()
list_rates <- list()

# Loop through each country folder
for (folder in country_folders) {
  # Extract the clean country name from the folder path
  raw_name <- basename(folder)
  # Handle country naming conventions
  if (raw_name == "uk pop") {
    country_name <- "UK (civil population only)"
  } else if (raw_name == "newzealand") {
    country_name <- "New Zealand"
  } else if (raw_name == "west germany") {
    country_name <- "West Germany"
  } else {
    # Capitalize the first letter for all other countries
    country_name <- paste0(toupper(substr(raw_name, 1, 1)), substr(raw_name, 2, nchar(raw_name)))
  }
  
  # Locate the specific files within the folder using pattern matching
  path_d <- list.files(folder, pattern = "Deaths_1x1.txt", full.names = TRUE)
  path_e <- list.files(folder, pattern = "Exposures_1x1.txt", full.names = TRUE)
  path_m <- list.files(folder, pattern = "Mx_1x1.txt", full.names = TRUE)
  
  # Helper function to read txt files (skips first 2 description lines)
  read_txt <- function(file_path, country) {
    if (length(file_path) == 0) return(NULL)
    df <- read.table(file_path, header = TRUE, na.strings = ".", skip = 2, stringsAsFactors = FALSE)
    df$Country <- country_name
    return(df)
  }
  
  # Read the txt files and store in lists
  list_deaths[[country_name]] <- read_txt(path_d, country_name)
  list_exposures[[country_name]] <- read_txt(path_e, country_name)
  list_rates[[country_name]] <- read_txt(path_m, country_name)
}

# Combine all country lists into single dataframes
all_deaths <- do.call(rbind, list_deaths)
all_exposures <- do.call(rbind, list_exposures)
all_rates <- do.call(rbind, list_rates)

# Clean up environment
rm(list_deaths, list_exposures, list_rates, country_folders, folder, path_d, path_e, path_m)

##### Creating merged database #####

# Create temporary copies to rename columns without affecting originals
temp_deaths <- all_deaths
temp_exposures <- all_exposures
temp_rates <- all_rates
# Define the value columns that need suffixes
val_cols <- c("Female", "Male", "Total")
# Rename columns in the temporary dataframes
names(temp_deaths)[names(temp_deaths) %in% val_cols] <- paste0(val_cols, "_deaths")
names(temp_exposures)[names(temp_exposures) %in% val_cols] <- paste0(val_cols, "_exposures")
names(temp_rates)[names(temp_rates) %in% val_cols] <- paste0(val_cols, "_rates")

# Merge the dataframes on Country, Year, and Age
# all = TRUE ensures that if a row is missing in one file, we don't lose the data
step1 <- merge(temp_deaths, temp_exposures, by = c("Country", "Year", "Age"), all = TRUE)
full_data <- merge(step1, temp_rates, by = c("Country", "Year", "Age"), all = TRUE)
# Clean up temporary objects
rm(temp_deaths, temp_exposures, temp_rates, step1)

##### Grouping 90+ #####
# Because there are a lot of small and zero numbers over 90, and since
# tiny populations in some ages create wild random noise that overpowers
# real mortality trends, ruining the forecast accuracy for the entire population,
# we create a new dataset 'grouped_data' where all ages >= 90 are merged

grouped_data <- full_data %>%
  # Filter Years, we'll need 1960 onwards later anyway
  filter(Year >= 1960) %>%
  mutate(Age_Numeric = as.numeric(gsub("\\+", "", Age))) %>%
  # Define the Group: 90 and above get ID 90
  mutate(Age_Temp = ifelse(Age_Numeric >= 90, 90, Age_Numeric)) %>%
  # Aggregate by Country, Year, and the numeric ID
  group_by(Country, Year, Age_Temp) %>%
  summarise(
    # Summing Deaths and Exposures
    Female_deaths    = sum(Female_deaths, na.rm = TRUE),
    Male_deaths      = sum(Male_deaths, na.rm = TRUE),
    Total_deaths     = sum(Total_deaths, na.rm = TRUE),
    Female_exposures = sum(Female_exposures, na.rm = TRUE),
    Male_exposures   = sum(Male_exposures, na.rm = TRUE),
    Total_exposures  = sum(Total_exposures, na.rm = TRUE),
    # Recalculating Rates - weighted average preserves the Deaths=Exposure*Rate relationship
    Female_rates = ifelse(sum(Female_exposures) > 0, sum(Female_deaths)/sum(Female_exposures), 0),
    Male_rates   = ifelse(sum(Male_exposures) > 0, sum(Male_deaths)/sum(Male_exposures), 0),
    Total_rates  = ifelse(sum(Total_exposures) > 0, sum(Total_deaths)/sum(Total_exposures), 0),
    .groups = "drop" # Ungroup after summarizing
  ) %>%
  # SORTING IS CRITICAL HERE (before converting back to categorical)
  arrange(Country, Year, Age_Temp) %>%
  # Apply the 90+ label
  mutate(Age = ifelse(Age_Temp == 90, "90+", as.character(Age_Temp))) %>%
  # Clean up: remove the temporary numeric column and put Age first
  select(Country, Year, Age, everything(), -Age_Temp)

rm(all_deaths, all_exposures, all_rates, full_data)

##### Train–test split #####

# Exclude COVID years (2020+)
grouped_data <- grouped_data %>% filter(Year <= 2019)
# Training set: 1960–2004
train_data <- grouped_data %>% filter(Year >= 1960, Year <= 2004)
# Test set: 2005–2019
test_data <- grouped_data %>% filter(Year >= 2005, Year <= 2019)

##### 2. HELPER FUNCTIONS #####

# Impute series helper (linear interp + forward/backward fill)
impute_series <- function(v) {
  # if everything is NA, leave as is for now (we'll handle it later)
  if (all(is.na(v))) return(v)
  # linear interpolation for internal NAs
  v <- na.approx(v, na.rm = FALSE)
  # carry first/last observed values outward
  v <- na.locf(v, na.rm = FALSE) # forward
  v <- na.locf(v, fromLast = TRUE, na.rm = FALSE) # backward
  v
}

# Helper: build age × year matrix and impute NAs
make_mx_imputed <- function(data, country, rate_col) {
  sub <- data %>%
    filter(Country == country) %>%
    mutate(Age_Start = as.numeric(gsub("\\+", "", Age))) %>%
    arrange(Age_Start, Year) %>%
    select(Age_Start, Year, rate = !!sym(rate_col))
  if (nrow(sub) == 0) return(NULL) # no data at all for this country/sex
  
  # Wide: rows = ages, cols = years
  wide <- sub %>% tidyr::pivot_wider(names_from = Year, values_from = rate) %>% arrange(Age_Start)
  ages <- wide$Age_Start
  mx <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mx) <- ages
  
  # Impute within each age over years
  mx_imp <- t(apply(mx, 1, impute_series))
  
  # Ages with all NA even after that (no data at all) -> impute from other ages
  all_na_rows <- apply(mx_imp, 1, function(z) all(is.na(z)))
  if (any(all_na_rows)) {
    col_means <- apply(mx_imp[!all_na_rows, , drop = FALSE], 2, function(z) {
      if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
    })
    for (i in which(all_na_rows)) mx_imp[i, ] <- col_means
  }
  # If there are still NAs (e.g. entire column had no data), fill with global mean
  if (any(is.na(mx_imp))) {
    overall_mean <- mean(mx_imp, na.rm = TRUE)
    if (is.na(overall_mean)) overall_mean <- 1e-8
    mx_imp[is.na(mx_imp)] <- overall_mean
  }
  years_train <- sort(unique(sub$Year))
  list(mx = mx_imp, years = years_train)
}

##### Poisson Negative Log-Likelihood Function #####
# Used for ALL methods (BFGS, Nelder-Mead, SANN)

lc_poisson_nll_safe <- function(params, D, E) {
  nA <- nrow(D); nY <- ncol(D)
  ax <- params[1:nA]
  bx <- params[(nA + 1):(2 * nA)]
  kt <- params[(2 * nA + 1):(2 * nA + nY)]
  
  log_m <- outer(bx, kt)
  log_m <- sweep(log_m, 1, ax, "+")
  
  # Guard: overflow protection (this is crucial for SANN)
  if (any(!is.finite(log_m))) return(1e50)
  log_m <- pmin(log_m, 700)
  log_m <- pmax(log_m, -700)
  
  lambda <- E * exp(log_m)
  
  # Guard: invalid lambda
  if (any(!is.finite(lambda))) return(1e50)
  if (any(lambda < 0, na.rm = TRUE)) return(1e50)
  
  # Log-likelihood with small eps for stability
  # We ignore the log(D!) term because it is constant with respect to parameters
  log_L <- sum(D * log(lambda + 1e-10) - lambda, na.rm = TRUE)
  
  if (!is.finite(log_L)) return(1e50)
  
  return(-log_L)
}

##### 3. ESTIMATION LOOP (SVD + 3 MLE METHODS) #####

countries <- sort(unique(grouped_data$Country))
sexes <- c("Female", "Male")
forecast_list <- list()

# Define the MLE methods we want to run
mle_methods_to_run <- c("BFGS", "Nelder-Mead", "SANN")

for (sex in sexes) {
  rate_col <- paste0(sex, "_rates")
  d_col <- paste0(sex, "_deaths")
  e_col <- paste0(sex, "_exposures")
  
  for (cty in countries) {
    # message("Processing: ", cty, " - ", sex)
    
    # get imputed matrix (crucial for SVD)
    mm <- make_mx_imputed(train_data, cty, rate_col)
    if (is.null(mm)) next
    mx_train <- mm$mx  # This has no NAs
    years_train <- mm$years
    
    # get raw matrices (D and E) for MLE
    sub_data <- train_data %>% 
      filter(Country == cty) %>% 
      mutate(Age_Start = as.numeric(gsub("\\+", "", Age))) %>% 
      arrange(Age_Start, Year)
    
    wide_d <- sub_data %>% select(Age_Start, Year, val = !!sym(d_col)) %>% 
      pivot_wider(names_from = Year, values_from = val) %>% arrange(Age_Start)
    D_mat <- as.matrix(wide_d[, -1])
    
    wide_e <- sub_data %>% select(Age_Start, Year, val = !!sym(e_col)) %>% 
      pivot_wider(names_from = Year, values_from = val) %>% arrange(Age_Start)
    E_mat <- as.matrix(wide_e[, -1])
    
    if (any(dim(D_mat) != dim(E_mat)) || nrow(D_mat) < 3) next
    
    # 1. SVD Estimation
    # Use the IMPUTED mx_train, but apply the 1e-6 floor logic 
    # to match our original previously successful SVD runs
    mx_clean <- mx_train
    mx_clean[mx_clean <= 0 | is.na(mx_clean)] <- 1e-6
    
    # Run SVD
    lc_svd <- try(MortCast::leecarter.estimate(mx_clean), silent = TRUE)
    if (inherits(lc_svd, "try-error")) next
    
    # Helper to forecast and store results
    forecast_and_store <- function(ax, bx, kt, label) {
      min_year <- min(years_train)
      max_year <- max(years_train)
      # Test years for this country AFTER last training year
      test_years <- sort(unique(test_data$Year[test_data$Country == cty]))
      test_years <- test_years[test_years > max_year]
      if (length(test_years) == 0) return(NULL)
      
      # Forecast k_t as a time series
      kt_ts <- ts(kt, start = min_year, frequency = 1)
      fit_kt <- forecast::auto.arima(kt_ts)
      fc_kt  <- forecast::forecast(fit_kt, h = length(test_years))
      kt_future <- as.numeric(fc_kt$mean)
      
      # Reconstruct rates: log m_x,t = a_x + b_x * k_t
      log_m_hat <- outer(bx, kt_future)
      log_m_hat <- sweep(log_m_hat, 1, ax, "+")
      dimnames(log_m_hat) <- list(Age_Start = as.numeric(names(ax)), Year = test_years)
      
      df_hat <- as.data.frame(as.table(log_m_hat))
      names(df_hat) <- c("Age_Start", "Year", "log_rate_hat")
      df_hat$Age_Start <- as.numeric(as.character(df_hat$Age_Start))
      df_hat$Year      <- as.numeric(as.character(df_hat$Year))
      df_hat$Country <- cty
      df_hat$Sex     <- sex
      df_hat$Method  <- label
      return(df_hat)
    }
    
    # Store SVD result
    forecast_list[[paste(cty, sex, "SVD", sep = "_")]] <- forecast_and_store(lc_svd$ax, lc_svd$bx, lc_svd$kt, "SVD")
    
    # 2. MLE Estimation (loop through all methods)
    
    # Get initial values from SVD (we think this is a good starting point)
    par_init <- c(lc_svd$ax, lc_svd$bx, lc_svd$kt)
    nA <- nrow(D_mat); nY <- ncol(D_mat)
    
    for (method_name in mle_methods_to_run) {
      
      # Setup control params based on the method
      if (method_name == "SANN") {
        # SANN specific: temp controls exploration, tmax controls cooling
        ctrl_use <- list(maxit = 1200, temp = 2, tmax = 10)
      } else if (method_name == "Nelder-Mead") {
        ctrl_use <- list(maxit = 3000, reltol = 1e-8)
      } else {
        # BFGS default
        ctrl_use <- list(maxit = 1000)
      }
      
      # We use a safe function for all methods to avoid BFGS crashes on NaN/Inf
      opt_res <- try(optim(par = par_init, fn = lc_poisson_nll_safe, D = D_mat, E = E_mat, 
                           method = method_name, control = ctrl_use), silent = TRUE)
      
      if (!inherits(opt_res, "try-error") && is.finite(opt_res$value)) {
        # Extract parameters
        est <- opt_res$par
        ax_m <- est[1:nA]
        bx_m <- est[(nA + 1):(2 * nA)]
        kt_m <- est[(2 * nA + 1):(2 * nA + nY)]
        
        # Normalize: sum(bx) should be 1, sum(kt) should be 0
        s_bx <- sum(bx_m)
        bx_norm <- bx_m / s_bx
        kt_temp <- kt_m * s_bx
        m_kt <- mean(kt_temp)
        kt_norm <- kt_temp - m_kt
        ax_norm <- ax_m + bx_norm * m_kt
        
        label <- paste0("MLE_", method_name)
        forecast_list[[paste(cty, sex, label, sep = "_")]] <- forecast_and_store(ax_norm, bx_norm, kt_norm, label)
      }
    }
  }
}

all_results <- bind_rows(forecast_list)

##### 4. EVALUATION #####

# Observed test data in long format (both sexes)
obs_df <- test_data %>%
  mutate(Age_Start = as.numeric(gsub("\\+", "", Age))) %>%
  select(Country, Year, Age_Start, Female = Female_rates, Male = Male_rates) %>%
  pivot_longer(cols = c(Female, Male), names_to = "Sex", values_to = "rate_obs")

# Join forecasts with observations
final_results <- all_results %>%
  inner_join(obs_df, by = c("Country", "Year", "Age_Start", "Sex")) %>%
  mutate(
    log_rate_obs = log(pmax(rate_obs, 1e-10)),
    log_rate_hat = log_rate_hat
  )

# RMSE on log death rates
rmse_table <- final_results %>%
  group_by(Country, Sex, Method) %>%
  summarise(RMSE = sqrt(mean((log_rate_obs - log_rate_hat)^2, na.rm = TRUE)), .groups = "drop")

rmse_wide <- rmse_table %>%
  pivot_wider(names_from = Method, values_from = RMSE)

# Find the best method for each country-sex pair
rmse_wide$Best_Method <- colnames(rmse_wide)[3:ncol(rmse_wide)][apply(rmse_wide[, 3:ncol(rmse_wide)], 1, which.min)]
print(table(rmse_wide$Best_Method))

# Print summary
print(head(rmse_wide))

##### 5. PLOTS #####

# Prepare data for plot: Select SVD and MLE_BFGS for visualization
comparison <- rmse_wide %>%
  select(Country, Sex, SVD, MLE = MLE_BFGS) # Renaming MLE_BFGS to MLE for the plot

##### SVD vs MLE (BFGS) RMSE scatter (global picture) #####
ggplot(comparison, aes(x = SVD, y = MLE, color = Sex, label = Country)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 3, alpha = 0.8) +
  # Labeling some extreme cases where difference is large
  ggrepel::geom_text_repel(
    data = subset(comparison, abs(SVD - MLE) > 0.1),
    size = 3, max.overlaps = 20, show.legend = FALSE
  ) +
  coord_equal() +
  labs(
    x = "RMSE (log mortality) – SVD",
    y = "RMSE (log mortality) – MLE (BFGS)", 
    color = "Sex",
    title = "Out-of-sample performance: SVD vs Poisson MLE (BFGS) Lee–Carter",
    subtitle = "Each point = country–sex; points above the line favour SVD"
  ) +
  theme_minimal(base_size = 12)

##### Map visualization #####

# Define regions (needed to create comparison_region later)
mediterranean <- c("Italy", "Portugal", "Spain")
post_socialist <- c("Bulgaria", "Czechia", "Estonia", "Hungary", 
                    "Latvia", "Lithuania", "Poland", "Slovakia")

# Prepare the data: Add Region to the existing 'comparison' object
comparison_region <- comparison %>%
  mutate(
    Region = case_when(
      Country %in% mediterranean ~ "Mediterranean",
      Country %in% post_socialist ~ "Post-Socialist",
      TRUE ~ "Other"
    )
  )

# Aggregate to country-level winner (ignoring sex)
country_win <- comparison_region %>%
  group_by(Country, Region) %>%
  summarise(
    SVD_better = mean(MLE > SVD, na.rm = TRUE),  # share of cases where SVD better
    .groups = "drop"
  ) %>%
  mutate(
    Winner_overall = ifelse(SVD_better >= 0.5, "SVD", "MLE")
  )

# Fix some country name typos to match rnaturalearth
country_win <- country_win %>%
  mutate(
    Country_clean = dplyr::recode(
      Country,
      "Usa"        = "United States of America",
      "Netherland" = "Netherlands",
      "New zealand"= "New Zealand",
      "Latbia"     = "Latvia",
      "Sovakia"    = "Slovakia",
      "UK (civil population only)" = "United Kingdom",
      "Czechia" = "Czech Republic",
      "West Germany" = "Germany",
      .default = Country
    )
  )

# Get World Map Data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Join data
world_join <- world %>%
  left_join(country_win, by = c("name_long" = "Country_clean"))

# Plot the map
ggplot(world_join) +
  geom_sf(aes(fill = Winner_overall), color = "grey40", size = 0.1) +
  scale_fill_manual(
    values = c("SVD" = "#1b9e77", "MLE" = "#d95f02"),
    na.value = "grey90"
  ) +
  coord_sf(xlim = c(-20, 40), ylim = c(35, 70)) +  # zoom to Europe
  labs(
    fill = "Better method",
    title = "Which method wins overall by country?",
    subtitle = "Winner decided by majority of sex-specific RMSE comparisons"
  ) +
  theme_void(base_size = 12)


##### 6. COMPARISON TABLES #####

##### BFGS vs SVD #####

rmse_long <- rmse_table # not to distort our table

bfgs_vs_svd <- rmse_long %>%
  filter(Method %in% c("SVD", "MLE_BFGS")) %>% # filtering for the two method
  pivot_wider(names_from = Method, values_from = RMSE) %>%
  mutate(
    Winner = case_when(
      MLE_BFGS < SVD ~ "BFGS",
      SVD < MLE_BFGS ~ "SVD",
      TRUE ~ "Tie"
    )
  ) %>% # account even for the ties if they occur
  count(Winner, name = "Wins")

bfgs_vs_svd # let' see the comparison

##### NELLDER-MEAD vs SVD #####

nelder_vs_svd <- rmse_long %>%
  filter(Method %in% c("SVD", "MLE_Nelder-Mead")) %>%
  pivot_wider(names_from = Method, values_from = RMSE) %>%
  mutate(
    Winner = case_when(
      `MLE_Nelder-Mead` < SVD ~ "Nelder-Mead",
      SVD < `MLE_Nelder-Mead` ~ "SVD",
      TRUE ~ "Tie"
    )
  ) %>%
  count(Winner, name = "Wins")

nelder_vs_svd

#### SANN vs SVD

sann_vs_svd <- rmse_long %>%
  filter(Method %in% c("SVD", "MLE_SANN")) %>%
  pivot_wider(names_from = Method, values_from = RMSE) %>%
  mutate(
    Winner = case_when(
      MLE_SANN < SVD ~ "SANN",
      SVD < MLE_SANN ~ "SVD",
      TRUE ~ "Tie"
    )
  ) %>%
  count(Winner, name = "Wins")

sann_vs_svd

##### Sex #####

bfgs_svd_by_sex <- rmse_table %>%
  filter(Method %in% c("SVD", "MLE_BFGS")) %>% # filter for the methods
  pivot_wider(names_from = Method, values_from = RMSE) %>%
  mutate(
    Winner = case_when(
      MLE_BFGS < SVD ~ "BFGS", # comparison of the two model for the win
      SVD < MLE_BFGS ~ "SVD",
      TRUE ~ "Tie" # in case of a tie
    )
  ) %>%
  count(Sex, Winner, name = "Wins") %>%
  tidyr::complete(Sex, Winner = c("BFGS", "SVD", "Tie"), fill = list(Wins = 0)) %>%
  pivot_wider(names_from = Winner, values_from = Wins, values_fill = 0) %>%
  mutate(
    Total = BFGS + SVD + Tie,
    Share_BFGS = ifelse(Total > 0, BFGS / Total, NA_real_) # get a share of the winners for the ratios
  )

bfgs_svd_by_sex

##### Mediterranean vs Post-Socialist #####

bfgs_svd_by_region <- rmse_table %>%
  filter(Method %in% c("SVD", "MLE_BFGS")) %>%
  pivot_wider(names_from = Method, values_from = RMSE) %>%
  mutate(
    Region = case_when(
      Country %in% mediterranean ~ "Mediterranean",
      Country %in% post_socialist ~ "Post-Socialist",
      TRUE ~ "Other"
    ),
    Winner = case_when(
      MLE_BFGS < SVD ~ "BFGS",
      SVD < MLE_BFGS ~ "SVD",
      TRUE ~ "Tie"
    )
  ) %>%
  count(Region, Winner, name = "Wins") %>%
  tidyr::complete(Region, Winner = c("BFGS", "SVD", "Tie"), fill = list(Wins = 0)) %>%
  pivot_wider(names_from = Winner, values_from = Wins, values_fill = 0) %>%
  mutate(
    Total = BFGS + SVD + Tie,
    Share_BFGS = ifelse(Total > 0, BFGS / Total, NA_real_)
  ) %>%
  arrange(match(Region, c("Mediterranean", "Post-Socialist", "Other")))

bfgs_svd_by_region

