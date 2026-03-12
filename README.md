# Mortality Forecasting: SVD vs. Poisson MLE Analysis

This project was developed for the **Computational Statistics** course at Corvinus University of Budapest. It explores and compares different estimation methods for the **Lee-Carter mortality model**, a standard framework for forecasting death rates.

## Team members
*   Bendegúz Birkmayer
*   Bojta Rácz
*   Kristóf Légrádi
*   Roland Tuboly

## Project Overview
The goal is to evaluate the out-of-sample performance of the Lee-Carter model using two primary estimation approaches across a dataset of 30 countries/regions:
1.  **Singular Value Decomposition (SVD):** The classic approach that assumes normally distributed log-mortality rates.
2.  **Poisson Maximum Likelihood Estimation (MLE):** A more robust approach that models death counts as Poisson-distributed. We compare three optimization algorithms for the MLE:
    *   **BFGS** (Broyden–Fletcher–Goldfarb–Shanno)
    *   **Nelder-Mead**
    *   **SANN** (Simulated Annealing)

### Key Features
*   **Data Preprocessing:** Aggregation of HMD-style data, age-grouping (90+), and linear imputation for missing values.
*   **Train-Test Split:** Models are trained on data from **1960–2004** and tested on **2005–2019**.
*   **Evaluation:** Comparative analysis using **RMSE** (Root Mean Square Error) on log-death rates.
*   **Visualizations:** performance scatter plots and geographical maps showing which method "wins" by country.

## Repository Structure
```text
.
├── data/                            # Mortality data for 30 countries (Deaths, Exposures, Rates)
├── mortality_analysis.R             # Main R analysis script
├── Mortality_Forecasting_Report.pdf # Final report
└── README.md                        # Project documentation
```

## How to Run
1.  Ensure you have **R** and **RStudio** installed.
2.  Install the required packages:
    ```r
    install.packages(c("dplyr", "tidyr", "MortCast", "forecast", "zoo", 
                       "ggplot2", "ggrepel", "sf", "rnaturalearth", 
                       "rnaturalearthdata"))
    ```
3.  Open `mortality_analysis.R` and run the script. It uses relative paths, so as long as your working directory is the root of this folder, it will work automatically.

## Results Summary
The analysis generally shows that while SVD is computationally simpler, **Poisson MLE (particularly via BFGS)** often provides superior forecasting accuracy for specific countries or sex-specific data where the Poisson assumption better captures the variance of death counts.
