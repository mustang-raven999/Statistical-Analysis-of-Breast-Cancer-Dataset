
# Statistical Analysis of Breast Cancer Dataset

This repository contains a statistical analysis of the BreastCancer dataset, part of the `mlbench` package, used to classify breast cancer tissue samples as benign or malignant based on nine cytological characteristics.

## Project Overview
Using various statistical and machine learning techniques, we aim to identify the most significant features that differentiate benign from malignant samples and create classifiers for accurate predictions.

## Dataset
- **Source**: Fine Needle Aspiration Cytology (FNAC) samples of breast tissues from 699 women.
- **Variables**: Nine cytological characteristics, such as cell size, shape, and mitosis, scored from 1 to 10.

## Methods
We employed several classifiers:
- **Logistic Regression** (with AIC and BIC for feature selection)
- **Lasso Regression** for regularization
- **Linear Discriminant Analysis (LDA)** for dimensionality reduction

## Results
- The best predictive model is **Logistic Regression** with BIC feature selection (5-variable model), offering low complexity and low error rate.
- **Cl.thickness**, **Bare.nuclei**, and **Bl.cromatin** are key predictors.

## Files
- **BreastCancer.csv**: The dataset used.
- **BreastCancerScript.R**: R script for the analysis.
- **notebookbc.Rmd**: R Markdown file documenting the analysis.
- **report.pdf**: Detailed report on the analysis and results.

## Libraries Used
- `mlbench`, `dplyr`, `ggplot2`, `GGally`, `bestglm`, `glmnet`, `MASS`, `caret`

## How to Run
1. Clone the repository.
2. Open and run the R scripts for analysis.
3. Review the findings in the **report.pdf** file.
