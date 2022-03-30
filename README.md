# Identification-of-ideotypes

Supporting code associated with *"Bayesian optimisation over trait-relationship space for the quantitative identification of crop ideotypes in Brassica"*, with example input/output data. The purpose of the repo is to document the data processing & modelling carried out.

Each analysis stage is carried out by a single **R script**, with *supporting scripts*, *input files*, and *output files* as specified below.


### data preprocessing
- **impute_data.R**: impute values for missing datapoints.

	- *input*:
		- `/data/final_data_normed_cleaned.rds` - data table of normalised data for all accessions. Transformations which have been already applied to normalise the data are shown in supplemental table 2.

	- *output*:
		- `/data/final_data_normed_cleaned_imputed.rds` - missing values in input have been imputed using the **`mice_3.11.0`** R package.


### identification of trait relationship structure
- **id_trait_structure.R**: Identification of trait-trait relationship structure directed acyclic graph (DAG). Uses **`bnlearn_4.5`** R package to search DAG structure space to find linear regression model structures (Gaussian Bayesian Networks) which best explain the data according to the BIC criterion. A separate model is estimated for each of 5-fold data splits, and an average consensus model computed from these.

	- *supporting script files*:
	  - `id_trait_structure_functions.R`

	- *input*:
		- `/data/final_data_normed_cleaned_imputed.rds` - data table of phenotype data used for model fitting.
		- `/data/final_traits_blacklist.txt` - specifies a hierarchy of traits such that traits on lower lines cannot be parents of traits on higher lines.
		- `/data/final_traits_whitelist.txt` - trait-trait links which must be included in any considered model.

	- *output*:
		`trait_structure_output/data/`
		- `bnlearn_data.rds` : filtered phenotype data for spring ecotype varieties only. Which was used to learn spring OSR model structures.
		- `bnlearn_results.rds` : list of results relating model structures inferred for spring OSR. Includes identified model structure for each of the k-fold models identified by bnlearn.
			- `results$str.df.list[[i]]` - dataframe of edges and edge consistency found during bootstrap in ith fold. 6th element is average over all 5 folds, by weighted sum of edge weights in each, where weight is proportional to BIC score.
			- `results$avg.dag.list[[i]]` - averaged network from bootstrapping in ith fold. 6th element is average over all 5 folds, by weighted sum of edge weights in each, where weight is proportional to BIC score.
			- `results$scores[[i]]` - BIC scores calculated for each of the averaged networks
			- `results$model.strings[[i]]` - averaged DAG model structure in string format, for use in building Stan model structure.
			- `results$xval.df` - dataframe of predicted phenotype values for witheld data using ith averaged network model.

		`trait_structure_output/graphs/`
		- `trait.value_vs_predictions.pdf` : plots of predicted values vs observed values for withheld data for each of the k-fold models (plot of `results$xval.df`).
		- `k-fold_DAGs.pdf` : DAGs showing each of the k-fold `results$avg.dag.list` models identified.
		- `avg_DAG_strength.pdf` : consensus DAG identified from `results$avg.dag.list` based on connections seen frequently & strongly in the k-fold identified structures. (this DAG corresponds to the 6th element of `results$avg.dag.list`).


### Predicting the yield consequence of modifying traits
- **trait_modification_prediction.R**
	- *supporting script files*:
		- `trait_modification_functions.R` - R functions used
		- `GP_ARD_functions.stan` - Stan functions defining kernel used for GP prior.

	- *input*:
		- `id_trait_structure_output/data/bnlearn_data.rds`
		- `id_trait_structure_output/data/bnlearn_results.rds`
		- (if running with **TESTING=TRUE**):
		`data/example_trait_modification_data/ki/height_fitted.rds` - precomputed predictions for trait values when height is modified).

	- *output*:
		`trait_modification_output/model_files/`
		- `{TRAIT}_NULL.stan` - Stan source code for model to sample   from normal distribution with observed mean and sd.
		- `{TRAIT}_NULL.rds` - compiled rStan model of the above
		- `{TRAIT}_NULL_fitted.rds` - sampled values from the above model. NB this file >100Mb, and so is not included in the github, but can be generated from the provided scripts and input data.
		- `k{i}-fold/{TRAIT}/{TRAIT}_opt.stan` - stan source code for model to use for max likelihood hyperparameter estimation, using perturbed trait TRAIT, and k-fold model structure i.
		- `k{i}-fold/{TRAIT}/{TRAIT}_opt.rds` - rstan model for the above
		- `/k{i}-fold/{TRAIT}/{TRAIT}_pred.stan` - stan source code for model to make predictions for values of other traits, when trait TRAIT is modified.
		- `/k{i}-fold/{TRAIT}/{TRAIT}_pred.rds` - rstan model for the above.
		- `/k{i}-fold/{TRAIT}/{TRAIT}_fitted.rds` - sampled values from the posterior defined in the above model. If running with TESTING=TRUE, this is just copied from the input `/data/example_trait_modification_data/k{i}/height_fitted.rds`.

		`trait_modification_output/graphs/`
		- `{TRAIT}_k-avg_predicted_perturbation_effect.pdf` - plot of CIs of trait values when perturbed trait is fixed to x-axis values.
		- `{TRAIT}_k-avg_predicted_perturbation_effect_back.pdf` - as above, after data is back transformed to undo transformations to make more normally distributed.
		- `k{i}-fold/{TRAIT}/`
			- `optimal_hyperparameters.csv` - values of alpha, sigma, rho hyper parameters estimated by regularised maximum likelihood.
			- `fit_summary.csv` - all values sampled from posterior (same data as `TRAIT_fitted.rds`, but in .csv format).
		- `k-averaged/{TRAIT}/`
			- `avg_pred_vs_measured.pdf` - plot of observed values vs mean predictions from sampled posterior.
			- `avg_pred_vs_meaured.data.rds` - data used in making above plot
			- `avg_pred_vs_measured_back.pdf` - plot of observed values vs mean predictions from sampled posterior, after normalisation transformation back transformed.
			- `avg_pred_vs_meaured.data_back.rds` - the data used in making the above plot.


### Bayesian optimisation for ideotype identification
- **bayesian_optimisation.R**
	- *supporting script files*:
		- `BO_functions.R` - functions used in bayesian_optimisation.R
		- `BO_plot_results.R` - script to plot q-EI points, and observations
	- *input*:
		`data/BO_data/`
		- `pred_model.stan` - stan source code used to sample posterior - prediction for next proposed q-point, conditioned on observations and hyperparamteres.
		- `pred_model.rds` - rstan of the above model.
		- `spring_trait_data.rds` - data for spring accessions only (as produced by id_trait_structure.R, renamed bnlearn_data.rds).
		- `winter_trait_data.rds` - same as above, but winter OSR accessions.
	- *output*:
		`BO_output/`
		- `{ECOTYPE}_hyperparam_model/`
			- `opt_model.stan` - stan source for model used to infer max likelihood values for GP prior hyperparams
			-	`opt_model.rds` - rstan compiled model of the above
			-	`opt_hyperparams.csv` - optimal hyperparameter values
		- `{ECOTYPE}_PCA_{SUFFIX}.csv` - summary table of variance explained by principle components.
		- `{ECOTYPE}_stan.data_{SUFFIX}.RDS` - the data used to the q-values for exploration
		- `{ECOTYPE}_q.values_{SUFFIX}.RDS` - the next q identified values for exploration.

### Simulation of causal SNP identification
- **simSNP_inference.R**
	- *supporting script files*:
	  - `simSNP_inference_functions.R`
