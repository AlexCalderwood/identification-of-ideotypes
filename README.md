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
- **trait_modification_prediction.R** :
    - Takes the model structures inferred by **id_trait_structure.R**, uses Gaussian Processes (GP) to model the potentially non-linear relationships between the connected nodes, and predict the effects of perturbing each yield trait.
    - MCMC sampling is carried out using Stan. This script parses model structures encoded as strings in `bnlearn$model.strings` to produce `.stan` files describing the GP modelled relationships between each node and its parent traits. It also calls Stan to compile and run no-U-turn MCMC sampling on these models.
    <br>
	- *supporting script files*:
		- `trait_modification_functions.R` : helper R functions
		- `GP_ARD_functions.stan` : Stan functions used to define the kernel used for the GP prior.

	- *input*:
		- `id_trait_structure_output/data/bnlearn_data.rds` : generated output of  **id_trait_structure.R** (described above).
		- `id_trait_structure_output/data/bnlearn_results.rds` : generated output of  **id_trait_structure.R** (described above).
		- (if running with `TESTING=TRUE`:
		`data/example_trait_modification_data/k{i}/height_fitted.rds` : These precomputed predictions (for other trait values when "height" is modified), are normally generated as output of this script. They are loaded for convenience if "testing", as they are otherwise slow to compute).

	- *output*:
		`trait_modification_output/model_files/`<br>
    `_NULL` files : traits without parent nodes in at least one, but not all of the k-fold DAGs are not modelled by a GP in models where it has no parents. Instead they are modelled as following a normal distribution, estimated from the data. These traits are called _NULL traits. Their sampling is required so that the predicted response to other perturbations for these traits can be averaged over all five-fold estimated DAG structures, not just those where they have parents.
		- `{TRAIT}_NULL.stan` : Stan source code for model used to sample from normal distribution of  all the _NULL traits with empirical mean and sd. {TRAIT} is the perturbed trait, not the name of the _NULL trait.
		- `{TRAIT}_NULL.rds` : compiled rStan model of the above `.stan` file
		- `{TRAIT}_NULL_fitted.rds` : sampled values for _NULL traits from the above model. NB that this file >100Mb, and so is not included in the github, but can be generated from the provided scripts and input data.
    <br>

      `_opt` files : Full Bayesian inference of the GP hyperparameters $\alpha$, $\rho$ & $\sigma$, by sampling at the same time as rest of the GP model is too computationally demanding. Instead here they are estimated by Regularised Maximum Marginal Likelihood, and point estimates are passed as parameters the rest of the model (as described [here](https://betanalpha.github.io/assets/case_studies/gp_part2/part2.html#3_regularized_maximum_marginal_likelihood)). Priors for the hyperparameters for each trait ($t$) are:

      $\alpha_{t} \sim Normal(0,1) $
      $\rho_{t} \sim InvGamma(5,5)$
      $\sigma_{t} \sim HalfNormal(0,1)$
      <br>

		 - `k{i}-fold/{TRAIT}/{TRAIT}_opt.stan` : stan source code for model used for maximum likelihood hyperparameter estimation. For use with perturbed trait {TRAIT}, and k-fold model structure {i}.
		  - `k{i}-fold/{TRAIT}/{TRAIT}_opt.rds` : compiled rstan model for the above.
      <br>

      GOT UP TO HERE

      `_pred` files.
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
