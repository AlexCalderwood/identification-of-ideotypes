# Identification-of-ideotypes

Supporting code associated with *"Bayesian optimisation over trait-relationship space for the quantitative identification of crop ideotypes in Brassica"*, with example input/output data. Scripts are provided to document the modelling carried out.

### data preprocessing
- **impute_data.R**: impute values for missing datapoints.

	- *input*: `/data/final_data_normed_cleaned.rds` - normalised data for all accessions. Transformations which have been already applied to normalise the data are shown in supplemental table 2.
	- *output*: `/data/final_data_normed_cleaned_imputed.rds` - missing values in input have been imputed using the **`mice_3.11.0`** R package.


### identification of trait relationship structure
- **id_trait_structure.R**: Identification of trait-trait relationship structure
directed acyclic graph (DAG).
	supporting script files: id_trait_structure_functions.R
	input: final_data_normed_cleaned_imputed.rds
	output:
		trait_structure_output/data/
			bnlearn_data.rds : data for spring ecotype accessions only. Used for learning model structures
			bnlearn_results.rds : includes identified model structure for each of the k-fold models identified by bnlearn.
		trait_structure_output/graphs/
			trait.value_vs_predictions.pdf : plots of predicted values vs observed values for withheld data for each of the k-fold models.
			k-fold_DAGs.pdf : DAGs showing each of the k-fold models identified.
			avg_DAG_strength.pdf : consensus DAG identified based on connections seen frequently & strongly in the k-fold identified structures.


### Predicting the yield consequence of modifying traits
- trait_modification_prediction.R
	supporting script files:
		trait_modification_functions.R - R functions used
		GP_ARD_functions.stan - Stan functions defining kernel used for GP prior.
	input:
		id_trait_structure_output/data/bnlearn_data.rds
		id_trait_structure_output/data/bnlearn_results.rds
		(if running with TESTING=TRUE:
		data/example_trait_modification_data/ki/height_fitted.rds - precomputed predictions for trait values when height is modified).
	output:
		trait_modification_output/model_files/
			height_NULL.stan - Stan source code for model to sample   from normal distribution with observed mean and sd.
			height_NULL.rds - compiled rStan model of the above
			height_NULL_fitted.rds - sampled values from the above model.

			/ki-fold/TRAIT/TRAIT_opt.stan - stan source code for model to use for max likelihood hyperparameter estimation, using perturbed trait TRAIT, and k-fold model structure i.
			/ki-fold/TRAIT/TRAIT_opt.rds - rstan model for the above

			/ki-fold/TRAIT/TRAIT_pred.stan - stan source code for model to make predictions for values of other traits, when trait TRAIT is modified.
			/ki-fold/TRAIT/TRAIT_pred.rds - rstan model for the above.
			/ki-fold/TRAIT/TRAIT_fitted.rds - sampled values from the posterior defined in the above model. If running with TESTING=TRUE, this is just copied from "data/example_trait_modification_data/ki/height_fitted.rds".

		trait_modification_output/graphs/TRAIT_k-avg_predicted_perturbation_effect.pdf - plot of CIs of trait values when perturbed trait is fixed to x-axis values.
		trait_modification_output/graphs/TRAIT_k-avg_predicted_perturbation_effect_back.pdf - as above, after data is back transformed to undo transformations to make more normally distributed.

		ki-fold/TRAIT/
			optimal_hyperparameters.csv - values of alpha, sigma, rho hyper parameters estimated by regularised maximum likelihood.
			fit_summary.csv - all values sampled from posterior (as TRAIT_fitter.rds, but in csv format).

		k-averaged/TRAIT/
			avg_pred_vs_measured.pdf - plot of observed values vs mean predictions from sampled posterior.
			avg_pred_vs_meaured.data.rds - data used in making above plot
			avg_pred_vs_measured_back.pdf - plot of observed values vs mean predictions from sampled posterior, after normalisation transformation back transformed.
			avg_pred_vs_meaured.data_back.rds - the data used in making the above plot.


### Bayesian optimisation for ideotype identification
- bayesian_optimisation.R
	supporting script files:
		BO_functions.R - functions used in bayesian_optimisation.R
		BO_plot_results.R - script to plot q-EI points, and observations
	input:
		data/BO_data/
			pred_model.stan - stan source code used to sample posterior prediction for next proposed q-point, conditioned on observations and hyperparamteres.
			pred_model.rds - rstan of the above model.
			spring_trait_data.rds - data for spring accessions only (as produced by id_trait_structure.R, renamed bnlearn_data.rds).
			witner_trait_data.rds - same as above, but winter OSR accessions.
	output:
		BO_output/
			ECOTYPE_hyperparam_model/
				opt_model.stan - stan source for model used to infer max likelihood values for GP prior hyperparams
				opt_model.rds - rstan compiled model of the above
				opt_hyperparams.csv - optimal hyperparameter values
		ECOTYPE_PCA_SUFFIX.csv - summary table of variance explained by principle components.
		ECOTYPE_stan.data_SIFFUX.RDS - the data used to the q-values for exploration
		ECOTYPE_q.values_SUFFIXX.RDS - the next q identified values for exploration.

### Simulation of causal SNP identification

-simSNP_inference.R
	supporting script files:
		simSNP_inference_functions.R
	output:
		simSNP_inference_out.csv - table of number of parent, child, & non-causal SNPS associated with child trait.
