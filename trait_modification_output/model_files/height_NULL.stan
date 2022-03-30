data {
	int<lower=0> N;
	vector[N] timeToFlower;
	vector[N] ovaryLen;
	vector[N] ovuleArea;
	vector[N] ovuleACov;
	vector[N] numFlowers;
	vector[N] ovuleNum;
	vector[N] gynLen;
	vector[N] percAbortedMain;
	vector[N] numPodsMain;
	vector[N] podLen;
	vector[N] percAborted2ary;
	vector[N] beakLen;
	vector[N] numPods2ary;
	vector[N] timeToMature;
	vector[N] tenPodNumber;
	vector[N] tenPodArea;
	vector[N] tenPodWeight;
	vector[N] totSeedArea;
	vector[N] totSeedNum;
	vector[N] seedACov;
	vector[N] totCompactness;
	vector[N] TGW;
	vector[N] tenPodCompactness;
	vector[N] totSeedW;
	vector[N] oilContent;
	vector[N] styleLen;
}

parameters {
	real mu_timeToFlower;
	real<lower=0> sigma_timeToFlower;

	real mu_ovaryLen;
	real<lower=0> sigma_ovaryLen;

	real mu_ovuleArea;
	real<lower=0> sigma_ovuleArea;

	real mu_ovuleACov;
	real<lower=0> sigma_ovuleACov;

	real mu_numFlowers;
	real<lower=0> sigma_numFlowers;

	real mu_ovuleNum;
	real<lower=0> sigma_ovuleNum;

	real mu_gynLen;
	real<lower=0> sigma_gynLen;

	real mu_percAbortedMain;
	real<lower=0> sigma_percAbortedMain;

	real mu_numPodsMain;
	real<lower=0> sigma_numPodsMain;

	real mu_podLen;
	real<lower=0> sigma_podLen;

	real mu_percAborted2ary;
	real<lower=0> sigma_percAborted2ary;

	real mu_beakLen;
	real<lower=0> sigma_beakLen;

	real mu_numPods2ary;
	real<lower=0> sigma_numPods2ary;

	real mu_timeToMature;
	real<lower=0> sigma_timeToMature;

	real mu_tenPodNumber;
	real<lower=0> sigma_tenPodNumber;

	real mu_tenPodArea;
	real<lower=0> sigma_tenPodArea;

	real mu_tenPodWeight;
	real<lower=0> sigma_tenPodWeight;

	real mu_totSeedArea;
	real<lower=0> sigma_totSeedArea;

	real mu_totSeedNum;
	real<lower=0> sigma_totSeedNum;

	real mu_seedACov;
	real<lower=0> sigma_seedACov;

	real mu_totCompactness;
	real<lower=0> sigma_totCompactness;

	real mu_TGW;
	real<lower=0> sigma_TGW;

	real mu_tenPodCompactness;
	real<lower=0> sigma_tenPodCompactness;

	real mu_totSeedW;
	real<lower=0> sigma_totSeedW;

	real mu_oilContent;
	real<lower=0> sigma_oilContent;

	real mu_styleLen;
	real<lower=0> sigma_styleLen;

}

model {
	timeToFlower ~ normal(mu_timeToFlower, sigma_timeToFlower);
	ovaryLen ~ normal(mu_ovaryLen, sigma_ovaryLen);
	ovuleArea ~ normal(mu_ovuleArea, sigma_ovuleArea);
	ovuleACov ~ normal(mu_ovuleACov, sigma_ovuleACov);
	numFlowers ~ normal(mu_numFlowers, sigma_numFlowers);
	ovuleNum ~ normal(mu_ovuleNum, sigma_ovuleNum);
	gynLen ~ normal(mu_gynLen, sigma_gynLen);
	percAbortedMain ~ normal(mu_percAbortedMain, sigma_percAbortedMain);
	numPodsMain ~ normal(mu_numPodsMain, sigma_numPodsMain);
	podLen ~ normal(mu_podLen, sigma_podLen);
	percAborted2ary ~ normal(mu_percAborted2ary, sigma_percAborted2ary);
	beakLen ~ normal(mu_beakLen, sigma_beakLen);
	numPods2ary ~ normal(mu_numPods2ary, sigma_numPods2ary);
	timeToMature ~ normal(mu_timeToMature, sigma_timeToMature);
	tenPodNumber ~ normal(mu_tenPodNumber, sigma_tenPodNumber);
	tenPodArea ~ normal(mu_tenPodArea, sigma_tenPodArea);
	tenPodWeight ~ normal(mu_tenPodWeight, sigma_tenPodWeight);
	totSeedArea ~ normal(mu_totSeedArea, sigma_totSeedArea);
	totSeedNum ~ normal(mu_totSeedNum, sigma_totSeedNum);
	seedACov ~ normal(mu_seedACov, sigma_seedACov);
	totCompactness ~ normal(mu_totCompactness, sigma_totCompactness);
	TGW ~ normal(mu_TGW, sigma_TGW);
	tenPodCompactness ~ normal(mu_tenPodCompactness, sigma_tenPodCompactness);
	totSeedW ~ normal(mu_totSeedW, sigma_totSeedW);
	oilContent ~ normal(mu_oilContent, sigma_oilContent);
	styleLen ~ normal(mu_styleLen, sigma_styleLen);
}

generated quantities {
	vector[N] timeToFlower_rep;
	vector[N] ovaryLen_rep;
	vector[N] ovuleArea_rep;
	vector[N] ovuleACov_rep;
	vector[N] numFlowers_rep;
	vector[N] ovuleNum_rep;
	vector[N] gynLen_rep;
	vector[N] percAbortedMain_rep;
	vector[N] numPodsMain_rep;
	vector[N] podLen_rep;
	vector[N] percAborted2ary_rep;
	vector[N] beakLen_rep;
	vector[N] numPods2ary_rep;
	vector[N] timeToMature_rep;
	vector[N] tenPodNumber_rep;
	vector[N] tenPodArea_rep;
	vector[N] tenPodWeight_rep;
	vector[N] totSeedArea_rep;
	vector[N] totSeedNum_rep;
	vector[N] seedACov_rep;
	vector[N] totCompactness_rep;
	vector[N] TGW_rep;
	vector[N] tenPodCompactness_rep;
	vector[N] totSeedW_rep;
	vector[N] oilContent_rep;
	vector[N] styleLen_rep;

	for (n in 1:N) {
		timeToFlower_rep[n] = normal_rng(mu_timeToFlower, sigma_timeToFlower);
		ovaryLen_rep[n] = normal_rng(mu_ovaryLen, sigma_ovaryLen);
		ovuleArea_rep[n] = normal_rng(mu_ovuleArea, sigma_ovuleArea);
		ovuleACov_rep[n] = normal_rng(mu_ovuleACov, sigma_ovuleACov);
		numFlowers_rep[n] = normal_rng(mu_numFlowers, sigma_numFlowers);
		ovuleNum_rep[n] = normal_rng(mu_ovuleNum, sigma_ovuleNum);
		gynLen_rep[n] = normal_rng(mu_gynLen, sigma_gynLen);
		percAbortedMain_rep[n] = normal_rng(mu_percAbortedMain, sigma_percAbortedMain);
		numPodsMain_rep[n] = normal_rng(mu_numPodsMain, sigma_numPodsMain);
		podLen_rep[n] = normal_rng(mu_podLen, sigma_podLen);
		percAborted2ary_rep[n] = normal_rng(mu_percAborted2ary, sigma_percAborted2ary);
		beakLen_rep[n] = normal_rng(mu_beakLen, sigma_beakLen);
		numPods2ary_rep[n] = normal_rng(mu_numPods2ary, sigma_numPods2ary);
		timeToMature_rep[n] = normal_rng(mu_timeToMature, sigma_timeToMature);
		tenPodNumber_rep[n] = normal_rng(mu_tenPodNumber, sigma_tenPodNumber);
		tenPodArea_rep[n] = normal_rng(mu_tenPodArea, sigma_tenPodArea);
		tenPodWeight_rep[n] = normal_rng(mu_tenPodWeight, sigma_tenPodWeight);
		totSeedArea_rep[n] = normal_rng(mu_totSeedArea, sigma_totSeedArea);
		totSeedNum_rep[n] = normal_rng(mu_totSeedNum, sigma_totSeedNum);
		seedACov_rep[n] = normal_rng(mu_seedACov, sigma_seedACov);
		totCompactness_rep[n] = normal_rng(mu_totCompactness, sigma_totCompactness);
		TGW_rep[n] = normal_rng(mu_TGW, sigma_TGW);
		tenPodCompactness_rep[n] = normal_rng(mu_tenPodCompactness, sigma_tenPodCompactness);
		totSeedW_rep[n] = normal_rng(mu_totSeedW, sigma_totSeedW);
		oilContent_rep[n] = normal_rng(mu_oilContent, sigma_oilContent);
		styleLen_rep[n] = normal_rng(mu_styleLen, sigma_styleLen);
	}
}
