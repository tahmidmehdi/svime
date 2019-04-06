/*
 * svime.h
 * Stochastic Variational Inference for Motif Elicitation
 *  Updated on: Feb 21, 2019
 *      Author: Tahmid Mehdi
 *
Copyright 2019 Tahmid Mehdi
This file is part of svime.

svime is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

svime is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with svime.  If not, see <http://www.gnu.org/licenses/>.

 *  References:
 *  [1] Dunson, D.B. & Xing, C. (2012) Nonparametric Bayes Modeling of
 *  Multivariate Categorical Data. J Am Stat Assoc., 104(487), 1042-1051.
 *
 *  [2] Hoffman, M.D., Blei, D.M., Wang, C., Paisley, J. (2013) Stochastic
 *  variational inference. Journal of Machine Learning Research, 14(1), 1303-1347.
 */

#ifndef SRC_SVIME_H_
#define SRC_SVIME_H_
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <Eigen/Dense>

/* This class implements a Dirichlet Process Mixture of Product-Multinomials
(DPMPM) for motif discovery as described in [1] and fits it with stochastic
variational inference (SVI) [2] */
class svime {
private:
	// alpha parameter for the Dirichlet Process. Controls how finely it searches for clusters.
	float alpha;
	float tol; // SVI stops when the gain in ELBO is less than tol in 2 consecutive epochs
	float step_pars[2]; // parameters for the step size function described in [2]
	int epochs; // maximum number of epochs
	int max_clusters; // maximum number of motifs
	int n_jobs; // number of cores
	int random_state; // random seed
	int batch_size; // number of w-mers in each batch
	int window; // number of bases in each window
public:
	/* Position score matrix with a vector for each base (A, C, G & T). Each
	vector should be window-length and the pth element of a vector represents
	its score for the corresponding base in position p. */
	struct psm {
		Eigen::VectorXd a, c, g, t;
	};
	/* Variational distribution learned by SVI. concentrations contain the
	parameters for the Dirichlets that generate base probabilities.
	e_ln_base_probs contain expected log base probabilities. sb contains
	stick-breaking parameters and e_ln_v contains expected logs of the
	stick-breaking parameters. */
	struct variationalDist {
		std::vector<psm> concentrations, e_ln_base_probs;
		std::vector<std::pair<double, double>> sb, e_ln_v;
		double elbo;
	};
	// Constructs a svime object and sets its private member
	svime(const int w, const float a, const int e, const int max, const float* step, const int bs = 1000, const float t = 0.001, const int jobs = 1, const int rs = 42);
	// Destructor
	virtual ~svime();
	// Take expectations of log base probabilities
	psm moments_base_probs(psm* conc);
	// Take expectations of log stick-breaking variables
	std::pair<double, double> moments_v(std::pair<double, double> sbPars);
	// Calculate cluster probabilities for cluster associated with moments
	Eigen::VectorXd clust_probs(Eigen::MatrixXd& Xa, Eigen::MatrixXd& Xc, Eigen::MatrixXd& Xg, Eigen::MatrixXd& Xt, psm* moments);
	// Produces binary matrix that indicates positions in X which are base
	Eigen::MatrixXd binarize_seq(const std::vector<std::string>& X, const char base);
	// Calculates cluster probability matrix where rows are data points and
	// columns are clusters
	Eigen::MatrixXd update_probs(Eigen::MatrixXd& Xa, Eigen::MatrixXd& Xc, Eigen::MatrixXd& Xg, Eigen::MatrixXd& Xt, variationalDist& q);
	// Calculates intermediate values of concentration parameters
	psm intermediate_conc(Eigen::MatrixXd& Xa, Eigen::MatrixXd& Xc, Eigen::MatrixXd& Xg, Eigen::MatrixXd& Xt, psm* hyperparameters, Eigen::MatrixXd& Ez, int k, float multiplier);
	// Calculates intermediate values of stick-breaking parameters
	std::pair<double, double> intermediate_sb(Eigen::MatrixXd& Ez, int k, float multiplier);
	// Calculates the ELBO
	double calculate_elbo(std::vector<std::string>& seqs, variationalDist& q, psm* hyperparameters, Eigen::MatrixXd& Ez, Eigen::VectorXd& lnPc);
	// Fits the DPMPM with SVI
	variationalDist fit_predict(std::string outDir, std::map<std::string, int> chrSizes, psm* hyperparameters = NULL);

};

#endif /* SRC_SVIME_H_ */
