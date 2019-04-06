/*
 * svime.cpp
 *
 *  Updated on: Feb 23, 2019
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
 */

#include "svime.h"
#include "distribution.h"
#include "util.h"
#include "asa103.hpp"
#include "processFasta.h"
#include <omp.h>
#include <boost/foreach.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <Eigen/StdVector>

svime::svime(const int w, const float a, const int e, const int max,
		const float* step, const int bs, const float t, const int jobs, const int rs) {
	window = w;
	alpha = a;
	epochs = e;
	max_clusters = max;
	batch_size = bs;
	step_pars[0] = step[0]; step_pars[1] = step[1];
	tol = t;
	n_jobs = jobs;
	random_state = rs;
}

svime::~svime() {
	// TODO Auto-generated destructor stub
}

Eigen::MatrixXd svime::binarize_seq(const std::vector<std::string>& X,
		const char base) {
	Eigen::MatrixXd binaryMatrix = Eigen::MatrixXd::Zero(X.size(), window);
	int n = X.size();
	for ( int i = 0; i < n; ++i ) { // iterate over w-mers
		for ( int p = 0; p < window; ++p ) { // iterate over positions in each w-mer
			// if position p in w-mer i is base then set the corresponding
			// element in binaryMatrix to 1
			if (X[i][p] == base) {
				binaryMatrix(i, p) = 1;
			}
		}
	}
	return binaryMatrix;
}

svime::psm svime::moments_base_probs(svime::psm* conc) {
	int ifault; // required for asa103's digamma implementation
	svime::psm moments = svime::psm(); // initialize return psm
	moments.a = Eigen::VectorXd::Zero(window);
	moments.c = Eigen::VectorXd::Zero(window);
	moments.g = Eigen::VectorXd::Zero(window);
	moments.t = Eigen::VectorXd::Zero(window);
	double concSum;
	for ( int p = 0; p < window; ++p ) {
		concSum = conc->a(p)+conc->c(p)+conc->g(p)+conc->t(p);
		// moments of log base probabilities
		moments.a(p) = digamma(conc->a(p), &ifault)-digamma(concSum, &ifault);
		moments.c(p) = digamma(conc->c(p), &ifault)-digamma(concSum, &ifault);
		moments.g(p) = digamma(conc->g(p), &ifault)-digamma(concSum, &ifault);
		moments.t(p) = digamma(conc->t(p), &ifault)-digamma(concSum, &ifault);
	}
	return moments;
}

std::pair<double, double> svime::moments_v(std::pair<double, double> sbPars) {
	int ifault;
	double sbSum = sbPars.first+sbPars.second;
	// moments of log stick-breaking variables
	std::pair<double, double> moments(digamma(sbPars.first, &ifault)-digamma(sbSum, &ifault), digamma(sbPars.second, &ifault)-digamma(sbSum, &ifault));
	return moments;
}

Eigen::VectorXd svime::clust_probs(Eigen::MatrixXd& Xa, Eigen::MatrixXd& Xc,
		Eigen::MatrixXd& Xg, Eigen::MatrixXd& Xt, svime::psm* moments) {
	Eigen::VectorXd Ezt;
	Ezt = Xa*moments->a + Xc*moments->c + Xg*moments->g + Xt*moments->t;
	return Ezt;
}

Eigen::MatrixXd svime::update_probs(Eigen::MatrixXd& Xa, Eigen::MatrixXd& Xc,
		Eigen::MatrixXd& Xg, Eigen::MatrixXd& Xt, svime::variationalDist& q) {
	int n = Xa.rows();
	Eigen::MatrixXd batchEz = Eigen::MatrixXd::Zero(n, max_clusters);
	omp_set_num_threads(n_jobs);
	// calculate cluster probabilities for each cluster on different threads
	#pragma omp parallel for
	for ( int k = 0; k < max_clusters; ++k ) {
		batchEz.col(k) = svime::clust_probs(Xa, Xc, Xg, Xt, &q.e_ln_base_probs[k]);
	}
	// log probabilities for cluster assignments
	std::vector<float> lnivs(max_clusters, 0);
	Eigen::RowVectorXd lnPc = Eigen::RowVectorXd::Zero(max_clusters);
	for ( int k = 0; k < max_clusters; ++k ) {
		lnivs[k] = q.e_ln_v[k].second;
		lnPc(k) = q.e_ln_v[k].first + std::accumulate(lnivs.begin(), lnivs.begin()+k, 0);
	}
	batchEz = batchEz.rowwise() + lnPc;
	// exp-normalize trick to avoid underflow
	batchEz = batchEz.colwise() - batchEz.rowwise().maxCoeff();
	batchEz = batchEz.array().exp().matrix();
	// normalize rows
	Eigen::VectorXd rowSums = batchEz.rowwise().sum();
	for ( int i = 0; i < n; ++i ) {
		batchEz.row(i) = batchEz.row(i)/rowSums(i);
	}
	return batchEz;
}

svime::psm svime::intermediate_conc(Eigen::MatrixXd& Xa, Eigen::MatrixXd& Xc,
		Eigen::MatrixXd& Xg, Eigen::MatrixXd& Xt, svime::psm* hyperparameters,
		Eigen::MatrixXd& Ez, int k, float multiplier) {
	svime::psm inter_conc = svime::psm();
	inter_conc.a = Eigen::VectorXd::Zero(window);
	inter_conc.c = Eigen::VectorXd::Zero(window);
	inter_conc.g = Eigen::VectorXd::Zero(window);
	inter_conc.t = Eigen::VectorXd::Zero(window);
	// calculate sum terms for natural gradients of concentration parameters
	Eigen::VectorXd XatEz = Xa.transpose()*Ez.col(k);
	Eigen::VectorXd XctEz = Xc.transpose()*Ez.col(k);
	Eigen::VectorXd XgtEz = Xg.transpose()*Ez.col(k);
	Eigen::VectorXd XttEz = Xt.transpose()*Ez.col(k);
	for ( int p = 0; p < window; ++p ) {
		// natural gradients of concentration parameters
		inter_conc.a(p) = hyperparameters->a(p) + multiplier*XatEz(p);
		inter_conc.c(p) = hyperparameters->c(p) + multiplier*XctEz(p);
		inter_conc.g(p) = hyperparameters->g(p) + multiplier*XgtEz(p);
		inter_conc.t(p) = hyperparameters->t(p) + multiplier*XttEz(p);
	}
	return inter_conc;
}

std::pair<double, double> svime::intermediate_sb(Eigen::MatrixXd& Ez, int k,
		float multiplier) {
	std::pair<double, double> inter_sb;
	// natural gradients of stick-breaking parameters
	inter_sb.first = 1 + multiplier*Ez.col(k).sum();
	if (k == max_clusters-1) {
		inter_sb.second = alpha;
	} else {
		inter_sb.second = alpha + multiplier*Ez.block(0, k+1, Ez.rows(), max_clusters-k-1).sum();
	}
	return inter_sb;
}

double svime::calculate_elbo(std::vector<std::string>& seqs,
		svime::variationalDist& q, svime::psm* hyperparameters,
		Eigen::MatrixXd& Ez, Eigen::VectorXd& lnPc) {
	int n = seqs.size();
	/* creates 3D vectors for moments and parameters for base probabilities for
	easy access. lnp[k][p](0) contains the expected log A probability for motif
	k at position p. lnp[k][p](1) has the C probability, lnp[k][p](2) has the G
	probability and lnp[k][p](3) has the T probability. conc is similar but for
	concentration parameters */
	std::vector<std::vector<Eigen::Vector4d,Eigen::aligned_allocator<Eigen::Vector4d>>> lnp(max_clusters, std::vector<Eigen::Vector4d,Eigen::aligned_allocator<Eigen::Vector4d>>(window, Eigen::Vector4d::Zero()));
	std::vector<std::vector<Eigen::Vector4d,Eigen::aligned_allocator<Eigen::Vector4d>>> conc(max_clusters, std::vector<Eigen::Vector4d,Eigen::aligned_allocator<Eigen::Vector4d>>(window, Eigen::Vector4d::Zero()));
	/* creates 2D vector for hyperparameter for base probabilities. hypers[p](0)
	contains the prior concentration parameter for A at position p. lnp[k][p](1)
	has the C parameter, lnp[k][p](2) has the G parameter and lnp[k][p](3) has
	the T parameter. */
	std::vector<Eigen::Vector4d,Eigen::aligned_allocator<Eigen::Vector4d>> hypers(window, Eigen::Vector4d::Zero());
	omp_set_num_threads(n_jobs);
	#pragma omp parallel for
	for ( int k = 0; k < max_clusters; ++k ) {
		for ( int p = 0; p < window; ++p ) {
			lnp[k][p](0) = q.e_ln_base_probs[k].a(p);
			lnp[k][p](1) = q.e_ln_base_probs[k].c(p);
			lnp[k][p](2) = q.e_ln_base_probs[k].g(p);
			lnp[k][p](3) = q.e_ln_base_probs[k].t(p);
			conc[k][p](0) = q.concentrations[k].a(p);
			conc[k][p](1) = q.concentrations[k].c(p);
			conc[k][p](2) = q.concentrations[k].g(p);
			conc[k][p](3) = q.concentrations[k].t(p);
		}
	}
	for ( int p = 0; p < window; ++p ) {
		hypers[p](0) = hyperparameters->a(p);
		hypers[p](1) = hyperparameters->c(p);
		hypers[p](2) = hyperparameters->g(p);
		hypers[p](3) = hyperparameters->t(p);
	}
	// vector of log likelihoods for each cluster
	std::vector<double> loglikSums(max_clusters, 0);
	#pragma omp parallel for
	for ( int k = 0; k < max_clusters; ++k ) {
		for ( int i = 0; i < n; ++i ) {
			double loglikLocal = 0;
			for ( int p = 0; p < window; ++p ) {
				loglikLocal += Dist::ln_multinomial(seqs[i][p], lnp[k][p]);
			}
			loglikSums[k] += Ez(i, k)*loglikLocal;
		}
	}
	// sum log likelihoods for each cluster
	double loglik = std::accumulate(loglikSums.begin(), loglikSums.end(), 0);
	// other terms for ELBO
	double lnPrior = 0, lnPv = 0, localSum = 0;
	#pragma omp parallel for
	for ( int k = 0; k < max_clusters; ++k ) {
		for ( int p = 0; p < window; ++p ) {
			lnPrior += Dist::ln_dirichlet(lnp[k][p], hypers[p]);
			localSum += Dist::ln_dirichlet(lnp[k][p], conc[k][p]);
		}
		lnPv += Dist::ln_beta(q.e_ln_v[k].first, q.e_ln_v[k].second, 1, alpha);
		localSum += Dist::ln_beta(q.e_ln_v[k].first, q.e_ln_v[k].second, q.sb[k].first, q.sb[k].second);
	}
	// log probabilities of cluster assignments
	double lnPz = Ez.colwise().sum().dot(lnPc);
	// variational entropy
	double lnq = Ez.array().pow(Ez.array()).log().sum() + localSum;
	double e = loglik+lnPrior+lnPv+lnPz-lnq; //ELBO
	return e;
}

svime::variationalDist svime::fit_predict(std::string outDir,
		std::map<std::string, int> chrSizes, svime::psm* hyperparameters) {
	svime::psm hyper = svime::psm();
	if (!hyperparameters) {
		std::cout << "hyperparameters are NULL. Setting default hyperparameters" << std::endl;
		hyperparameters = &hyper;
		hyperparameters->a = Eigen::VectorXd::Ones(window);
		hyperparameters->c = Eigen::VectorXd::Ones(window);
		hyperparameters->g = Eigen::VectorXd::Ones(window);
		hyperparameters->t = Eigen::VectorXd::Ones(window);
	}

	srand(random_state);
	std::pair<std::string, int> item;
	std::vector<std::string> chromosomes;
	std::string chr, probFile;
	int size, clustIdx, last;
	int n = 0; // total number of sequences
	Eigen::MatrixXd Ez;
	// make files to store cluster probability matrices for each chromosomes
	BOOST_FOREACH(item, chrSizes) {
		chr = item.first;
		size = item.second;
		chromosomes.push_back(chr);
		n += size;
		Ez = Eigen::MatrixXd::Zero(size, max_clusters);
		for ( int i = 0; i < size; ++i ) {
			clustIdx = rand() % max_clusters;
			Ez(i, clustIdx) = 1;
		}
		probFile = outDir+"/"+chr+"_prob.dat";
		Util::write_binary(probFile.c_str(), Ez);
	}
	// initialize variational distribution
	svime::variationalDist q = svime::variationalDist();
	for ( int k = 0; k < max_clusters; ++k ) {
		q.concentrations.push_back(*hyperparameters);
		q.sb.push_back(std::make_pair(1, alpha));
		q.e_ln_base_probs.push_back(svime::moments_base_probs(&q.concentrations[k]));
		q.e_ln_v.push_back(svime::moments_v(q.sb[k]));
	}
	q.elbo = std::numeric_limits<double>::lowest();
	// choose random test chromosome
	int testChrIdx = rand() % chromosomes.size();
	std::string testChr = chromosomes[testChrIdx];
	std::vector<std::string> testSeqs, chrSeqs, batchSeqs, regions;
	Eigen::MatrixXd testProbs, chrProbs, testA, testC, testG, testT,
	batchA, batchC, batchG, batchT;
	Util::loadSeqs(outDir+"/"+testChr+"_seq.txt", testSeqs);
	std::string testProbFile = outDir+"/"+testChr+"_prob.dat";
	testA = svime::binarize_seq(testSeqs, 'A');
	testC = svime::binarize_seq(testSeqs, 'C');
	testG = svime::binarize_seq(testSeqs, 'G');
	testT = svime::binarize_seq(testSeqs, 'T');
	float multiplier, step_size;

	int iteration = 1;
	std::vector<float> lnivs(max_clusters, 0);
	Eigen::VectorXd lnPc = Eigen::VectorXd::Zero(max_clusters);
	bool ELBOdescending = false; // track whether test ELBO is decreasing
	double prevELBO, ELBOgain;
	omp_set_num_threads(n_jobs);
	for ( int epoch = 0; epoch < epochs; ++epoch ) {
		std::cout << "Running epoch " << epoch+1 << std::endl;
		// permute chromosomes
		std::random_shuffle(chromosomes.begin(), chromosomes.end());
		for (auto const& batchChr: chromosomes) {
			std::cout << "Updating variational parameters with sequences from " << batchChr << std::endl;
			// load batch chromosome
			Util::loadSeqs(outDir+"/"+batchChr+"_seq.txt", chrSeqs);
			probFile = outDir+"/"+batchChr+"_prob.dat";
			Util::read_binary(probFile.c_str(), chrProbs);
			// iterate over batches
			for ( int b = 0; b < chrSizes[batchChr]; b += batch_size ) {
				// pick last w-mer in batch
				if (b+batch_size > chrSizes[batchChr]) {
					last = chrSizes[batchChr];
				} else {
					last = b+batch_size;
				}
				Ez.resize(last-b, max_clusters);
				batchSeqs = Util::slice(chrSeqs, b, last); // extract batch w-mers
				batchA = svime::binarize_seq(batchSeqs, 'A');
				batchC = svime::binarize_seq(batchSeqs, 'C');
				batchG = svime::binarize_seq(batchSeqs, 'G');
				batchT = svime::binarize_seq(batchSeqs, 'T');
				// cluster probabilities for batch
				Ez = svime::update_probs(batchA, batchC, batchG, batchT, q);
				// re-assign batch probabilities
				chrProbs.block(b, 0, last-b, max_clusters) = Ez;
				// multiplier for intermediate values
				multiplier = (float)n / (float)(last-b);
				// step size function
				step_size = std::pow((float)(iteration+step_pars[0]), -step_pars[1]);
				// Update parameters and moments on different threads
				#pragma omp parallel for
				for ( int k = 0; k < max_clusters; ++k ) {
					svime::psm inter_conc = svime::intermediate_conc(batchA,
							batchC, batchG, batchT, hyperparameters, Ez, k, multiplier);
					std::pair<double, double> inter_sb = svime::intermediate_sb(Ez, k, multiplier);
					// parameters are weighted averages of their current &
					// intermediate values
					q.concentrations[k].a = (1-step_size)*q.concentrations[k].a + step_size*inter_conc.a;
					q.concentrations[k].c = (1-step_size)*q.concentrations[k].c + step_size*inter_conc.c;
					q.concentrations[k].g = (1-step_size)*q.concentrations[k].g + step_size*inter_conc.g;
					q.concentrations[k].t = (1-step_size)*q.concentrations[k].t + step_size*inter_conc.t;
					q.sb[k].first = (1-step_size)*q.sb[k].first + step_size*inter_sb.first;
					q.sb[k].second = (1-step_size)*q.sb[k].second + step_size*inter_sb.second;
					q.e_ln_base_probs[k] = svime::moments_base_probs(&q.concentrations[k]);
					q.e_ln_v[k] = svime::moments_v(q.sb[k]);
				}
				iteration++;
			}
			// write cluster probabilities to file
			Util::write_binary(probFile.c_str(), chrProbs);
		}
		// Calculate test ELBO
		testProbs = svime::update_probs(testA, testC, testG, testT, q);
		Util::write_binary(testProbFile.c_str(), testProbs);
		for ( int k = 0; k < max_clusters; ++k ) {
			lnivs[k] = q.e_ln_v[k].second;
			lnPc(k) = q.e_ln_v[k].first + std::accumulate(lnivs.begin(), lnivs.begin()+k, 0);
		}
		prevELBO = q.elbo;
		q.elbo = svime::calculate_elbo(testSeqs, q, hyperparameters, testProbs, lnPc);
		ELBOgain = q.elbo-prevELBO;
		std::cout << "Sample ELBO: " << q.elbo << " ... gained " << ELBOgain << std::endl;
		// check convergence
		if ( (ELBOgain < tol) && ELBOdescending ) {
			std::cout << "Sample ELBO converged!" << std::endl;
			break;
		} else if (ELBOgain < tol) {
			ELBOdescending = true;
		} else {
			ELBOdescending = false;
		}
	}
	// write cluster assignments to csv
	std::cout << "Writing clusters to " << outDir << "/results/clusters.csv" << std::endl;
	std::vector<Eigen::MatrixXd::Index> batchCluster;
	std::ofstream clustersCSV;
	clustersCSV.open(outDir+"/results/clusters.csv");
	BOOST_FOREACH(item, chrSizes) {
		chr = item.first;
		size = item.second;
		Util::loadSeqs(outDir+"/"+chr+"_seq.txt", chrSeqs);
		probFile = outDir+"/"+chr+"_prob.dat";
		Util::read_binary(probFile.c_str(), chrProbs);
		Util::loadSeqs(outDir+"/"+chr+"_coords.txt", regions);
		for ( int b = 0; b < size; b += batch_size ) {
			if (b+batch_size > size) {
				last = size;
			} else {
				last = b+batch_size;
			}
			Ez.resize(last-b, max_clusters);
			batchSeqs = Util::slice(chrSeqs, b, last);
			batchA = svime::binarize_seq(batchSeqs, 'A');
			batchC = svime::binarize_seq(batchSeqs, 'C');
			batchG = svime::binarize_seq(batchSeqs, 'G');
			batchT = svime::binarize_seq(batchSeqs, 'T');
			Ez = svime::update_probs(batchA, batchC, batchG, batchT, q);
			chrProbs.block(b, 0, last-b, max_clusters) = Ez;
			batchCluster = Util::max_index(Ez);
			for ( int i = 0; i < last-b; ++i ) {
				clustersCSV << regions[b+i] << "," << batchCluster[i] << "\n";
			}
		}
		Util::write_binary(probFile.c_str(), chrProbs);
	}
	clustersCSV.close();
	// write concentration parameters to text file
	std::cout << "Writing variational position weight distributions to " << outDir << "/results/variationalPWD_motif*.txt" << std::endl;
	std::ofstream variationalPWDs;
	for ( int k = 0; k < max_clusters; ++k ) {
		variationalPWDs.open(outDir+"/results/variationalPWD_motif"+std::to_string(k)+".txt");
		variationalPWDs << "Motif " << k << "\n";
		for ( int p = 0; p < window; ++p ) {
			variationalPWDs << q.concentrations[k].a(p) << "\t";
		}
		variationalPWDs << "\n";
		for ( int p = 0; p < window; ++p ) {
			variationalPWDs << q.concentrations[k].c(p) << "\t";
		}
		variationalPWDs << "\n";
		for ( int p = 0; p < window; ++p ) {
			variationalPWDs << q.concentrations[k].g(p) << "\t";
		}
		variationalPWDs << "\n";
		for ( int p = 0; p < window; ++p ) {
			variationalPWDs << q.concentrations[k].t(p) << "\t";
		}
		variationalPWDs.close();
	}
	std::cout << "Done!\n";
	return q;
}
