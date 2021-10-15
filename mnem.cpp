#include <iostream>
#include <fstream>
#include "functions.h"
#include <string>
#include <cmath>
#include <cstring>
using namespace std;

int main(int argc, char** argv)
{
	srand((unsigned)time(NULL));
	int n = 0;
	int m = 0;
	int l = 0;
	string file = "none";
	string result = "mnem_result.txt";
	int simulate = 0;
	int K = 2;
	int type = 2;
	double eps = 0.00001;
	int maxIter = 100;
	int init = 1;
	int null = 1;
	// simulation variables:
	string sfile = "mnem_gtn.txt";
	double mu = 0;
	double sigma = 0.1;
	double p = 0.8;
	double* pi;
	double** P;
	double*** data;
	double** phi1;
	double** theta1;
	double*** phi;
	double*** theta;
	int* r;
	for (int i = 0; i < argc; i++) {
		if (std::strcmp(argv[i], "-h") == 0 || argc == 1) {
			cout << "Mixture Nested Effects Models (version 0.9.9)\n";
			cout << "\n";
			cout << "-n		number of S-genes" << endl;
			cout << "-m		number of E-genes" << endl;
			cout << "-l		number of samples" << endl;
			cout << "-K		number of components (default: 2)" << endl;
			cout << "-type		exhaustive (1) or greedy (2, default) search" << endl;
			cout << "-eps		log-likelihood convergence margin (default: 0.00001)" << endl;
			cout << "-maxIter	number of maximum iterations for EM algorithm (default: 100)" << endl;
			cout << "-in		filename of the input data; must be tab delimited and ordered S-gene names in first line;\n";
			cout << "		S-genes must be named numerically from 0 to n-1" << endl;
			cout << "-out		filename of output information (default: mnem_result.txt)\n";
			cout << "-s		set to 1 for data simulation or to 0 (default) for inference\n";
			cout << "-sout		filename for the ground truth networks of the simulated data (default: mnem_gtn.txt)\n";
			cout << "-mu		mean distance to 1 respectivley -1 for effects (>0) and no effects (<0) (default: 0),\n";
			cout << "		e.g., mu = 1 means both distribution have mean 0 and are indistinguishable\n";
			cout << "-sigma		standard deviation of the effects (default: 0.1)\n";
			cout << "-p		edge probability (default: 0.8)\n";
			cout << "-init		initialise with random networks (1, default) or random posteriors (2)\n";
			cout << "-null		include (1, default) null model (no (self-)edges) in the inference or not (0)\n";
			cout << "\n";
			cout << "Reference:\n";
			cout << "Pirkl M, Beerenwinkel N (2018). Single cell network analysis with a mixture\n";
			cout << "of Nested Effects Models. Bioinformatics, 34(258202).https://doi.org/10.1093/bioinformatics/bty602.";
		}
		if (std::strcmp(argv[i], "-n") == 0) {
			n = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-m") == 0) {
			m = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-l") == 0) {
			l = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-K") == 0) {
			K = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-type") == 0) {
			type = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-eps") == 0) {
			eps = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-maxIter") == 0) {
			maxIter = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-in") == 0) {
			file = argv[i + 1];
		}
		else if (std::strcmp(argv[i], "-out") == 0) {
			result = argv[i + 1];
		}
		else if (std::strcmp(argv[i], "-s") == 0) {
			simulate = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-sout") == 0) {
			sfile = argv[i + 1];
		}
		else if (std::strcmp(argv[i], "-mu") == 0) {
			mu = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-sigma") == 0) {
			sigma = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-p") == 0) {
			p = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-init") == 0) {
			init = atof(argv[i + 1]);
		}
		else if (std::strcmp(argv[i], "-null") == 0) {
			null = atof(argv[i + 1]);
		}
	}
	if (((file == "none" && simulate == 0) || n == 0 || m == 0 || l == 0) && argc <= 1) {
		cout << "\n";
	}
	else if (((file == "none" && simulate == 0) || n == 0 || m == 0 || l == 0) && argc > 1) {
		cout << "-n,-m,-l and -in are essential, but not provided. Use -h for more information." << endl;
	}
	else {
		if (simulate == 1) {
			pi = new double[K];
			double piSum = 0;
			pi = runif(K);
			for (int i = 0; i < K; i++) {
				//pi[i] = 1 / K;
				piSum += pi[i];
			}
			for (int i = 0; i < K; i++) {
				pi[i] = pi[i] / piSum;
			}
			r = new int[n];
			int** r2;
			if (K <= 1) {
				phi1 = randPhi(n, n, p);
				phi1 = transClose(phi1, n, n);
				theta1 = randTheta(n, m);
				theta1 = transpose(theta1, m, n);
				for (int i = 0; i < n; i++) {
					r[i] = 1;
				}
				data = nemSim(phi1, theta1, n, m, r, mu, sigma);
				printMat(phi1, n, n);
			}
			else {
				phi = new double** [K];
				theta = new double** [K];
				r2 = new int* [n];
				for (int i = 0; i < n; i++) {
					r2[i] = new int[K];
					r[i] = 0;
					for (int k = 0; k < K; k++) {
						r2[i][k] = round((pi[k] * l) / n);
						r[i] = r[i] + r2[i][k];
					}
				}
				for (int k = 0; k < K; k++) {
					phi[k] = randPhi(n, n, p);
					phi[k] = transClose(phi[k], n, n);
					theta[k] = randTheta(n, m);
					theta[k] = transpose(theta[k], m, n);
				}
				data = mnemSim(phi, theta, pi, n, m, K, l, mu, sigma);
				phi1 = phi[1];
				theta1 = theta[1];
				writePhis(phi, pi, n, K, sfile);
			}
			writeData(data, n, m, r, result);
		}
		else {
			fstream dataFile;
			double*** data;
			data = new double** [n];
			double Num;
			int ncount;
			int* r;
			r = new int[n];
			int* colnames;
			colnames = new int[l];
			for (int j = 0; j < n; j++) {
				data[j] = new double* [m];
				ncount = 0;
				dataFile.open(file, ios::in);
				for (int k = 0; k <= m; k++) {
					if (k > 0) {
						data[j][k - 1] = new double[r[j]];
					}
					ncount = 0;
					for (int i = 0; i < l; i++) {
						dataFile >> Num;
						if (k == 0) {
							colnames[i] = Num;
							if (Num == j) {
								ncount++;
							}
						}
						else {
							if (colnames[i] == j && k > 0) {
								data[j][k - 1][ncount] = Num;
								ncount++;
							}
						}
					}
					if (k == 0) {
						r[j] = ncount;
					}
				}
				dataFile.close();
			}
			if (null == 0 && K == 1) {
				double** phis;
				phis = nem(data, n, m, r, type, result);
			}
			else {
				double*** phis;
				phis = mnem(data, n, m, r, l, K, type, eps, maxIter, init, p, result, null);
			}
		}
	}
	return 0;
}