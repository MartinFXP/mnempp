#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstring>
using namespace std;

float acc(double** phi1, double** phi2, int r, int c, int type = 1) {
	double acc = 0;
	double tp = 0;
	double fp = 0;
	double fn = 0;
	if (type == 1) {
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				if (phi1[i][j] == phi2[i][j]) {
					acc++;
				}
			}
		}
		acc = acc / (r * c);
	}
	else if (type == 2 || type == 3) {
		for (int i = 0; i < r; i++) {
			for (int j = 0; j < c; j++) {
				if (phi1[i][j] == 1 && phi2[i][j] == 1) {
					tp++;
				}
				if (phi1[i][j] == 0 && phi2[i][j] == 1) {
					fp++;
				}
				if (phi1[i][j] == 1 && phi2[i][j] == 0) {
					fn++;
				}
			}
		}
	}
	if (type == 2) {
		acc = tp / (tp + fp);
	}
	if (type == 3) {
		acc = tp / (tp + fn);
	}
	return acc;
}

double* runif(int N, double l = 0, double u = 1) {
	double* data;
	data = new double[N];
	double a = 0;
	double b = 0;
	for (int i = 0; i < N; i++) {
		a = rand();
		b = rand();
		if (a >= b) {
			data[i] = (b / a) * abs(u - l) + l;
		}
		else {
			data[i] = (a / b)* abs(u - l) + l;
		}
	}
	return data;
}

double pi(int iter = 100) {
	double piEst = 0;
	double l = 0;
	double* a;
	a = new double[2];
	for (int i = 0; i < iter; i++) {
		a = runif(2);
		l = sqrt(a[0] * a[0] + a[1] * a[1]);
		if (l <= 1) {
			piEst++;
		}
	}
	piEst = (piEst / iter) * 4;
	return piEst;
}

double pdf(double x, double mu, double sigma)
{
	double f;
	f = (1 / (sigma * sqrt(2 * (21.991 / 7)))) * exp(-0.5 * (((x - mu) / sigma) * ((x - mu) / sigma)));
	return f;
}

double* rnorm(int N, double mu, double sigma, int ilen = 10, int burn = 100, int type = 1) {
	int count = 0;
	double a;
	double b;
	double cand;
	double* data = new double[N];
	if (type == 1) {
		for (int i = 0; i < N; i++) {
			data[i] = 0;
		}
		double ilenD = ilen;
		double frac;
		double sample;
		double dOld;
		double dNew;
		sample = mu;
		dOld = pdf(sample, mu, sigma);
		for (int j = 0; j < N * ilen + burn + 1; j++) {
			a = rand();
			b = rand();
			if (a >= b) {
				cand = sample - (b / a) * sigma;
			}
			else {
				cand = sample + (a / b) * sigma;
			}
			dNew = pdf(cand, mu, sigma);
			a = rand();
			b = rand();
			if (a >= b) {
				a = b / a;
			}
			else {
				a = a / b;
			}
			frac = j / ilenD;
			if (a < dNew / dOld) {
				dOld = dNew;
				sample = cand;
			}
			if (round(frac) == frac && j > burn) {
				data[count] = cand;
				count++;
			}
		}
	}
	else if (type == 2) {
		double muNew = 0;
		double uniSum = 0;
		for (int i = 0; i < N; i++) {
			uniSum = 0;
			for (int j = 0; j < 100; j++) {
				uniSum = uniSum + rand();
			}
			data[i] = uniSum;
			muNew = muNew + (uniSum / N);
		}
		double sigmaNew = 0;
		for (int i = 0; i < N; i++) {
			sigmaNew += ((data[i] - muNew) * (data[i] - muNew)) / N;
		}
		sigmaNew = sqrt(sigmaNew);
		for (int i = 0; i < N; i++) {
			data[i] = data[i] / sigmaNew - muNew / sigmaNew;
			data[i] = data[i] * sigma + mu;
		}
	}
	else if (type == 3) {
		int reject = 1;
		double dens = 0;
		for (int i = 0; i < N; i++) {
			reject = 1;
			while (reject == 1) {
				a = rand();
				b = rand();
				if (a >= b) {
					a = b / a;
				}
				else {
					a = a / b;
				}
				cand = a * sigma * 10 - 5 * sigma + mu;
				dens = pdf(cand, mu, sigma);
				dens = dens / (10 * (1 / (sigma * 10)));
				a = rand();
				b = rand();
				if (a >= b) {
					a = b / a;
				}
				else {
					a = a / b;
				}
				if (dens > a) {
					reject = 0;
					data[i] = cand;
				}
				//cout << dens << " " << a << endl;
			}
		}
	}
	return data;
}

double* rnormMix(int N, int K, double* pi, double* mu, double* sigma, int ilen = 100, int burn = 1000, int type = 1)
{
	double* data = new double[N];
	for (int i = 0; i < N; i++) {
		data[i] = 0;
	}
	int* Nk = new int[K];
	double* data_tmp;
	int count = 0;
	for (int i = 0; i < K; i++) {
		Nk[i] = round(pi[i] * N);
		data_tmp = rnorm(Nk[i], mu[i], sigma[i], ilen, burn, type);
		for (int j = 0; j < Nk[i]; j++) {
			data[count] = data_tmp[j];
			count++;
		}
	}
	ofstream myfile("temp.txt");
	if (myfile.is_open())
	{
		for (int i = 0; i < N; i++) {
			myfile << data[i] << endl;
		}
		myfile.close();
	}
	else cout << "Unable to open file";
	return data;
}

double* gmm(int N, int K, double* data, double* pi, double* mu, double* sigma, int maxIter) {
	double** post = 0;
	post = new double* [N];
	for (int i = 0; i < N; i++) {
		post[i] = new double[K];
	}
	double* postSum = new double[K];
	for (int i = 0; i < K; i++) {
		postSum[i] = 0;
	}
	for (int i = 0; i < N; i++) {
		double llSum = 0;
		for (int j = 0; j < K; j++) {
			post[i][j] = pdf(data[i], mu[j], sigma[j]);
			llSum += post[i][j] * pi[j];
		}
		for (int j = 0; j < K; j++) {
			post[i][j] = (post[i][j] * pi[j]) / llSum;
			postSum[j] += post[i][j];
		}
	}
	for (int i = 0; i < K; i++) {
		pi[i] = 0;
		for (int j = 0; j < N; j++) {
			pi[i] += post[j][i] / N;
		}
	}
	for (int i = 0; i < K; i++) {
		mu[i] = 0;
		for (int j = 0; j < N; j++) {
			mu[i] += (data[j] * post[j][i]) / postSum[i];
		}
		sigma[i] = 0;
		for (int j = 0; j < N; j++) {
			sigma[i] += (post[j][i] * (data[j] - mu[i]) * (data[j] - mu[i])) / postSum[i];
		}
		sigma[i] = sqrt(sigma[i]);
	}
	int count = 0;
	for (int i = 0; i < K; i++) {
		cout << mu[i] << " " << sigma[i] << " " << pi[i] << endl;
	}
	while (count < maxIter)
	{
		count++;
		//cout << "iteration " << count << endl;
		for (int i = 0; i < K; i++) {
			postSum[i] = 0;
		}
		for (int i = 0; i < N; i++) {
			double llSum = 0;
			for (int j = 0; j < K; j++) {
				post[i][j] = pdf(data[i], mu[j], sigma[j]);
				llSum += post[i][j] * pi[j];
			}
			for (int j = 0; j < K; j++) {
				post[i][j] = (post[i][j] * pi[j]) / llSum;
				postSum[j] += post[i][j];
			}
		}
		for (int i = 0; i < K; i++) {
			pi[i] = 0;
			for (int j = 0; j < N; j++) {
				pi[i] += post[j][i] / N;
			}
		}
		for (int i = 0; i < K; i++) {
			mu[i] = 0;
			for (int j = 0; j < N; j++) {
				mu[i] += (data[j] * post[j][i]) / postSum[i];
			}
			sigma[i] = 0;
			for (int j = 0; j < N; j++) {
				sigma[i] += (post[j][i] * (data[j] - mu[i]) * (data[j] - mu[i])) / postSum[i];
			}
			sigma[i] = sqrt(sigma[i]);
		}
		/*for (int i = 0; i < K; i++) {
			cout << mu[i] << " " << sigma[i] << " " << pi[i] << endl;
		}*/
	}
	for (int i = 0; i < K; i++) {
		cout << mu[i] << " " << sigma[i] << " " << pi[i] << endl;
	}
	return 0;
}

void printMat(double** x, int r, int c) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			std::cout << x[i][j] << " ";
		}
		std::cout << endl;
	}
}

void printData(double*** data, int n, int m, int* r, int rnd = 2) {
	for (int k = 0; k < m; k++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < r[i]; j++) {
				std::cout << round(pow(10, rnd) * data[i][k][j]) / pow(10, rnd) << " ";
			}
		}
		std::cout << endl;
	}
}

void writeData(double*** data, int n, int m, int* r, string file = "mnem_data.txt") {
	ofstream myfile(file);
	if (myfile.is_open())
	{
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < r[j]; k++) {
				myfile << j;
				myfile << "\t";
			}
		}
		myfile << endl;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < r[j]; k++) {
					myfile << data[j][i][k];
					myfile << "\t";
				}
			}
			myfile << endl;
		}
		myfile.close();
	}
	else cout << "Unable to open file";
}

void writePhis(double*** phis, double* pi, int n, int K, string file = "mnem_phis.txt") {
	ofstream myfile(file);
	if (myfile.is_open())
	{
		for (int i = 0; i < K; i++) {
			myfile << "COMPONENT " << i << endl;
			myfile << "mixture weight:" << endl;
			myfile << pi[i] << endl;
			myfile << "phi:" << endl;
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					myfile << phis[i][j][k];
					myfile << "\t";
				}
				myfile << endl;
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file";
}

double** transClose(double** x, int r, int c) {
	for (int k = 0; k < r; k++) {
		x[k][k] = 1;
		for (int i = 0; i < c; i++) {
			for (int j = 0; j < r; j++) {
				x[j][i] = (x[j][i] || (x[k][i] && x[j][k]));
			}
		}
	}
	return x;
}

double** copyPhi(double** phi, int r, int c) {
	static double** phi2;
	phi2 = new double* [r];
	for (int i = 0; i < r; i++) {
		phi2[i] = new double[c];
		for (int j = 0; j < c; j++) {
			phi2[i][j] = phi[i][j];
		}
	}
	return phi2;
}

double** randPhi(int r = 5, int c = 5, double p = 0.75)
{
	double** M;
	M = new double*[r];
	double a;
	double b;
	for (int i = 0; i < r; i++) {
		M[i] = new double[c];
		for (int j = 0; j < c; j++) {
			a = rand();
			b = rand();
			if (a >= b) {
				a = b / a;
			}
			else {
				a = a / b;
			}
			if (a > p && j > i) {
				M[i][j] = 1;
			} else {
				M[i][j] = 0;
			}
		}
	}
	return M;
}

int* hex2bin(double h, int l) {
	int* b;
	b = new int[l];
	for (int i = l-1; i >= 0; i -= 1) {
		if (h - pow(2, i) >= 0) {
			b[i] = 1;
			h = h - pow(2, i);
		}
		else {
			b[i] = 0;
		}
	}
	return b;
}

double** idxPhi(int r, int c, double idx) {
	int* b;
	b = hex2bin(idx, r * c);
	double** M = 0;
	M = new double* [r];
	int count = 0;
	for (int i = 0; i < r; i++) {
		M[i] = new double[c];
		for (int j = 0; j < c; j++) {
			M[i][j] = b[count];
			count++;
		}
	}
	delete[] b;
	b = NULL;
	M = transClose(M, r, c);
	return M;
}

double** transpose(double** x, int r, int c) {
	double** M;
	M = new double* [c];
	for (int i = 0; i < c; i++) {
		M[i] = new double[r];
		for (int j = 0; j < r; j++) {
			M[i][j] = x[j][i];
		}
	}
	return M;
}

double** randTheta(int r = 5, int c = 10) {
	double** M = 0;
	M = new double* [c];
	double a;
	double b;
	for (int i = 0; i < c; i++) {
		M[i] = new double[r];
		a = rand();
		b = rand();
		if (a >= b) {
			a = (b / a) * r;
		}
		else {
			a = (a / b) * r;
		}
		for (int j = 0; j < r; j++) {
			M[i][j] = 0;
			if (a > j && a <= j + 1) {
				M[i][j] = 1;
			}
		}
	}
	return M;
}

double** matMult(double** A, double** B, int r, int c, int c2)
{
	double** C;
	C = new double* [r];
	for (int i = 0; i < r; i++) {
		C[i] = new double[c2];
		for (int j = 0; j < c2; j++) {
			C[i][j] = 0;
			for (int k = 0; k < c; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return C;
}

double nemLlr(double** P, int n, int m) {
	double llr = 0;
	double ellr = 0;
	for (int i = 0; i < m; i++) {
		ellr = 0;
		for (int j = 0; j < n; j++) {
			if (ellr < P[i][j]) {
				ellr = P[i][j];
			}
		}
		llr = llr + ellr;
	}
	return llr;
}

double*** nemSim(double** phi, double** theta, int n, int m, int* r, double mu = 0, double sigma = 1) {
	double*** data;
	data = new double** [n];
	double** F;
	F = matMult(phi, theta, n, n, m);
	F = transpose(F, n, m);
	double* noise;
	int rSum = 0;
	for (int i = 0; i < n; i++) {
		rSum += r[i];
	}
	noise = rnorm(n * m * rSum, mu, sigma);
	int count = 0;
	for (int j = 0; j < n; j++) {
		data[j] = new double* [m];
		for (int i = 0; i < m; i++) {
			data[j][i] = new double[r[j]];
			for (int k = 0; k < r[j]; k++) {
				if (F[i][j] == 0) {
					 data[j][i][k] = -1 + mu + noise[count];
				}
				else {
					data[j][i][k] = 1 - mu + noise[count];
				}
				count++;
			}
		}
	}
	return data;
}

double** nemMean(double*** data, int n, int m, int* r) {
	double llrM;
	double** dataM;
	dataM = new double* [m];
	for (int i = 0; i < m; i++) {
		dataM[i] = new double[n];
		for (int j = 0; j < n; j++) {
			llrM = 0;
			for (int k = 0; k < r[j]; k++) {
				llrM += data[j][i][k];
			}
			dataM[i][j] = llrM / r[j];
		}
	}
	return dataM;
}

double** nem(double*** data, int n, int m, int* r, int type = 2, string file = "none") {
	double** dataM;
	dataM = nemMean(data, n, m, r);
	double llr = 0;
	double llrTmp = -1;
	double** phi;
	phi = randPhi(n, n, 1);
	if (type == 1) {
		double maxIdx = pow(2, n * n);
		for (double i = 0; i < maxIdx; i = i + 1) {
			double** phiTmp;
			phiTmp = idxPhi(n, n, i);
			double** P;
			P = matMult(dataM, phiTmp, m, n, n);
			llrTmp = nemLlr(P, n, m);
			for (int j = 0; j < m; j++) {
				delete[] P[j];
			}
			delete[] P;
			P = NULL;
			for (int j = 0; j < n; j++) {
				delete[] phiTmp[j];
			}
			delete[] phiTmp;
			phiTmp = NULL;
			if (llrTmp > llr || i == 0) {
				llr = llrTmp;
				phi = idxPhi(n, n, i);
			}
		}
	}
	else if (type == 2) {
		int iBest = 0;
		int jBest = 0;
		double llrTmp2 = 0;
		llrTmp = 1;
		while(llrTmp > llr) {
			llr = llrTmp;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					double** phiTmp;
					phiTmp = copyPhi(phi, n, n);
					phiTmp[i][j] = 1 - phiTmp[i][j];
					if (phiTmp[j][i] == phiTmp[i][j]) {
						phiTmp[j][i] = 1 - phiTmp[i][j];
					}
					phiTmp = transClose(phiTmp, n, n);
					double** P;
					P = matMult(dataM, phiTmp, m, n, n);
					for (int i = 0; i < n; i++) {
						delete[] phiTmp[i];
					}
					delete[] phiTmp;
					phiTmp = NULL;
					llrTmp2 = nemLlr(P, n, m);
					for (int j = 0; j < m; j++) {
						delete[] P[j];
					}
					delete[] P;
					P = NULL;
					if (llrTmp2 > llrTmp) {
						llrTmp = llrTmp2;
						iBest = i;
						jBest = j;
					}
				}
			}
			if (llrTmp > llr) {
				phi[iBest][jBest] = 1 - phi[iBest][jBest];
				if (phi[jBest][iBest] == phi[iBest][jBest]) {
					phi[jBest][iBest] = 1 - phi[iBest][jBest];
				}
				phi = transClose(phi, n, n);
			}
		}
	}
	for (int i = 0; i < m; i++) {
		delete[] dataM[i];
	}
	delete[] dataM;
	dataM = NULL;
	if (file == "none") {
	}
	else {
		ofstream myfile(file);
		if (myfile.is_open()) {
			myfile << "log-likelihood ratio: " << endl;
			myfile << llr << endl;
			myfile << "phi:" << endl;
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					myfile << phi[j][k] << "\t";
				}
				myfile << endl;
			}
			myfile.close();
		}
	}
	return phi;
}

double*** mnemSim(double*** phi, double*** theta, double* pi, int n, int m, int K = 2, int l = 1000, double mu = 0, double sigma = 1) {
	int** r;
	r = new int* [n];
	int count;
	double** F;
	double*** data;
	double*** noise;
	int* rSum;
	data = new double** [n];
	noise = new double** [n];
	rSum = new int[n];
	for (int i = 0; i < n; i++) {
		r[i] = new int[K];
		rSum[i] = 0;
		for (int j = 0; j < K; j++) {
			r[i][j] = round((pi[j] * l) / n);
			rSum[i] = rSum[i] + r[i][j];
		}
	}
	for (int i = 0; i < n; i++) {
		noise[i] = new double* [m];
		data[i] = new double* [m];
		for (int j = 0; j < m; j++) {
			noise[i][j] = new double[rSum[i]];
			data[i][j] = new double[rSum[i]];
			noise[i][j] = rnorm(rSum[i], mu, sigma);
			count = 0;
			for (int k = 0; k < K; k++) {
				F = matMult(phi[k], theta[k], n, n, m);
				F = transpose(F, n, m);
				for (int s = 0; s < r[i][k]; s++) {
					if (F[j][i] == 0) {
						data[i][j][count] = -1 + mu + noise[i][j][count];
					}
					else {
						data[i][j][count] = 1 - mu + noise[i][j][count];
					}
					count++;
				}
			}
		}
	}
	return data;
}

double*** mnem(double*** data, int n, int m, int* r, int l, int K = 1, int type = 2, double eps = 0.00001, int maxIter = 100, int init = 1, double p = 0.75, string file = "mnem_result.txt", int null = 1) {
	if (null == 1) {
		K++;
	}
	double** post;
	double a;
	double b;
	double postSum = 0;
	post = new double* [l];
	double* llrVec;
	llrVec = new double[l];
	int tmpIdx = 0;
	for (int i = 0; i < l; i++) {
		post[i] = new double[K];
		llrVec[i] = 0;
		a = rand();
		b = rand();
		if (a >= b) {
			a = b / a;
		}
		else {
			a = a / b;
		}
		a = a * K;
		for (int j = 0; j < K; j++) {
			tmpIdx = j + 1;
			if (a >= j && a < tmpIdx) {
				post[i][j] = 1;
			}
			else {
				post[i][j] = 0;
			}
		}
	}
	double* pi;
	pi = new double[K];
	for (int i = 0; i < K; i++) {
		pi[i] = 0;
	}
	double*** phis;
	phis = new double** [K];
	for (int i = 0; i < K; i++) {
		phis[i] = new double* [n];
		for (int j = 0; j < n; j++) {
			phis[i][j] = new double[n];
			for (int k = 0; k < n; k++) {
				phis[i][j][k] = 0;
			}
		}
	}
	double** P;
	P = new double* [m];
	double** dataM;
	dataM = new double* [m];
	for (int j = 0; j < m; j++) {
		P[j] = new double[n];
		dataM[j] = new double[n];
		for (int i = 0; i < n; i++) {
			P[j][i] = 0;
			dataM[j][i] = 0;
		}
	}
	double** F;
	F = new double* [n];
	double** theta;
	theta = new double* [n];
	double*** dataW;
	dataW = new double** [n];
	for (int i = 0; i < n; i++) {
		F[i] = new double[m];
		theta[i] = new double[m];
		dataW[i] = new double* [m];
		for (int j = 0; j < m; j++) {
			F[i][j] = 0;
			theta[i][j] = 0;
			dataW[i][j] = new double[r[i]];
			for (int k = 0; k < r[i]; k++) {
				dataW[i][j][k] = 0;
			}
		}
	}
	int count = 0;
	double tmp = 0;
	double llr = -2;
	double llrTmp = -1;
	int iter = 0;
	while (llrTmp > llr && llrTmp - llr > eps && iter < maxIter) {
		iter++;
		llr = llrTmp;
		for (int i = 0; i < K; i++) {
			pi[i] = 0;
			for (int j = 0; j < l; j++) {
				pi[i] += post[j][i] / l;
			}
			count = 0;
			for (int s = 0; s < n; s++) {
				for (int k = 0; k < r[s]; k++) {
					for (int j = 0; j < m; j++) {
						dataW[s][j][k] = data[s][j][k] * post[count][i];
					}
					count++;
				}
			}
			if (init == 1 && iter == 1) {
				phis[i] = randPhi(n, n, p);
			}
			else {
				phis[i] = nem(dataW, n, m, r, type);
			}
			if (null == 1 && i == K - 1) {
				for (int j = 0; j < n; j++) {
					for (int k = 0; k < n; k++) {
						phis[i][j][k] = 0;
					}
				}
			}
			dataM = nemMean(dataW, n, m, r);
			P = matMult(dataM, phis[i], m, n, n);
			for (int s = 0; s < m; s++) {
				tmp = 0;
				for (int j = 0; j < n; j++) {
					if (P[s][j] > tmp) {
						for (int k = 0; k < n; k++) {
							theta[k][s] = 0;
						}
						theta[j][s] = 1;
						tmp = P[s][j];
					}
				}
			}
			F = matMult(phis[i], theta, n, n, m);
			count = 0;
			for (int j = 0; j < n; j++) {
				for (int s = 0; s < r[j]; s++) {
					tmp = 0;
					post[count][i] = 0;
					for (int k = 0; k < m; k++) {
						tmp += F[j][k] * data[j][k][s];
					}
					post[count][i] = exp(tmp) * pi[i];
					count++;
				}
			}
		}
		llrTmp = 0;
		for (int i = 0; i < l; i++) {
			llrVec[i] = 0;
			for (int j = 0; j < K; j++) {
				llrVec[i] += post[i][j];
			}
			llrVec[i] = log(llrVec[i]);
			llrTmp += llrVec[i];
		}
		for (int i = 0; i < l; i++) {
			postSum = 0;
			for (int j = 0; j < K; j++) {
				postSum += post[i][j];
			}
			for (int j = 0; j < K; j++) {
				post[i][j] = post[i][j] / postSum;
			}
		}
		std::cout << llrTmp << endl;
	}
	ofstream myfile(file);
	if (myfile.is_open()) {
		myfile << "log-likelihood ratio: " << endl;
		myfile << llrTmp << endl;
		for (int i = 0; i < K; i++) {
			if (null == 0 || i < K - 1) {
				myfile << "COMPONENT " << i << endl;
				myfile << "mixture weight:" << endl;
				myfile << pi[i] << endl;
				myfile << "phi:" << endl;
				for (int j = 0; j < n; j++) {
					for (int k = 0; k < n; k++) {
						myfile << phis[i][j][k] << "\t";
					}
					myfile << endl;
				}
			}
		}
		for (int i = 0; i < K; i++) {
			if (null == 0 || i < K - 1) {
				myfile << "COMPONENT " << i << endl;
				myfile << "posterior: " << endl;
				for (int j = 0; j < l; j++) {
					myfile << post[j][i] << "\t";
				}
				myfile << endl;
			}
		}
		myfile.close();
	}
	else cout << "Unable to open file";
	return phis;
}
