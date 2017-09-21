#define _USE_MATH_DEFINES
#include <thread>
#include <tuple>
#include <iostream>
#include <random>
#include <cmath>

double estimate_pi(int nsamples) {
	std::default_random_engine rnd(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_real_distribution<double> dist(-1.0, 1.0);

	double numHit = 0;
	double x, y;
	// generate 100 random doubles
	for (int i = 0; i<nsamples; i++) {
		x = dist(rnd);
		y = dist(rnd);

		//std::cout << x*x + y*y << std::endl;
		
		if (x*x + y*y <= 1.0) {
			numHit++;
		}
	}

	return numHit / nsamples * 4;
}

// generates a random sample and sets the value of `inside`
// to true if within the unit circle
void pi_sampler(std::vector<bool>& hits, int idx) {

	// single instance of random engine and distribution
	static std::default_random_engine rnd;
	static std::uniform_real_distribution<double> dist(-1.0, 1.0);

	// YOUR CODE HERE
	double x, y;
	x = dist(rnd);
	y = dist(rnd);
	hits[idx] = (x*x + y*y <= 1.0);
}


// naively uses multithreading to try to speed up computations
double estimate_pi_multithread_naive(int nsamples) {
	// stores whether each sample falls within circle
	std::vector<bool> hits(nsamples, false);

	// create and store all threads
	std::vector<std::thread> threads;
	for (int i = 0; i<nsamples; ++i) {
		threads.push_back(std::thread(pi_sampler, std::ref(hits), i));
	}

	// wait for all threads to complete
	for (int i = 0; i<nsamples; ++i) {
		threads[i].join();
	}

	// estimate our value of PI
	double pi = 0;
	for (int i = 0; i<nsamples; ++i) {
		std::cout << hits[i] << std::endl;
		if (hits[i]) {
			pi = pi + 1;
		}
	}
	pi = pi / nsamples * 4;

	return pi;
}

 //count number of hits using nsamples, populates hits[idx]
void pi_hits(std::vector<int>& hits, int idx, int nsamples) {

	// single instance of random engine and distribution
	static std::default_random_engine rnd;
	static std::uniform_real_distribution<double> dist(-1.0, 1.0);

	// YOUR CODE HERE
	double x;
	double y;
	for (int i = 0; i < nsamples; i++) {
		x = dist(rnd);
		y = dist(rnd);

		if (x*x + y*y <= 1.0) {
			hits[idx]++;
		}
	}
}

// divides work among threads intelligently
double estimate_pi_multithread(int nsamples) {

	// number of available cores
	int nthreads = std::thread::hardware_concurrency();
	std::cout << nthreads << std::endl;

	// hit counts
	std::vector<int> hits(nthreads, 0);

	// create and store threads
	std::vector<std::thread> threads;
	int msamples = 0; // samples accounted for
	for (int i = 0; i<nthreads - 1; ++i) {
		threads.push_back(
			std::thread(pi_hits, std::ref(hits), i, nsamples / nthreads));
		msamples += nsamples / nthreads;
	}
	// remaining samples
	threads.push_back(
		std::thread(pi_hits, std::ref(hits), nthreads - 1, nsamples - msamples));

	// wait for threads to finish
	for (int i = 0; i<nthreads; ++i) {
		threads[i].join();
	}

	// estimate pi
	double pi = 0;
	for (int i = 0; i<nthreads; ++i) {
		pi += hits[i];
	}
	pi = pi / nsamples * 4;

	return pi;
}

double SumOfSquare(double x, double y, double z) {
	return pow(x, 2) + pow(y, 2) + pow(z, 2);
}

//double estimateIntegral(double minX, double maxX, double(*func)(double, double, double), int nSample) {
//
//	std::default_random_engine rnd(std::chrono::system_clock::now().time_since_epoch().count());
//	static std::uniform_real_distribution<double> dist(minX, maxX);
//
//	double tempSum = 0;
//	double Vol = maxX - minX;
//	double div = Vol / (double)nSample;
//	double x = 1, y = 1, z = 1;
//
//
//	for (int i = 0; i < nSample; i++) {
//		do {
//			x = dist(rnd);
//			y = dist(rnd);
//			z = dist(rnd);
//		} while (sqrt(SumOfSquare(x, y, z))>1);
//
//
//		tempSum += Vol*(*func)(x, y, z);
//	}
//
//	double avgIntegral = tempSum / (double)nSample;
//
//	return avgIntegral;
//}

std::tuple<double, double, double> estimateIntegral(double radius, double(*func)(double, double, double), int nSample) {

	std::default_random_engine rnd(std::chrono::system_clock::now().time_since_epoch().count());
	static std::uniform_real_distribution<double> dist(-radius, radius);
	
	double Vol = 4.0 / 3.0*M_PI*pow(radius, 3);

	std::vector<double> tempSum(3);
	std::vector<double> x(3);
	std::vector<double> mass(3);


	for (int i = 0; i < nSample; i++) {
		do {
			x[0] = dist(rnd);
			x[1] = dist(rnd);
			x[2] = dist(rnd);
		} while (sqrt(SumOfSquare(x[0], x[1], x[2]))>radius);

		for (int j = 0; j < x.size(); j++) {
			tempSum[j] += Vol*(*func)(x[0], x[1], x[2]);
		}
	}

	for (int i = 0; i < x.size(); i++) {
		mass[i] = tempSum[i] / (double)nSample;
	}

	return std::make_tuple(mass[0], mass[1], mass[2]);
}

std::tuple<double, double, double> estimateMoment(double radius, double(*func)(double, double, double), int nSample) {

	std::default_random_engine rnd(std::chrono::system_clock::now().time_since_epoch().count());
	static std::uniform_real_distribution<double> dist(-radius, radius);

	double Vol = 4.0 / 3.0*M_PI*pow(radius, 3);

	std::vector<double> tempSum(3);	
	double pointMom;
	std::vector<double> x(3);
	std::vector<double> moment(3);

	for (int i = 0; i < nSample; i++) {
		do {
			x[0]= dist(rnd);
			x[1] = dist(rnd);
			x[2] = dist(rnd);
		} while (sqrt(SumOfSquare(x[0], x[1], x[2]))>radius);

		for (int j = 0; j < x.size(); j++) {
			pointMom = x[j]*(*func)(x[0], x[1], x[2]);
			tempSum[j] += Vol*pointMom;
		}
	}

	for (int i = 0; i < x.size(); i++) {
		moment[i] = tempSum[i] / (double)nSample;
	}

	return std::make_tuple(moment[0], moment[1], moment[2]);
}

std::tuple<double, double, double> calcCenterOfMass(double radius, double(*func)(double, double, double), int nSample) {
	double cx, cy, cz;
	double tempMoment, tempMass;
	//double Vol = 4.0 / 3.0*M_PI*pow(radius, 3);

	auto moments = estimateMoment(radius, (*func), nSample);
	auto mass = estimateIntegral(radius, (*func), nSample);

	tempMoment = std::get<0>(moments);
	tempMass = std::get<0>(mass);
	cx = tempMoment / tempMass;

	tempMoment = std::get<1>(moments);
	tempMass = std::get<1>(mass);
	cy = tempMoment / tempMass;

	tempMoment = std::get<2>(moments);
	tempMass = std::get<2>(mass);
	cz = tempMoment / tempMass;

	return std::make_tuple(cx, cy, cz);
}

double DensityFunction0(double x, double y, double z) {
	if (z > 0)
		return 1.0;
	else
		return 0.01;
}

double DensityFunction1(double x, double y, double z) {
	return exp(-(pow(abs(x), 2)));

}

double DensityFunction2(double x, double y, double z) {
	return abs(x + y + z);
}

double DensityFunction3(double x, double y, double z) {
	return pow(x - 1, 2) + pow(y - 2, 2) + pow(z - 3, 2);
}

void estimateCenterMultithread(double radius, std::vector<std::vector<double>>& mSum, 
								std::vector<std::vector<double>>& iSum, 
								int idx, double(*func)(double, double, double), int nSample) {
	auto moments = estimateMoment(radius, (*func), nSample);
	auto inertia = estimateIntegral(radius, (*func), nSample);

	mSum[0][idx] = std::get<0>(moments)*nSample;
	mSum[1][idx] = std::get<1>(moments)*nSample;
	mSum[2][idx] = std::get<2>(moments)*nSample;
	
	iSum[0][idx] = std::get<0>(inertia)*nSample;
	iSum[1][idx] = std::get<1>(inertia)*nSample;
	iSum[2][idx] = std::get<2>(inertia)*nSample;
}

std::tuple<double, double, double> calcCenterOfMassMultithread(double radius, double(*func)(double, double, double), int nSamples) {
	std::vector<double> centers(3);
	std::vector<double> moments(3,0);
	std::vector<double> inertia(3,0);
	double tempMoment, tempInertia;
	// number of available cores
	int nthreads = std::thread::hardware_concurrency();
	std::vector<std::vector<double>> mSum(3, std::vector<double>(nthreads, 0));
	std::vector<std::vector<double>> iSum(3, std::vector<double>(nthreads, 0));

	std::vector<std::thread> threads;
	int msamples = 0; // samples accounted for
	for (int i = 0; i<nthreads - 1; ++i) {
		threads.push_back(std::thread(estimateCenterMultithread, radius, std::ref(mSum), std::ref(iSum), i, (*func), nSamples/nthreads));
		msamples += nSamples / nthreads;
	}
	// remaining samples
	threads.push_back(std::thread(estimateCenterMultithread, radius, std::ref(mSum), std::ref(iSum), nthreads-1, (*func), nSamples - msamples));

	// wait for threads to finish
	for (int i = 0; i<nthreads; ++i) {
		threads[i].join();
	}

	for (int i = 0; i < nthreads; i++) {
		for (int j = 0; j < 3; j++) {
			moments[j] += mSum[j][i];
			inertia[j] += iSum[j][i];
		}
	}

	for (int i = 0; i < 3; i++) {
		tempMoment = moments[i] / nSamples;
		tempInertia = inertia[i] / nSamples;
		centers[i] = tempMoment / tempInertia;
	}

	return std::make_tuple(centers[0], centers[1], centers[2]); // 0 is x, 1 is y, 2 is z
}

int main() {
	/*double pi = estimate_pi_multithread(10000000);
	std::cout << "My estimate of PI is: " << pi << std::endl;*/

	/*double intgrl = estimateIntegral(0, 1, [](double x) {return x; }, 10000);
	std::cout << intgrl << std::endl;*/
	
	//auto com = calcCenterOfMass(1, DensityFunction2, 10000);

	auto com = calcCenterOfMassMultithread(1, DensityFunction0, 1000000);

	std::cout << std::fixed << std::get<0>(com) << std::endl;
	std::cout << std::fixed << std::get<1>(com) << std::endl;
	std::cout << std::fixed << std::get<2>(com) << std::endl;

	system("pause");
	return 0;
}