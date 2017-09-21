#include <iostream>
#include <random>
#include <thread>

int main() {
	// specify the engine and distribution.
	std::default_random_engine rnd(	
		std::chrono::system_clock::now().time_since_epoch().count());	
	std::uniform_real_distribution<double> dist(-1.0, 1.0);

	// generate 100 random doubles
	for (int i = 0; i<100; ++i) {
		std::cout << dist(rnd) << std::endl;
	}

	system("pause");
	return 0;
}