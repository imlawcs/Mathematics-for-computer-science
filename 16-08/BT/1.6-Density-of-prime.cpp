#include <iostream>
#include <cmath>
using namespace std;

void densityOfPrimes(int n) {
	int count = 0;
	count = n / log(n);
	cout << "Number of primes between 1 and " << n << ": " << count;
}

int main() {
	int n = 50;
	densityOfPrimes(n);
	return 0;
}
