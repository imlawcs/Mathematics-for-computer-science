#include <iostream>
using namespace std;

void sumOfFactor(int n) {
	int sum = 0;
	cout << "Factor of " << n << ": ";
	for(int i = 1; i <= n; i++) {
		if(n % i == 0) {
			sum += i;
			cout << i << " ";
		}
	}
	cout << "\nSum of factor: " << sum;
}

int main() {
	int n = 50;
	sumOfFactor(n);
	return 0;
}
