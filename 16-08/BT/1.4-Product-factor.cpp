#include <iostream>
using namespace std;

void productOfFactor(int n) {
	int pro = 1;
	cout << "Factor of " << n << ": ";
	for(int i = 1; i <= n; i++) {
		if(n % i == 0) {
			pro *= i;
			cout << i << " ";
		}
	}
	cout << "\nProduct of factor: " << pro;
}

int main() {
	int n = 50;
	productOfFactor(n);
	return 0;
}
