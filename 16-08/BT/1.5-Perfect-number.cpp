#include <iostream>
using namespace std;

void perfectNumber(int n) {
	int sum = 0;
	cout << "Factor of " << n << ": ";
	for(int i = 1; i <= n/2; i++) {
		if(n % i == 0) {
			sum += i;
			cout << i << " ";
		}
	}
	cout << endl;
	if(sum == n) {
		cout << n << " is perfect number.";
	}
	else cout << n << " is not perfect number!";
}

int main() {
	int n = 50;
	perfectNumber(n);
	return 0;
}
