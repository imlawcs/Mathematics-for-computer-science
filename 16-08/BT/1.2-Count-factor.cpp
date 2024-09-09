#include <iostream>
using namespace std;

void countOfFactor(int n) {
	int count = 0;
	cout << "Factor of " << n << ": ";
	for(int i = 1; i < n; i++) {
		if(n % i == 0) {
			count++;
			cout << i << " ";
		}
	}
	cout << n; 
	cout << "\nNumber of factor: " << ++count;
}

int main() {
	int n = 50;
	countOfFactor(n);
	return 0;
}
