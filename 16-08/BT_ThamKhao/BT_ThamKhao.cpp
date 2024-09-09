#include "bits/stdc++.h"
using namespace std;

vector<int> p;
vector<int> anpha;

bool prime (int n) {
	if (n < 2) return false;
	for (int i = 2; i <= sqrt(n); i++) {
		if (n%i == 0) return false;
	}
	return true;
}

vector<int> factors (int n) {
	vector<int> f;
	for (int i = 2; i < sqrt(n); i++) {
		while (n%i == 0) {
			f.push_back(i);
			n /= i;
		}
	}
	if (n > 1) f.push_back (n);
	return f;
}

void countFreq (vector<int> arr, int n) {
    vector<bool> visited(n, false);
    for (int i = 0; i < n; i++) {
        if (visited[i] == true)
            continue;
        int count = 1;
        for (int j = i + 1; j < n; j++) {
            if (arr[i] == arr[j]) {
                visited[j] = true;
                count++;
            }
        }
        p.push_back(arr[i]);
        anpha.push_back(count);
    }
}

void print_Bai1 (int n) {
	cout << "\nExercise 01: Prime factorization for " << n << " is: " << n << " = ";
	vector<int> printVector = factors(n);
	countFreq (printVector, printVector.size());
	for (int i = 0; i < p.size(); i++) {
		if (i != p.size() - 1)
			cout << p[i] << "^" << anpha[i] << " * ";
		else 
			cout << p[i] << "^" << anpha[i] << endl;
	}
}

void print_Bai2_C1 (int n) {
	cout << "\nExercise 02 _ Method 1:" << endl;
	cout << "- Factors of " << n << ": " ;
	int count = 0;
	for (int i = 1; i <= n; i++) {
		if (n%i == 0) {
			count++;
			cout << i << " ";
		}
	}
	cout << "\n- Number of factors of " << n << ": " << count << endl;
}

void print_Bai2_C2 (int n) {
	cout << "\nExercise 02 _ Method 2:" << endl;
	cout << "- Number of factors of " << n << " (using formular): " ;
	int product = 1;
	for (int i = 0; i < anpha.size(); i++) {
		product *= (anpha[i] + 1);
	}
	cout << product << endl;
}

void print_Bai3_C1 (int n) {
	cout << "\nExercise 03 _ Method 1:" << endl;
	cout << "- Sum of the factors of " << n <<" is: " ;
	int sum = 0;
	for (int i = 1; i <= n; i++) {
		if (n%i == 0) {
			sum += i;
		}
	}
	cout << sum << endl;
}

void print_Bai3_C2 (int n) {
	cout << "\nExercise 03 _ Method 2:" << endl;
	cout << "- Sum of the factors of " << n << " is (using formular): " ;
	long long sum = 1;
	for (int i = 0; i < anpha.size(); i++) {
		sum *= 1ll*((pow(p[i], anpha[i] + 1) - 1) / (p[i] - 1));
	}
	cout << sum << endl;
}

void print_Bai4_C1 (int n) {
	cout << "\nExercise 04 _ Method 1:" << endl;
	cout << "- Product of the factors of " << n << " is: " ;
	long long product = 1;
	for (int i = 1; i <= n; i++) {
		if (n%i == 0) {
			product *= i;
		}
	}
	cout << product << endl;
}

void print_Bai4_C2 (int n) {
	cout << "\nExercise 04 _ Method 2:" << endl;
	cout << "- Product of the factors of " << n << " is (using formular): " ;
	long long t = 1;
	for (int i = 0; i < anpha.size(); i++) {
		t *= (anpha[i] + 1);
	}
	long long product = pow (n, t/2);
	cout << product << endl;
}

void print_Bai5_C1 (int n) {
	cout << "\nExercise 05 _ Method 1:" << endl;
	cout << "- Is " << n << " a perfect number? " ;
	int sum = 0;
	for (int i = 1; i <= n/2; i++) {
		if (n%i == 0) {
			sum += i;
		}
	}
	string result = (sum == n) ? "Yes" : "No";
	cout << result << endl;
}

void print_Bai5_C2 (int n) {
	cout << "\nExercise 05 _ Method 2:" << endl;
	cout << "- Is " << n << " a perfect number (using formular)? " ;
	long long sum = 1;
	for (int i = 0; i < anpha.size(); i++) {
		sum *= 1ll*((pow(p[i], anpha[i] + 1) - 1) / (p[i] - 1));
	}
	string result = (sum == 2*n) ? "Yes" : "No";
	cout << result << endl;
}

void print_Bai6 (int n) {
	cout << "\nExercise 06: Density of primes from 1 to " << n << ": " << n/log(n) << " ~ " << ceil(n/log(n)) << endl;
}

void print_Bai7_C1 (int n) {
	cout << "\nExercise 07 \n- Part 1: Sieve of Eratosthenes: ";
	int sieve[n + 1];
	for (int x = 0; x <= n; x++) {
		sieve[x] = 0;
	}
	for (int x = 2; x <= n; x++) {
		if (sieve[x]) continue;
		for (int u = 2*x; u <= n; u += x) {
			sieve[u] = x;
		}
	}
	for (int i = 2; i <= n; i++) {
		if (sieve[i] == 0) {
			cout << i << "  ";
		} 
	}
}

void print_Bai7_C2 (int n) {
	cout << "- Part 2:";
	int a;
	cout << "\n  Enter another positive integer: "; cin >> a;
	cout << "  => Primes from " << a << " to " << n << " are: ";
	for (int i = a; i <= n; i++) {
		if (prime(i)) cout << i << "  ";
	}
}

int main() {
	int n;
	cout << "Enter a positive integer: " ;
	cin >> n;
	
	print_Bai1(n);
	
	print_Bai2_C1(n);
	print_Bai2_C2(n);
	
	print_Bai3_C1(n);
	print_Bai3_C2(n);
	
	print_Bai4_C1(n);
	print_Bai4_C2(n);
	
	print_Bai5_C1(n);
	print_Bai5_C2(n);
	
	print_Bai6(n);
	
	print_Bai7_C1(n);
	print_Bai7_C2(n);
}

