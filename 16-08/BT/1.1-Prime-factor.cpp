#include <iostream>
#include <cmath>
#include <bits/stdc++.h>
using namespace std;

vector<int> p;
vector<int> anpha;

void primeFactors(int n) 
{ 
	cout << "Exercise 01: Prime factorization for " << n << " is: " << n << " = ";

    while (n % 2 == 0) 
    { 
        std::cout << 2 << " "; 
        n = n / 2; 
    } 
    for (int i = 3; i <= std::sqrt(n); i = i + 2) 
    { 
        while (n % i == 0) 
        { 
            std::cout << i << " "; 
            n = n / i; 
        } 
    } 
    if (n > 2) 
        std::cout << n << " "; 
} 

void countOfFactor(int n) {
	cout << "\n- Factors of " << n << ": " ;

	int count = 0;
	for(int i = 1; i < n; i++) {
		if(n % i == 0) {
			count++;
			cout << i << " ";
		}
	}
	cout << n; 
	cout << "\nExercise 02:" << endl;
	cout << "- Number of factor: " << ++count;
}

//void countOfFactor2(int n) {
//	cout << "\nExercise 02 _ Method 2:" << endl;
//	cout << "- Number of factors of " << n << " (using formular): " ;
//	int product = 1;
//	for (int i = 0; i < anpha.size(); i++) {
//		product *= (anpha[i] + 1);
//	}
//	cout << product << endl;
//}

void sumOfFactor(int n) {
	cout << "\nExercise 03:" << endl;
	cout << "- Sum of the factors of " << n <<" is: " ;
	
	int sum = 0;
	for(int i = 1; i <= n; i++) {
		if(n % i == 0) {
			sum += i;
//			cout << i << " ";
		}
	}
	cout << sum;
}

//void sumOfFactor2 (int n) {
//	cout << "\nExercise 03 _ Method 2:" << endl;
//	cout << "- Sum of the factors of " << n << " is (using formular): " ;
//	long long sum = 1;
//	for (int i = 0; i < anpha.size(); i++) {
//		sum *= 1ll*((pow(p[i], anpha[i] + 1) - 1) / (p[i] - 1));
//	}
//	cout << sum << endl;
//}

void productOfFactor(int n) {
	cout << "\nExercise 04:" << endl;
	cout << "- Product of the factors of " << n << " is: " ;
	
	int pro = 1;
	for(int i = 1; i <= n; i++) {
		if(n % i == 0) {
			pro *= i;
//			cout << i << " ";
		}
	}

	cout << pro;
}

//void productOfFactor2 (int n) {
//	cout << "\nExercise 04 _ Method 2:" << endl;
//	cout << "- Product of the factors of " << n << " is (using formular): " ;
//	long long t = 1;
//	for (int i = 0; i < anpha.size(); i++) {
//		t *= (anpha[i] + 1);
//	}
//	long long product = pow (n, t/2);
//	cout << product << endl;
//}

void perfectNumber(int n) {
	int sum = 0;
	cout << "\nExercise 05:" << endl;
	cout << "- Is " << n << " a perfect number? " << endl ;

//	cout << "Factor of " << n << ": ";
	for(int i = 1; i <= n/2; i++) {
		if(n % i == 0) {
			sum += i;
//			cout << i << " ";
		}
	}
//	cout << endl;
	if(sum == n) {
		cout << n << " is perfect number.";
	}
	else cout << n << " is not perfect number!";
}

//void perfectNumber2 (int n) {
//	cout << "\nExercise 05 _ Method 2:" << endl;
//	cout << "- Is " << n << " a perfect number (using formular)? " ;
//	long long sum = 1;
//	for (int i = 0; i < anpha.size(); i++) {
//		sum *= 1ll*((pow(p[i], anpha[i] + 1) - 1) / (p[i] - 1));
//	}
//	string result = (sum == 2*n) ? "Yes" : "No";
//	cout << result << endl;
//}

void densityOfPrimes(int n) {
	int count = 0;
	count = n / log(n);
	cout << "\nExercise 06: Density of primes from 1 to " << n << ": ";
	cout << count;
}

void SieveOfEratosthenes(int n)
{
	cout << "\nExercise 07 \n- Part 1: Sieve of Eratosthenes: ";
    cout << endl;
    cout << "Following are the prime numbers smaller "
		<< "than or equal to " << n << endl;
	bool prime[n + 1];
	memset(prime, true, sizeof(prime));

	for (int p = 2; p * p <= n; p++) {
		if (prime[p] == true) {
			for (int i = p * p; i <= n; i += p)
				prime[i] = false;
		}
	}

	for (int p = 2; p <= n; p++)
		if (prime[p])
			cout << p << " ";
}

bool prime (int n) {
	if (n < 2) return false;
	for (int i = 2; i <= sqrt(n); i++) {
		if (n%i == 0) return false;
	}
	return true;
}

void SieveOfEratosthenes2 (int n) {
	cout << "\n- Part 2:";
	int a;
	cout << "\nEnter another positive integer: "; cin >> a;
	cout << "=> Primes from " << a << " to " << n << " are: ";
	for (int i = a; i <= n; i++) {
		if (prime(i)) cout << i << " ";
	}
}

int main() 
{ 
    int n = 315; 
    primeFactors(n); 
    countOfFactor(n);
//    countOfFactor2(n);
    sumOfFactor(n);
//    sumOfFactor2(n);
    productOfFactor(n);
//    productOfFactor2(n);
    perfectNumber(n);
//    perfectNumber2(n);
    densityOfPrimes(n);
	SieveOfEratosthenes(n);
	SieveOfEratosthenes2(n);
    return 0; 
}

