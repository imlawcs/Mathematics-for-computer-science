#include "bits/stdc++.h"
using namespace std;

vector<int> p;
vector<int> anpha;

int Cau1_UCLN_C1(int a, int b) {
	while ( a != b) {
        if (a > b)
            a = a - b;
        else
            b = b - a;
    }
	return a;
}

int Cau1_UCLN_C2 (int a, int b) {
	if (b == 0) return a;
	return Cau1_UCLN_C2(b, a%b);
}

int Cau1_BCNN (int a, int b) {
	return a * b / Cau1_UCLN_C1(a, b);
}

// ham phan tich TSNT
vector<int> factors (int n) {
	vector<int> f;
	for (int i = 2; i <= sqrt(n); i++) {
		while (n%i == 0) {
			f.push_back(i);
			n /= i;
		}
	}
	if (n > 1) f.push_back (n);
	return f;
}

// Thiet lap vecto p, anpha
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

void Cau2_Coprime () {
	int n;
	cin >> n;
	vector<int> printVector = factors(n);
	countFreq (printVector, printVector.size());
	int count = 1;
	for (int i = 0; i < p.size(); i++) {
		count *= (pow(p[i], (anpha[i] - 1)) * (p[i] - 1));
	}
	cout << "+ So luong cac so nguyen to cung nhau cua " << n << ": " << count << endl;
	cout << "+ Liet ke: ";
	for (int i = 1; i < n; i++) {
		if (Cau1_UCLN_C1(i, n) == 1) {
			if (i == n - 1) cout << i;
			else cout << i << ", ";
		}
	}
}

void Cau3_CheckFormula () {
	int x, y, m;
	cout << "+ Nhap gia tri cua x: x = "; cin >> x;
	cout << "+ Nhap gia tri cua y: y = "; cin >> y;
	cout << "+ Nhap gia tri cua m: m = "; cin >> m;
	
	string result;
	bool check;
	
	check = ((x + y) % m == ((x%m + y%m)) % m);
	result = (check == true) ? "Dung" : "Sai";
	cout << "\n   (" << x << " + " << y << ") mod " << m << " = (" << x << " mod " << m << " + " << y << " mod " << m << ") mod " << m << " => " << result;
	
	check = ((x - y) % m == ((x%m - y%m)) % m);
	result = (check == true) ? "Dung" : "Sai";
	cout << "\n   (" << x << " - " << y << ") mod " << m << " = (" << x << " mod " << m << " - " << y << " mod " << m << ") mod " << m << " => " << result;
	
	check = ((x * y) % m == ((x%m * y%m)) % m);
	result = (check == true) ? "Dung" : "Sai";
	cout << "\n   (" << x << " * " << y << ") mod " << m << " = (" << x << " mod " << m << " * " << y << " mod " << m << ") mod " << m << " => " << result << endl;
	
	check = ((int)pow(x, y) % m == ((int)pow(x%m, y)) % m);
	result = (check == true) ? "Dung" : "Sai";
	cout << "   " << x << "^" << y << " mod " << m << " = (" << x << " mod " << m << ")^" << y << " mod " << m << " => " << result;
}

bool checkBack (int y, int x, int k, int m) {
	if (k * m + 1 == x * y) return true;
	else return false;
}

int Cau4_ModuleInverse(int x, int m) {
	int count = 0;
	while (true) {
		int temp = (int)(1.0 * (m*count + 1)/x);
		if (checkBack(temp, x, count, m)) {
			return temp;
		}
		else count++;
	}
}

vector<int> arrayA;
vector<int> arrayM;

void init() {
	int level;
	cout << "   Nhap so phuong trinh ban muon giai: so phuong trinh = "; cin >> level;
	int count = 0;
	while (count < level) {
		int a1, m1;
		cout << "- Nhap he so phuong trinh " << (count + 1) << ":" << endl;
		cout << "  a" << (count + 1) << " = "; cin >> a1;
		cout << "  m" << (count + 1) << " = "; cin >> m1;
		arrayA.push_back(a1);
		arrayM.push_back(m1);
		count++;
	}
}

bool checkCondition() {
	for (int i = 0; i < arrayM.size(); i++) {
		int temp = arrayM[i];
		for (int j = 0; j < arrayM.size(); j++) {
			if (i == j) continue;
			if (Cau1_UCLN_C1(temp, arrayM[j]) == 1) continue;
			else return false;
		}
	}
	return true;
}

void Cau5_ChineseRemainder() {
	if (checkCondition()) {
		int M = 1;
		int final = 0;
		for (int i = 0; i < arrayM.size(); i++) {
			int mul = 1;
			for (int j = 0; j < arrayM.size(); j++) {
				if (i != j) mul *= arrayM[j];
			}
			int modInverse = Cau4_ModuleInverse(mul, arrayM[i]);
            final += (arrayA[i] * mul * modInverse);
			M *= arrayM[i];
		}
		cout << "Mot nghiem cua x la: x = " << final << endl;
		cout << "Ho nghiem cua x la:  x = " << (final%M) << " + k * " << M << " (k thuoc N)" << endl; 
	}
	else {
		cout << "He phuong trinh khong hop le!";
	}
}

int main() {
	
	int a, b;
	cout << "(*)Cau 1: Nhap 2 so nguyen duong a va b de tim UCLN, BCNN:" << endl;
	cout << "   a = "; cin >> a;
	cout << "   b = "; cin >> b;
	
	cout << "UCLN (" << a << ", " << b << ") = " << Cau1_UCLN_C1(a, b) << " (cach 1)\n";
	cout << "UCLN (" << a << ", " << b << ") = " << Cau1_UCLN_C2(a, b) << " (cach 2)\n";
	cout << "BCNN (" << a << ", " << b << ") = " << Cau1_BCNN(a, b);
	
	cout << "\n\n(*)Cau 2: Nhap 1 so nguyen duong n de tim so nguyen to cung nhau:" << endl;
	cout << "   n = ";
	Cau2_Coprime();
	
	cout << "\n\n(*)Cau 3: Kiem tra tinh dung dan cua cac bieu thuc:" << endl;
	Cau3_CheckFormula();
	
	cout << "\n\n(*)Cau 4: Tinh module dao nguoc:\n";
	int x, m;
	cout << "   Nhap x: x = "; cin >> x;
	cout << "   Nhap m: m = "; cin >> m;
	cout << "-> Module dao nguoc: " << x << "*x^(-1) mod " << m << " = 1 -> x^(-1) = ";
	cout << Cau4_ModuleInverse(x, m);
	
	cout << "\n\n(*)Cau 5: Giai he phuong trinh dong du:\n";
	init();
	Cau5_ChineseRemainder();
}

