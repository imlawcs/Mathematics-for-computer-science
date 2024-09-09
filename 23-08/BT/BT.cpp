#include <bits/stdc++.h>
using namespace std;

int Cau1_UCLN (int a, int b) {
	if (b == 0) return a;
	return Cau1_UCLN(b, a%b);
}

int Cau1_BCNN (int a, int b) {
	return a * b / Cau1_UCLN(a, b);
}

bool soNguyenTo(int n) {
    for(int i = 2; i <= n/2; i++) {
        if(n % i == 0) {
            return false;
        }
    }
    return true;
}

void Cau2() {
	int n;
	cout << "Nhap n tu ban phim: ";
	cin >> n;
	if(soNguyenTo(n)) {
		cout << "So luong cac so nguyen to cung nhau voi " << n << " la: " << n - 1 << endl;
		for(int i = 1; i <= n - 1; i++) cout << i << " ";
		cout << endl;
		return;
	}
	else {
		int count = 0;
		cout << "Cac so nguyen to cung nhau voi " << n << " la: " << endl;
		for(int i = 1; i <= n-1; i++) {
			if(Cau1_UCLN(i, n) == 1) {
				count++;
				cout << i << " ";
			}
		}
		cout << endl << "So luong cac so nguyen to cung nhau voi " << n << " la: " << count << endl;	
	} 
	
	map<int, int> m;
	cout << "So luong cac so nguyen to cung nhau voi " << n << " (Cach 2) la: ";
    for(int i = 2; i <= n; i++) {
        while(n % i == 0) {
            m[i]++;
            n /= i;
        }
    }
    int sl = 1;
    for(auto i : m) {
    	sl *= pow(i.first, i.second - 1) * (i.first - 1);
    }
    cout << sl << endl;
}

void modular1(int x, int y, int m) {
	string result;
	bool check;
	
	check = ((x + y) % m == ((x%m + y%m)) % m);
	result = (check == true) ? "Dung" : "Sai";
	cout << "\n   (" << x << " + " << y << ") mod " << m << " = (" << x << " mod " << m << " + " << y << " mod " << m << ") mod " << m << " => " << result;
}

void modular2(int x, int y, int m) {
	string result;
	bool check;
	
	check = ((x - y) % m == ((x%m - y%m)) % m);
	result = (check == true) ? "Dung" : "Sai";
	cout << "\n   (" << x << " - " << y << ") mod " << m << " = (" << x << " mod " << m << " - " << y << " mod " << m << ") mod " << m << " => " << result;
}

void modular3(int x, int y, int m) {
	string result;
	bool check; 
	
	check = ((x * y) % m == ((x%m * y%m)) % m);
	result = (check == true) ? "Dung" : "Sai";
	cout << "\n   (" << x << " * " << y << ") mod " << m << " = (" << x << " mod " << m << " * " << y << " mod " << m << ") mod " << m << " => " << result << endl;
}

void modular4(int x, int y, int m){
	string result;
	bool check; 
	
	check = ((int)pow(x, y) % m == ((int)pow(x%m, y)) % m);
	result = (check == true) ? "Dung" : "Sai";
	cout << "   " << x << "^" << y << " mod " << m << " = (" << x << " mod " << m << ")^" << y << " mod " << m << " => " << result;
}

void Cau3(){
	int x, y, m;
	cout << "+ Nhap gia tri cua x: x = "; cin >> x;
	cout << "+ Nhap gia tri cua y: y = "; cin >> y;
	cout << "+ Nhap gia tri cua m: m = "; cin >> m;
	modular1(x, y, m);
	modular2(x, y, m);
	modular3(x, y, m);
	modular4(x, y, m);
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
			if (Cau1_UCLN(temp, arrayM[j]) == 1) continue;
			else return false;
		}
	}
	return true;
}

void Cau5() {
	vector<int> res_nhan;
	vector<int> mod;
	if (checkCondition()) {
		int M = 1;
		int final = 0;
		for (int i = 0; i < arrayM.size(); i++) {
			int mul = 1;
			for (int j = 0; j < arrayM.size(); j++) {
				if (i != j) mul *= arrayM[j];
			}
			res_nhan.push_back(mul);
			mod.push_back(Cau4_ModuleInverse(mul, arrayM[i]));
			final += (arrayA[i]*mul*mod[i]);
			M *= arrayM[i];
		}
		cout << "Mot nghiem cua x la: x = " << final << endl;
		cout << "Ho nghiem cua x la:  x = " << (final%M) << " + k*" << M << " (k thuoc N)" << endl; 
	}
}

int main(){
	int a, b;
	cout << "(*)Cau 1: Nhap 2 so nguyen duong a va b de tim UCLN, BCNN:" << endl;
	cout << "   a = "; cin >> a;
	cout << "   b = "; cin >> b;
	
	cout << "UCLN (" << a << ", " << b << ") = " << Cau1_UCLN(a, b) << endl;
	cout << "BCNN (" << a << ", " << b << ") = " << Cau1_BCNN(a, b);
	
	cout << "\n\n(*)Cau 2: Nhap 1 so nguyen duong n de tim so nguyen to cung nhau:" << endl;
	Cau2();
	
	cout << "\n\n(*)Cau 3: Kiem tra tinh dung dan cua cac bieu thuc:" << endl;
	Cau3();
	
	cout << "\n\n(*)Cau 4: Tinh module dao nguoc:\n";
	int x, m;
	cout << "   Nhap x: x = "; cin >> x;
	cout << "   Nhap m: m = "; cin >> m;
	cout << "-> Module dao nguoc: " << x << "*x^(-1) mod " << m << " = 1 -> x^(-1) = ";
	cout << Cau4_ModuleInverse(x, m);
	
    cout << "\n\n(*)Cau 5: Giai he phuong trinh dong du:\n";
    init();
	Cau5();
}

