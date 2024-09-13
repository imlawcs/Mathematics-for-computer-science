//#include <iostream>
//#include <cmath>
//#include <algorithm>
//
//using namespace std;
//
//// Hàm tính ð?nh th?c c?a ma tr?n n x n
//double determinant(int n, float a[][100]) {
//    double det = 1.0;
//    for (int i = 0; i < n; i++) {
//        // T?m hàng v?i ph?n t? l?n nh?t ? c?t hi?n t?i
//        int pivot = i;
//        for (int j = i + 1; j < n; j++) {
//            if (abs(a[j][i]) > abs(a[pivot][i])) {
//                pivot = j;
//            }
//        }
//
//        // Hoán ð?i hàng ð? ðýa hàng v?i ph?n t? l?n nh?t lên trên
//        if (pivot != i) {
//            swap(a[i], a[pivot]);
//            det *= -1;  // Thay ð?i d?u c?a ð?nh th?c khi hoán ð?i hàng
//        }
//
//        // N?u ph?n t? chéo chính b?ng 0 th? ð?nh th?c là 0
//        if (a[i][i] == 0) {
//            return 0;
//        }
//
//        // Nhân ð?nh th?c v?i ph?n t? chéo chính
//        det *= a[i][i];
//
//        // Th?c hi?n lo?i Gauss ð? ðýa ma tr?n v? d?ng tam giác trên
//        for (int j = i + 1; j < n; j++) {
//            double factor = a[j][i] / a[i][i];
//            for (int k = i + 1; k < n; k++) {
//                a[j][k] -= factor * a[i][k];
//            }
//        }
//    }
//    return det;
//}
//
//// Ví d? s? d?ng
//int main() {
//    int n = 3;
//    float matrix[100][100] = {
//        {2, 1, 0},
//        {1, 3, 1},
//        {0, 1, 2}
//    };
//
//    double det = determinant(n, matrix);
//    cout << "Determinant: " << det << endl;
//
//    return 0;
//}

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

// Hàm tính ð?nh th?c c?a ma tr?n vuông
double determinant(const vector<vector<double>>& matrix, int n) {
    if (n == 1) {
        return matrix[0][0];
    }

    if (n == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }

    double det = 0;
    vector<vector<double>> submatrix(n - 1, vector<double>(n - 1));
    int sign = 1;

    for (int x = 0; x < n; x++) {
        int subi = 0;
        for (int i = 1; i < n; i++) {
            int subj = 0;
            for (int j = 0; j < n; j++) {
                if (j == x) continue;
                submatrix[subi][subj] = matrix[i][j];
                subj++;
            }
            subi++;
        }
        det += sign * matrix[0][x] * determinant(submatrix, n - 1);
        sign = -sign;
    }
    return det;
}

// Ví d? s? d?ng
int main() {
    vector<vector<double>> matrix = {
        {2, 1, 0},
        {1, 3, 1},
        {0, 1, 2}
    };

    int n = matrix.size();
    double det = determinant(matrix, n);
    cout << "Determinant: " << det << endl;

    return 0;
}

