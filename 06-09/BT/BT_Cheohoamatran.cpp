#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Hàm tính định thức của ma trận vuông (bằng phương pháp đệ quy)
double determinant(const vector<vector<double>>& matrix, int n) {
    double det = 0;
    if (n == 1) {
        return matrix[0][0];
    }

    vector<vector<double>> submatrix(n - 1, vector<double>(n - 1));
    int sign = 1;

    for (int x = 0; x < n; x++) {
        int subi = 0;
        for (int i = 1; i < n; i++) {
            int subj = 0;
            for (int j = 0; j < n; j++) {
                if (j == x) {
                    continue;
                }
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

// Hàm tính ma trận nghịch đảo 2x2
vector<vector<double>> inverse2x2(const vector<vector<double>>& matrix) {
    vector<vector<double>> inverse(2, vector<double>(2));
    double det = determinant(matrix, 2);

    if (det == 0) {
        cout << "Matrix is singular and cannot be inverted." << endl;
        exit(0);
    }

    inverse[0][0] = matrix[1][1] / det;
    inverse[0][1] = -matrix[0][1] / det;
    inverse[1][0] = -matrix[1][0] / det;
    inverse[1][1] = matrix[0][0] / det;

    return inverse;
}

// Hàm nhân hai ma trận
vector<vector<double>> multiplyMatrices(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    vector<vector<double>> result(A.size(), vector<double>(B[0].size(), 0));
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < B[0].size(); j++) {
            for (int k = 0; k < A[0].size(); k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// Hàm tìm các giá trị riêng của ma trận (chỉ giải cho ma trận 2x2)
vector<double> findEigenvalues(const vector<vector<double>>& matrix) {
    vector<double> eigenvalues(2);

    double a = 1;
    double b = -(matrix[0][0] + matrix[1][1]);
    double c = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    double discriminant = b * b - 4 * a * c;
    if (discriminant >= 0) {
        eigenvalues[0] = (-b + sqrt(discriminant)) / (2 * a);
        eigenvalues[1] = (-b - sqrt(discriminant)) / (2 * a);
    } else {
        cout << "Eigenvalues are complex numbers, diagonalization is not possible." << endl;
        exit(0);
    }

    return eigenvalues;
}

// Hàm tìm vector riêng của ma trận (chỉ giải cho ma trận 2x2)
vector<vector<double>> findEigenvectors(const vector<vector<double>>& matrix, const vector<double>& eigenvalues) {
    vector<vector<double>> eigenvectors(2, vector<double>(2));

    for (int i = 0; i < 2; i++) {
        double lambda = eigenvalues[i];
        eigenvectors[i][0] = matrix[0][1];
        eigenvectors[i][1] = lambda - matrix[0][0];
    }

    return eigenvectors;
}

// Hàm in ma trận
void printMatrix(const vector<vector<double>>& matrix) {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

int main() {
    vector<vector<double>> matrix = {
        {4, 1},
        {2, 3}
    };

    cout << "Original matrix:" << endl;
    printMatrix(matrix);

    // Bước 1: Tìm các giá trị riêng
    vector<double> eigenvalues = findEigenvalues(matrix);

    cout << "Eigenvalues: " << eigenvalues[0] << ", " << eigenvalues[1] << endl;

    // Bước 2: Tìm vector riêng
    vector<vector<double>> eigenvectors = findEigenvectors(matrix, eigenvalues);

    cout << "Eigenvectors:" << endl;
    printMatrix(eigenvectors);

    // Bước 3: Xây dựng ma trận P và D
    vector<vector<double>> P = eigenvectors;
    vector<vector<double>> D = {
        {eigenvalues[0], 0},
        {0, eigenvalues[1]}
    };

    cout << "Matrix P (eigenvectors):" << endl;
    printMatrix(P);

    cout << "Matrix D (diagonal matrix):" << endl;
    printMatrix(D);

    // Bước 4: Tính ma trận P^-1
    vector<vector<double>> P_inv = inverse2x2(P);

    // Bước 5: Kiểm tra lại P^-1 * A * P = D
    vector<vector<double>> P_inv_AP = multiplyMatrices(P_inv, multiplyMatrices(matrix, P));

    cout << "Matrix P^-1 * A * P:" << endl;
    printMatrix(P_inv_AP);

    cout << "Matrix D (diagonal matrix):" << endl;
    printMatrix(D);

    // So sánh P_inv_AP với D
    bool isEqual = true;
    for (int i = 0; i < D.size(); i++) {
        for (int j = 0; j < D[i].size(); j++) {
            if (abs(P_inv_AP[i][j] - D[i][j]) > 1e-9) {
                isEqual = false;
                break;
            }
        }
        if (!isEqual) break;
    }

    cout << "Verification of P^-1 * A * P = D: " << (isEqual ? "Correct" : "Incorrect") << endl;

    return 0;
}
