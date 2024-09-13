#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
using namespace std;

// Cấp phát ma trận
vector<vector<double>> allocated_matrix(int row, int col) {
    return vector<vector<double>>(row, vector<double>(col, 0.0));
}

// Hàm in ma trận
void print_matrix(const vector<vector<double>>& matrix) {
    int r = matrix.size();
    int c = matrix[0].size();
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            if (fabs(matrix[i][j]) < 1e-6) cout << setw(5) << 0 << " ";
            else cout << setw(5) << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// Hàm nhân ma trận
vector<vector<double>> matrix_multiplication(const vector<vector<double>>& matrixA, const vector<vector<double>>& matrixB) {
    int n = matrixA.size();
    int m = matrixB[0].size();
    vector<vector<double>> resultMatrix = allocated_matrix(n, m);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            resultMatrix[i][j] = 0;
            for (int k = 0; k < matrixB.size(); ++k) {
                resultMatrix[i][j] += matrixA[i][k] * matrixB[k][j];
            }
            if (fabs(resultMatrix[i][j]) < 1e-6) resultMatrix[i][j] = 0;
        }
    }
    return resultMatrix;
}

// Hàm tính giá trị đa thức tại x
double evaluate_polynomial(const vector<double>& coeffs, double x) {
    double result = 0;
    int n = coeffs.size();
    for (int i = 0; i < n; ++i) {
        result += coeffs[i] * pow(x, n - 1 - i);
    }
    return result;
}

// Hàm tính đạo hàm của đa thức tại x
double evaluate_derivative(const vector<double>& coeffs, double x) {
    double result = 0;
    int n = coeffs.size();
    for (int i = 0; i < n - 1; ++i) {
        result += (n - i - 1) * coeffs[i] * pow(x, n - 2 - i);
    }
    return result;
}

double newton_raphson(const vector<double>& coeffs, double initGuessValue) {
    double x = initGuessValue;
    for (int i = 0; i < 1000; ++i) {
        double f = evaluate_polynomial(coeffs, x);
        double f_derivative = evaluate_derivative(coeffs, x);
        double x_next = x - f / f_derivative;
        if (fabs(x_next - x) < 1e-6) return x_next;
        x = x_next;
    }
    return x;
}

vector<double> divide_polynomial(const vector<double>& coeffs, double root) {
    int n = coeffs.size();
    vector<double> newCoeffs(n - 1);
    newCoeffs[0] = coeffs[0];
    for (int i = 1; i < n - 1; ++i) {
        newCoeffs[i] = coeffs[i] + root * newCoeffs[i - 1];
    }
    return newCoeffs;
}

vector<double> find_all_roots(const vector<double>& coeffs, double initGuesses) {
    vector<double> roots;
    vector<double> currentCoeffs = coeffs;
    while (currentCoeffs.size() > 1) {
        double root = newton_raphson(currentCoeffs, initGuesses);
        currentCoeffs = divide_polynomial(currentCoeffs, root);
        roots.push_back(root);
    }
    return roots;
}

vector<double> danilevsky(const vector<vector<double>>& matrix, vector<vector<double>>& matrixProductM) {
    int n = matrix.size();
    vector<vector<double>> A = matrix;
    vector<vector<double>> M = allocated_matrix(n, n);
    vector<vector<double>> M_minus1 = allocated_matrix(n, n);
    matrixProductM = allocated_matrix(n, n);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            M[i][j] = (i == j) ? 1 : 0;
            M_minus1[i][j] = 0;
            matrixProductM[i][j] = (i == j) ? 1 : 0;
        }
    }

    for (int k = n - 2; k >= 0; --k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != k) {
                    if (i == j) {
                        M[i][j] = 1;
                        M_minus1[i][j] = 1;
                    } else {
                        M[i][j] = 0;
                        M_minus1[i][j] = 0;
                    }
                } else {
                    M_minus1[i][j] = A[k + 1][j];
                    M[i][j] = (j == k) ? 1 / A[k + 1][k] : -A[k + 1][j] / A[k + 1][k];
                }
            }
        }
        vector<vector<double>> B = matrix_multiplication(A, M);
        A = matrix_multiplication(M_minus1, B);
        matrixProductM = matrix_multiplication(matrixProductM, M);
    }

    vector<double> coeffs(n + 1);
    coeffs[0] = 1;
    for (int i = 0; i < n; ++i) {
        coeffs[i + 1] = -A[0][i];
    }
    return coeffs;
}

// Tìm ma trận P
vector<vector<double>> find_matrix_P(const vector<vector<double>>& matrixA, vector<vector<double>>& matrixD) {
    int n = matrixA.size();
    vector<vector<double>> matrixProductM;
    vector<vector<double>> matrixP = allocated_matrix(n, n);
    vector<double> eigenValues = find_all_roots(danilevsky(matrixA, matrixProductM), -100);
    
    for (int i = 0; i < n; ++i) {
        double eigenValue = eigenValues[i];
        for (int j = 0; j < n; ++j) {
            matrixP[j][i] = pow(eigenValue, n - 1 - j);
            matrixD[i][j] = 0;
        }
        matrixD[i][i] = eigenValue;
    }
    
    matrixP = matrix_multiplication(matrixProductM, matrixP);
    return matrixP;
}

// Ma trận nghịch đảo P ^ -1
vector<vector<double>> find_inverse_matrix(const vector<vector<double>>& matrix) {
    int n = matrix.size();
    vector<vector<double>> P_minus1 = allocated_matrix(n, 2 * n);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            P_minus1[i][j] = matrix[i][j];
        }
        for (int j = n; j < 2 * n; ++j) {
            P_minus1[i][j] = (j == n + i) ? 1 : 0;
        }
    }
    
    for (int i = 0; i < n - 1; ++i) {
        if (P_minus1[i][i] == 0) {
            int cnt = (i < n) ? i : n;
            while (P_minus1[cnt][i] == 0 && cnt <= n) ++cnt;
            if (cnt == n + 1) continue;
            swap(P_minus1[i], P_minus1[cnt]);
        }
        for (int j = i + 1; j < n; ++j) {
            double ratio = P_minus1[j][i] / P_minus1[i][i];
            for (int k = 0; k < 2 * n; ++k) {
                P_minus1[j][k] -= ratio * P_minus1[i][k];
            }
        }
    }
    
    double det = 1;
    for (int i = 0; i < n; ++i) det *= P_minus1[i][i];
    if (det == 0) {
        cout << "Matrix is singular.\n";
        return vector<vector<double>>(n, vector<double>(n, 0));
    }
    
    for (int i = n - 2; i >= 0; --i) {
        for (int j = i + 1; j < n; ++j) {
            double ratio = P_minus1[i][j] / P_minus1[j][j];
            for (int k = 0; k < 2 * n; ++k) {
                P_minus1[i][k] -= ratio * P_minus1[j][k];
            }
        }
    }
    
    for (int i = 0; i < n; ++i) {
        double k = P_minus1[i][i];
        for (int j = 0; j < 2 * n; ++j) {
            P_minus1[i][j] /= k;
        }
    }
    
    vector<vector<double>> result = allocated_matrix(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            result[i][j] = P_minus1[i][j + n];
            if (fabs(result[i][j]) < 1e-6) result[i][j] = 0;
        }
    }
    return result;
}

int main() {
    int n;
    cout << "Enter the size of the matrix n = "; cin >> n;
    vector<vector<double>> matrixA = allocated_matrix(n, n);
    
    cout << "Enter matrix A:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> matrixA[i][j];
        }
    }
    
    vector<vector<double>> matrixD = allocated_matrix(n, n);
    vector<vector<double>> matrixP = find_matrix_P(matrixA, matrixD);
    
    cout << "Matrix P:\n";
    print_matrix(matrixP);
    
    cout << "Matrix D:\n";
    print_matrix(matrixD);
    
    cout << "Matrix P^(-1):\n";
    vector<vector<double>> matrixP_inv = find_inverse_matrix(matrixP);
    print_matrix(matrixP_inv);

    // Kiểm tra P * D * P^-1
    vector<vector<double>> PD = matrix_multiplication(matrixP, matrixD);
    vector<vector<double>> P_D_P_inv = matrix_multiplication(PD, matrixP_inv);

    cout << "Matrix P * D * P^(-1):\n";
    print_matrix(P_D_P_inv);

    // So sánh P * D * P^(-1) với A
    bool isEqual = true;
    for (int i = 0; i < matrixA.size(); i++) {
        for (int j = 0; j < matrixA[i].size(); j++) {
            if (fabs(P_D_P_inv[i][j] - matrixA[i][j]) > 1e-9) {
                isEqual = false;
                break;
            }
        }
        if (!isEqual) break;
    }

    cout << "Verification of P * D * P^(-1) = A: " << (isEqual ? "Correct" : "Incorrect") << endl;
    
    return 0;
}
