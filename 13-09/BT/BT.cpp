#include <iostream>
#include <E:\SCHOOL\CODE\HK5\TOAN UNG DUNG CNTT\eigen-3.4.0\eigen-3.4.0\Eigen\Dense>
#include <vector>
#include <iomanip>
#include <math.h>

using namespace Eigen;
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

// Ma trận nghịch đảo A ^ -1
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

// Hàm tính giá trị riêng
vector<double> computeEigenValues(const vector<vector<double>>& matrix) {
    int n = matrix.size();
    MatrixXd A(n, n);
    
    // Chuyển đổi ma trận từ vector sang MatrixXd của Eigen
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i, j) = matrix[i][j];
        }
    }

    // Sử dụng EigenSolver để tính giá trị riêng
    EigenSolver<MatrixXd> solver(A);
    VectorXd eigenValues = solver.eigenvalues().real(); // Chỉ lấy phần thực

    // Chuyển đổi từ VectorXd sang vector
    vector<double> eigenVals(n);
    for (int i = 0; i < n; ++i) {
        eigenVals[i] = eigenValues(i);
    }

    return eigenVals;
}

// Hàm tính vector riêng
vector<vector<double>> computeEigenVectors(const vector<vector<double>>& matrix) {
    int n = matrix.size();
    MatrixXd A(n, n);
    
    // Chuyển đổi ma trận từ vector sang MatrixXd của Eigen
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i, j) = matrix[i][j];
        }
    }

    // Sử dụng EigenSolver để tính vector riêng
    EigenSolver<MatrixXd> solver(A);
    MatrixXd eigenVectors = solver.eigenvectors().real(); // Chỉ lấy phần thực

    // Chuyển đổi từ MatrixXd sang vector
    vector<vector<double>> eigenVecs(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            eigenVecs[i][j] = eigenVectors(i, j);
        }
    }

    return eigenVecs;
}

int main () {
    int n;
    cout << "Nhập bậc của ma trận cần tính:";
    cin >> n;
    vector<vector<double>> A = allocated_matrix(n, n);

    cout << "Enter matrix A:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> A[i][j];
        }
    }

    vector<vector<double>> AT = allocated_matrix(n, n);
    vector<vector<double>> S = allocated_matrix(n, n);

    S = matrix_multiplication(A.transpose(), A);

    print_matrix(S);

    vector<double> gtr;
    gtr = computeEigenValues(S);
    for(int i = 0; i <= n; i++) {
        cout << gtr[i] << " ";
    }
    return 0;
}

