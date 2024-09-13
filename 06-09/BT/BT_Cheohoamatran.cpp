#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

using namespace std;

vector<vector<double>> allocated_matrix(int rows, int cols) {
    return vector<vector<double>>(rows, vector<double>(cols, 0.0));
}

vector<vector<double>> matrix_multiplication(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    int rowsA = A.size();
    int colsA = A[0].size();
    int colsB = B[0].size();

    vector<vector<double>> result(rowsA, vector<double>(colsB, 0.0));

    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < colsB; ++j) {
            for (int k = 0; k < colsA; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

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

// Hàm tính ma tr?n ngh?ch ð?o 3x3
vector<vector<double>> inverse3x3(const vector<vector<double>>& matrix) {
    vector<vector<double>> inverse(3, vector<double>(3));
    double det = determinant(matrix, 3);

    if (det == 0) {
        throw runtime_error("Matrix is singular and cannot be inverted.");
    }

    inverse[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / det;
    inverse[0][1] = -(matrix[0][1] * matrix[2][2] - matrix[0][2] * matrix[2][1]) / det;
    inverse[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / det;
    inverse[1][0] = -(matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) / det;
    inverse[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / det;
    inverse[1][2] = -(matrix[0][0] * matrix[1][2] - matrix[0][2] * matrix[1][0]) / det;
    inverse[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / det;
    inverse[2][1] = -(matrix[0][0] * matrix[2][1] - matrix[0][1] * matrix[2][0]) / det;
    inverse[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / det;

    return inverse;
}

// Hàm nhân hai ma tr?n
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

// Hàm ki?m tra tính ð?c l?p tuy?n tính c?a vector
bool areVectorsLinearlyIndependent(const vector<vector<double>>& vectors) {
    // S? d?ng ð?nh th?c c?a ma tr?n t?o t? các vector
    double det = determinant(vectors, 3);
    return fabs(det) > 1e-9;
}

vector<double> divide_polynomial(const vector<double>& coeffs, double root) {
    vector<double> result(coeffs.size() - 1);
    result[0] = coeffs[0];

    for (size_t i = 1; i < coeffs.size() - 1; ++i) {
        result[i] = result[i - 1] * root + coeffs[i];
    }

    return result;
}

double newton_raphson(const vector<double>& coeffs, double guess) {
    double tolerance = 1e-9;
    double max_iter = 1000;
    double x = guess;
    int iter = 0;

    while (iter < max_iter) {
        double fx = 0;
        double dfx = 0;
        
        // Tính giá tr? c?a ða th?c t?i x
        for (size_t i = 0; i < coeffs.size(); ++i) {
            fx += coeffs[i] * pow(x, coeffs.size() - 1 - i);
        }

        // Tính ð?o hàm c?a ða th?c t?i x
        for (size_t i = 0; i < coeffs.size() - 1; ++i) {
            dfx += (coeffs.size() - 1 - i) * coeffs[i] * pow(x, coeffs.size() - 2 - i);
        }

        if (fabs(dfx) < tolerance) {
            throw runtime_error("Derivative near zero. Newton-Raphson method failed.");
        }

        double x_new = x - fx / dfx;

        if (fabs(x_new - x) < tolerance) {
            return x_new;
        }

        x = x_new;
        ++iter;
    }

    throw runtime_error("Newton-Raphson method did not converge.");
}


// Hàm t?m t?t c? các nghi?m c?a ða th?c
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
    vector<vector<double>> B = allocated_matrix(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            M[i][j] = 0;
            M_minus1[i][j] = 0;
            matrixProductM[i][j] = (i == j) ? 1 : 0;
        }
    }

    for (int k = n - 2; k >= 0; --k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != k) {
                    M[i][j] = (i == j) ? 1 : 0;
                    M_minus1[i][j] = 0;
                } else {
                    M_minus1[i][j] = A[k + 1][j];
                    M[i][j] = (j == k) ? 1 / A[k + 1][k] : -A[k + 1][j] / A[k + 1][k];
                }
            }
        }
        B = matrix_multiplication(A, M);
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


// Hàm t?m ma tr?n P và D
vector<vector<double>> find_matrix_P(const vector<vector<double>>& matrixA, vector<vector<double>>& matrixD) {
    int n = matrixA.size();
    vector<vector<double>> matrixProductM = allocated_matrix(n, n);
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

vector<double> findEigenvalues(const vector<vector<double>>& matrix) {
    int n = matrix.size();
    vector<vector<double>> matrixProductM(n, vector<double>(n));
    vector<double> coeffs = danilevsky(matrix, matrixProductM);
    return find_all_roots(coeffs, 0.0);
}


// Hàm gi?i h? phýõng tr?nh tuy?n tính
vector<double> solveLinearSystem(const vector<vector<double>>& A, const vector<double>& b) {
    // S? d?ng phýõng pháp gi?i h? phýõng tr?nh, ví d?: phýõng pháp Gauss
    int n = A.size();
    vector<vector<double>> augmented(n, vector<double>(n + 1));
    vector<double> x(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][n] = b[i];
    }

    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (fabs(augmented[k][i]) > fabs(augmented[maxRow][i])) {
                maxRow = k;
            }
        }

        swap(augmented[i], augmented[maxRow]);

        for (int k = i + 1; k < n; ++k) {
            double factor = augmented[k][i] / augmented[i][i];
            for (int j = i; j <= n; ++j) {
                augmented[k][j] -= factor * augmented[i][j];
            }
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        x[i] = augmented[i][n] / augmented[i][i];
        for (int k = i - 1; k >= 0; --k) {
            augmented[k][n] -= augmented[k][i] * x[i];
        }
    }

    return x;
}

vector<vector<double>> findEigenvectors(const vector<vector<double>>& matrix, const vector<double>& eigenvalues) {
    int n = matrix.size();
    vector<vector<double>> eigenvectors(n, vector<double>(n));

    for (size_t i = 0; i < eigenvalues.size(); ++i) {
        double lambda = eigenvalues[i];
        vector<vector<double>> A = matrix;
        
        // Tr? ma tr?n ?I
        for (int j = 0; j < n; ++j) {
            A[j][j] -= lambda;
        }

        // Gi?i h? phýõng tr?nh (A - ?I)x = 0 ð? t?m vector eigen
        vector<double> b(n, 0.0);
        b[i] = 1; // Ch?n vector b[i] làm vector không nông
        vector<double> x = solveLinearSystem(A, b);

        for (int j = 0; j < n; ++j) {
            eigenvectors[j][i] = x[j];
        }
    }

    return eigenvectors;
}

int main() {
    vector<vector<double>> matrix = {
        {2, 1, 0},
        {1, 3, 1},
        {0, 1, 2}
    };

    cout << "Original matrix:" << endl;
    for (const auto& row : matrix) {
        for (double value : row) {
            cout << value << " ";
        }
        cout << endl;
    }

    // Bý?c 1: T?m các giá tr? riêng
    vector<double> eigenvalues = findEigenvalues(matrix);

    cout << "Eigenvalues: " << eigenvalues[0] << ", " << eigenvalues[1] << ", " << eigenvalues[2] << endl;

    // Bý?c 2: T?m vector riêng
    vector<vector<double>> eigenvectors = findEigenvectors(matrix, eigenvalues);

    cout << "Eigenvectors:" << endl;
    for (const auto& row : eigenvectors) {
        for (double value : row) {
            cout << value << " ";
        }
        cout << endl;
    }

    // Bý?c 3: Xây d?ng ma tr?n P và D
    vector<vector<double>> P = eigenvectors;
    vector<vector<double>> D = {
        {eigenvalues[0], 0, 0},
        {0, eigenvalues[1], 0},
        {0, 0, eigenvalues[2]}
    };

    cout << "Matrix P (eigenvectors):" << endl;
    for (const auto& row : P) {
        for (double value : row) {
            cout << value << " ";
        }
        cout << endl;
    }

    cout << "Matrix D (diagonal matrix):" << endl;
    for (const auto& row : D) {
        for (double value : row) {
            cout << value << " ";
        }
        cout << endl;
    }

    try {
        // Bý?c 4: Tính ma tr?n P^-1
        vector<vector<double>> P_inv = inverse3x3(P);

        // Bý?c 5: Ki?m tra l?i P^-1 * A * P = D
        vector<vector<double>> P_inv_AP = multiplyMatrices(P_inv, multiplyMatrices(matrix, P));

        cout << "Matrix P^-1 * A * P:" << endl;
        for (const auto& row : P_inv_AP) {
            for (double value : row) {
                cout << value << " ";
            }
            cout << endl;
        }

        cout << "Matrix D (diagonal matrix):" << endl;
        for (const auto& row : D) {
            for (double value : row) {
                cout << value << " ";
            }
            cout << endl;
        }

        // So sánh P_inv_AP v?i D
        bool isEqual = true;
        for (int i = 0; i < D.size(); i++) {
            for (int j = 0; j < D[i].size(); j++) {
                if (fabs(P_inv_AP[i][j] - D[i][j]) > 1e-9) {
                    isEqual = false;
                    break;
                }
            }
            if (!isEqual) break;
        }

        cout << "Verification of P^-1 * A * P = D: " << (isEqual ? "Correct" : "Incorrect") << endl;
    } catch (const runtime_error& e) {
        cout << e.what() << endl;
    }

    return 0;
}

