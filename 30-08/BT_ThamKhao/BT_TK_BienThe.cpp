#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

bool CholeskyDecomposition(vector<vector<double>>& A, vector<vector<double>>& L) {
    int n = A.size();

    // Initialize L with 0s
    L.assign(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;

            for (int k = 0; k < j; ++k)
                sum += L[i][k] * L[j][k];

            if (i == j)
                L[i][j] = sqrt(A[i][i] - sum);
            else
                L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
        }

        // If the matrix isn't positive definite, the diagonal will be <= 0
        if (L[i][i] <= 0)
            return false;
    }

    return true;
}

int main() {
    vector<vector<double>> A = {
        {4, 12, -16},
        {12, 37, -43},
        {-16, -43, 98}
    };

    vector<vector<double>> L;

    if (CholeskyDecomposition(A, L)) {
        cout << "Matrix L (Cholesky factor):\n";
        for (const auto& row : L) {
            for (double val : row)
                cout << val << " ";
            cout << "\n";
        }
    } else {
        cout << "Matrix is not positive definite.\n";
    }

    return 0;
}
