#include <bits/stdc++.h>


using namespace std;


int rows = 3, columns = 3;
float matrix[100][100];
float oldCholesky[100][100];
float k[100][100], d[100][100];


void inputValue() {
    cout << "   Nhap so hang cho ma tran: "; cin >> rows;
    cout << "   Nhap so cot cho ma tran: "; cin >> columns;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            cout << "- Nhap phan tu matrix[" << i + 1 << "][" << j + 1 << "] = "; cin >> matrix[i][j];
            oldCholesky[i][j] = 0;
        }
    }
}


void outputValue(float matrixx[][100]) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            cout << matrixx[i][j] << "\t";
        }
        cout << "\n";
    }
}


bool checkSquare() {
    return (rows == columns);
}


bool checkReflection(float matrixx[][100]) {
    if (checkSquare()) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                if (i != j && matrixx[i][j] != matrixx[j][i])
                    return false;
            }
        }
        return true;
    }
    else return false;
}


double determinant(int n, float a[][100]) {
    double det = 1.0;
    for (int i = 0; i < n; i++) {
        int pivot = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(a[j][i]) > abs(a[pivot][i])) {
                pivot = j;
            }
        }
        if (pivot != i) {
            swap(a[i], a[pivot]);
            det *= -1;
        }
        if (a[i][i] == 0) {
            return 0;
        }
        det *= a[i][i];
        for (int j = i + 1; j < n; j++) {
            double factor = a[j][i] / a[i][i];
            for (int k = i + 1; k < n; k++) {
                a[j][k] -= factor * a[i][k];
            }
        }
    }
    return det;
}


bool checkPositiveDefinite(float matrixx[][100]) {
    bool check = true;
    for (int i = 0; i < rows; i++) {
        if (matrixx[i][i] <= 0) check = false;
    }
    if (checkReflection(matrixx) == true) {
        for (int i = 0; i < rows; i++) {
            float temp[100][100];
            for (int k = 0; k <= i; k++) {
                for (int j = 0; j <= i; j++) {
                    temp[k][j] = matrixx[k][j];
                }
            }
            if (determinant(i + 1, temp) < 0) check = false;
        }
    }
    return check;
}


void multiplyMatrix (float matrix1[][100], float matrix2[][100], float matrix3[][100]) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            matrix3[i][j] = 0;
            for (int k = 0; k < rows; k++) {
                matrix3[i][j] += (1.0) * matrix1[i][k] * matrix2[k][j];
            }
        }
    }
}


void calCholeskyBienThe() {
    if (checkPositiveDefinite(matrix)) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j <= i; j++) {
                int sum = 0;
                if (j == i)
                {
                    for (int k = 0; k < j; k++)
                        sum += pow(oldCholesky[j][k], 2);
                        oldCholesky[j][j] = sqrt(matrix[j][j] - sum);
                } else {
                    for (int k = 0; k < j; k++)
                        sum += (oldCholesky[i][k] * oldCholesky[j][k]);
                        oldCholesky[i][j] = (matrix[i][j] - sum) / oldCholesky[j][j];
                }
            }
        }
        float sMuTru1[100][100];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                if (i == j) sMuTru1[i][j] = 1.0 / oldCholesky[i][j];
                else sMuTru1[i][j] = 0;
            }
        }
        multiplyMatrix(oldCholesky, sMuTru1, k);
        cout << "(*)Ma tran K = L*S^(-1):" << endl;
        outputValue(k);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                if (i == j) d[i][j] = oldCholesky[i][j] * oldCholesky[i][j];
                else d[i][j] = 0;
            }
        }
        cout << "(*)Ma tran D = S*S:" << endl;
        outputValue(d);
        float kMuTru1[100][100];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                kMuTru1[i][j] = k[j][i];
            }
        }
        cout << "(*)Ma tran K^(-1): " << endl;
        outputValue(kMuTru1);
        cout << "=> A = K*D*K^(-1)." << endl;
    }
}


int main() {
    inputValue();
    cout << "Cau 2: Phan ra Cholesky bien the: " << endl;
    calCholeskyBienThe();
}

