#include<iostream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<algorithm>
using namespace std;


#define EPSILON 1e-6
#define TOLERANCE 1000
//* Cấp phát động ma trận
double **allocated_matrix(int row, int col){
  double **matrix = new double *[row];
  for(int i = 0; i < row; ++i){
    matrix[i] = new double[col];
  }
  return matrix;
}


//* Xóa ma trận
void free_matrix(double **matrix, int row){
  for(int i = 0; i < row; ++i){
    delete[] matrix[i];
  }
  delete[] matrix;
}


//* Hàm in ma trận
void print_matrix(double **matrix, int r, int c){
  for(int i = 0; i < r; ++i){
    for(int j = 0; j < c; ++j){
      if(fabs(matrix[i][j]) < EPSILON) cout << setw(5) << 0 << " ";
      else cout << setw(5) << matrix[i][j] << " ";
    }
    cout << endl;
  }
}


//* Hàm nhân ma trận
double **matrix_multiplication(double **matrixA, double **matrixB, int n, int m, int p){
  double **resultMatrix = allocated_matrix(n, p);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < p; ++j){
      resultMatrix[i][j] = 0;
      for(int k = 0; k < m; ++k){
        resultMatrix[i][j] += matrixA[i][k] * matrixB[k][j]; 
      }
      if(fabs(resultMatrix[i][j]) < EPSILON) resultMatrix[i][j] = 0;
    }
  }
  return resultMatrix;  
}


//* Hàm tính giá trị đa thức tại x
double evaluate_polynomial(const vector<double> coeffs, double x){
  double result = 0;
  int n = coeffs.size();
  for(int i = 0; i < n; ++i){
    result += coeffs[i] * pow(x, n - 1 - i);
  }
  return result;
}


//* Hàm tính đạo hàm của đa thức tại x
double evaluate_derivative(const vector<double> coeffs, double x) {
  double result = 0;
  int n = coeffs.size();
  for (int i = 0; i < n - 1; ++i) {
    result += (n - i - 1) * coeffs[i] * pow(x, n - 2 - i);
  }
  return result;
}


double newton_raphson(const vector<double> coeffs, double initGuessValue){
  double x = initGuessValue;
  for(int i = 0; i < TOLERANCE; ++i){
    double f = evaluate_polynomial(coeffs, x);
    double f_derivative = evaluate_derivative(coeffs, x);
    double x_next = x - f / f_derivative;
    if(fabs(x_next - x) < EPSILON) return x_next;
    x = x_next;
  }
  return x;
}


vector<double> divide_polynomial(const vector<double> coeffs, double root){
  int n = coeffs.size();
  vector<double> newCoeffs(n - 1);
  newCoeffs[0] = coeffs[0];
  for(int i = 1; i < n - 1; ++i){
    newCoeffs[i] = coeffs[i] + root * newCoeffs[i - 1];
  }
  return newCoeffs;
}


vector<double> find_all_roots(const vector<double> &coeffs, const double initGuesses){
  vector<double> roots;
  vector<double> currentCoeffs = coeffs;
  while(currentCoeffs.size() > 1){
    double root = newton_raphson(currentCoeffs, initGuesses);
    currentCoeffs = divide_polynomial(currentCoeffs, root);
    roots.push_back(root);
  }
  return roots;
}


vector<double> danilevsky(double **matrix, double **&matrixProductM,int n){
  double **A = allocated_matrix(n,n);
  double **M = allocated_matrix(n,n);
  double **M_minus1 = allocated_matrix(n,n);
  double **B = allocated_matrix(n,n);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      A[i][j] = matrix[i][j];
      M[i][j] = 0;
      M_minus1[i][j] = 0;
      if(i != j) matrixProductM[i][j] = 0;
      else matrixProductM[i][j] = 1;
    }
  }
  for(int k = n-2; k >= 0; --k){
    for(int i = 0; i < n; ++i){
      for(int j = 0; j < n; ++j){
        if(i != k){
          if(i == j){
            M[i][j] = 1;
            M_minus1[i][j] = 1;
          }
          else{
            M[i][j] = 0;
            M_minus1[i][j] = 0;
          }
        }
        else{
          M_minus1[i][j] = A[k+1][j];
          if(j == k) M[i][j] = 1/A[k+1][k];
          else M[i][j] = -A[k+1][j] / A[k+1][k];
        }
      }
    }
    double **B = matrix_multiplication(A,M,n,n,n);
    A = matrix_multiplication(M_minus1,B,n,n,n);
    matrixProductM = matrix_multiplication(matrixProductM,M,n,n,n);
  }
  vector<double> coeffs(n+1);
  coeffs[0] = 1;
  for(int i = 0; i < n; ++i){
    coeffs[i+1] = -A[0][i];
  }
  return coeffs;
}


//* Tìm trị riêng, vector riêng
double **find_matrix_P(double **matrixA, double** &matrixD, int n){
  double **matrixProductM = allocated_matrix(n, n);
  double **matrixP = allocated_matrix(n,n);
  vector<double> eigenValues = find_all_roots(danilevsky(matrixA,matrixProductM,n), -100);
  for(int i = 0; i < n; ++i){
    double eigenValue = eigenValues[i];
    for(int j = 0; j < n; ++j){
      matrixP[j][i] = pow(eigenValue, n-1-j);
      matrixD[i][j] = 0;
    }
    matrixD[i][i] = eigenValue;
  }
  matrixP = matrix_multiplication(matrixProductM,matrixP,n,n,n);
  return matrixP;
}


//* Ma trận nghịch đảo P ^ -1
void swap_rows(double **matrix, int r1, int r2, int c){
  for(int i = 0; i < c; ++i){
    double temp = matrix[r1][i];
    matrix[r1][i] = matrix[r2][i];
    matrix[r2][i] = temp;
  }
}


double** find_inverse_matrix(double **matrix, int n){
  double **P_minus1 = allocated_matrix(n,2 * n);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      P_minus1[i][j] = matrix[i][j];
    }
    for(int j = n; j < 2 * n; ++j){
      P_minus1[i][j] = 0;
    }
    P_minus1[i][n+i] = 1;
  }
  for(int i = 0; i < n - 1; ++i){
    if(P_minus1[i][i] == 0){
      int cnt = (i < n) ? i : n;
      while(P_minus1[cnt][i] == 0 && cnt <= n) ++cnt;
      if(cnt == n+1) continue;
      else{
        swap_rows(P_minus1, i, cnt, 2 * n);
      }
    }
    for(int j = i + 1; j < n; ++j){
      double ratio = P_minus1[j][i] / P_minus1[i][i];
      for(int k = 0; k < 2 * n; ++k){
        P_minus1[j][k] -= ratio * P_minus1[i][k];
      }
    }
  }
  double det = 1;
  for(int i = 0; i < n; ++i) det *= P_minus1[i][i];
  if(det == 0){
    cout << "Ma tran khong kha nghich.\n";
    return NULL;
  }
  for(int i = n - 2; i >= 0; --i){
    for(int j = i + 1; j < n; ++j){
      double ratio = P_minus1[i][j] / P_minus1[j][j];
      for(int k = 0; k < 2 * n; ++k){
        P_minus1[i][k] -= ratio * P_minus1[j][k];
      }
    }
  }
  for(int i = 0; i < n; ++i){
    double k = P_minus1[i][i];
    for(int j = 0; j < 2 * n; ++j){
      P_minus1[i][j] /= k;
    }
  }
  double **result = allocated_matrix(n,n);
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      result[i][j] = P_minus1[i][j + n];
      if(fabs(result[i][j]) < EPSILON) result[i][j] = 0;
    }
  }
  return result;
}


int main(){
  int n;
  cout << "Nhap kich thuoc cua ma tran vuong n = "; cin >> n;
  double **matrixA = allocated_matrix(n,n);
  cout << "Nhap ma tran A:\n";
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      cin >> matrixA[i][j];
    }
  }
  double **matrixD = allocated_matrix(n,n);
  double **matrixP = find_matrix_P(matrixA,matrixD,n);
  cout << "Ma tran P:\n";
  print_matrix(matrixP,n,n);
  cout << "Ma tran D:\n";
  print_matrix(matrixD,n,n);
  cout << "Ma tran P^(-1):\n";
  double **matrixP_inv = find_inverse_matrix(matrixP,n);
  print_matrix(matrixP_inv,n,n);


  cout << "Ma tran P * D * P^(-1) = A:\n";
  double **checkA = matrix_multiplication(matrixP, matrixD,n,n,n);
  checkA = matrix_multiplication(checkA, matrixP_inv,n,n,n);
  print_matrix(checkA,n,n);


  free_matrix(matrixA,n);
  free_matrix(matrixP,n);
  free_matrix(matrixP_inv,n);
  return 0;
}