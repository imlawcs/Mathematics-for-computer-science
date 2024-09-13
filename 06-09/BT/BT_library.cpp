#include <iostream>
#include "E:\SCHOOL\CODE\HK5\TOAN UNG DUNG CNTT\eigen-3.4.0\eigen-3.4.0\Eigen\Dense"

using namespace std;
using namespace Eigen;

int main() {
    // Khởi tạo ma trận vuông
    MatrixXd A(3, 3);
    A << 4, 12, -16,
         12, 37, -43,
         -16, -43, 98;

    // Tìm giá trị riêng và vector riêng
    Eigen::EigenSolver<MatrixXd> solver(A);

    // Ma trận đường chéo
    MatrixXd D = solver.pseudoEigenvalueMatrix();

    // Ma trận vector riêng
    MatrixXd P = solver.eigenvectors().real();

    // Ma trận nghịch đảo của ma trận vector riêng
    MatrixXd P_inv = P.inverse();

    // Kiểm tra P_inv * A * P = D
    MatrixXd A_reconstructed = P_inv * A * P;
    
    cout << "Original Matrix A:" << endl << A << endl << endl;
    cout << "Diagonal Matrix D:" << endl << D << endl << endl;
    cout << "Matrix P (Eigenvectors):" << endl << P << endl << endl;
    cout << "Reconstructed Matrix A from P_inv * A * P:" << endl << A_reconstructed << endl << endl;

    return 0;
}

// #include <iostream>
// #include <Eigen/Dense>

// using namespace std;
// using namespace Eigen;

// int main() {
//     int n;
//     cout << "Enter the size of the matrix (n): ";
//     cin >> n;

//     // Khởi tạo ma trận vuông A kích thước n x n
//     MatrixXd A(n, n);

//     // Nhập các phần tử cho ma trận A
//     cout << "Enter the elements of the matrix:" << endl;
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
//             cin >> A(i, j);
//         }
//     }

//     cout << "Original matrix A:" << endl;
//     cout << A << endl;

//     // Tìm trị riêng và vector riêng sử dụng Eigen
//     EigenSolver<MatrixXd> es(A);

//     // Ma trận D là ma trận chéo (diagonal matrix) chứa các trị riêng
//     MatrixXd D = es.eigenvalues().asDiagonal();

//     // Ma trận P là ma trận các vector riêng
//     MatrixXd P = es.eigenvectors().real();

//     // Ma trận P^-1 (nghịch đảo của ma trận P)
//     MatrixXd P_inv = P.inverse();

//     cout << "\nEigenvalues (Diagonal matrix D):" << endl;
//     cout << D << endl;

//     cout << "\nEigenvectors (Matrix P):" << endl;
//     cout << P << endl;

//     cout << "\nMatrix P^-1 (Inverse of P):" << endl;
//     cout << P_inv << endl;

//     // Kiểm tra P^-1 * A * P = D
//     MatrixXd P_inv_AP = P_inv * A * P;

//     cout << "\nP^-1 * A * P:" << endl;
//     cout << P_inv_AP << endl;

//     // So sánh P^-1 * A * P với D
//     bool isEqual = (P_inv_AP.isApprox(D, 1e-9));
//     cout << "\nVerification of P^-1 * A * P = D: " << (isEqual ? "Correct" : "Incorrect") << endl;

//     return 0;
// }

