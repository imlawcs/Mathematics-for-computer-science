#include <iostream>
#include <Eigen/Dense>

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
