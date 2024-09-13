#include <iostream>
#include <E:\SCHOOL\CODE\HK5\TOAN UNG DUNG CNTT\eigen-3.4.0\eigen-3.4.0\Eigen\Dense>

using namespace Eigen;
using namespace std;

void computeEigenValuesAndVectors(const MatrixXd& matrix) {
    // Kiểm tra ma trận có phải vuông không
    if (matrix.rows() != matrix.cols()) {
        cout << "Ma trận không phải vuông." << endl;
        return;
    }

    EigenSolver<MatrixXd> solver(matrix);

    // Lấy các giá trị riêng (eigenvalues)
    cout << "Giá trị riêng:" << endl;
    cout << solver.eigenvalues() << endl;

    // Lấy các vector riêng (eigenvectors)
    cout << "Vector riêng:" << endl;
    cout << solver.eigenvectors() << endl;
}

int main() {
    // Ma trận ví dụ 3x3
    MatrixXd matrix(3, 3);
    matrix << 1, 2, 3,
              0, 4, 5,
              0, 0, 6;

    cout << "Ma trận ban đầu:" << endl;
    cout << matrix << endl << endl;

    computeEigenValuesAndVectors(matrix);

    return 0;
}
