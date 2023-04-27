#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <memory>
#include <string>

using namespace std;

class Matrix {
public:
    Matrix(int rows, int cols) : data(rows, vector<double>(cols, 0)) {}

    int numRows() const { return data.size(); }
    int numCols() const { return data[0].size(); }

    double& operator()(int row, int col) { return data[row][col]; }
    const double& operator()(int row, int col) const { return data[row][col]; }

    Matrix transpose() const;
    Matrix inverse() const;
    void print(const string& label) ;

private:
    vector<vector<double>> data;
};

Matrix operator*(const Matrix& A, const Matrix& B);
Matrix operator*(const Matrix& A, const vector<double>& b);

Matrix Matrix::transpose() const {
    Matrix AT(numCols(), numRows());
    for (int i = 0; i < numRows(); i++) {
        for (int j = 0; j < numCols(); j++) {
            AT(j, i) = (*this)(i, j);
        }
    }
    return AT;
}

Matrix operator*(const Matrix& A, const Matrix& B) {
    if (A.numCols() != B.numRows()) {
        throw runtime_error("Matrix dimensions do not match for multiplication.");
    }

    Matrix C(A.numRows(), B.numCols());
    for (int i = 0; i < A.numRows(); i++) {
        for (int j = 0; j < B.numCols(); j++) {
            for (int k = 0; k < A.numCols(); k++) {
                C(i, j) += A(i, k) * B(k, j);
            }
        }
    }
    return C;
}

Matrix operator*(const Matrix& A, const vector<double>& b) {
    if (A.numCols() != b.size()) {
        throw runtime_error("Matrix and vector dimensions do not match for multiplication.");
    }

    Matrix C(A.numRows(), 1);
    for (int i = 0; i < A.numRows(); i++) {
        for (int j = 0; j < A.numCols(); j++) {
            C(i, 0) += A(i, j) * b[j];
        }
    }
    return C;
}

Matrix Matrix::inverse() const {
    if (numRows() != numCols()) {
        throw runtime_error("Matrix must be square to compute its inverse.");
    }

    int n = numRows();
    Matrix B(n, n);
    for (int i = 0; i < n; i++) {
        B(i, i) = 1;
    }

    Matrix C = *this;
    for (int i = 0; i < n; i++) {
        double pivot = C(i, i);
        /*if (abs(pivot) < 1e-9) {
            throw runtime_error("Matrix is singular and cannot be inverted.");
        }*/

        for (int j = 0; j < n; j++) {
            C(i, j) /= pivot;
            B(i, j) /= pivot;
        }

        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            double factor = C(j, i);
            for (int k = 0; k < n; k++) {
                C(j, k) -= factor * C(i, k);
                B(j, k) -= factor * B(i, k);
            }
        }
    }
    return B;
}

void Matrix::print(const string& label)  {
    cout << label << endl;
    for (auto& row : data) {
        for (int i = 0; i < row.size(); i++) {
            if ((row[i] < 0.0000005)&&(row[i] > -0.0000005))
                row[i] = 0.0000;

            cout << fixed << setprecision(4) << row[i] << " ";
        }
        cout << endl;
    }
}

void print_vector(const vector<double>& v) {
    for (int i = 0; i < v.size(); i++) {
        if(i<v.size()-1){
            cout << fixed << setprecision(4) << v[i] << endl;
        }else{
            cout << fixed << setprecision(4) << v[i]<<endl;
        }
    }
}
#ifdef WIN32
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif
int main() {
#ifdef WIN32
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE* pipe = popen(GNUPLOT_NAME, "w");
#endif
    int m = 10;
    vector<double> t(m), b(m);

    t = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    b = {0, 1, 2, 1, 2, 6, 3, 2, 7, 5};

    int n = 3;
    Matrix A(10, 4);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j <= n; j++) {
            A(i, j) = pow(t[i], j);
        }
    }


    Matrix AT = A.transpose();
    Matrix ATA = AT * A;
    Matrix ATA_inv = ATA.inverse();
    Matrix ATb = AT * b;
    Matrix x = ATA_inv * ATb;

    //A.print("A:");
    //ATA.print("A_T*A:");
    //ATA_inv.print("(A_T*A)^-1:");
    //ATb.print("A_T*b:");
    //cout<<"x~:"<<endl;
    /*for(int i=0;i<=n;i++){
        print_vector({x(i, 0)});
    }*/
    fprintf(pipe, "plot '-' using 1:2 title 'My Graph'\n");
    double ca = x(3,0);
    double cb = x(2, 0);
    double cc = x(1,0);
    double cd = x(0, 0);
    for(int i = 0; i < 150; i++){
        double x = 0 + 0.066*i;
        double y = ca*pow(x,3) + cb*pow(x,2) + cc*pow(x,1) + cd*pow(x,0);
        fprintf(pipe, "%f\t%f\n", x, y);
    }
#ifdef WIN32
    _pclose(pipe);
#else
    pclose(pipe);
#endif
    return 0;
}
