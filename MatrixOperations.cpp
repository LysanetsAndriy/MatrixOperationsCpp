// MatrixOperations.cpp

#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <chrono>
#include <fstream>

#include "MatrixOperations.h"

using namespace std;

namespace MatOp
{
    int Matrix::getRow() const {return m;}
    int Matrix::getColumn() const {return n;}
    double** Matrix::getMatrix() {return M;}
    double Matrix::getElement(int i, int j) const {return M[i][j];}
    double Matrix::getDet() const {return dt;}
    bool Matrix::getValid() const {return valid;}



    Matrix::Matrix(){}

    Matrix::Matrix(const Matrix& mx)
    {
        m = mx.getRow();
        n = mx.getColumn();
        valid = mx.getValid();
        dt =  mx.getDet();

        M = new double*[m];
        for(int i = 0; i < m; i++) {
            M[i] = new double[n];
            for(int j = 0; j < n; j++) {
                M[i][j] = mx.getElement(i, j);
            }
        }
    }

    Matrix::Matrix(int a, int b)
    {
        if (a > 0 && b > 0)
        {
            m = a; n = b;
            M = new double*[m];
            for(int i = 0; i < m; i++)
            {
                M[i] = new double[n];
                for(int j = 0; j < n; j++)
                {
                    M[i][j] = 0;
                }
            }
        }
        else
        {
            cout << "Wrong size!\n\n";
            valid = false;
        }
    }

    Matrix::~Matrix()
    {
        for (int i = 0; i < this->m; i++)
        {
            delete[] this->M[i];
        }
        delete[] this->M;
    }

    void Matrix::setM()
    {
        cin >> m >> n;
        if (m <= 0 || n <= 0)
        {
            cout << "Wrong size!\n\n";
            valid = false;
        }
        else
        {
            this->M = new double*[m];
            for(int i = 0; i < m; i++)
                this->M[i] = new double[n];

            for(int i = 0; i < m; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    cin >> this->M[i][j];
                }
            }
        }
        cout << '\n';
    }

    void Matrix::setSize(int a, int b)
    {
        if (a <= 0 || b <= 0)
        {
            cout << "Wrong size!\n\n";
            valid = false;
        }
        else
        {
            m = a;
            n = b;

            this->M = new double*[m];
            for(int i = 0; i < m; i++)
                M[i] = new double[n];
        }

    }

    void Matrix::setM(double** mat)
    {
        for(int i = 0; i < m; i++)
        {
            for(int j = 0; j < n; j++)
            {
                M[i][j] = mat[i][j];
            }
        }

    }

    void Matrix::setDet(double d)
    {
        dt = d;
    }

    void Matrix::setValid(bool b)
    {
        valid = b;
    }

    void Matrix::printM()
    {
        if (this->valid)
        {
            for(int i = 0; i < this->m; i++)
            {
                for(int j = 0; j < this->n; j++)
                {
                    cerr <<  this->M[i][j] << ' ';
                }
                cerr << '\n';
            }
            cerr << '\n';
        }
        else
            cerr << "This matrix is empty!\n\n";
    }

    void Matrix::printM(string s)
    {
        ofstream fout(s);
        if (this->valid)
        {
            for(int i = 0; i < this->m; i++)
            {
                for(int j = 0; j < this->n; j++)
                {
                    fout <<  this->M[i][j] << ' ';
                }
                fout << '\n';
            }
            fout << '\n';
        }
        else
            fout << "This matrix is empty!\n\n";
    }

    void Matrix::printM(ofstream& fout)
    {
        if (this->valid)
        {
            for(int i = 0; i < this->m; i++)
            {
                for(int j = 0; j < this->n; j++)
                {
                    fout <<  this->M[i][j] << ' ';
                }
                fout << '\n';
            }
            fout << '\n';
        }
        else
            fout << "This matrix is empty!\n\n";
    }

    Matrix Matrix::operator - ()
    {
        if (!this->valid)
        {
            cout << "Matrices are not valid!\n\n";
            return Matrix();
        }
        Matrix tmp(this->m, this->n);
        for(int i = 0; i < this->m; i++)
        {
            for(int j = 0; j < this->n; j++)
            {
                tmp.M[i][j] = -(this->M[i][j]);
            }
        }
        return tmp;
    }

    Matrix Matrix::operator + (const Matrix& mx) const
    {
        if (!valid || !mx.valid)
        {
            cout << "Matrices are not valid!\n\n";
            return Matrix();
        }

        if (mx.getRow() == m && mx.getColumn() == 1)
        {
            Matrix tmp(m, n);
            for(int i = 0; i < m; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    tmp.M[i][j] = M[i][j] + mx.M[i][0];
                }
            }
            return tmp;
        }
        else
        {
            if (mx.getRow() != m || mx.getColumn() != n)
            {
                cout << "Can't add! The dimensions of the matrices are not the same! ";
                cout << "(" << m << ", " << n << ") != ";
                cout << "(" << mx.getRow() << ", " << mx.getColumn() << ")";
            }
            else
            {
                Matrix tmp(m, n);
                for(int i = 0; i < m; i++)
                {
                    for(int j = 0; j < n; j++)
                    {
                        tmp.M[i][j] = M[i][j] + mx.M[i][j];
                    }
                }
                return tmp;
            }
        }

    }

    Matrix Matrix::operator - (const Matrix& mx) const
    {
        if (!valid || !mx.valid)
        {
            cout << "Matrices are not valid!\n\n";
            return Matrix();
        }
        if (mx.getRow() == m && mx.getColumn() == 1)
        {
             Matrix tmp(m, n);
             for(int i = 0; i < m; i++)
             {
                 for(int j = 0; j < n; j++)
                 {
                     tmp.M[i][j] = M[i][j] - mx.M[i][0];
                 }
             }
             return tmp;
        }
        else
            if (mx.getRow() != m || mx.getColumn() != n)
            {
                cout << "Can't substruct! The dimensions of the matrices are not the same! ";
                cout << "(" << m << ", " << n << ") != ";
                cout << "(" << mx.getRow() << ", " << mx.getColumn() << ")";
            }
            else
            {
                Matrix tmp(m, n);
                for(int i = 0; i < m; i++)
                {
                    for(int j = 0; j < n; j++)
                    {
                        tmp.M[i][j] = M[i][j] - mx.M[i][j];
                    }
                }
                return tmp;
            }
    }

    Matrix Matrix::operator * (const Matrix& mx) const
    {
        if (!valid || !mx.valid)
        {
            cout << "Matrices are not valid!\n\n";
            return Matrix();
        }
        if (mx.getRow() != m || mx.getColumn() != n)
        {
            cout << "Can't multiply! The dimensions of the matrices are not the same! ";
            cout << "(" << m << ", " << n << ") != ";
            cout << "(" << mx.getRow() << ", " << mx.getColumn() << ")";
        }
        else
        {
            Matrix tmp(m, n);
            for(int i = 0; i < m; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    tmp.M[i][j] = M[i][j] * mx.M[i][j];
                }
            }
            return tmp;
        }
    }

    Matrix Matrix::operator / (const Matrix& mx) const
    {
        if (!valid || !mx.valid)
        {
            cout << "Matrices are not valid!\n\n";
            return Matrix();
        }
        if (mx.getRow() != m || mx.getColumn() != n)
        {
            cout << "Can't divide! The dimensions of the matrices are not the same! ";
            cout << "(" << m << ", " << n << ") != ";
            cout << "(" << mx.getRow() << ", " << mx.getColumn() << ")";
        }
        else
        {
            Matrix tmp(m, n);
            try
            {
                for(int i = 0; i < m; i++)
                {
                    for(int j = 0; j < n; j++)
                    {
                        if (mx.M[i][j] != 0)
                            tmp.M[i][j] = M[i][j] / mx.M[i][j];
                        else
                            throw 0;
                    }
                }
                return tmp;
            }
            catch (int)
            {
                cout << "Division by zero!\n";
                tmp.valid = false;
                return tmp;
            }
        }
    }

    bool Matrix::operator == (const Matrix& mx) const
    {
        if (!valid || !mx.valid)
        {
            cout << "Matrices are not valid!\n\n";
            return false;
        }
        if (mx.getRow() != m || mx.getColumn() != n)
        {
            cout << "The dimensions of the matrices are not the same! ";
            cout << "(" << m << ", " << n << ") != ";
            cout << "(" << mx.getRow() << ", " << mx.getColumn() << ")";
        }
        else
        {
            for(int i = 0; i < m; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    if (M[i][j] != mx.M[i][j])
                        return false;
                }
            }
            return true;
        }
    }

    Matrix Matrix::operator + (double d) // M + d
    {
        if (!valid)
        {
            cout << "Matrix is not valid!\n\n";
            return Matrix();
        }
        Matrix tmp(m, n);
        for(int i = 0; i < m; i++)
        {
            for(int j = 0; j < n; j++)
            {
                tmp.M[i][j] = M[i][j] + d;
            }
        }
        return tmp;
    }

    Matrix Matrix::operator - (double d) // M - d
    {
        if (!valid)
        {
            cout << "Matrix is not valid!\n\n";
            return Matrix();
        }
        Matrix tmp(m, n);
        for(int i = 0; i < m; i++)
        {
            for(int j = 0; j < n; j++)
            {
                tmp.M[i][j] = M[i][j] - d;
            }
        }
        return tmp;
    }

    Matrix Matrix::operator * (double d) // M * d
    {
        if (!valid)
        {
            cout << "Matrix is not valid!\n\n";
            return Matrix();
        }
        Matrix tmp(m, n);
        for(int i = 0; i < m; i++)
        {
            for(int j = 0; j < n; j++)
            {
                tmp.M[i][j] = M[i][j] * d;
            }
        }
        return tmp;
    }

    Matrix Matrix::operator / (double d) // M / d
    {
        if (!valid)
        {
            cout << "Matrix is not valid!\n\n";
            return Matrix();
        }
        Matrix tmp(m, n);
        if ((d == 0))
        {
            cout << "Division by zero!\n\n";
            tmp.setM(M);
            return tmp;
        }
        else
        {
            for(int i = 0; i < m; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    tmp.M[i][j] = M[i][j] / d;
                }
            }
            return tmp;
        }
    }

    Matrix& Matrix::operator = (const Matrix& mx)
    {
        if (valid)
        {
                dt = mx.getDet();
                m = mx.getRow();
                n = mx.getColumn();
                M = new double*[m];
                for (int i = 0; i < m; i++)
                {
                    M[i] = new double[n];
                    for (int j = 0; j < n; j++)
                    {
                        M[i][j] = mx.getElement(i, j);
                    }
                }
                return *this;
        }
        else
            cout << "This matrix is empty!\n\n";

    }

    Matrix Matrix::operator ^ (const Matrix& mx) const // Matrix multiplication (dot product)
    {
        if (!valid || !mx.valid)
        {
            cerr << "Matrices are not valid!\n\n";
            return Matrix();
        }
        // first matrix size: (m,n); second: (k, h)
        int k, h;
        k = mx.getRow(); h = mx.getColumn();
        Matrix tmp(m, h);
        if ((n != k))
        {
            cerr << "Matrices can't be multiplied! ";
            cerr << "(" << m << ", " << n << "); ";
            cerr << "(" << k << ", " << h << "); ";
            cerr << "(" << n << " != " << k << ")\n\n";
            tmp.valid = false;
            return tmp;
        }
        else
        {
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < h; j++)
                {
                    tmp.M[i][j] = 0;
                    for (int p = 0; p < n; p++)
                    {
                        tmp.M[i][j] = tmp.M[i][j] + M[i][p] * mx.M[p][j];
                    }
                }
            }
            return tmp;

        }
    }

    Matrix Matrix::log()
    {
        if (!this->valid)
        {
            cout << "Matrices are not valid!\n\n";
            return Matrix();
        }
        Matrix tmp(this->m, this->n);
        for(int i = 0; i < this->m; i++)
        {
            for(int j = 0; j < this->n; j++)
            {
                tmp.M[i][j] = std::log(this->M[i][j]);
            }
        }
        return tmp;
    }

    Matrix Matrix::exp()
    {
        if (!this->valid)
        {
            cout << "Matrices are not valid!\n\n";
            return Matrix();
        }
        Matrix tmp(this->m, this->n);
        for(int i = 0; i < this->m; i++)
        {
            for(int j = 0; j < this->n; j++)
            {
                tmp.M[i][j] = std::exp(this->M[i][j]);
            }
        }
        return tmp;
    }

    Matrix Matrix::T() // Transpose matrix
    {
        if (valid)
        {
            Matrix tmp(n, m);
            for(int i = 0; i < n; i++)
            {
                for(int j = 0; j < m; j++)
                {
                    tmp.M[i][j] = M[j][i];
                }
            }
            return tmp;
        }
        else
            cout << "This matrix is empty!\n\n";
    }

    Matrix Matrix::minor(int a, int b)
    {
        if (valid)
        {
            if(a+1 > m || b+1 > n)
            {
                cout << "Wrong value of a or b!\n";
                return *this;
            }
            Matrix tmp(m-1, n-1);
            int x = 0, y = 0;
            for(int i = 0; i < m; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    if (i != a && j != b)
                    {
                        tmp.M[x][y] = M[i][j];
                        y++;
                        if (y == n-1)
                        {
                            y = 0;
                            x++;
                        }
                    }
                }
            }
            return tmp;
        }
        else
            cout << "This matrix is empty!\n\n";
    }

    double Matrix::det()
    {
        if (valid)
        {
            if(m == n)
            {
                double d = getDet();
                double **A = getMatrix();
                if (n == 1)
                    return A[0][0];
                for(int i = 0; i < m; i++)
                {
                    if (i % 2 == 0)
                        d += (A[0][i]*minor(0,i).det());
                    else
                        d += (-A[0][i]*minor(0,i).det());
                }
                return d;
            }
            else
            {
                cout << "Wrong matrix size!\n";
                return 0;
            }
        }
        else
            cout << "This matrix is empty!\n\n";
    }

    Matrix Matrix::invMat()
    {
        if (valid)
        {
            if (this->det() != 0)
            {
                Matrix tmp(m, n);
                double **A = tmp.getMatrix();
                for(int i = 0; i < m; i++)
                {
                    for(int j = 0; j < n; j++)
                    {
                        if ((i + j + 2) % 2 == 0)
                            A[i][j] = this->minor(i,j).det();
                        else
                            A[i][j] = -this->minor(i,j).det();
                    }
                }
                tmp.setM(A);
                tmp = tmp / this->det();
                return tmp.T();
            }
            else
            {
                cout << "This matrix doesn't have inverse matrix(det = 0)!";
            }
        }
        else
        {
            cout << "This matrix is empty!\n\n";
        }
    }

    double Matrix::sum()
    {
        int m = this->getRow();
        int n = this->getColumn();
        double **A = this->getMatrix();
        double sm = 0;
        for(int i = 0; i < m; i++)
        {
            for(int j = 0; j < n; j++)
            {
                sm += A[i][j];
            }
        }
        return sm;
    }

    Matrix Matrix::sum(int axis)
    {
        int m = this->getRow();
        int n = this->getColumn();
        double **A = this->getMatrix();

        if (axis == 0)
        {
            Matrix sm = Matrix::zeros(1, n);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    sm.M[0][j] += A[i][j];
                }
            }
            return sm;
        }
        if (axis == 1)
        {
            Matrix sm = Matrix::zeros(1, m);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    sm.M[0][i] += A[i][j];
                }
            }
            return sm;
        }
    }

    Matrix Matrix::argmax(int axis)
    {
        int m = this->getRow();
        int n = this->getColumn();
        double **A = this->getMatrix();

        if (axis == 0)
        {
            Matrix mx = Matrix::zeros(1, n) - 100000;
            Matrix mxi = Matrix::zeros(1, n);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (mx.M[0][j] < A[i][j])
                    {
                        mx.M[0][j] = A[i][j];
                        mxi.M[0][j] = i;
                    }
                }
            }
            return mxi;
        }
        if (axis == 1)
        {
            Matrix mx = Matrix::zeros(1, m) - 1000000;
            Matrix mxi = Matrix::zeros(1, m);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (mx.M[0][i] < A[i][j])
                    {
                        mx.M[0][i] = A[i][j];
                        mxi.M[0][i] = j;
                    }
                }
            }
            return mxi;
        }
    }

    Matrix Matrix::argmin(int axis)
    {
        int m = this->getRow();
        int n = this->getColumn();
        double **A = this->getMatrix();

        if (axis == 0)
        {
            Matrix mn = Matrix::zeros(1, n) + 100000;
            Matrix mni = Matrix::zeros(1, n);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (mn.M[0][j] > A[i][j])
                    {
                        mn.M[0][j] = A[i][j];
                        mni.M[0][j] = i;
                    }
                }
            }
            return mni;
        }
        if (axis == 1)
        {
            Matrix mn = Matrix::zeros(1, m) + 1000000;
            Matrix mni = Matrix::zeros(1, m);
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                     if (mn.M[0][i] > A[i][j])
                    {
                        mn.M[0][i] = A[i][j];
                        mni.M[0][i] = j;
                    }
                }
            }
            return mni;
        }
    }

    Matrix Matrix::eye(int a)
    {
        if (a > 0)
        {
            Matrix tmp(a,a);
            double **A = tmp.getMatrix();
            for(int i = 0; i < a; i++)
            {
                for(int j = 0; j < a; j++)
                {
                    if(i == j)
                        A[i][j] = 1;
                    else
                        A[i][j] = 0;
                }
            }
            tmp.setM(A);
            return tmp;
        }
        else
        {
            cout << "Wrong size!\n\n";
        }
    }

    Matrix Matrix::eye()
    {
        int a;
        cin >> a;
        if (a > 0)
        {
            Matrix tmp(a,a);
            double **A = tmp.getMatrix();
            for(int i = 0; i < a; i++)
            {
                for(int j = 0; j < a; j++)
                {
                    if(i == j)
                        A[i][j] = 1;
                    else
                        A[i][j] = 0;
                }
            }
            tmp.setM(A);
            return tmp;
        }
        else
        {
            cout << "Wrong size!\n\n";
        }
    }

    Matrix Matrix::zeros(int a, int b)
    {
        if (a > 0 && b > 0)
        {
            Matrix tmp(a,b);
            //double **A = tmp.getMatrix();
            for(int i = 0; i < a; i++)
            {
                for(int j = 0; j < b; j++)
                {
                    tmp.M[i][j] = 0;
                }
            }
            //tmp.setM(A);
            return tmp;
        }
        else
        {
            cout << "Wrong size!\n\n";
        }
    }

    Matrix Matrix::zeros()
    {
        int a, b;
        cin >> a >> b;
        if (a > 0 && b > 0)
        {
            Matrix tmp(a,b);
            double **A = tmp.getMatrix();
            for(int i = 0; i < a; i++)
            {
                for(int j = 0; j < b; j++)
                {
                    A[i][j] = 0;
                }
            }
            tmp.setM(A);
            return tmp;
        }
        else
        {
            cout << "Wrong size!\n\n";
        }
    }

    Matrix Matrix::ones(int a, int b)
    {
        if (a > 0 && b > 0)
        {
            Matrix tmp(a,b);
            double **A = tmp.getMatrix();
            for(int i = 0; i < a; i++)
            {
                for(int j = 0; j < b; j++)
                {
                    A[i][j] = 1;
                }
            }
            tmp.setM(A);
            return tmp;
        }
        else
        {
            cout << "Wrong size!\n\n";
        }
    }

    Matrix Matrix::ones()
    {
        int a, b;
        cin >> a >> b;
        if (a > 0 && b > 0)
        {
            Matrix tmp(a,b);
            double **A = tmp.getMatrix();
            for(int i = 0; i < a; i++)
            {
                for(int j = 0; j < b; j++)
                {
                    A[i][j] = 1;
                }
            }
            tmp.setM(A);
            return tmp;
        }
        else
        {
            cout << "Wrong size!\n\n";
        }
    }

    Matrix Matrix::norm(double mu, double sigma, int a, int b)
    {
        if (a > 0 && b > 0)
        {
            static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            static std::default_random_engine gen(seed);
            std::normal_distribution<double> dist(mu, sigma);
            Matrix tmp(a,b);
            double **A = tmp.getMatrix();
            for(int i = 0; i < a; i++)
            {
                for(int j = 0; j < b; j++)
                {

                    A[i][j] = dist(gen);
                }
            }
            tmp.setM(A);
            return tmp;
        }
        else
        {
            cout << "Wrong size!\n\n";
        }

    }

    Matrix operator + ( double d, Matrix mx) // d + M
    {
        return mx + d;
    }

    Matrix operator - ( double d, Matrix mx) // d - M
    {
        return ((-mx) + d);
    }

    Matrix operator * ( double d, Matrix mx) // d * M
    {
        return mx * d;
    }

    Matrix operator / ( double d, Matrix mx) // d / M
    {
        if (!(mx.getValid()))
        {
            cout << "Matrix is not valid!\n\n";
            return Matrix();
        }
        int r = mx.getRow();
        int c = mx.getColumn();
        Matrix tmp(r, c);
        double **A = mx.getMatrix();
        for(int i = 0; i < r; i++)
        {
            for(int j = 0; j < c; j++)
            {
                A[i][j] = d / A[i][j];
            }
        }
        tmp.setM(A);
        return tmp;
    }

}
