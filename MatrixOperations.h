// MatrixOperations.h

#ifndef MATRIXOPERATIONS_H_INCLUDED
#define MATRIXOPERATIONS_H_INCLUDED

#include <cmath>
#include <fstream>

using namespace std;

namespace MatOp
{
    class Matrix
    {
    private:
        bool valid = true;
        int m, n;
        double **M;
        double dt=0;

    public:
        //Getters
        int getRow() const;
        int getColumn() const;
        double** getMatrix() ;
        double getElement(int i, int j) const;
        double getDet() const;
        bool getValid() const;

        //Constructor and destructor
        Matrix();
        Matrix(int a, int b);
        Matrix(const Matrix& mx);
        ~Matrix();

        //Setters
        void setM();
        void setM(double** mat);
        void setSize(int a, int b);
        void setDet(double d);
        void printM();
        void printM(string s);
        void printM(ofstream& fout);
        void setValid(bool v);

        // Operations

        Matrix operator - ();
        Matrix operator + (const Matrix& mx) const;
        Matrix operator - (const Matrix& mx) const;
        Matrix operator * (const Matrix& mx) const;
        Matrix operator / (const Matrix& mx) const;
        bool operator == (const Matrix& mx) const;
        Matrix operator + (double d);
        Matrix operator - (double d);
        Matrix operator * (double d);
        Matrix operator / (double d);
        Matrix& operator=(const Matrix& mx);
        Matrix operator ^ (const Matrix& mx) const; // Matrix multiplication

        // Functions
        Matrix exp();
        Matrix log();
        Matrix T(); // Transpose matrix
        Matrix minor(int a, int b); // Matrix minor
        double det(); // determinant
        Matrix invMat(); // inverse matrix

        // Additional functions
        double sum();
        Matrix sum(int axis); // result (1, m) == row
        Matrix argmax(int axis); // indices of the maximum values along an axis
        Matrix argmin(int axis); // indices of the minimum values along an axis

        // Different types of matrices
        static Matrix eye(int a); // Diagonal matrix
        static Matrix eye(); // Diagonal matrix
        static Matrix zeros(int a, int b); // Matrix with zeros
        static Matrix zeros(); // Matrix with zeros
        static Matrix ones(int a, int b); // Matrix with ones
        static Matrix ones(); // Matrix with ones
        static Matrix norm(double mean, double stddev, int a, int b); // Matrix of random numbers with normal distribution

    };

    Matrix operator + ( double d, Matrix mx); // d + M
    Matrix operator - ( double d, Matrix mx); // d - M
    Matrix operator * ( double d, Matrix mx); // d * M
    Matrix operator / ( double d, Matrix mx); // d / M

}

#endif // MATRIXOPERATIONS_H_INCLUDED
