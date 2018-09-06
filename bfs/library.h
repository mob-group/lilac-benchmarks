#ifndef _LIBRARY_
#define _LIBRARY_
#include <iostream>
#include <vector>

extern "C" void spmv_harness_(double* ov, double* a, double* iv,
                   int* rowstr, int* colidx, int* rows);

class CSRMatrix
{
public:
    int                 rows;
    int                 columns;
    std::vector<int>    rowstr;
    std::vector<int>    colidx;
    std::vector<double> values;

    void print() const;

    int get_cols() const { return columns; }
};

class COOMatrix
{
public:
    struct Value { int x; int y; };

    static COOMatrix read(std::istream&);

    int                rows;
    int                columns;
    std::vector<Value> values;

    void print() const;

    void normalize();

    CSRMatrix convert();
};

class Vector
{
public:
    Vector(size_t N=0) : values(N) { }

    static Vector rand(size_t N);

    double l1norm() const;
    double l2norm() const;
    double average() const;

    std::vector<double> values;
};

Vector operator*(double a, const Vector& b);

Vector operator+(const Vector& a, const Vector& b);

Vector operator+(const Vector& a, double b);

Vector operator-(const Vector& a, const Vector& b);

CSRMatrix operator*(double a, const CSRMatrix& b);

Vector operator*(const CSRMatrix& m, const Vector& v);

#endif
