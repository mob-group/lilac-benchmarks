#include "library.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

void COOMatrix::normalize()
{
    std::vector<int> column_sums(columns);

    for(size_t i = 0; i < values.size(); i++)
    {
        if(values[i].x != values[i].y)
        {
            column_sums[values[i].y-1] += 1;
        }
        else
        {
            values[i] = values.back();
            values.pop_back();
            i--;
        }
    }

    for(int i = 0; i < column_sums.size(); i++)
    {
        if(column_sums[i] == 0)
        {
            for(int j = 0; j < rows; j++)
            {
                if(i != j)
                    values.push_back(Value{j+1, i+1});
            }
        }
    }
}

CSRMatrix COOMatrix::convert()
{
    std::vector<double> column_sums(columns);

    for(int i = 0; i < values.size(); i++)
    {
        column_sums[values[i].y-1] += 1;
    }

    std::vector<std::vector<int>> vectors(rows);

    for(int i = 0; i < values.size(); i++)
    {
        vectors[values[i].x-1].push_back(values[i].y-1);
    }

    for(int i = 0; i < rows; i++)
    {
        std::sort(vectors[i].begin(), vectors[i].end());
    }

    CSRMatrix result;
    result.rows = rows;
    result.columns = columns;
    result.rowstr.push_back(0);

    for(int i = 0; i < rows; i++)
    {
        if(vectors[i].size() == 0)
        {
            for(int j = 0; j < columns; j++)
            {
                result.colidx.push_back(j);
                result.values.push_back(1.0);
            }
        }
        else
        {
            for(int j = 0; j < vectors[i].size(); j++)
            {
                if(j == 0 || vectors[i][j] != vectors[i][j-1])
                {
                    result.colidx.push_back(vectors[i][j]);
                    result.values.push_back(1.0);
                }
                else
                {
                    result.values.back() += 1.0;
                }
            }
        }

        result.rowstr.push_back(result.values.size());
    }

    for(int i = 0; i < result.values.size(); i++)
    {
        result.values[i] /= column_sums[result.colidx[i]];
    }

    return result;
}


void COOMatrix::print() const
{
    std::vector<std::vector<double>> matrix(rows,
        std::vector<double>(columns, 0.0));

    for(int i = 0; i < values.size(); i++)
    {
        matrix[values[i].x-1][values[i].y-1] += 1.0;
    }

    std::cout.precision(1);

    std::cout<<"Begin COO Matrix\n";
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < columns; j++)
        {
            std::cout<<" "<<std::fixed<<matrix[i][j];
        }
        std::cout<<"\n";
    }
    std::cout<<"End COO Matrix\n";
}
void CSRMatrix::print() const
{
    std::vector<std::vector<double>> matrix(rows,
        std::vector<double>(columns, 0.0));

    for(int i = 0; i < rows; i++)
    {
        for(int j = rowstr[i]; j < rowstr[i+1]; j++)
            matrix[i][colidx[j]] = values[j];
    }

    std::cout.precision(3);

    std::cout<<"Begin CSR Matrix\n";
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < columns; j++)
        {
            std::cout<<" "<<std::fixed<<matrix[i][j];
        }
        std::cout<<"\n";
    }
    std::cout<<"End CSR Matrix\n";
}

COOMatrix COOMatrix::read(std::istream& istr)
{
    size_t size;
    COOMatrix result;

    istr>>result.rows>>result.columns>>size;

    result.values.resize(size);

    for(int i = 0; i < size; i++)
    {
        istr>>result.values[i].x>>result.values[i].y;
    }

    return result;
}


Vector Vector::rand(size_t N)
{
    Vector result(N);

    for(int i = 0; i < N; i++)
        result.values[i] = (double)::rand() / (double)RAND_MAX;

    return result;
};

double Vector::l2norm() const {
    double value = 0.0;
    for(int i = 0; i < values.size(); i++)
        value += values[i]*values[i];
    return sqrt(value);
}

double Vector::l1norm() const {
    double value = 0.0;
    for(int i = 0; i < values.size(); i++)
        value += values[i];
    return value;
}

double Vector::average() const {
    double value = 0.0;
    for(int i = 0; i < values.size(); i++)
        value += values[i];
    return value / (double)values.size();
}

Vector operator*(double a, const Vector& b)
{
    Vector result(b.values.size());

    for(int i = 0; i < b.values.size(); i++)
        result.values[i] = a*b.values[i];

    return result;
}

Vector operator+(const Vector& a, const Vector& b)
{
    Vector result(a.values.size());

    for(int i = 0; i < a.values.size(); i++)
        result.values[i] = a.values[i] + b.values[i];

    return result;
}

Vector operator+(const Vector& a, double b)
{
    Vector result(a.values.size());

    for(int i = 0; i < a.values.size(); i++)
        result.values[i] = a.values[i] + b;

    return result;
}

Vector operator-(const Vector& a, const Vector& b)
{
    Vector result(a.values.size());

    for(int i = 0; i < a.values.size(); i++)
        result.values[i] = a.values[i] - b.values[i];

    return result;
}

CSRMatrix operator*(double a, const CSRMatrix& b)
{
    CSRMatrix result = b;

    for(int i = 0; i < b.values.size(); i++)
        result.values[i] *= a;

    return result;
}

Vector operator*(const CSRMatrix& m, const Vector& v)
{
    Vector out(m.rows);

    spmv_harness(&out.values[0], &m.values[0], &v.values[0],
                 &m.rowstr[0], &m.colidx[0], m.rows);

    return out;
}