#include "library.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>

Vector PR(const CSRMatrix& M, double d, int iterations)
{
    size_t N = M.rows;
    Vector v = Vector::rand(N);
    v = 1.0/v.l1norm() * v;

    Vector last_v;
    volatile double error;
    for(int i = 0; i < iterations; i++)
    {
        last_v = v;
        v = (d * M) * v + ((1 - d)) * v.average();

        error = (v - last_v).l2norm();
    }

    return v;
}

int main(int args, char** argv)
{
    COOMatrix coo_mat = COOMatrix::read(std::cin);

//    coo_mat.print();
    coo_mat.normalize();
//    coo_mat.print();

    CSRMatrix csr_mat = coo_mat.convert();

//    csr_mat.print();

    clock_t start_time = clock();
    Vector result = PR(csr_mat, 0.85, 100000);
    clock_t end_time = clock();

    double lapsed_secs = (double)(end_time-start_time)/(double)CLOCKS_PER_SEC;

    std::cout<<"time passed: "<<lapsed_secs<<"s\n";
/*
    std::cout<<"Page rank result:\n";

    for(int i = 0; i < result.values.size(); i++)
    {
        std::cout<<i+1<<":\t"<<result.values[i]<<"\n";
    }
*/

}
