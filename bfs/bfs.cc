#include "library.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <random>

Vector operator&(const Vector& A, const Vector& B)
{
    Vector result(A.values.size());

    for(int i = 0; i < result.values.size(); i++)
    {
        result.values[i] = (A.values[i] != 0) && (B.values[i] != 0);
    }

    return result;
}

Vector operator~(const Vector& A)
{
    Vector result(A.values.size());

    for(int i = 0; i < result.values.size(); i++)
    {
        result.values[i] = (A.values[i] == 0);
    }

    return result;
}

// https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7046157
Vector BFS(const CSRMatrix& M, int s)
{
    size_t N = M.rows;
    Vector front(N);
    Vector distances(N);

    for(int i = 0; i < N; i++)
    {
        front.values[i] = 0;
        distances.values[i] = 0;
    }

    front.values[s] = 1;
    distances.values[s] = 1;

    bool cont = true;
    for(int i = 1; cont; i++)
    {
        front = M * front & ~distances;

        cont = false;
        for(int j = 0; j < front.values.size(); j++)
        {
            if(front.values[j])
            {
                cont = true;
                distances.values[j] = i+1;
            }
        }
    }

    return distances;
}

int main(int args, char** argv)
{
    COOMatrix coo_mat = COOMatrix::read(std::cin);

//    coo_mat.print();

    CSRMatrix csr_mat = coo_mat.convert();

//    csr_mat.print();

    std::random_device rd{};
    std::uniform_int_distribution<> dis(0, csr_mat.columns - 1);

    auto start = std::chrono::system_clock::now();

    auto iters = 16;
    for(auto i = 0; i < iters; ++i) {
      auto s = dis(rd);
      Vector result = BFS(csr_mat, s);
      volatile double result_val = result.values.at(0);
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> time = end - start;

    std::cout<<"time passed: "<< time.count() <<"s\n";
/*
    std::cout<<"Page rank result:\n";

    for(int i = 0; i < result.values.size(); i++)
    {
        std::cout<<i+1<<":\t"<<result.values[i]<<"\n";
    }
*/

}
