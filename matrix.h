#pragma once
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

typedef  double db;
typedef vector<db> dvector;
typedef vector<dvector> dmatrix;


// Вывести вектор
ostream& operator<<(ostream& out, const dvector& v)
{
    for (double el : v) {
        out << right << setw(8) << setprecision(3) << el << ' ';
    }
    return out;
}

// Сложение векторов
dvector operator+(const dvector& a, const dvector& b)
{
    dvector res = a;
    for (int i = 0; i < a.size(); i++)
        res[i] += b[i];
    return res;
}

// Умножение вектора на число
dvector operator*(const dvector& a, const db num)
{
    dvector res = a;
    for (int i = 0; i < a.size(); i++)
        res[i] *= num;
    return res;
}

// Умножение матрицы на вектор
dvector operator*(const dmatrix& a, const dvector& x)
{
    int m = a.size();
    dvector r(m, 0);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            r[i] += a[i][j] * x[j];
    return r;
}

// Умножение матрицы на другую матрицу
dmatrix operator*(const dmatrix& a, const dmatrix& b)
{
    dmatrix res(a.size(), dvector(b[0].size(), 0));
    for (int i = 0; i < a.size(); i++)
        for (int j = 0; j < b[0].size(); j++)
            for (int k = 0; k < b.size(); k++)
                res[i][j] += a[i][k] * b[k][j];
    return res;
}
