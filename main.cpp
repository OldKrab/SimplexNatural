#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <cfloat>
#include <cassert>
#include "matrix.h"

using namespace std;

#define FILE_IN R"(D:\source\clion\igri7\in.txt)"

const db EPS = 1e-5;

void fillVector(vector<int> &v, int from, int to) {
    for (int i = from; i < to; i++)
        v[i - from] = i;
}

class NaturalSimplex {
public:
    /// Ввод исходных данных
    void Input() {
        ifstream fin(FILE_IN);
        isMin = InputOptimum(fin);
        c = InputTarget(fin);
        InputConds(fin);

        db temp;
        fin >> temp;
        if (!fin.fail())
            throw invalid_argument("В исходных данных лишние данные!");
    }

    /// Решить задачу
    void Solve() {
        CreateTable();
        SolveSystemWithEps();
        /// Удалим фи-функцию и поставим на ее место исходную
        table[0] = table[m + 1];
        table.pop_back();
        SolveSystem();
        result = table[0][0] * (-isMin);
        AddX();

        HandleZeroCoefs();
    }

    /// Вывести результаты
    void PrintResults() {
        cout << "F " << (isMin == 1 ? "min" : "max") << " = " << result << endl;
        for (int i = 0; i < Xs.size(); i++) {
            cout << "X_" << i + 1 << " = (" << Xs[i][0];
            for (int j = 1; j < Xs[i].size(); j++)
                cout << ", " << Xs[i][j];
            cout << ")\n";
        }
        cout << "X = a_1*X_1";
        for (int i = 1; i < Xs.size(); i++) {
            cout << " + a_" << i + 1 << "*X_" << i + 1;
        }
        cout << ",\n0 <= a_i <= 1, sum(a_i) = 1, i = 1, " << (Xs.size() > 2 ? "..., " : "") << Xs.size() << ".";
        cout << "\nУпрощенная форма:\n";
        for (int i = 0; i < Xs.size(); i++) {
            cout << "X_" << i + 1 << " = (" << Xs[i][0];
            for (int j = 1; j < mainVarsCnt; j++)
                cout << ", " << Xs[i][j];
            cout << ")\n";
        }
    }

private:
    /// Обработать нулевые коэфы в целевой функции
    void HandleZeroCoefs() {
        auto tableCopy = table;
        auto basisCopy = basisVars, nonBasisCopy = nonBasisVars;
        for (int i = 1; i < table[0].size(); i++)
            if (abs(table[0][i]) < EPS) {
                int r = GetPermissiveRow(i);
                if (r != -1) {
                    SwapVars(i, r);
                    AddX();
                }
                table = tableCopy;
                basisVars = basisCopy;
                nonBasisVars = nonBasisCopy;
            }
    }

    /// Добавить результат их текущей таблицы
    void AddX() {
        dvector x(varsCnt, 0);
        for (int i = 0; i < basisVars.size(); i++)
            x[basisVars[i]] = table[i + 1][0];
        Xs.push_back(x);
    }

    /// Решить исходную систему
    void SolveSystem() {
        int k = GetPermissiveCol();
        int r = GetPermissiveRow(k);
        while (k != -1) {
            if (r == -1)
                throw runtime_error("Система не имеет конечного оптимума!");
            SwapVars(k, r);
            k = GetPermissiveCol();
            r = GetPermissiveRow(k);
        }
    }

    /// Решить систему для фи-функции
    void SolveSystemWithEps() {
        int k = GetPermissiveCol();
        int r = GetPermissiveRow(k);
        while (k != -1 && epsCnt) {
            if (r == -1)
                throw runtime_error("Система не имеет конечного оптимума!");
            SwapVars(k, r);
            if (nonBasisVars[k - 1] < 0) {
                DeleteCol(k);
                epsCnt--;
            }
            k = GetPermissiveCol();
            r = GetPermissiveRow(k);
        }
        if (epsCnt || (abs(table[0][0]) > EPS)) {
            throw runtime_error("Система не имеет допустимого значения!");
        }
    }

    /// Поменять базисную переменную со свободной
    /// \param k-ый столбец разрешающий
    /// \param r-ая строка разрешающая
    void SwapVars(int k, int r) {
        // Преобразуем разрешающий элемент
        table[r][k] = 1. / table[r][k];
        // Преобразуем разрешающий столбец
        for (int i = 0; i < table.size(); i++)
            if (i != r)
                table[i][k] *= -table[r][k];
        // Преобразуем остальные элементы
        for (int i = 0; i < table.size(); i++)
            for (int j = 0; j < table[0].size(); j++)
                if (i != r && j != k)
                    table[i][j] += table[r][j] * table[i][k];
        // Преобразуем разрешающую строку
        for (int j = 0; j < table[r].size(); j++)
            if (j != k)
                table[r][j] *= table[r][k];

        swap(basisVars[r - 1], nonBasisVars[k - 1]);
    }

    /// Удалить стоблец из таблицы
    void DeleteCol(int col) {
        for (auto &row : table)
            row.erase(row.begin() + col);
        nonBasisVars.erase(nonBasisVars.begin() + col - 1);
    }

    /// Получить разрешающий столбец
    int GetPermissiveCol() {
        auto L = table[0];
        int col = -1;
        db mn = -EPS;
        for (int i = 1; i < L.size(); i++)
            if (L[i] < mn) {
                mn = L[i];
                col = i;
            }
        return col;
    }

    /// Получить разрешающую строку
    int GetPermissiveRow(int k) {
        int row = -1;
        db min_value = DBL_MAX;
        for (int i = 1; i < m + 1; i++)
            if (table[i][k] > EPS && table[i][0] / table[i][k] < min_value) {
                row = i;
                min_value = table[i][0] / table[i][k];
            }
        return row;
    }

    /// Создать симплекс таблицу
    void CreateTable() {
        c = c * isMin;
        FixNegativeB();

        basisVars.resize(m);
        nonBasisVars.resize(n);
        fillVector(nonBasisVars, 0, n);
        varsCnt = n;

        for (int i = 0; i < m; i++)
            HandleCond(i, signs[i]);

        table.resize(m + 2);
        table[m + 1] = c;
        for (int i = 1; i <= m; i++) {
            table[i].resize(n + 1);
            table[i][0] = b[i - 1];
            for (int j = 1; j <= n; j++)
                table[i][j] = a[i - 1][j - 1];
        }
        table[0].resize(n + 1, 0);
        for (int i = 0; i < m; i++)
            if (basisVars[i] < 0)
                table[0] = table[0] + table[i + 1];
        table[0] = table[0] * -1;
    }

    /// Обработать условие: неравенство привести к равенству,
    /// добавить добавочные и искусственные переменные
    void HandleCond(int row, const string &sign) {
        if (sign == ">=") {
            AddNonBasisVar(row);
            basisVars[row] = -epsCnt - 1;
            epsCnt++;
        } else if (sign == "=") {
            basisVars[row] = -epsCnt - 1;
            epsCnt++;
        } else
            basisVars[row] = varsCnt++;
    }

    /// Добавить свободную переменную
    void AddNonBasisVar(int row) {
        for (int i = 0; i < m; i++)
            if (i == row)
                a[i].push_back(-1);
            else
                a[i].push_back(0);
        nonBasisVars.push_back(varsCnt++);
        c.push_back(0);
        n++;
    }

    /// Исправить отрицательные свободные члены в условиях
    void FixNegativeB() {
        for (int i = 0; i < m; i++)
            if (b[i] < 0) {
                a[i] = a[i] * -1;
                b[i] *= -1;
                if (signs[i] == "<=")
                    signs[i] = ">=";
                else if (signs[i] == ">=")
                    signs[i] = "<=";
            }
    }

    /// Ввод оптимума
    /// \return isMin
    static int InputOptimum(ifstream &fin) {
        string type;
        fin >> type;
        transform(type.begin(), type.end(), type.begin(), ::tolower);
        int isMin = 1;
        if (type == "max")
            isMin = -1;
        else if ((type) != "min")
            throw invalid_argument("Оптимум (max/min) задан неверно!");
        return isMin;
    }

    /// Ввод целевой функции
    dvector InputTarget(ifstream &fin) {
        fin >> n;
        mainVarsCnt = n;
        c.resize(n + 1);
        for (int i = 0; i < n + 1; i++)
            fin >> c[i];
        if (fin.fail())
            throw invalid_argument("Функция цели задана неверно!");
        return c;
    }

    /// Ввод условий
    void InputConds(ifstream &fin) {
        vector<string> allowsSigns = {"=", "<=", ">="};
        fin >> m;
        b.resize(m);
        a.resize(m);
        signs.resize(m);
        for (int i = 0; i < m; i++) {
            a[i].resize(n);
            for (int j = 0; j < n; j++)
                fin >> a[i][j];
            fin >> signs[i];
            fin >> b[i];
            if (fin.fail())
                throw invalid_argument("Коэффициенты условий заданы неверно!");
            if (find(allowsSigns.begin(), allowsSigns.end(), signs[i]) == allowsSigns.end())
                throw invalid_argument("Знаки условий заданы неверно!");
        }
    }

    int n = 0, m = 0, varsCnt = 0, epsCnt = 0, mainVarsCnt = 0;
    dmatrix table;
    vector<int> basisVars, nonBasisVars;

    dmatrix a;
    dvector b, c;
    int isMin = 1;
    vector<string> signs;

    dmatrix Xs;
    db result = 0;
};


int main() {
    setlocale(LC_ALL, "rus");
    NaturalSimplex ns;
    try {
        ns.Input();
        ns.Solve();
    } catch (exception &e) {
        cout << e.what();
        return -1;
    }
    ns.PrintResults();
    return 0;
}
