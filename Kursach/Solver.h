#pragma once
#include <vector>
using namespace std;
struct MatrixProf
{
    int size;
    vector<double> DI;
    vector<double> AL;
    vector<double> AU;
    vector<int> IA;
    vector<int> JA;
};

struct RawMatrix
{
public:
    int N;
    vector<double> DI;
    vector<double> AL;
    vector<int>IA;
    vector<int>JA;
};

class MSG
{
public:
    MSG();
    ~MSG();

    int MaxIterCount = 100000;
    int IterCount = 0;
    double Eps = 1e-15;
    double Difference = 0.0;

    int N = 0;
    vector<double> Ax;
    vector<double> r;
    vector<double> z;
    vector<double> xPrev;

    double DotProduct(vector<double> a, vector<double> b)
    {
        double result = 0.0;
        for (int i = 0; i < a.size(); i++)
            result += a[i] * b[i];

        return result;
    }

    double Norm(vector<double> a)
    {
        double result = 0.0;

        for (size_t i = 0; i < a.size(); i++)
        {
            result += a[i] * a[i];
        }

        return sqrt(result);
    }

    double Error(vector<double> a, vector<double> b)
    {
        double result = 0.0;
        int N = a.size();

        for (int i = 0; i < N; i++)
            result += (a[i] - b[i]) * (a[i] - b[i]);

        return sqrt(result);
    }

    RawMatrix LLT;

    vector <double> MVMultiply(vector<double> x, MatrixProf A)
    {
        vector<double> result(x.size());


        for (int i = 0; i < A.size; i++)
        {
            result[i] = x[i] * A.DI[i];
            for (int k = A.IA[i]; k < A.IA[i + 1]; k++)
            {
                int j = A.JA[k];
                result[i] += A.AL[k] * x[j];
                result[j] += A.AU[k] * x[i];
            }
        }
        return result;
    }

    vector<double> Solve(MatrixProf matrix, std::vector<double> B)
    {
        N = matrix.size;
        InitAuxVectors(N);
        for (size_t i = 0; i < matrix.IA.size(); i++)
        {
            matrix.IA[i] -= 1;
        }
        LLTFactorization(matrix);
        vector<double> x(N);
        Ax = MVMultiply(x, matrix);
        for (int i = 0; i < N; i++)
            r[i] = B[i] - Ax[i];

        Forward(LLT, z, r);
        Backward(LLT, z, z);

        Difference = Norm(r) / Norm(B);

        double dot1 = 0;
        double dot2 = 0;

        dot1 = DotProduct(z, r);

        while (IterCount < MaxIterCount && Difference >= Eps && Error(x, xPrev) > 1.0e-10)
        {
            Ax = MVMultiply(z, matrix);

            double a = dot1 / DotProduct(Ax, z);

            for (int i = 0; i < N; i++)
            {
                xPrev[i] = x[i];
                x[i] += a * z[i];
                r[i] -= a * Ax[i];
            }

            Forward(LLT, Ax, r);
            Backward(LLT, Ax, Ax);

            dot2 = DotProduct(Ax, r);
            double b = dot2 / dot1;
            dot1 = dot2;

            for (int i = 0; i < N; i++)
                z[i] = Ax[i] + b * z[i];

            Difference = Norm(r) / Norm(B);
            IterCount++;
        }

        return x;
    }

    void LLTFactorization(MatrixProf matrix)
    {
        RawMatrix LLT;
        LLT.N = matrix.size;
        LLT.IA.resize(LLT.N + 1);

        for (int i = 0; i < matrix.size + 1; i++)
            LLT.IA[i] = matrix.IA[i];

        LLT.AL.resize(LLT.IA[LLT.N]);
        LLT.JA.resize(LLT.IA[LLT.N]);
        LLT.DI.resize(LLT.N);

        for (int i = 0; i < matrix.IA[matrix.size]; i++)
            LLT.JA[i] = matrix.JA[i];

        for (int i = 0; i < matrix.size; i++)
        {
            double sumD = 0;
            int i0 = matrix.IA[i], i1 = matrix.IA[i + 1];

            for (int k = i0; k < i1; k++)
            {
                double sumL = 0, sumU = 0;
                int j = matrix.JA[k];

                // Calculate L[i][j], U[j][i]
                int j0 = matrix.IA[j], j1 = matrix.IA[j + 1];

                int kl = i0, ku = j0;

                for (; kl < i1 && ku < j1;)
                {
                    int j_kl = matrix.JA[kl];
                    int j_ku = matrix.JA[ku];

                    if (j_kl == j_ku)
                    {
                        sumL += LLT.AL[kl] * LLT.AL[ku];
                        kl++;
                        ku++;
                    }
                    if (j_kl > j_ku)
                        ku++;
                    if (j_kl < j_ku)
                        kl++;
                }

                LLT.AL[k] = (matrix.AL[k] - sumL) / LLT.DI[j];

                // Calculate sum for DI[i]
                sumD += LLT.AL[k] * LLT.AL[k];
            }

            // Calculate DI[i]
            LLT.DI[i] = sqrt(matrix.DI[i] - sumD);
        }

        this->LLT = LLT;
    }

    void InitAuxVectors(int N)
    {
        Ax.resize(N);
        r.resize(N);
        z.resize(N);
        xPrev.resize(N);

        for (int i = 0; i < N; i++)
            xPrev[i] = 1.0;
    }

    void Forward(RawMatrix A, vector<double> &x, vector<double> b)
    {
        vector<double> di = A.DI;
        vector<double> al = A.AL;
        vector<int> ia = A.IA;
        vector<int> ja = A.JA;
        int N = A.N;


        for (int i = 0; i < N; i++)
        {
            double sum = 0;
            int i0 = ia[i], i1 = ia[i + 1];
            for (int k = i0; k < i1; k++)
            {
                int j = ja[k];
                sum += al[k] * x[j];
            }
            x[i] = (b[i] - sum) / di[i];
        }
    }

    void Backward(RawMatrix A, vector<double> &x, vector<double> b)
    {
        vector<double> di = A.DI;
        vector<double> al = A.AL;
        vector<int> ia = A.IA;
        vector<int> ja = A.JA;
        int N = A.N;

        for (int i = 0; i < N; i++)
            x[i] = b[i];

        for (int i = N - 1; i >= 0; i--)
        {
            int i0 = ia[i], i1 = ia[i + 1];
            x[i] /= di[i];
            for (int k = i0; k < i1; k++)
            {
                int j = ja[k];
                x[j] -= al[k] * x[i];

            }
        }
    }

private:

};

MSG::MSG()
{
}

MSG::~MSG()
{
}