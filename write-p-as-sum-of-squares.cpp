#include <bits/stdc++.h>
using namespace std;

template <typename T, size_t N, size_t M>
struct Matrix
{
    T x[N][M];

    T *operator[](size_t i) { return x[i]; }

    template <size_t L>
    Matrix<T, N, L> operator*(Matrix<T, M, L> mat)
    {
        Matrix<T, N, L> z;
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < L; ++j)
                for (size_t k = z[i][j] = 0; k < M; ++k)
                    z[i][j] = z[i][j] + x[i][k] * mat[k][j];
        return z;
    }

    Matrix<T, N, M> operator+(Matrix<T, N, M> mat)
    {
        Matrix<T, N, M> z;
        for (size_t i = 0; i < N; i++)
            for (size_t j = 0; j < M; j++)
                z[i][j] = x[i][j] + mat[i][j];
        return z;
    }

    array<T, N> operator*(array<T, M> vec)
    {
        array<T, N> z;
        for (size_t i = 0; i < N; i++)
            for (size_t j = z[i] = 0; j < M; j++)
                z[i] = z[i] + x[i][j] * vec[j];
        return z;
    }
};

template <typename T, size_t N>
Matrix<T, N, N> make_identity()
{
    Matrix<T, N, N> z;
    memset(z.x, 0, sizeof z.x);
    for (size_t i = 0; i < N; i++)
        z[i][i] = 1;
    return z;
}

template <typename T>
T binary_exp(T x, T y, T mod)
{
    T z = 1;

    while (y)
    {
        if (y & 1)
            z = (z * x) % mod;
        x = (x * x) % mod;
        y >>= 1;
    }

    return z;
}

template <typename T>
T legendre(T n, T p) { return binary_exp(n, (p - 1) >> 1, p); }

int main()
{
    int64_t p, m = 2;
    cin >> p;
    assert(!((p - 1) & 3));
    while (1)
    {
        if (legendre(m, p) == p - 1)
            break;
        ++m;
    }
    m = binary_exp(m, (p - 1) >> 2, p);

    int64_t a = p, b = 2 * m, c = (m * m + 1) / p;
    assert(a > 0 && b * b - ((a * c) << 2) < 0);
    Matrix<int64_t, 2, 2> M = make_identity<int64_t, 2>(), N = make_identity<int64_t, 2>();

    Matrix<int64_t, 2, 2> shear, shear_inv, rotation, rotation_inv;
    shear[0][0] = shear[0][1] = shear[1][1] = 1, shear[1][0] = 0;
    shear_inv[0][0] = shear_inv[1][1] = 1, shear_inv[0][1] = -1, shear_inv[1][0] = 0;
    rotation[0][0] = rotation[1][1] = 0, rotation[0][1] = 1, rotation[1][0] = -1;
    rotation_inv[0][0] = rotation_inv[1][1] = 0, rotation_inv[0][1] = -1, rotation_inv[1][0] = 1;

    cout << "Reducing binary quadratic form (" << a << ", " << b << ", " << c << ")\n";

    while (!((-a < b && b <= a && a < c) || (0 <= b && b <= a && a == c)))
    {
        int64_t a_, b_, c_;
        if (a > c || (a == c && b < 0))
        {
            a_ = c, b_ = -b, c_ = a;
            M = M * rotation;
            N = rotation_inv * N;
        }
        else
        {
            if (a < c)
            {
                if (b <= -a)
                {
                    a_ = a, b_ = b + 2 * a, c_ = c + b + a;
                    M = M * shear;
                    N = shear_inv * N;
                }
                else
                {
                    a_ = a, b_ = b - 2 * a, c_ = c - b + a;
                    M = M * shear_inv;
                    N = shear * N;
                }
            }
            else
            {
                a_ = a, b_ = b - 2 * a, c = c - b + a;
                M = M * shear_inv;
                N = shear * N;
            }
        }

        swap(a, a_);
        swap(b, b_);
        swap(c, c_);

        cout << a << ' ' << b << ' ' << c << '\n';
    }

    cout << "matrix:\n"
         << M[0][0] << ' ' << M[0][1] << '\n'
         << M[1][0] << ' ' << M[1][1] << '\n'
         << "matrix inverse:\n"
         << N[0][0] << ' ' << N[0][1] << '\n'
         << N[1][0] << ' ' << N[1][1] << '\n'
         << p << " = " << abs(N[0][0]) << "^2 + " << abs(N[1][0]) << "^2\n";
}