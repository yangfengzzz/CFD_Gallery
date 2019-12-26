#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>

//const double Ea = 14.;
//const double Q = 14.;
//const double f = 1.;
//const double pgamma = 1.4;
//const double A = 10.;

//const double Ea = 150.;
//const double Q = 50.;
//const double f = 1.8;
//const double pgamma = 1.2;
//const double A = 4000000.;

//const double Ea = 50.;
//const double Q = 50.;
//const double f = 1.8;
//const double pgamma = 1.2;
//double A = 200.;
//const double L_half = 1.;

const double Ea = 20.;
const double Q = 2.;
const double f = 1.1;
const double pgamma = 1.2;
double A = 1.13437e+6;
const double L_half = 1.;

/// below parameters are from Detonatino book.
//const double Ea = 27.0;
//const double Q = 50.;
//const double f = 1.; /// for CJ detonation
//const double pgamma = 1.2;
//double A = 55.;
//const double L_half = 1.;

/// below parameters are from Detonatino book.
//const double Ea = 150.0;
//const double Q = 50.;
//const double f = 1.8;
//const double pgamma = 1.2;
//double A = 1.17021e+07;//60500.;//55.;
//const double L_half = 1.;

/// below parameters are from Detonatino book.
//const double Ea = 27.8;
//const double Q = 50.;
//const double f = 1.; /// for CJ detonation
//const double pgamma = 1.2;
//double A = 56.2272;//60500.;//55.;
//const double L_half = 1.;

const double Pu = 1.;
const double rhou = 1.;
const double uu = 0.;
const double Vu = 1. / rhou;

const double mCJ = sqrt(pgamma * Pu / Vu + (pgamma * pgamma - 1.) * Q / (Vu * Vu) * (1. + sqrt(1. + 2. * pgamma * Pu * Vu / (Q * (pgamma * pgamma - 1.)))));
const double sCJ = (rhou * uu + mCJ) / rhou;
const double s = sqrt(f) * sCJ;
const double m = rhou * (s - uu);

double beta(const double &Y)
{
    return sqrt((m * m * Vu - pgamma * Pu) * (m * m * Vu - pgamma * Pu) - 2. * (pgamma * pgamma - 1.) * m * m * Q * (1. - Y));
}

double Pressure(const double &Y)
{
    return (m * m * Vu + Pu) / (pgamma + 1.) + beta(Y) / (pgamma + 1.);
}

double Vol(const double &Y)
{
    return pgamma * (m * m * Vu + Pu) / (m * m * (pgamma + 1.)) - beta(Y) / (m * m * (pgamma + 1.));
}

double velocity(const double &Y)
{
    return s - m * Vol(Y);
}

double omega(const double &rho, const double &P, const double Y)
{
    double T = P / rho;
    return -A * exp(-Ea / T) * rho * Y;
}

void updateState(const double &Y, double &P, double &u, double &V, double &rho)
{
    P = Pressure(Y);
    V = Vol(Y);
    u = velocity(Y);
    rho = 1. / V;
}

int main(int argc, char *argv[])
{

    double beg_x = 0.;
    double end_x = -1100.;
    double dx = -0.0001;
    const int len = int((end_x - beg_x) / dx) + 1;
    std::vector<double> rho(len);
    std::vector<double> P(len);
    std::vector<double> Y(len);
    std::vector<double> u(len);
    std::vector<double> V(len);
    std::vector<double> the_x(len);

    Y[0] = 1.; /// initial condition
    P[0] = Pu;
    rho[0] = rhou;
    u[0] = uu;
    V[0] = Vu;
    the_x[0] = beg_x;

    //#pragma omp parallel for
    for (int i = 0; i < len; ++i)
    {
        Y[i + 1] = Y[i] - dx * omega(rho[i], P[i], Y[i]) / m;
        updateState(Y[i + 1], P[i + 1], u[i + 1], V[i + 1], rho[i + 1]);
        the_x[i + 1] = the_x[i] + dx;
    }

    /// for calculating A to unify L1/2
    double tmpA = 0.;
    for (int i = 0; i < len - 1; ++i)
    {
        if (Y[i] > 0.5)
            tmpA += L_half * m * 0.5 * (-A / omega(rho[i], P[i], Y[i]) - A / omega(rho[i + 1], P[i + 1], Y[i + 1])) * (Y[i] - Y[i + 1]);
    }

    A = tmpA;

    /// output the results, need to recalculate all solutions, with the updated A, preexponential factor.
    std::fstream file, file1;
    file.open("ZND_configuration_solution.dat", std::ios::in | std::ios::out | std::ios::trunc);

    file << Ea << " " << Q << " " << f << " " << pgamma << " " << A << " " << L_half << " " << s << std::endl;

    file << len << std::endl; /// lengh of the following records
                              //#pragma omp parallel for
    for (int i = 0; i < len; ++i)
    {
        Y[i + 1] = Y[i] - dx * omega(rho[i], P[i], Y[i]) / m;
        updateState(Y[i + 1], P[i + 1], u[i + 1], V[i + 1], rho[i + 1]);
        the_x[i + 1] = the_x[i] + dx;
    }

    for (int i = 0; i < len; ++i)
    {
        file << std::setprecision(10)
             << std::fixed
             << the_x[i] << " "
             << 1. / V[i] << " "
             << u[i] << " "
             << P[i] << " "
             << Y[i]
             << std::endl;
    }

    file.close();

    file1.open("ZND_solution_output.dat", std::ios::in | std::ios::out | std::ios::trunc);

    //#pragma omp parallel for
    for (int i = 0; i < len; ++i)
    {
        Y[i + 1] = Y[i] - dx * omega(rho[i], P[i], Y[i]) / m;
        updateState(Y[i + 1], P[i + 1], u[i + 1], V[i + 1], rho[i + 1]);
        the_x[i + 1] = the_x[i] + dx;
    }

    for (int i = 0; i < len; ++i)
    {
        file1 << std::setprecision(10) << the_x[i] << " " << Y[i] << " " << P[i] << " " << P[i] / rho[i] << " " << omega(rho[i], P[i], Y[i]) << std::endl;
    }

    file1.close();
}
