#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
using namespace std;

//break when any float is NaN


#define pi 3.14159265358979323846
#define nmax 2000
#define mmax 200

int n;
double xme, xmh, xm, nbeta, xn;
double omega, dam;
double xnue, xnuh;
vector<double> s(nmax), w(nmax);
double xieh;
vector<double> selfe(nmax), selfh(nmax), gammb(nmax);

double funci(double ideh, double xnuy)
{
    double xmi = (1 - ideh) * xme + ideh * xmh, temp = 0;
    for (int i = 1; i <= n; i++)
    {
        double selfi = (1 - ideh) * selfe[i] + ideh * selfh[i];
        double ek = nbeta * (s[i] * s[i] * xm / xmi + selfi - xnuy);
        //double ek = nbeta * (s[i] * pow(-1,ieh) * xm / xmi + selfi - xnuy);
        double sum;
        if (ek > 50)
            sum = 0;
        else
            sum = w[i] * 1. / (exp(ek) + 1);
        temp += sum;
    }
    return xn - 1. / pi * temp;
}

double rtbis(double func, double x1, double x2, double xacc)
{
    int jmax = 40;
    double ans, dx, xmid;
    double fmid = funci(func, x2);
    double f = funci(func, x1);
    if (f * fmid < 0)
    {
        if (f < 0.0)
        {
            ans = x1;
            dx = x2 - x1;
        }
        else
        {
            ans = x2;
            dx = x1 - x2;
        }
        for (int j = 1; j <= jmax; j++)
        {
            dx = dx * .5;
            xmid = ans + dx;
            fmid = funci(func, xmid);
            if (fmid <= 0.)
                ans = xmid;
            if ((fabs(dx) < xacc) || (fmid == 0.))
                return ans;
        }
    }
    cout << "root must be bracketed in rtbis" << endl;
    exit(1);
}

void gauleg(double x1, double x2, vector<double> &t, vector<double> &w)
{
    int m = (n + 1) / 2;
    double xm = 0.5 * (x2 + x1);
    double xl = 0.5 * (x2 - x1);
    double eps = 3.0e-14;
    for (int i = 1; i <= m; i++)
    {
        double pp, z1, z = cos(pi * (i - .25) / (n + .5));
        do
        {
            double p1 = 1.0;
            double p2 = 0.0;
            for (int j = 1; j <= n; j++)
            {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (abs(z - z1) <= eps);
        t[i] = xm - xl * z;
        t[n + 1 - i] = xm + xl * z;
        w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[n + 1 - i] = w[i];
    }
}

void chgvar()
{
    vector<double> z(n), v(n);
    gauleg(-1., 1., z, v); //gauleg(0, 1, z, v);
    double alpha = 0.5;
    for (int i = 1; i <= n; i++)
    {
        s[i] = alpha * z[i] / (1. - z[i] * z[i]);
        w[i] = v[i] * alpha * (z[i] * z[i] + 1.) / (z[i] * z[i] - 1.) / (z[i] * z[i] - 1.);
    }
}

double fermi(double x, double ieh)
{
    double emu = nbeta * (x - xnue * (1 - ieh) - xnuh * ieh);
    //double emu = nbeta * (x* pow(-1,ieh) - xnue * (1 - ieh) - xnuh * ieh);
    if (emu > 50)
        return 0;
    else
        return 1. / (exp(emu) + 1.);
}

void xeh()
{
    double xie = 0, xih = 0;
    for (int i = 1; i <= n; i++)
    {
        double fei = fermi(s[i] * s[i] * xm / xme, 0);
        double fhi = fermi(s[i] * s[i] * xm / xmh, 1);
        xie = xie + w[i] * fei * (1. - fei);
        xih = xih + w[i] * fhi * (1. - fhi);
    }
    xie = 1 / pi * xie;
    xih = 1 / pi * xih;
    xieh = xie + xih;
}

double vjj(double q)
{
    double x = q * q / (2. * omega);
    return 2. * cyl_bessel_k(0, x);
}

double exe(double p, double q)
{
    double fe = fermi((p - q) * (p - q) * xm / xme, 0);
    return -vjj(q) * fe;
}

double veh(double q)
{
    return -vjj(q);
}

double wpl(double q)
{
    return sqrt(2. * xn * q * q * vjj(q));
}

double wq(double q)
{
    return sqrt(pow(wpl(q), 2) + 2. * xn * q * q / (nbeta * xieh) + pow(q, 4));
}

double g(double x)
{
    double amu = nbeta * x;
    if (amu > 50)
        return 0;
    else
        return 1. / (exp(amu) - 1.);
}

double core(double p, double q)
{
    double eps = dam;
    double ei = p * p * xm / xme;
    double ej = (p - q) * (p - q) * xm / xme;
    double fe = fermi(ei, 0);
    double wwq = wq(q);
    double gg = g(wwq);
    return vjj(q) * pow(wpl(q), 2) / (2. * wwq) * ((1. + gg - fe) / (pow(ei - ej - wwq, 2) + eps * eps) * (ei - ej - wwq) + (gg + fe) / (pow(ei - ej + wwq, 2) + eps * eps) * (ei - ej + wwq));
}

double exh(double p, double q)
{
    double fh = fermi((p - q) * (p - q) * xm / xmh, 1);
    return -vjj(q) * fh;
}

double corh(double p, double q)
{
    double eps = dam;
    double ei = p * p * xm / xmh;
    double ej = (p - q) * (p - q) * xm / xmh;
    double fh = fermi(ei, 1);
    double wwq = wq(q);
    double gg = g(wwq);
    return vjj(q) * pow(wpl(q), 2) / (2. * wwq) * ((1. + gg - fh) / (pow(ei - ej - wwq, 2) + eps * eps) * (ei - ej - wwq) + (gg + fh) / (pow(ei - ej + wwq, 2) + eps * eps) * (ei - ej + wwq));
}

double xime(double p, double q)
{
    double eps = dam;
    double ei = p * p * xm / xme;
    double ej = (p - q) * (p - q) * xm / xme;
    double fe = fermi(ej, 0);
    double wwq = wq(q);
    double gg = g(wwq);
    return vjj(q) * pow(wpl(q), 2) / (2. * wwq) * ((1. + gg - fe) / (pow(ei - ej - wwq, 2) + eps * eps) * (ei - ej - wwq) + (gg + fe) / (pow(ei - ej + wwq, 2) + eps * eps) * (ei - ej + wwq));
}

double ximh(double p, double q)
{
    double eps = dam;
    double ei = p * p * xm / xmh;
    double ej = (p - q) * (p - q) * xm / xmh;
    double fh = fermi(ej, 1);
    double wwq = wq(q);
    double gg = g(wwq);
    return vjj(q) * pow(wpl(q), 2) / (2. * wwq) * ((1. + gg - fh) / (pow(ei - ej - wwq, 2) + eps * eps) * (ei - ej - wwq) + (gg + fh) / (pow(ei - ej + wwq, 2) + eps * eps) * (ei - ej + wwq));
}

double core0(double q)
{
    double ee = q * q * xm / xme;
    double fe = fermi(ee, 0);
    double wwq = wq(q);
    double gg = g(wwq);
    return vjj(q) * pow(wpl(q), 2) / (2. * wwq) * ((1. + gg - fe) / (pow(ee - wwq, 2) + dam * dam) * (ee - wwq) + (gg + fe) / (pow(ee + wwq, 2) + dam * dam) * (ee + wwq));
}

double corh0(double q)
{
    double ee = q * q * xm / xmh;
    double fh = fermi(ee, 1);
    double wwq = wq(q);
    double gg = g(wwq);
    return vjj(q) * pow(wpl(q), 2) / (2. * wwq) * ((1. + gg - fh) / (pow(ee - wwq, 2) + dam * dam) * (ee - wwq) + (gg + fh) / (pow(ee + wwq, 2) + dam * dam) * (ee + wwq));
}

void se()
{
    for (int i = 1; i <= n; i++)
    {
        double sum1 = 0, sum2 = 0, sum3 = 0;
        for (int j = 1; j <= n; j++)
        {
            sum1 += w[j] * (exe(s[i], s[j]) + core(s[i], s[j]));
            sum2 += w[j] * (exh(s[i], s[j]) + corh(s[i], s[j]));
            sum3 += w[j] * (xime(s[i], s[j]) + ximh(s[i], s[j]));
        }
        selfe[i] = 1. / (2 * pi) * sum1;
        selfh[i] = 1. / (2 * pi) * sum2;
        gammb[i] = 1. / (2 * pi) * sum3;
    }
    double sum1 = 0, sum2 = 0, sum3 = 0;
    for (int j = 1; j <= n; j++)
    {
        sum1 += w[j] * (exe(0.0, s[j]) + core0(s[j]));
        sum2 += sum2 + w[j] * (exh(0.0, s[j]) + corh0(s[j]));
    }
}

void susz(double &hw, vector<complex<double>> &sz)
{
    double idx = 1, eei, ehi, fei, fhi, resz, aisz;
    if (idx == 0)
        for (int i = 1; i <= n; i++)
        {
            eei = s[i] * s[i] * xm / xme + selfe[i];
            ehi = s[i] * s[i] * xm / xmh + selfh[i];
            fei = fermi(eei, 0);
            fhi = fermi(ehi, 1);
            resz = (hw - eei - ehi) * (fei + fhi - 1.) / ((hw - eei - ehi) * (hw - eei - ehi) + dam * dam);
            aisz = (1. - fei - fhi) * dam / ((hw - eei - ehi) * (hw - eei - ehi) + dam * dam);
            sz[i].real(resz);
            sz[i].imag(aisz);
        }
    else if (idx == 1)
        for (int i = 1; i <= n; i++)
        {
            eei = s[i] * s[i] * xm / xme + selfe[i];
            ehi = s[i] * s[i] * xm / xmh + selfh[i];
            fei = fermi(eei, 0);
            fhi = fermi(ehi, 1);
            resz = (hw - eei - ehi) * (fei + fhi - 1.) / ((hw - eei - ehi) * (hw - eei - ehi) + dam * dam);
            aisz = (1. - fei - fhi) * dam / ((hw - eei - ehi) * (hw - eei - ehi) + dam * dam);
            sz[i].real(resz);
            sz[i].imag(aisz);
        }
}

void vsm(vector<vector<complex<double>>> &xv, double &hw)
{
    complex<double> ss;
    vector<double> t(8, 0), ares(8, 0), aims(8, 0);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
        {
            if (j == i)
                xv[i][j] = {0.0, 0.0};
            else
            {
                double eei = s[i] * s[i] * xm / xme + selfe[i];
                double ehi = s[i] * s[i] * xm / xmh + selfh[i];
                double eej = s[j] * s[j] * xm / xme + selfe[j];
                double ehj = s[j] * s[j] * xm / xmh + selfh[j];
                // double eei = s[i] * xm / xme + selfe[i];
                // double ehi = s[i] * xm / xmh + selfh[i];
                // double eej = s[j] * xm / xme + selfe[j];
                // double ehj = s[j] * xm / xmh + selfh[j];
                double fei = fermi(eei, 0);
                double fej = fermi(eej, 0);
                double fhi = fermi(ehi, 1);
                double fhj = fermi(ehj, 1);
                double anehi = 1.0 - fei - fhi;
                double anehj = 1.0 - fej - fhj;
                double sij = s[i] - s[j];
                double vveh = veh(sij);
                double wwpl = wpl(sij);
                double wwq = wq(sij);
                double gg = g(wwq);
                double gam = 0.5;
                t[1] = -(gg + fej) * (fei - fermi(eej - wwq, 0)) / (eei - eej + wwq);
                t[2] = -(1 + gg - fej) * (fei - fermi(eej + wwq, 0)) / (eei - eej - wwq);
                t[3] = -(gg + fhj) * (fhi - fermi(ehj - wwq, 1)) / (ehi - ehj + wwq);
                t[4] = -(1 + gg - fhj) * (fhi - fermi(ehj + wwq, 1)) / (ehi - ehj - wwq);
                t[5] = -(gg + fej) * (fermi(eej - wwq, 0) + fhi - 1);
                t[5] /= ((hw - eej - ehi + wwq) * (hw - eej - ehi + wwq) + gam * gam);
                ares[5] = t[5] * (hw - eej - ehi + wwq);
                aims[5] = -t[5] * gam;
                t[6] = -(1 + gg - fej) * (fermi(eej + wwq, 0) + fhi - 1);
                t[6] /= ((hw - eej - ehi - wwq) * (hw - eej - ehi - wwq) + gam * gam);
                ares[6] = t[6] * (hw - eej - ehi - wwq);
                aims[6] = -t[6] * gam;
                t[7] = -(gg + fhj) * (fei + fermi(ehj - wwq, 1) - 1);
                t[7] = t[7] / ((hw - eei - ehj + wwq) * (hw - eei - ehj + wwq) + gam * gam);
                ares[7] = t[7] * (hw - eei - ehj + wwq);
                aims[7] = -t[7] * gam;
                t[8] = -(1 + gg - fhj) * (fei + fermi(ehj + wwq, 1) - 1) / ((hw - eei - ehj - wwq) * (hw - eei - ehj - wwq) + gam * gam);
                ares[8] = t[8] * (hw - eei - ehj - wwq);
                aims[8] = -t[8] * gam;
                double aress = 1. / (anehi * anehj) * wwpl * wwpl / (2 * wwq) * (t[1] + t[2] + t[3] + t[4] + ares[5] + ares[6] + ares[7] + ares[8]);
                double aimss = 1. / (anehi * anehj) * wwpl * wwpl / (2 * wwq) * (aims[5] + aims[6] + aims[7] + aims[8]);
                ss = {aress, 0 * aimss};
                xv[i][j] = vveh * (ss + 1.);
            }
        }
    ofstream vsmfile("vsm.txt");
    for (int j = 1; j <= n; j++)
        if (j != n / 2)
            vsmfile << s[j] << " " << vjj(s[n / 2] - s[j]) << " " << -real(xv[n / 2][j]) << endl;
}

void xxx(vector<complex<double>> &sz, vector<vector<complex<double>>> &xv, vector<vector<complex<double>>> &xa, vector<complex<double>> &y)
{
    for (int i = 1; i <= n; i++)
    {
        y[i] = sz[i];
        for (int j = 1; j <= n; j++)
            if (i == j)
                xa[i][j] = {1, 0};
            else
                xa[i][j] = complex<double>(0, 0) + 1. / (2 * pi) * sz[i] * w[j] * xv[i][j];
    }
}

void lud(vector<vector<complex<double>>> &ax)
{
    complex<double> sum;
    for (int j = 1; j <= n; j++)
    {
        for (int i = 1; i <= j; i++)
        {
            sum = ax[i][j];
            for (int k = 1; k <= i - 1; k++)
                sum = sum - ax[i][k] * ax[k][j];
            ax[i][j] = sum;
        }
        for (int i = j + 1; i <= n; i++)
        {
            sum = ax[i][j];
            for (int k = 1; k <= j - 1; k++)
                sum = sum - ax[i][k] * ax[k][j];
            sum = sum / ax[j][j];
            ax[i][j] = sum;
        }
    }
}

void luge(vector<vector<complex<double>>> &ax, vector<complex<double>> &b)
{
    complex<double> sum;
    for (int i = 1; i <= n; i++)
    {
        sum = b[i - 1];
        for (int j = 1; j <= i - 1; j++)
            sum = sum - ax[i - 1][j - 1] * b[j - 1];
        b[i - 1] = sum;
    }
    for (int i = n; i >= 1; i--)
    {
        sum = b[i - 1];
        for (int j = i + 1; j <= n; j++)
            sum = sum - ax[i - 1][j - 1] * b[j - 1];
        b[i - 1] = sum / ax[i - 1][i - 1];
    }
}

void run(double &x1, double &x2, int &m, vector<double> &xw, vector<double> &aimsus, vector<double> &aresus)
{
    complex<double> sum, susc;
    vector<complex<double>> sz(nmax), vec(nmax);
    vector<vector<complex<double>>> xv(nmax, vector<complex<double>>(nmax)), xa(nmax, vector<complex<double>>(nmax));
    double hw = x1, step = (x2 - x1) / m;
    for (int k = 1; k <= m + 1; k++)
    {
        xw[k] = hw;
        susz(hw, sz);
        vsm(xv, hw);
        xxx(sz, xv, xa, vec);
        lud(xa);
        luge(xa, vec);
        sum = 0;
        for (int i = 1; i <= n; i++)
            sum += w[i] * vec[i];
        susc = 1. / (2 * pi) * sum;
        aimsus[k] = imag(susc);
        aresus[k] = real(susc);
        // cout << xw[k] << " " << aimsus[k] << endl;
        hw += step;
    }
}

int main()
{
    double tem, x1, x2;
    int m;
    vector<double> hw(mmax), ab(mmax), re(mmax);
    string outfile;
    // cout << "Harmonic frequency (in Rydberg):";
    // cin >> omega;
    omega = 5;
    // cout << "Damping rate (in Rydberg):";
    // cin >> dam;
    dam = 0.5;
    // cout << "Plasma density (in 1 / Bohr radius):";
    // cin >> xn;
    xn = 2;
    // cout << "Temperature (K):";
    // cin >> tem;
    tem = 300;
    // cout << "Gaussian point number (max2000):";
    // cin >> n;
    n = 30;
    // cout << "Output filename:";
    // cin >> outfile;
    //int nstar = 1
    //double B = 2 or 6 or 12
    outfile = "test.txt";
    // cout << "hw-Eg: min  max  step_number(max200)";
    // cin >> x1 >> x2 >> m;
    x1 = -10.0, x2 = 10.0, m = 1;
    xme = 0.067, xmh = 0.457; //xme = 1, xmh = -1, xm = sqrt(abs(nstar)*B); <Change Here>
    xm = xme * xmh / (xme + xmh);
    double e0 = 4.6;
    nbeta = e0 / (0.0862 * tem);
    for (int i = 1; i <= n; i++)
    {
        selfe[i] = 0;
        selfh[i] = 0;
    }
    chgvar();
    // Chemical potential:
    xnue = rtbis(0, -1.0e+02, 1.0e+02, 1.0e-6);
    xnuh = rtbis(1, -1.0e+02, 1.0e+02, 1.0e-6);
    cout << "Chemical potential: " << xnue + xnuh << endl;
    // Renormalized chemical potential:
    xeh();
    se();
    xnue = rtbis(0, -1.0e+02, 1.0e+02, 1.0e-6);
    xnuh = rtbis(1, -1.0e+02, 1.0e+02, 1.0e-6);
    double xnu = xnue + xnuh;
    cout << "Renormalized chemical potential: " << xnu << endl;
    cout << "BGR:" << selfe[n / 2] + selfh[n / 2] << endl;
    cout << "gamma:" << gammb[n / 2] << endl;
    run(x1, x2, m, hw, ab, re);
    ofstream ofile(outfile);
    ofile << "omega=" << omega << "; dam=" << dam << "; tem=" << tem << "; xn=" << xn << endl;
    for (int i = 1; i <= m + 1; i++)
    {
        double em = g(hw[i] - xnu) * ab[i];
        ofile << hw[i] << " " << re[i] << " " << ab[i] << " " << em << endl;
    }
    ofstream sefile("se.txt");
    for (int i = 1; i <= n; i++)
        sefile << s[i] << " " << selfe[i] + selfh[i] << " " << gammb[i] << endl;
}