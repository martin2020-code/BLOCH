#include <lapack.h>
#include "ula/include/ula/matrixtypes.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PI 3.141549265358979323846

using namespace ula;
using namespace std;

Complex I(0., 1.);
int N, N_av, Ic;
double Wc = 1., Wa = 1., Om = 1.;
double tf = 100., dt = 1.;

double Norm_(ComplexVector phi)
{
  double Norm = 0.;
  for (int i = 0; i < 2 * N; i++)
    Norm += real(phi(i) * conj(phi(i)));
  return Norm;
}

double delta_(int i, int j)
{
  int k = 0;
  if (i == j)
    k = 1;
  return k;
}

ComplexVector IniState(Real theta, Real beta)
{
  ComplexVector phi(2 * N);
  double Norm = 0.;
  boost::math::poisson_distribution<> p(N_av);
  for (int n = 0; n < N; n++)
  {
    phi(2 * n + 0) = pdf(p, n) * cos(theta / 2.);
    phi(2 * n + 1) = pdf(p, n) * sin(theta / 2.) * exp(I * beta);
    // phi(2*n+0)=delta_(2,n)*cos(theta/2.);
    // phi(2*n+1)=delta_(2,n)*sin(theta/2.)*exp(I*beta);
  }
  // phi(6)=1.;
  Norm = Norm_(phi);
  for (int i = 0; i < 2 * N; i++)
    phi(i) = phi(i) / sqrt(Norm);
  //  cout<<"phi = "<<phi<<endl;
  return phi;
}

ComplexMatrix Tr_E(ComplexVector phi)
{
  ComplexMatrix rho(2, 2);
  for (int n = 0; n < N; n++)
  {
    rho(0, 0) += phi(2 * n + 0) * conj(phi(2 * n + 0));
    rho(1, 1) += phi(2 * n + 1) * conj(phi(2 * n + 1));
    rho(0, 1) += phi(2 * n + 0) * conj(phi(2 * n + 1));
  }
  rho(1, 0) = conj(rho(0, 1));
  return rho;
}

ComplexVector Evo_t(ComplexVector phi_0, Real t)
{
  ComplexVector phi_t(2 * N);
  Real Delta = Wc - Wa;
  for (int n = 1; n < N; n++)
  {
    Real Theta_n = atan2(sqrt(n) * Om, -Delta);
    Real CThn = cos(Theta_n / 2.);
    Real SThn = sin(Theta_n / 2.);
    RealVector Epm(2);
    Epm(0) = n * Wc - Delta / 2. + sqrt(Delta * Delta + n * Om * Om) / 2.;
    Epm(1) = n * Wc - Delta / 2. - sqrt(Delta * Delta + n * Om * Om) / 2.;
    ComplexVector phipm(2);
    phipm(0) = CThn * phi_0(2 * n - 1) + SThn * phi_0(2 * n);
    phipm(1) = -SThn * phi_0(2 * n - 1) + CThn * phi_0(2 * n);

    phi_t(2 * n - 1) = phipm(0) * CThn * exp(-I * Epm(0) * t) - SThn * phipm(1) * exp(-I * Epm(1) * t);
    phi_t(2 * n) = phipm(0) * SThn * exp(-I * Epm(0) * t) + CThn * phipm(1) * exp(-I * Epm(1) * t);
  }
  phi_t(0) = exp(-I * (-Delta / 2.) * t) * phi_0(0);
  phi_t(2 * N - 1) = exp(-I * (N * Wc + Delta / 2.) * t) * phi_0(2 * N - 1);
  return phi_t;
}

int main()
{
  N_av = 25;
  N = 100;
  int N_Calc_T = 0;
  for (int i = 0; i <=16; i++)
  {
    double theta = i * PI / 16.;
    for (double beta = 0; beta <= 32. * sin(theta); beta += PI / 16.)
    {
      N_Calc_T++;
    }
  }
  cout << "Total number of files: " << N_Calc_T << endl;
  int d2 = tf / dt;
  cout << "d2 = " << d2 << " N_Calc_T = " << N_Calc_T << endl;
  double M[3][d2][N_Calc_T];
  cout << "Ok" << endl;
  int N_File = 0;
  for (int i = 0; i <= 16; i++)
  {
    double theta = i * PI / 16.;
    for (double beta = 0; beta <= 32. * sin(theta); beta += PI / 16.)
    {
      FILE *dskw1;
      int Err;
      char FileName[200];
      ComplexVector phi_0(2 * N);
      //Err = sprintf(FileName, "/home/arturo/Scratch/JayC/N00025/denmat-theta_%07.6g-beta_%07.6g.dat", theta, beta);
      Err = sprintf(FileName, "./data/N00025/denmat-theta_%07.6g-beta_%07.6g.dat", theta, beta);
      Err++;
      cout << FileName << endl;
      dskw1 = fopen(FileName, "w+");
      phi_0 = IniState(theta, beta);
      ComplexVector phi_t(2 * N);
      ComplexMatrix rho_s(2, 2);
      for (int Tindex = 0; Tindex < d2; Tindex++)
      {
        double t = Tindex * dt;
        phi_t = Evo_t(phi_0, t);
        //   cout<<"phi_t = "<<phi_t<<endl;
        rho_s = Tr_E(phi_t);
        cout << "Tindex = " << Tindex << " N_File = " << N_File << endl;
        M[0][Tindex][N_File] = 2. * real(rho_s(0, 1));
        M[1][Tindex][N_File] = 2. * imag(rho_s(1, 0));
        M[2][Tindex][N_File] = real(rho_s(0, 0) - rho_s(1, 1));
        fprintf(dskw1, "%g %g %g %g %g\n", t,
                2. * real(rho_s(0, 1)),
                2. * imag(rho_s(1, 0)),
                real(rho_s(0, 0) - rho_s(1, 1)),
                real(rho_s(1, 1)));
        //  fprintf(dskw1,"%g %g %g %g %g\n",t,real(rho_s(0,0)),imag(rho_s(0,0)),real(rho_s(1,1)),imag(rho_s(1,1)));
        fflush(dskw1);
      }
      N_File++;
      fclose(dskw1);
    }
  }
  for (int N = 0; N < d2; N++)
  {
    FILE *dskw2;
    int Err;
    char FileName[200];
    // Err=sprintf(FileName,"/home/arturo/Scratch/JayC/VideoData/snapshot_%04i.dat",N);
    Err = sprintf(FileName, "./data/VideoData/snapshot_%04i.dat", N);
    Err++;
    cout << FileName << endl;
    dskw2 = fopen(FileName, "w+");
    for (int N_ = 0; N_ < N_Calc_T; N_++)
      fprintf(dskw2, "%g %g %g\n", M[0][N][N_], M[1][N][N_], M[2][N][N_]);
    fclose(dskw2);
  }
  return 0;
}

/*

clear
make Jaynes_Cummings_v2.20912
rm /home/arturo/Scratch/JayC/VideoData/*.dat
time ./Jaynes_Cummings_v2.20912
ls /home/arturo/Scratch/JayC/VideoData/*.dat > arg-file
rm /home/arturo/Scratch/JayC/VideoData/*.png
time parallel -v -j+0 A={} /home/arturo/Scratch/JayC/VideoData/parallel_plotter :::: arg-file
rm /home/arturo/Scratch/JayC/bloch.mp4
ffmpeg -framerate 10 -i /home/arturo/Scratch/JayC/VideoData/snapshot_%04d.dat.png /home/arturo/Scratch/JayC/bloch.mp4
 */
