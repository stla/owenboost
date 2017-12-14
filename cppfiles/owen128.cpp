#include <boost/multiprecision/float128.hpp>
//#include <boost/math/special_functions/owens_t.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <cmath>
#include <vector>

namespace mp = boost::multiprecision;
namespace m = boost::math;

const double one_div_root_two_pi = 3.989422804014326779399460599343818684e-01;
const double one_div_root_two = 7.071067811865475244008443621048490392e-01;
const double root_two_pi = 2.506628274631000502415765284811045253;
const double one_div_two_pi = 0.159154943091895335768883763372514362;
const mp::float128 root_two_pi128 = 2.5066282746310005024157652848110452530069867406099Q;
const mp::float128 root_two128 = 1.4142135623730950488016887242096980785696718753769Q;

double dnorm(double x){
  return exp(-x*x/2) * one_div_root_two_pi;
}

double pnorm(double q){
  return m::erfc(-q * one_div_root_two)/2.0;
}

mp::float128 dnorm128(mp::float128 x){
  return mp::exp(-x*x/2)/root_two_pi128;
}

mp::float128 pnorm128(mp::float128 q){
  return m::erfc(-q / root_two128)/2;
}

double* testvectorout(double* in, size_t n){
  double* out = new double[n];
  int i;
  for(i=0; i<n; i++){
    out[i] = pnorm(in[i]);
  }
  return out;
}

//********* Owen T-function **************************************************//
//****** http://people.sc.fsu.edu/~jburkardt/cpp_src/owens/owens.html ********//
double znorm1(double x){
  return 0.5 * m::erf ( x * one_div_root_two );
}

double znorm2(double x){
  return 0.5 * m::erfc ( x * one_div_root_two );
}

double tfun ( double h, double a, double ah ){
  double ai;
  double arange[7] = {0.025, 0.09, 0.15, 0.36, 0.5, 0.9, 0.99999};
  double as;
  double c2[21] = {0.99999999999999987510,
   -0.99999999999988796462,      0.99999999998290743652,
   -0.99999999896282500134,      0.99999996660459362918,
   -0.99999933986272476760,      0.99999125611136965852,
   -0.99991777624463387686,      0.99942835555870132569,
   -0.99697311720723000295,      0.98751448037275303682,
   -0.95915857980572882813,      0.89246305511006708555,
   -0.76893425990463999675,      0.58893528468484693250,
   -0.38380345160440256652,      0.20317601701045299653,
   -0.82813631607004984866E-01,  0.24167984735759576523E-01,
   -0.44676566663971825242E-02,  0.39141169402373836468E-03 };
  double hrange[14] = {
    0.02, 0.06, 0.09, 0.125, 0.26,
    0.4,  0.6,  1.6,  1.7,   2.33,
    2.4,  3.36, 3.4,  4.8 };
  double hs;
  int i;
  int iaint;
  int icode;
  int ifail;
  int ihint;
  int ii;
  int m;
  int maxii;
  int meth[18] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6 };
  double normh;
  int ord[18] = {2, 3, 4, 5, 7,10,12,18,10,20,30,20, 4, 7, 8,20,13, 0 };
  double pts[13] = {0.35082039676451715489E-02,
      0.31279042338030753740E-01,  0.85266826283219451090E-01,
      0.16245071730812277011,      0.25851196049125434828,
      0.36807553840697533536,      0.48501092905604697475,
      0.60277514152618576821,      0.71477884217753226516,
      0.81475510988760098605,      0.89711029755948965867,
      0.95723808085944261843,      0.99178832974629703586 };
  double r;
  int select[15*8] = {
    1, 1, 2,13,13,13,13,13,13,13,13,16,16,16, 9,
    1, 2, 2, 3, 3, 5, 5,14,14,15,15,16,16,16, 9,
    2, 2, 3, 3, 3, 5, 5,15,15,15,15,16,16,16,10,
    2, 2, 3, 5, 5, 5, 5, 7, 7,16,16,16,16,16,10,
    2, 3, 3, 5, 5, 6, 6, 8, 8,17,17,17,12,12,11,
    2, 3, 5, 5, 5, 6, 6, 8, 8,17,17,17,12,12,12,
    2, 3, 4, 4, 6, 6, 8, 8,17,17,17,17,17,12,12,
    2, 3, 4, 4, 6, 6,18,18,18,18,17,17,17,12,12 };
  double value;
  double vi;
  double wts[13] = {0.18831438115323502887E-01,
      0.18567086243977649478E-01,  0.18042093461223385584E-01,
      0.17263829606398753364E-01,  0.16243219975989856730E-01,
      0.14994592034116704829E-01,  0.13535474469662088392E-01,
      0.11886351605820165233E-01,  0.10070377242777431897E-01,
      0.81130545742299586629E-02,  0.60419009528470238773E-02,
      0.38862217010742057883E-02,  0.16793031084546090448E-02 };
  //double x;
  double y;
  double yi;
  double z;
  double zi;
  //
  //  Determine appropriate method from t1...t6
  //
  ihint = 15;

  for ( i = 1; i <= 14; i++ )
  {
    if ( h <= hrange[i-1] )
    {
      ihint = i;
      break;
    }
  }

  iaint = 8;

  for ( i = 1; i <= 7; i++ )
  {
    if ( a <= arange[i-1] )
    {
      iaint = i;
      break;
    }
  }

  icode = select[ihint-1+(iaint-1)*15];
  m = ord[icode-1];
  //
  //  t1(h, a, m) ; m = 2, 3, 4, 5, 7, 10, 12 or 18
  //  jj = 2j - 1 ; gj = exp(-h*h/2) * (-h*h/2)**j / j//
  //  aj = a**(2j-1) / (2*pi)
  //
  if ( meth[icode-1] == 1 )
  {
    hs = - 0.5 * h * h;
    double dhs = exp ( hs );
    as = a * a;
    int j = 1;
    int jj = 1;
    double aj = one_div_two_pi * a;
    value = one_div_two_pi * atan ( a );
    double dj = dhs - 1.0;
    double gj = hs * dhs;

    for ( ; ; )
    {
      value = value + dj * aj / ( double ) ( jj );

      if ( m <= j )
      {
        return value;
      }
      j = j + 1;
      jj = jj + 2;
      aj = aj * as;
      dj = gj - dj;
      gj = gj * hs / ( double ) ( j );
    }
  }
  //
  //  t2(h, a, m) ; m = 10, 20 or 30
  //  z = (-1)**(i-1) * zi ; ii = 2i - 1
  //  vi = (-1)**(i-1) * a**(2i-1) * exp[-(a*h)**2/2] / sqrt(2*pi)
  //
  else if ( meth[icode-1] == 2 )
  {
    maxii = m + m + 1;
    ii = 1;
    value = 0.0;
    hs = h * h;
    as = - a * a;
    vi = one_div_root_two_pi * a * exp ( - 0.5 * ah * ah );
    z = znorm1 ( ah ) / h;
    y = 1.0 / hs;

    for ( ; ; )
    {
      value = value + z;

      if ( maxii <= ii )
      {
        value = value * one_div_root_two_pi * exp ( - 0.5 * hs );
        return value;
      }
      z = y * ( vi - ( double ) ( ii ) * z );
      vi = as * vi;
      ii = ii + 2;
    }
  }
  //
  //  t3(h, a, m) ; m = 20
  //  ii = 2i - 1
  //  vi = a**(2i-1) * exp[-(a*h)**2/2] / sqrt(2*pi)
  //
  else if ( meth[icode-1] == 3 )
  {
    i = 1;
    ii = 1;
    value = 0.0;
    hs = h * h;
    as = a * a;
    vi = one_div_root_two_pi * a * exp ( - 0.5 * ah * ah );
    zi = znorm1 ( ah ) / h;
    y = 1.0 / hs;

    for ( ; ; )
    {
      value = value + zi * c2[i-1];

      if ( m < i )
      {
        value = value * one_div_root_two_pi * exp ( - 0.5 * hs );
        return value;
      }
      zi = y  * ( ( double ) ( ii ) * zi - vi );
      vi = as * vi;
      i = i + 1;
      ii = ii + 2;
    }
  }
  //
  //  t4(h, a, m) ; m = 4, 7, 8 or 20;  ii = 2i + 1
  //  ai = a * exp[-h*h*(1+a*a)/2] * (-a*a)**i / (2*pi)
  //
  else if ( meth[icode-1] == 4 )
  {
    maxii = m + m + 1;
    ii = 1;
    hs = h * h;
    as = - a * a;
    value = 0.0;
    ai = one_div_two_pi * a * exp ( - 0.5 * hs * ( 1.0 - as ) );
    yi = 1.0;

    for ( ; ; )
    {
      value = value + ai * yi;

      if ( maxii <= ii )
      {
        return value;
      }
      ii = ii + 2;
      yi = ( 1.0 - hs * yi ) / ( double ) ( ii );
      ai = ai * as;
    }
  }
  //
  //  t5(h, a, m) ; m = 13
  //  2m - point gaussian quadrature
  //
  else if ( meth[icode-1] == 5 )
  {
    value = 0.0;
    as = a * a;
    hs = - 0.5 * h * h;
    for ( i = 1; i <= m; i++ )
    {
      r = 1.0 + as * pts[i-1];
      value = value + wts[i-1] * exp ( hs * r ) / r;
    }
    value = a * value;
  }
  //
  //  t6(h, a);  approximation for a near 1, (a<=1)
  //
  else if ( meth[icode-1] == 6 )
  {
    normh = znorm2 ( h );
    value = 0.5 * normh * ( 1.0 - normh );
    y = 1.0 - a;
    r = atan ( y / ( 1.0 + a ) );

    if ( r != 0.0 )
    {
      value = value - one_div_two_pi * r * exp ( - 0.5 * y * h * h / r );
    }
  }
  return value;
}

double owent(double h, double a){
  double cut = 0.67;
  double normah;
  double normh;
  double value;

  double absh = fabs ( h );
  double absa = fabs ( a );
  double ah = absa * absh;

  if ( absa <= 1.0 )
  {
    value = tfun ( absh, absa, ah );
  }
  else if ( absh <= cut )
  {
    value = 0.25 - znorm1 ( absh ) * znorm1 ( ah )
      - tfun ( ah, 1.0 / absa, absh );
  }
  else
  {
    normh = znorm2 ( absh );
    normah = znorm2 ( ah );
    value = 0.5 * ( normh + normah ) - normh * normah
    - tfun ( ah, 1.0 / absa, absh );
  }

  if ( a < 0.0 )
  {
    value = - value;
  }

  return value;
}
//****************************************************************************//

int sign(double x){
//  return (x>0 ? 1 : -1);
  // if(x == 0.0){
  //   return 0;
  // }
  return (signbit(x) ? -1 : 1);
}

double* studentC(double q, int nu, double* delta, size_t J){
  const double a = sign(q)*sqrt(q*q/nu);
  const double sb = sqrt(nu/(nu+q*q));
  double* C = new double[J];
  int j;
  for(j=0; j<J; j++){
    C[j] = owent(delta[j] * sb, a) + pnorm(-delta[j]*sb); //+
  }
  return C;
}

double* studentCDF(double q, int nu, double* delta, size_t J){
  if(nu==1){
    return studentC(q, nu, delta, J);
  }
  const mp::float128 qq(q*q);
  const mp::float128 a = sign(q)*mp::sqrt(qq/nu);
  const mp::float128 b = nu/(nu+qq);
  const mp::float128 sb = mp::sqrt(b);
  std::vector<mp::float128> dsb(J);
  int j;
  for(j=0; j<J; j++){
    dsb[j] = delta[j] * sb;
  }
  //std::vector< std::vector<mp::float128> > M(nu-1, std::vector<mp::float128>(J));
  mp::float128 M[nu-1][J];
  for(j=0; j<J ; j++){
    M[0][j] = a * sb * dnorm128(dsb[j]) * pnorm128(a*dsb[j]);
  }
  if(nu>2){
    for(j=0; j<J; j++){
      M[1][j] = b * (delta[j] * a * M[0][j] + a * dnorm128(delta[j]) / root_two_pi128);
    }
    if(nu>3){
      std::vector<mp::float128> A(nu-3);
      A[0] = 1.0;
      int k;
      if(nu>4){
        for(k=1; k<nu-3; k++){
          A[k] = 1.0/k/A[k-1];
        }
      }
      for(k=2; k<nu-1; k++){
        for(j=0; j<J; j++){
          M[k][j] = (k-1) * b * (A[k-2] * delta[j] * a * M[k-1][j] + M[k-2][j]) / k;
        }
      }
    }
  }
  if(nu%2==1){
    double* C = studentC(q, nu, delta, J);
    std::vector<mp::float128> sum(J);
    int i;
    for(i=1; i<nu-1; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M[i][j];
      }
    }
    double* out = new double[J];
    for(j=0; j<J; j++){
      out[j] = C[j] + 2*sum[j].convert_to<double>();
    }
    delete[] C;
    return out;
  }
  int i;
  std::vector<mp::float128> sum(J);
  for(i=0; i<nu-1; i+=2){
    for(j=0; j<J; j++){
      sum[j] += M[i][j];
    }
  }
  double* out = new double[J];
  for(j=0; j<J; j++){
    out[j] = pnorm(-delta[j]) + (root_two_pi128*sum[j]).convert_to<double>();
  }
  return out;
}

double* owenC(int nu, double t, double* delta, double* R, size_t J){
  const double a = sign(t)*sqrt(t*t/nu);
  const double b = nu/(nu+t*t);
  const double sb = sqrt(b);
  double ab;
  if(fabs(t)>DBL_MAX){
    ab = 0;
  }else{
    ab = a*b;
  }
  double* C = new double[J];
  int i;
  for(i=0; i<J; i++){
    double C1 = owent(delta[i]*sb, a);
    double C2 = owent(R[i], a-delta[i]/R[i]);
    double C3 = owent(delta[i]*sb, (ab-R[i]/delta[i])/b);
    C[i] = pnorm(R[i]) - (delta[i] >= 0) + 2*(C1 - C2 - C3);
  }
  return C;
}

double* owenQ(int nu, double t, double* delta, double* R, size_t J){
  if(nu == 1){
    return owenC(nu, t, delta, R, J);
  }
  const double a = sign(t)*sqrt(t*t/nu);
  const double b = nu/(nu+t*t);
  const double sb = sqrt(b);
  double ab;
  double asb;
  if(fabs(t) > DBL_MAX){
    ab = 0;
    asb = sign(t);
  }else{
    ab = a*b;
    asb = sign(t)*sqrt(t*t/(nu+t*t));
  }
  std::vector<double> dsb(J);
  std::vector<double> dnormdsb(J);
  std::vector<double> dabminusRoversb(J);
  std::vector<double> dnormR(J);
  int j;
  for(j=0; j<J; j++){
    dsb[j] = delta[j] * sb;
    dnormdsb[j] = dnorm(dsb[j]);
    dabminusRoversb[j] = (delta[j]*ab - R[j])/sb;
    dnormR[j] = dnorm(R[j]);
  }
  const int n = nu-1;
  // plante si j'initialise les deux vecteurs !!
  //std::vector< std::vector<double> > M(n, std::vector<double>(J));
  double H[n][J];
  double M[n][J];
  //std::vector< std::vector<double> > H(n, std::vector<double>(J));
  for(j=0; j<J; j++){
    H[0][j] = -dnormR[j] * pnorm(a*R[j]-delta[j]);
    M[0][j] = asb * dnormdsb[j] * (pnorm(dsb[j]*a) - pnorm(dabminusRoversb[j]));
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      H[1][j] = R[j] * H[0][j];
      M[1][j] = delta[j]*ab*M[0][j] + ab * dnormdsb[j] * (dnorm(dsb[j]*a) - dnorm(dabminusRoversb[j]));
    }
    if(nu >= 4){
      double A[n]; //std::vector<double> A(n);
      //std::vector< std::vector<double> > L(n-2, std::vector<double>(J));
      double L[n-2][J];
      A[0] = 1;
      A[1] = 1;
      for(j=0; j<J; j++){
        L[0][j] = ab * R[j] * dnormR[j] * dnorm(a*R[j]-delta[j])/2;
      }
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          for(j=0; j<J; j++){
            L[k][j] = A[k+2] * R[j] * L[k-1][j];
          }
        }
      }
      for(k=2; k<n; k++){
        for(j=0; j<J; j++){
          H[k][j] = A[k] * R[j] * H[k-1][j];
          M[k][j] = (k-1.0)/k * (A[k-2] * delta[j] * ab * M[k-1][j] + b*M[k-2][j]) - L[k-2][j];
        }
      }
    }
  }
  if(nu % 2 == 0){
    int i;
    std::vector<double> sum(J);
    for(i=0; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M[i][j]+H[i][j];
      }
    }
    double* out = new double[J];
    for(j=0; j<J; j++){
      out[j] = pnorm(-delta[j]) + root_two_pi*sum[j];
    }
    return out;
  }else{
    std::vector<double> sum(J);
    int i;
    for(i=1; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M[i][j]+H[i][j];
      }
    }
    double* out = new double[J];
    double* C = owenC(nu, t, delta, R, J);
    for(j=0; j<J; j++){
      out[j] = C[j] + 2*sum[j];
    }
    delete[] C;
    return out;
  }
}

double* owenQ128(int nu, double t, double* delta, double* R, size_t J){
  if(nu == 1){
    return owenC(nu, t, delta, R, J);
  }
  const mp::float128 tt(t*t);
  const mp::float128 a = sign(t)*mp::sqrt(tt/nu);
  const mp::float128 b = nu/(nu+tt);
  const mp::float128 sb = mp::sqrt(b);
  mp::float128 ab;
  mp::float128 asb;
  if(fabs(t) > DBL_MAX){
    ab = 0;
    asb = sign(t);
  }else{
    ab = a*b;
    asb = sign(t)*mp::sqrt(tt/(nu+tt));
  }
  mp::float128 dsb[J];
  mp::float128 dnormdsb[J];
  mp::float128 dabminusRoversb[J];
  mp::float128 dnormR[J];
  int j;
  for(j=0; j<J; j++){
    dsb[j] = delta[j] * sb;
    dnormdsb[j] = dnorm128(dsb[j]);
    dabminusRoversb[j] = (delta[j]*ab - R[j])/sb;
    dnormR[j] = dnorm128(R[j]);
  }
  const int n = nu-1;
  mp::float128 H[n][J];
  mp::float128 M[n][J];
  for(j=0; j<J; j++){
    H[0][j] = -dnormR[j] * pnorm128(a*R[j]-delta[j]);
    M[0][j] = asb * dnormdsb[j] * (pnorm128(dsb[j]*a) - pnorm128(dabminusRoversb[j]));
  }
  if(nu >= 3){
    for(j=0; j<J; j++){
      H[1][j] = R[j] * H[0][j];
      M[1][j] = delta[j]*ab*M[0][j] + ab * dnormdsb[j] *
                 (dnorm128(dsb[j]*a) - dnorm128(dabminusRoversb[j]));
    }
    if(nu >= 4){
      mp::float128 A[n];
      mp::float128 L[n-2][J];
      A[0] = 1;
      A[1] = 1;
      for(j=0; j<J; j++){
        L[0][j] = ab * R[j] * dnormR[j] * dnorm128(a*R[j]-delta[j])/2;
      }
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          for(j=0; j<J; j++){
            L[k][j] = A[k+2] * R[j] * L[k-1][j];
          }
        }
      }
      for(k=2; k<n; k++){
        for(j=0; j<J; j++){
          H[k][j] = A[k] * R[j] * H[k-1][j];
          M[k][j] = (k-1.0)/k * (A[k-2] * delta[j] * ab * M[k-1][j] + b*M[k-2][j]) - L[k-2][j];
        }
      }
    }
  }
  if(nu % 2 == 0){
    int i;
    mp::float128 sum[J];
    for(i=0; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M[i][j]+H[i][j];
      }
    }
    double* out = new double[J];
    for(j=0; j<J; j++){
      out[j] = pnorm(-delta[j]) + root_two_pi*sum[j].convert_to<double>();
    }
    return out;
  }else{
    mp::float128 sum[J];
    int i;
    for(i=1; i<n; i+=2){
      for(j=0; j<J; j++){
        sum[j] += M[i][j]+H[i][j];
      }
    }
    double* out = new double[J];
    double* C = owenC(nu, t, delta, R, J);
    for(j=0; j<J; j++){
      out[j] = C[j] + 2*sum[j].convert_to<double>();
    }
    delete[] C;
    return out;
  }
}
