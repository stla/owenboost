#include <boost/multiprecision/float128.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <cmath>

namespace mp = boost::multiprecision;
namespace m = boost::math;

// double dnorm(double x){
//   return exp(-x*x/2)/sqrt(2*M_PI);
// }

double owent(double h, double a){
  return h*a;
}

double pnorm(double q){
  return m::erfc(-q/1.414213562373095)/2.0;
}

mp::float128 dnorm128(mp::float128 x){
  return mp::exp(-x*x/2)/2.5066282746310005024157652848110452530069867406099Q;
}

mp::float128 pnorm128(mp::float128 q){
  return m::erfc(-q / 1.4142135623730950488016887242096980785696718753769Q)/2;
}

double* testvectorout(double* in, size_t n){
  double* out = new double[n];
  int i;
  for(i=0; i<n; i++){
    out[i] = pnorm(in[i]);
  }
  return out;
}

int sign(double x){
//  return (x>0 ? 1 : -1);
  // if(x == 0.0){
  //   return 0;
  // }
  return (signbit(x) ? -1 : 1);
}

double* studentC(double q, int nu, double* delta, size_t J){
  double a = sqrt(q*q/nu);
  double sb = sqrt(nu/(nu+q*q));
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
  std::vector< std::vector<mp::float128> > M(nu-1, std::vector<mp::float128>(J));
  for(j=0; j<J ; j++){
    M[0][j] = a * sb * dnorm128(dsb[j]) * pnorm128(a*dsb[j]);
  }
  const mp::float128 sqrt2pi = 2.5066282746310005024157652848110452530069867406099Q;
  if(nu>2){
    for(j=0; j<J; j++){
      M[1][j] = b * (delta[j] * a * M[0][j] + a * dnorm128(delta[j]) / sqrt2pi);
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
    out[j] = pnorm(-delta[j]) + (sqrt2pi*sum[j]).convert_to<double>();
  }
  return out;
}
