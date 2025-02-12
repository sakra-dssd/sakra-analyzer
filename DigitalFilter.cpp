#include <complex>
#include <vector>
#include <algorithm>

using namespace std;

vector<complex<double>> FFT(vector<complex<double>> data, bool inverse=false){
  int n=data.size();
  int h=0;
  for(int i=0;(1<<i)<n;++i) ++h;

  for(int i=0;i<n;++i){
    int j=0;
    for(int k=0;k<n;++k) j |= ((i>>k) & 1) << (h-1-k);
    if(i<j) swap(data[i],data[j]);
  }

  for(int b=1;b<n;b*=2){
    for(int j=0;j<b;++j){
      complex<double> w=polar(1.0, (2*M_PI)/(2*b)*j*(inverse ? 1 : -1));
      for(int k=0;k<n;k+=b*2){
	complex<double> s=data[j+k];
	complex<double> t=data[j+k+b]*w;
	data[j+k]=s+t;
	data[j+k+b]=s-t;
      }
    }
  }

  if(inverse)
    for(int i=0;i<n;++i) data[i]/=n;
  
  return data;
}


vector<complex<double>> FFT(vector<uint16_t> a, bool inverse=false){
  vector<complex<double>> a_complex(a.size());
  for(int i=0;i<(int)a.size();++i) a_complex[i]=complex<double>(a[i],0);
  return FFT(a_complex,inverse);
}

vector<complex<double>> FFT(vector<double> a, bool inverse=false){
  vector<complex<double>> a_complex(a.size());
  for(int i=0;i<(int)a.size();++i) a_complex[i]=complex<double>(a[i],0);
  return FFT(a_complex,inverse);
}

vector<double> hanning_window(vector<uint16_t>a){
  int L=a.size();
  vector<double> temp(L);
  for(int i=0;i<L;++i) temp[i]=(0.5-0.5*cos(2*M_PI*i/L))*a[i];
  return temp;
}

vector<double> hamming_window(vector<uint16_t>a){
  int L=a.size();
  vector<double> temp(L);
  for(int i=0;i<L;++i) temp[i]=(0.54-0.46*cos(2*M_PI*i/L))*a[i];
  return temp;
}

vector<double> blackman_window(vector<uint16_t>a){
  int L=a.size();
  vector<double> temp(L);
  for(int i=0;i<L;++i) temp[i]=(0.42-0.5*cos(2*M_PI*i/L)+0.08*cos(4*M_PI*i/L))*a[i];
  return temp;
}

#include "filter_par/LPF_f10_N101.hpp"
#include "filter_par/LPF_f15_N401.hpp"
#include "filter_par/HPF_f15_N401.hpp"

vector<double> FIR_filter(vector<double> data, vector<double> par){
  int M=par.size()-1;
  vector<double> ret;
  vector<double> u(M+1,0);
  for(size_t i=0;i<(int)data.size();++i){
    int xn=data[i];
    for(int m=0;m<M;++m) u[m]=par[m]*xn+u[m+1];
    u[M]=par[M]*xn;
    ret.push_back(u[0]);
  }
  return ret;
}
/*
void FIR_filter_2(vector<double> &data, int &board, int &ch, vector<double> &par, int &div, int num){
  int M=par.size()-1;
  //vector<double> ret;
  vector<double> u(M+1,0);
  int N=data.size()/div;
  for(size_t i=max(0,N*num-(M+1));i<N*(num+1);++i){
    int xn=data[i];
    for(int m=0;m<M;++m) u[m]=par[m]*xn+u[m+1];
    u[M]=par[M]*xn;
    if(N*num<=i) data_fil[board][ch][i]=u[0];
    //ret.push_back(u[0]);
  }
  //return ret;
}
*/
vector<double> differential(vector<double> data, int N){
  vector<double> ret;
  for(int i=0;i<(int)data.size();++i){
    double diff=0;
    int n=N;
    if((int)data.size()<=i+N/2){
      diff+=data[data.size()-1];
      --n;
    }
    else diff+=data[i+N/2];
    if(i-N/2<0){
      diff-=data[0];
      --n;
    }
    else diff-=data[i-N/2];
    ret.push_back(100*diff/n);
  }
  return ret;
}

vector<double> shaping_filter(vector<double> data, int N){
  vector<double> ret;
  ret=differential(data,21);
  for(int i=0;i<N;++i) ret=FIR_filter(ret,LPF_f15_N401_par);
  return ret;
}
