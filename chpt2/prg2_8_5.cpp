#include "stdafx.h"
#include "itpp/base/random.h"
using std::cin;
using std::cout;
using std::endl;
using namespace itpp;

int main() {


  cout<<"randu()="<<randu()<<endl;
  cout<<"randu(3,4)="<<randu(3,4)<<endl;
  cout<<"randi(3,9)="<<randi(3,9)<<endl;
  cout<<"randi(4,3,9)="<<randi(4,3,9)<<endl;
  cout<<"randi(2,3,3,9)="<<randi(2,3,3,9)<<endl;
  cout<<"randray(3,1.0)="<<randray(3,1.0)<<endl;
  cout<<"randexp(3,1.0)="<<randexp(3,1.0)<<endl;
  cout<<"randn()="<<randn()<<endl;

  cout<<"randn(2)="<<randn(2)<<endl;
  cout<<"randn(2,3)="<<randn(2,3)<<endl;
  cout<<"randn_c()="<<randn_c()<<endl;
  cout<<"randn_c(2)="<<randn_c(2)<<endl;
  cout<<"randn_c(2,3)="<<randn_c(2,3)<<endl;

  return 0;
  }
