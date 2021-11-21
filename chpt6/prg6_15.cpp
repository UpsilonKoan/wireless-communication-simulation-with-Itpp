#include <itpp/itcomm.h>
#include <iomanip>

using namespace std;
using namespace itpp;

int main()
{  
  ND_UPSK chan;
  chan.set_M(1, 8); // scalar channel, 8-PSK (3 bits per symbol)
  cout << chan << endl;
  bvec b = randb(3);
  cvec x = chan.modulate_bits(b);
  cmat H = "1.0+1.0i";      // scalar channel
  double sigma2 = 0.01;
  cvec y= H*x + sqrt(sigma2)*randn_c(); // transmit vector x
  QLLRvec llr;
  QLLRvec llr_ap = zeros_i(3);
  chan.demodulate_soft_bits(y, H, sigma2, llr_ap, llr);
  cout << "True bits:" << b << endl;
  cout << "LLRs:" << chan.get_llrcalc().to_double(llr) << endl;
 }
