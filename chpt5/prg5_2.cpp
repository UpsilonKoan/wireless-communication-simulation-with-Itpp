#include <itpp/itoptim.h>
#include <iomanip>

using namespace std;
using namespace itpp;

double rosenbrock(const vec &x)
{
  double f1 = x(1) - sqr(x(0)), f2 = 1 - x(0);
  double F = 50 * sqr(f1) + 0.5 * sqr(f2) + 0.0;  //计算目标函数
  return F;
}

vec rosenbrock_gradient(const vec &x)
{
  double f1 = x(1) - sqr(x(0)), f2 = 1 - x(0);   //计算梯度
  vec g(2);
  g(0) = -200.0 * x(0) * f1 - f2;
  g(1) = 100.0 * f1;

  return g;
}

int main(void)
{
  cout << "=====================================" << endl;
  cout << "    Test of Numerical optimization " << endl;
  cout << "=====================================" << endl;

  cout << setprecision(6);
  {
    cout << "    Line_Search " << endl;
    cout << "--------------------" << endl;
    vec x0 = "-1.2 1";  //设定初始值
    double F0 = rosenbrock(x0);    //带入初始值获得的目标函数值
    vec g0 = rosenbrock_gradient(x0);//带入初始值获得的梯度
    vec h = -g0;

    Line_Search line;
    line.set_functions(rosenbrock, rosenbrock_gradient);

    cout << "x0 = " << x0 << endl;
    cout << "F0 = " << F0 << endl;
    cout << "g0 = " << g0 << endl;
    cout << "h = " << h << endl;

    double Fn;
    vec xn, gn;

    line.enable_trace();

    if (!line.search(x0, F0, g0, h, xn, Fn, gn))//soft线性搜索
      cout << "Line search failed" << endl;

    cout << "Soft: " << endl;
    cout << "xn = " << xn << endl;
    cout << "Fn = " << Fn << endl;
    cout << "gn = " << gn << endl;
    cout << "no_feval = " << line.get_no_function_evaluations() << endl;

    vec alpha_values, F_values, dF_values;
    line.get_trace(alpha_values, F_values, dF_values);//获得相关参数值

    cout << endl << "trace:" << endl;
    cout << "alpha = " << alpha_values << endl;
    cout << "F = " << F_values << endl;
    cout << "dF = " << dF_values << endl;

    line.set_method(Exact); //设定exact搜索方法
    if (!line.search(x0, F0, g0, h, xn, Fn, gn)) //exact 线性搜索
    cout << "Line search failed" << endl;

    cout << endl << "Exact: " << endl;
    cout << "xn = " << xn << endl;
    cout << "Fn = " << Fn << endl;
    cout << "gn = " << gn << endl;
    cout << "no_feval = " << line.get_no_function_evaluations() << endl;

    line.get_trace(alpha_values, F_values, dF_values);//获得相关参数值

    cout << endl << "trace:" << endl;
    cout << "alpha = " << alpha_values << endl;
    cout << "F = " << F_values << endl;
    cout << "dF = " << dF_values << endl;
  }

  {
    cout << endl << "    Newton_Search " << endl;
    cout << "--------------------" << endl;
    vec x0 = "-1.2 1";

    Newton_Search newton;
    newton.set_functions(rosenbrock, rosenbrock_gradient);
    newton.enable_trace();

    cout << "x0 = " << x0 << endl;

    vec xn, gn;

    if (!newton.search(x0, xn)) //Newton搜索法
      cout << "Newton search failed" << endl;

    cout << "xn = " << xn << endl;
    cout << "F = " << newton.get_function_value() << endl;
    cout << "norm(f') = " << newton.get_stop_1() << endl;
    cout << "norm(dx) = " << newton.get_stop_2() << endl;
    cout << "no_feval = " << newton.get_no_function_evaluations() << endl;
    cout << "no_iter = " << newton.get_no_iterations() << endl;

    Array<vec> xv;
    vec Fv, ngv, dv;
    newton.get_trace(xv, Fv, ngv, dv);//获得相关参数值

    cout << endl << "trace:" << endl;
    cout << "xv = " << xv << endl;
    cout << "Fv = " << Fv << endl;
    cout << "ngv = " << ngv << endl;
    cout << "dv = " << dv << endl;

    newton.disable_trace();

    cout << endl << "    fminunc " << endl;
    cout << "--------------------" << endl;
    xn = fminunc(rosenbrock, rosenbrock_gradient, x0);//使目标函数最小的最优解
    cout << "x0 = " << x0 << endl;
    cout << "xn = " << xn << endl;
  }

  return 0;

}
