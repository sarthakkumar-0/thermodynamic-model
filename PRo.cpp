#include<bits/stdc++.h>
#include<fstream>
using namespace std;
#define initialGuess  1.0
#define tolerance  0.00001
//right all
double A1(double tt){
    // double tt = t - 273.15;
    double ans = (3.2891 - (2.391 *pow(10,-3)*tt) +( 2.8446 * pow(10,-4)* pow(tt,2)) -(2.82 *pow(10,-6)*pow(tt,3))   + (8.477 * pow(10,-9)*pow(tt,4)));
return ans;}   
double A2 (double tt){
    // double tt = t-273.15;
    double ans = ((6.245 * pow(10,-5)) - (3.913 *pow(10,-6)* tt) - (3.499 * pow(10,-8)* pow(tt,2)) + (7.942 * pow(10,-10)* pow(tt,3))- (3.299 * pow(10,-12)* pow(tt,4)));
    return ans;
}
double B ( double tt){
    // double tt = t - 273.15;
    double ans = (19654.32 + (147.037*tt) - (2.2155* pow(tt,2)) + (1.0478* pow(10,-2)* pow(tt,3)) -( 2.2789* pow(10,-5)* pow(tt,4)));
    return ans ;
}
double V0 ( double tt){
    // double tt = t-273.15;
    double num = 1 + 18.1597*pow(10,-3)*tt;
    double den = 0.9998 + (18.2249*pow(10,-3)*tt )-( 7.9222 *pow(10,-6)* pow(tt,2) )- (55.4485 *pow(10,-9)* pow(tt,3)) + (149.7562 *pow(10,-12)* pow(tt,4)) - (393.2952 *pow(10,-15)* pow(tt,5));
    double ans = num/den;
    return ans;
}   
double den( double p , double t){
    double num = V0(t)*p;
    double deno = B(t) + A1(t)*p +A2(t)*pow(p,2);
    double k = V0(t)- num/deno;
    double ans = 1/k;
    return ans;
}

//right both
double Ps(double t){
    double tc = 647.096;
double a = 1-((t+273.15)/647.096);
double k = (tc/(t+273.15))* (-7.8595178*(a) + 1.8440825* pow(a,1.5) - 11.786649* pow(a,3)+22.680741* pow( a,3.5) - 15.9618719* pow ( a,4) + 1.8012250* pow( a,7.5) );
double ans = 220.6 * exp(k);
return ans;
}
double fh20 (double p,double t){
double num = 18.0152 * (p-Ps(t));
double deno = den( p,t) * 83.14 * (t+273.15);
double ans = (Ps(t)* exp(num/deno));
return ans ;
}
//all good
double aco2( double t){
    double m =  0.37464 + (1.54226 * 0.224) -( 0.26992 * 0.224 * 0.224);
     double tt = (t+273.15)/304.2;
     double d = (pow((83.14*304.2),2))/73.8;
     double ans = 0.457236 * d *pow((1.00+( m*( 1-pow(tt,0.5)))),2);
     return ans;
}
double Aco2 ( double t,double p){
    double c = pow((83.14*(t+273.15)),2);
    double ans = (aco2(t)* p)/c;
    return ans;
}
double Bco2 (double t,double p){
    double b = (0.077796 * 83.14*304.2)/73.8;
    double ans = (b * p)/(83.14*(t+273.15));
    return ans;
}

// all good 
double cubicFunction(double x, double a, double b, double c, double d) {
    return (a * x * x * x + b * x * x + c * x + d);
}
double cubicFunctionDerivative(double x, double a, double b, double c) {
    return (3 * a * x * x + 2 * b * x + c);
}
vector<double> newtonRaphson(double t,double p) {
   vector<double> roots;
double A = Aco2(t,p);
double B = Bco2(t,p);
double a =1;
double b = -( 1-B);
double c =  (A - 2*B -3*B*B );
double d = -(A*B - B*B - B*B*B);
    double x = initialGuess;
    double f_value, f_derivative_value;

    do {
        f_value = cubicFunction(x, a, b, c, d);
        f_derivative_value = cubicFunctionDerivative(x, a, b, c);
        x = x - f_value / f_derivative_value;
    } while (fabs(f_value) > tolerance);

    roots.push_back(x);

    // Now perform synthetic division to find the quadratic equation from which we can find remaining roots
    double new_a = a;
    double new_b = a * x + b;
    double new_c = new_b * x + c;
    double discriminant = new_b * new_b - 4 * new_a * new_c;

    if (discriminant >= 0) {
        double root1 = (-new_b + sqrt(discriminant)) / (2 * new_a);
        double root2 = (-new_b - sqrt(discriminant)) / (2 * new_a);

        roots.push_back(root1);
        roots.push_back(root2);
    }

    return roots;
}
double choosez (vector<double>v,double t,double p){
    if (v.size()==1){
        return v[0];
    }
    else{
   double Zg = v[2];
   double Zl = v[0];
   double c = (Zg - Zl + (log((Zl-Bco2(t,p))/(Zg-Bco2(t,p)))) - (Aco2(t,p)/(Bco2(t,p)*((-2)*pow(2,0.5))))* log(((Zl + (1+pow(2,0.5))*Bco2(t,p))/(Zl + (1-pow(2,0.5))*Bco2(t,p)))* ((Zg+ (1-pow(2,0.5))*Bco2(t,p))/(Zg+ (1+pow(2,0.5))*Bco2(t,p)))));
   if (c<0){
    return Zg;
   }
   else {
    return Zl;
   }
}}
// all good
double phi(double t,double p ,vector<double>v){
double z = choosez(v,t,p);
double k = (z-1) - log(z - Bco2(t,p)) + ((Aco2(t,p))/(Bco2(t,p)* 2*pow(2,0.5)))* log(((z+ (1-pow(2,0.5))*Bco2(t,p))/(z+ (1+pow(2,0.5))*Bco2(t,p))));
double ans = exp(k);
    return ans;
}

// right delb  hi right but little differ
double delB ( double t){

    double k = pow((1000/(t+273.15)),0.5);
    double ans = -5.279063 + 6.187967* k;
    return ans;
}
double hi(double t,double p){
    double n = -0.114535;
    double pp = den(p,t);
    double dd = ((83.14) * (t + 273.15)* pp)/(18.015);

    double ans = ((1-n)* log(fh20(p,t))) +  (n * log(dd)) + (2*pp*delB(t));
return exp(ans);
}

//right
double y (double t,double p){
    double tt = t+273.15;
    double c1 = -0.0652869;
    double c2= 1.6790636E-04;
    double c3 = 40.838951;
    double c6 = -3.9266518E-02;
    double c8 = 2.1157167E-02;
    double c9 = 6.5486487E-06;
    double ans = c1  + (c2 * tt) + (c3/tt )+ (c6* p/tt )+ (c8 *( p/(630 - tt))) + (c9 * tt* log(p));
   return ans;
   }
double e (double t,double p){
    double tt = t+273.15;
    double c1 = -1.144624E-02;
    double c2= 2.8274958E-05;
    double c6 = 1.3980876E-02;
    double c8 = -1.4349005E-02;
    double ans = c1  + (c2 * tt) +  (c6* p/tt )+ (c8 *( p/(630 - tt))) ;
   return ans;
   }
double gamai( double t,double p, double m){
    double l = 2* m * y(t,p) + 2*m*m*e(t,p);
    double ans = exp(l);
    return ans;
}
double kco2(double t,double p, double m,vector<double>v){
    double ans = (hi(t,p) *gamai(t,p,m))/(p * phi(t,p,v));
    return ans;
}

double koh20(double t){
    double ans = pow(10, -2.209 + (3.097 *pow(10,-2)*t )- (1.098 * pow(10,-4)*t*t) + (2.048* pow(10,-7)*t*t*t) );

    return ans;
}
double kh20 (double t, double p ){
    double d = ((koh20(t))/(fh20(p,t) * p))* exp (((p-1)*18.18)/(83.14 * (t+273.15)));
    return d ;
}

double yh20(double t, double p, double m , vector<double>v){
    double c = 1/kco2(t,p,m,v);
    double h = 1/kh20(t,p);
    double ans = (1- c)/(h-c);
    return ans;
}
double yco2(double t, double p, double m , vector<double>v){
    double ans= 1.0/(1+yh20(t,p,m,v));
    return ans;
}
double xco2(double t, double p, double m , vector<double>v){
    double ans = yco2(t,p,m,v)/kco2(t,p,m,v);
    return ans ;
}


int main(){
ifstream input("Code_Input.txt");
double t0,p,m;
input>>t0>>p>>m;

double t = t0-273.15;

vector<double>v = newtonRaphson(t,p);
sort(v.begin(),v.end());

ofstream codeout("Output_Code.txt");

// codeout <<"delB: "<< delB(t) <<endl;
// codeout <<"den: "<<den(p,t)<<endl;
// codeout <<"fh20: "<<fh20(p,t)<<endl;
// codeout <<"hco2: "<<<<endl;
// codeout <<"gamai: "<<gamai(t,p,m)<<endl;
// codeout <<"z: "<<choosez(v,t,p)<<endl;
// codeout <<"phi: "<<phi(t,p,v)<<endl;
// codeout <<"kco29999: "<<kco2(t,p,m,v)<<endl;
// codeout <<"kh20: "<<<<endl;
// codeout <<"yh20: "<<<<endl;
// codeout <<"yco2: "<<yco2(t,p,m,v)<<endl;
// codeout <<"kh20: "<<<<endl;
// codeout <<"ac02: "<<Aco2(t,p)<<endl;
// codeout <<"bc02: "<<Bco2(t,p)<<endl;
codeout<<"Fugacity Coefficient_CO2 = "<<phi(t,p,v)<<endl;
codeout<<"Henry's Constant_CO2 = "<<hi(t,p)<<endl;
codeout<<"Equilibrium Constant_H2O = "<<kh20(t,p)<<endl;
codeout<<"Vapour Phase Mole Fraction_H2O = "<<yh20(t,p,m,v)<<endl;
codeout<<"Liquid Phase Mole Fraction_CO2 = "<<xco2(t,p,m,v)<<endl;

codeout.close();

// cout<<"Fugacity coefficient of CO2 = "<<phi(t,p,v)<<endl;
// cout<<"Henry's Constant for CO2 = "<<hi(t,p)<<endl;
// cout<<"Equilibrium Constant for H2O = " <<kh20(t,p)<<endl;
// cout<<"Vapour Phase Mlle Fractrion for H20 = "<<yh20(t,p,m,v)<<endl;
// cout<<"Liquid Phase Mole Fraction for CO2 = "<< xco2(t,p,m,v)<<endl;

    return 0;
}
