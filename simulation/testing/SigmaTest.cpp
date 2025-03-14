#include <cmath>
#include <vector>
#include <iostream>
#include <string>

int GF= 1;

int mSigma =1200;
int ma= 940;
int m1= 940;
int m2= 600;

int vects= 1;
int Vus=1;
double D= 0.804;
double F= 0.463;
double ACRSigma= -D+F;
double kappaN= -1.913;
double kappaP= 1.793;
double fPi= 92.4;

std::string p1 = "p1";
std::string p2 = "p2";
std::string q1 = "q1";
std::string q2 = "q2";
std::string q3 = "q3";
std::string q = "q";

std::vector<double> Momentum(std::string type, int param)
{
  if(type == "p1")
  {
    return {950,3,-9,8};
  }
  else if(type == "p2")
  {
    return {70,8,88,-9};
  }
  else if(type == "q1")
  {
    return {330,8,-80,-7};
  }
  else if(type == "q2")
  {
    return {600,70,5,8};
  }
  else if(type == "q3")
  {
    return {40,-100,0,0};
  }
  else if(type == "q")
  {
    return {30,108,88,-9};
  }
  else
  {
    return {0};
  }
}

double Power(double arg, int exponent)
{
  return pow(arg, exponent);
}
double Power (std::string arg1, std::string arg2, int exponent)
{
  std::vector<double> vec1{};
  vec1= Momentum(arg1, vects);

  std::vector<double> vec2{};
  vec2= Momentum(arg2, vects);

  std::vector<double> vecTot{};
  vecTot[0]= vec1[0]-vec2[0];
  vecTot[1]= vec1[1]-vec2[1];
  vecTot[2]= vec1[2]-vec2[2];
  vecTot[3]= vec1[3]-vec2[3];

  return vecTot[0]*vecTot[0]-vecTot[1]*vecTot[1]-vecTot[2]*vecTot[2]-vecTot[3]*vecTot[3];
}

double Pair(std::vector<double> vec1, std::vector<double> vec2)
{
  return vec1[0]*vec2[0]-vec1[1]*vec2[1]-vec1[2]*vec2[2]-vec1[3]*vec2[3];
}

int main()
{
double CRSigma= (Power(ACRSigma,2)*Power(GF,2)*Power(Vus,2)*((32*(D - F)*((2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)))*
             (-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2) + 4*Pair(Momentum(p1,vects),Momentum(q2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               (2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                Pair(Momentum(q2,vects),Momentum(q2,vects))) + 2*((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                Pair(Momentum(p1,vects),Momentum(q2,vects)) + (ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
             (2*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/(Power(mSigma,2)-Power(p1,q2,2)) - 
       (32*Power(D - F,2)*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) + m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) - 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(2*ma*mSigma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
               Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*Pair(Momentum(p2,vects),Momentum(q3,vects))*((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                  ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) - 
               2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               Power(mSigma,2)*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/
               (Power(mSigma,2)-Power(p1,q2,2)) - (32*
          (4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*(Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
               m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               Power(mSigma,2)*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(-2*ma*mSigma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
               Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               m1*Pair(Momentum(p2,vects),Momentum(q3,vects))*((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                  ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
               2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))))/
        (Power(mSigma,2)-Power(p1,q2,2)) - (16*Power(D - F,2)*
          (-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             ((m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(q,vects),Momentum(q1,vects))*
                ((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
               Pair(Momentum(q,vects),Momentum(q,vects))*((-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                   Pair(Momentum(q1,vects),Momentum(q2,vects)) + (-(m1*ma) + m1*mSigma + 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(-2*Pair(Momentum(q,vects),Momentum(q1,vects))*
                (-2*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*(2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - 
                     Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
               Pair(Momentum(q,vects),Momentum(q,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                   (m1*ma + 2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p1,vects),Momentum(q1,vects))*(mSigma*(2*ma + mSigma) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  ma*(-2*mSigma*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))))))*
          (Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))/
        (Power(Power(m2,2) - Power(p2,q3,2),2)*(Power(mSigma,2)-Power(p1,q2,2))) + 
       (8*(2*kappaN + kappaP)*(-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects)) + 
               mSigma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               (m1*Pair(Momentum(p1,vects),Momentum(p2,vects)) + mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects)) + Pair(Momentum(p2,vects),Momentum(p2,vects))*
                ((-ma + mSigma)*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  m1*(-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
               (m1*(Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
                  (-2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects)) + (ma - mSigma)*Pair(Momentum(p2,vects),Momentum(q1,vects)))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
                (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   ((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) - 
                  2*m1*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) - 
               (-2*m1*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
                  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - 
                     Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p2,vects),Momentum(q1,vects))*((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                     ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/
        (ma*(Power(mSigma,2)-Power(p1,q2,2))) + 
       (8*(2*kappaN + kappaP)*(4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects)) + 
               mSigma*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               (m1*Pair(Momentum(p1,vects),Momentum(p2,vects)) + mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects)) + Pair(Momentum(p2,vects),Momentum(p2,vects))*
                ((-ma + mSigma)*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  m1*(-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
               (m1*(Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
                  (-2*m1*Pair(Momentum(p1,vects),Momentum(p2,vects)) + (ma - mSigma)*Pair(Momentum(p2,vects),Momentum(q1,vects)))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(-(m1*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*
                  (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))) - 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   ((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) - 
                  2*m1*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
               (-2*m1*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects)) + 
                  m1*Pair(Momentum(p1,vects),Momentum(p2,vects))*(2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - 
                     Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p2,vects),Momentum(q1,vects))*((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + 
                     ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))))*Pair(Momentum(q3,vects),Momentum(q3,vects)))))/
        (ma*(Power(mSigma,2)-Power(p1,q2,2))) + 
       (32*Power(D - F,2)*(-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) - 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*((m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q,vects),Momentum(q3,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
               ((m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
                ((-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q,vects),Momentum(q2,vects))*
                   Pair(Momentum(q1,vects),Momentum(q3,vects)) + 2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + 2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + Pair(Momentum(q,vects),Momentum(q3,vects))*
                   ((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
                     (m1*ma - m1*mSigma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
                  Power(mSigma,2)*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
                  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects))) + 
               ((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                   Pair(Momentum(q,vects),Momentum(q1,vects)) - Power(mSigma,2)*Pair(Momentum(p2,vects),Momentum(q,vects))*
                   Pair(Momentum(q1,vects),Momentum(q2,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q,vects))*
                   Pair(Momentum(q1,vects),Momentum(q2,vects)) - m1*ma*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  m1*mSigma*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  Pair(Momentum(p2,vects),Momentum(q1,vects))*((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                      Pair(Momentum(q,vects),Momentum(q2,vects)) - 2*Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))))*
                Pair(Momentum(q3,vects),Momentum(q3,vects))) + Pair(Momentum(q2,vects),Momentum(q2,vects))*
             (Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*
                (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(-(Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                     (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
                  Pair(Momentum(q,vects),Momentum(q3,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                      (m1*ma + 2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mSigma*(2*ma + mSigma) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     ma*(-2*mSigma*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects))))) + 
                  2*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                   (Pair(Momentum(q,vects),Momentum(q2,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
                     Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
               (-2*ma*mSigma*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) - 
                  Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
                  2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects)) + 
                  2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(q,vects),Momentum(q2,vects)) - 
                  Pair(Momentum(p1,vects),Momentum(q,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                   (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  Pair(Momentum(p2,vects),Momentum(q,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                      (m1*ma + 2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mSigma*(2*ma + mSigma) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     ma*(-2*mSigma*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects))))))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)))))/((Power(m2,2) - Power(p2,q3,2))*(Power(mSigma,2)-Power(p1,q2,2))) + 
       (4*Power(2*kappaN + kappaP,2)*(Pair(Momentum(q2,vects),Momentum(q2,vects))*
             (4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
               4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects)) - 2*m1*ma*Power(mSigma,2)*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               2*m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*m1*mSigma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*ma*mSigma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q3,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2)*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects)) + 4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q3,vects)) + 4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   ((ma + 2*mSigma)*Pair(Momentum(p1,vects),Momentum(p1,vects)) + ma*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))) + 
                  2*Pair(Momentum(q1,vects),Momentum(q3,vects))*(Pair(Momentum(p1,vects),Momentum(q3,vects))*
                      (2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) - 
                     2*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (2*ma*mSigma + Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
               4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Power(mSigma,2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*m1*mSigma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects))) - 4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (Pair(Momentum(p2,vects),Momentum(p2,vects))*(m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))) + 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + m1*mSigma*Pair(Momentum(q3,vects),Momentum(q3,vects))))\
             + 2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(Pair(Momentum(p2,vects),Momentum(p2,vects))*
                (m1*(-ma + mSigma)*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                  2*Pair(Momentum(q1,vects),Momentum(q3,vects))*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                     (-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q2,vects),Momentum(q3,vects)))) + 
               2*Pair(Momentum(p2,vects),Momentum(q1,vects))*((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*
                   Pair(Momentum(p2,vects),Momentum(q2,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(2*m1*(ma - mSigma)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) - 4*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + 2*Power(mSigma,2)*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - m1*ma*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
                  m1*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))))))/
        (Power(ma,2)*(Power(mSigma,2)-Power(p1,q2,2))) + 
       (4*Power(2*kappaN + kappaP,2)*(-4*Power(Pair(Momentum(p1,vects),Momentum(q2,vects)),2)*
             (2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (-Pair(Momentum(p2,vects),Momentum(q3,vects)) + Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*((2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(p2,vects),Momentum(q3,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
                  2*(m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(-2*m1*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) + 
                  (2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
            Pair(Momentum(q2,vects),Momentum(q2,vects))*(4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects)) + 4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects)) - 
               2*m1*ma*Power(mSigma,2)*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               2*m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*m1*mSigma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2) - 
               4*ma*mSigma*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q3,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects)) - 
               2*m1*ma*Power(Pair(Momentum(p2,vects),Momentum(q3,vects)),2)*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
               2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects)) + 4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q3,vects)) + 4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) - 
               4*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               4*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*m1*ma*Power(mSigma,2)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*m1*ma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 4*m1*mSigma*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*ma*mSigma*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + Power(mSigma,2)*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p1,vects),Momentum(q1,vects))*
                Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               2*ma*mSigma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) - 
               2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) + 2*m1*ma*Pair(Momentum(p2,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*
                Pair(Momentum(q3,vects),Momentum(q3,vects)) - Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   (mSigma*(2*ma + mSigma) + Pair(Momentum(p1,vects),Momentum(p1,vects)) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  Pair(Momentum(p2,vects),Momentum(q3,vects))*(Pair(Momentum(p1,vects),Momentum(p1,vects))*
                      (2*m1*ma + 4*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mSigma*(2*ma + mSigma) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     2*ma*(-(mSigma*Pair(Momentum(q1,vects),Momentum(q2,vects))) + m1*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects))))) - 
                  4*(ma*mSigma + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - 2*
                   (Pair(Momentum(p1,vects),Momentum(p1,vects))*(m1*ma + 2*m1*mSigma + Pair(Momentum(p1,vects),Momentum(q1,vects)) - 
                        2*Pair(Momentum(q1,vects),Momentum(q2,vects))) + 
                     Pair(Momentum(p1,vects),Momentum(q1,vects))*(mSigma*(2*ma + mSigma) - Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                     ma*(-2*mSigma*Pair(Momentum(q1,vects),Momentum(q2,vects)) + m1*(Power(mSigma,2) + Pair(Momentum(q2,vects),Momentum(q2,vects)))))*
                   Pair(Momentum(q3,vects),Momentum(q3,vects)))) + 
            2*Pair(Momentum(p1,vects),Momentum(q2,vects))*(2*Pair(Momentum(p2,vects),Momentum(q1,vects))*
                ((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(p2,vects),Momentum(q2,vects)) - 
                  2*Pair(Momentum(p1,vects),Momentum(p2,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                (Pair(Momentum(p2,vects),Momentum(q3,vects)) - Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(q3,vects))*(2*m1*(ma - mSigma)*Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) - 4*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q2,vects)) + 2*Power(mSigma,2)*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - 2*Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                   Pair(Momentum(q2,vects),Momentum(q3,vects)) - Power(mSigma,2)*Pair(Momentum(q1,vects),Momentum(q2,vects))*
                   Pair(Momentum(q3,vects),Momentum(q3,vects)) + Pair(Momentum(p1,vects),Momentum(p1,vects))*Pair(Momentum(q1,vects),Momentum(q2,vects))*
                   Pair(Momentum(q3,vects),Momentum(q3,vects)) - 2*m1*ma*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
                  2*m1*mSigma*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects)) + 
                  2*Pair(Momentum(p1,vects),Momentum(q1,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects))*Pair(Momentum(q3,vects),Momentum(q3,vects))) + 
               Pair(Momentum(p2,vects),Momentum(p2,vects))*(Pair(Momentum(p2,vects),Momentum(q3,vects))*
                   ((-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
                     2*(m1*(-ma + mSigma) + Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects))) + 
                  2*(2*Pair(Momentum(p1,vects),Momentum(q3,vects))*Pair(Momentum(q1,vects),Momentum(q3,vects))*Pair(Momentum(q2,vects),Momentum(q2,vects)) + 
                     (-Power(mSigma,2) + Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q3,vects))*
                      Pair(Momentum(q2,vects),Momentum(q3,vects)) + 
                     ((Power(mSigma,2) - Pair(Momentum(p1,vects),Momentum(p1,vects)))*Pair(Momentum(q1,vects),Momentum(q2,vects)) + 
                        (m1*ma - m1*mSigma - 2*Pair(Momentum(p1,vects),Momentum(q1,vects)))*Pair(Momentum(q2,vects),Momentum(q2,vects)))*
                      Pair(Momentum(q3,vects),Momentum(q3,vects)))))))/(Power(ma,2)*(Power(mSigma,2)-Power(p1,q2,2)))))/(8.*Power(fPi,2));
  std::cout<<CRSigma;
                      return 0;
}