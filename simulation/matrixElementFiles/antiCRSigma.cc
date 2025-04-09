int GF{1};
double numeratorConstants{(Power(ACRSigma,2)*Power(D - F,2)*Power(GF,2)*Power(Vus,2))};

double denominatorConstants{(4*Power(fPi,2)*Power(-Power(mSigma,2) + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar)+
   Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar)+Momentum(s,q,W,theta,thetaStar,phiStar)),2))};

// Matrix element body
double diagonals{16*Power(D - F,2)*Power(mSigma,2)*(-(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
 (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
 (-(Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
   2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
16*Power(mSigma,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
 Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
 (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
 (-(Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
   2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
16*Power(D - F,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
 (-(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
 (8*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
   4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
    (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
      2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) + 
16*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
 (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
 (8*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
   4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
    (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
      2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) + 
16*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
 Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
 Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
 (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
    (4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
      2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
 (4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
   2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
    (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
      2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) - 
16*Power(D - F,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
 Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
 Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
 (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
    (-4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
      2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
 (-4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
   2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
    (-(Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
      2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) + 
(8*Power(D - F,2)*Power(Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)),2)*
(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  (m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Power(q,2),2) + 
(8*Power(D - F,2)*(Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)),2)*
   (8*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   (-2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     (m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   (-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Power(q,2),2) - 
(8*Power(D - F,2)*Power(mSigma,2)*(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  (m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   (-4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Power(q,2),2) - 
(Power(2*kappaN + kappaP,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        8*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         (-4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
           2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
        4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-8*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
     (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         (4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) - 
(Power(2*kappaN + kappaP,2)*Power(mSigma,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        4*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   ((Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        8*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         (-4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
           2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
        4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
     (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         (4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*Power(mSigma,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        4*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        4*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   ((Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*(-8*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
           4*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-((Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*(-8*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
           4*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      ((Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2)};

double offDiagonals{2*(-16*mSigma*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
16*Power(D - F,2)*mSigma*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
16*Power(D - F,2)*mSigma*(-(m1*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) - 
16*mSigma*(m1*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
16*(D - F)*(Complex(0,-1)*Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
16*(D - F)*(Complex(0,1)*Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Complex(0,1)*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
16*(D - F)*mSigma*(2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
16*(D - F)*mSigma*(2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
16*(D - F)*Power(mSigma,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
16*(D - F)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
16*Power(D - F,2)*(4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   (-(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) + 
16*(4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) + 
16*(D - F)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      (-4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) + 
(Complex(0,2)*(D - F)*(2*kappaN + kappaP)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
  (-(Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
        Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
    2*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) + 
(Complex(0,2)*(D - F)*(2*kappaN + kappaP)*mSigma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
  (-(Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
        Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))) - 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(ma*(Power(m2,2) - Power(q,2))) - 
(Complex(0,2)*(D - F)*(2*kappaN + kappaP)*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
  (Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
    2*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) - 
(Complex(0,2)*(D - F)*(2*kappaN + kappaP)*mSigma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
  (Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(ma*(Power(m2,2) - Power(q,2))) + 
(Complex(0,4)*(D - F)*(2*kappaN + kappaP)*Power(mSigma,2)*(-(m1*
       Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
        Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    ma*Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
    m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(ma*(Power(m2,2) - Power(q,2))) + 
(Complex(0,4)*(D - F)*(2*kappaN + kappaP)*mSigma*(-(m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
        Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
    m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(ma*(Power(m2,2) - Power(q,2))) - 
(8*Power(D - F,2)*mSigma*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    (2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)),2) - 
       (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar)) + 
          ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Power(q,2),2) - 
(8*Power(D - F,2)*mSigma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
  (2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Power(q,2),2) + 
(8*Power(D - F,2)*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       (m1*ma - 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))))*
  (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Power(q,2),2) + 
(Complex(0,8)*(D - F)*ma*mSigma*(Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) + 
(Complex(0,8)*(D - F)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
  (Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) + 
(Complex(0,16)*(D - F)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
  (Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) - 
(4*(2*kappaN + kappaP)*Power(mSigma,2)*(2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(8*(D - F)*(2*kappaN + kappaP)*mSigma*(Complex(0,1)*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(8*(D - F)*(2*kappaN + kappaP)*Power(mSigma,2)*(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(4*(2*kappaN + kappaP)*mSigma*(-2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     (-(Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
          Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       (m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(4*(2*kappaN + kappaP)*mSigma*(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(4*(2*kappaN + kappaP)*mSigma*(-2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
       2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
       2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(4*(2*kappaN + kappaP)*(2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(4*(2*kappaN + kappaP)*(-4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma - 
(4*(D - F)*(2*kappaN + kappaP)*(4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(Complex(0,2)*(D - F)*(2*kappaN + kappaP)*mSigma*(Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*((Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    2*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma - 
(Complex(0,2)*(D - F)*(2*kappaN + kappaP)*mSigma*(-(Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
        Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    2*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar))*Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2) + 
    Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*(-2*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
    Complex(0,4)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(4*(D - F)*(2*kappaN + kappaP)*mSigma*(Complex(0,-2)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    Complex(0,1)*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(2*(D - F)*(2*kappaN + kappaP)*(Complex(0,-2)*m1*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Complex(0,2)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*ma*Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma - 
(8*Power(D - F,2)*mSigma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
  (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (-(ma*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
          Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       2*ma*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    (-2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       (-(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) + 
          ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) - 
(8*Power(D - F,2)*mSigma*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*(-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    (2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        (-4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) - 
(Power(2*kappaN + kappaP,2)*mSigma*(2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    4*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*mSigma*(2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    4*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(ma,2) + 
(4*(2*kappaN + kappaP)*Power(mSigma,2)*(-(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma + 
(8*Power(D - F,2)*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
  ((m1*ma - 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(Power(m2,2) - Power(q,2)) + 
(8*Power(D - F,2)*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
  (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(Power(m2,2) - Power(q,2)) + 
(8*Power(D - F,2)*Power(mSigma,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(Power(m2,2) - Power(q,2)) + 
(8*Power(D - F,2)*mSigma*(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (-2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       (-2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
          ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(Power(m2,2) - Power(q,2)) + 
(8*Power(D - F,2)*mSigma*(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
          ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(Power(m2,2) - Power(q,2)) + 
(Power(2*kappaN + kappaP,2)*Power(mSigma,2)*(2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*m1*ma*Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*(-4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    8*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    8*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*m1*ma*Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        (-4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    8*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) - 
(4*(2*kappaN + kappaP)*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (-2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma - 
(4*(2*kappaN + kappaP)*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (-2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma + 
(4*(2*kappaN + kappaP)*mSigma*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
       2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
       2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma + 
(4*(2*kappaN + kappaP)*mSigma*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma + 
(Power(2*kappaN + kappaP,2)*mSigma*(Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     (4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (m1*(-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
             (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
               2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
          Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
    m1*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (-2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*(-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*mSigma*(Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     (4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (m1*(-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
             Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
             2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
    m1*(-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (-2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*(-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*mSigma*(Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     (4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    (-2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
             2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (-2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       (-4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
             2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*mSigma*(Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     (4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    (2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*(-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
             Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
             Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (-2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
             2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (-2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       (-4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
             2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) - 
(8*(D - F)*(2*kappaN + kappaP)*mSigma*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (-((m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
          (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
            Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
          Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       (m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma - 
(8*(D - F)*(2*kappaN + kappaP)*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    4*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(2*Power(2*kappaN + kappaP,2)*mSigma*(-2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) - 
(4*(2*kappaN + kappaP)*(4*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma - 
(4*(2*kappaN + kappaP)*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(8*(D - F)*(2*kappaN + kappaP)*(-2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     ((-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(8*(D - F)*(2*kappaN + kappaP)*(-2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     ((-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma - 
(4*(2*kappaN + kappaP)*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(4*(2*kappaN + kappaP)*mSigma*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (-((Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
             (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
               Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
            m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
          Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        ((m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(4*(2*kappaN + kappaP)*mSigma*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (-((Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
             (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
               Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
            m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
          Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     (-(Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
          Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        ((m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(4*(2*kappaN + kappaP)*mSigma*(2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        ((m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          (m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma - 
(8*Power(D - F,2)*(m1*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)),2)*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          4*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       2*(-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
       (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           (-4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
             2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
       2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/(Power(m2,2) - Power(q,2)) + 
(8*Power(D - F,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
       2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
             Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (8*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
             Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/(Power(m2,2) - Power(q,2)) + 
(4*(2*kappaN + kappaP)*(4*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (-4*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/ma + 
(Power(2*kappaN + kappaP,2)*(-4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (2*m1*ma*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
             4*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (-((Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
                Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
               2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
                Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
             (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
               Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
          2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*(-4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
             4*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (-((Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
                Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
               2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
                Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
             (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
               Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
          2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*(-4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
             2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
             4*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        ((Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) - 
(Power(2*kappaN + kappaP,2)*(4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
             2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (-2*m1*ma*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             4*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           (-2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
             4*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (-((Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
                Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
               2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
                Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
             (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
               Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
          4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
(Power(2*kappaN + kappaP,2)*(4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
       4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     (-2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
          2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
          (2*m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        ((2*m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          2*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
              (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
                 Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
                2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
                 Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
             (m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
              Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2))};