int GF{1};
double numeratorConstants{(Power(ACRLambda,2)*Power(GF,2)*Power(Vus,2))};

double denominatorConstants{(16*Power(fPi,2)*Power(-Power(mLambda,2) + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar)-
   Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar)-Momentum(s,q2,W,theta,thetaStar,phiStar)),2))};

// Matrix element body
double diagonals{((32*Power(D + 3*F,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2))/9. + 
32*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
 Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
 Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2) + 
32*Power(mLambda,2)*(-((Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
 (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
(32*Power(D + 3*F,2)*Power(mLambda,2)*(-((Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
       m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/9. + 
32*(4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
 (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
 (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
   m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
 Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
 Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
 (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
(32*Power(D + 3*F,2)*(4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/9. + 
(16*Power(D + 3*F,2)*((m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Power(q,2),2)) + 
(16*Power(D + 3*F,2)*Power(mLambda,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  ((m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Power(q,2),2)) - 
(16*Power(D + 3*F,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   (4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   (-4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
     m1*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Power(q,2),2)) + 
(2*Power(kappaP,2)*Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*m1*ma*Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
(2*Power(kappaP,2)*Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*m1*ma*Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
(2*Power(kappaP,2)*Power(mLambda,2)*(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*m1*ma*Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
           2*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
(2*Power(kappaP,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*m1*ma*Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
           2*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
(2*Power(kappaP,2)*Power(mLambda,2)*(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*m1*ma*Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
(2*Power(kappaP,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*m1*ma*Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
  4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2))/2};

double offDiagonals{((-64*(D + 3*F)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2))/3. - 
64*mLambda*(2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
 Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
(ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
   m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
   ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
 Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
(64*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  (ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))))/9. + 
(64*Power(D + 3*F,2)*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/9. + 
64*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
 Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
 Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
 Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
(64*(D + 3*F)*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-(m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/3. - 
(64*(D + 3*F)*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  (-(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/3. + 
(64*(D + 3*F)*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/3. + 
(64*(D + 3*F)*mLambda*(ma*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/3. + 
(64*(D + 3*F)*mLambda*(ma*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/3. + 
(32*(D + 3*F)*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*m1*ma*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/3. + 
(64*(D + 3*F)*Power(mLambda,2)*((Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/3. - 
(64*(D + 3*F)*(4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   (-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/3. - 
(64*Power(D + 3*F,2)*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/9. - 
64*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
 (2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
   m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
   2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
 (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
    Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
(16*kappaP*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(16*kappaP*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(16*kappaP*mLambda*(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-(m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) - 
  2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(16*kappaP*(-4*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(16*kappaP*Power(mLambda,2)*(2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(16*kappaP*Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma + 
(16*kappaP*Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma - 
(Complex(0,2.6666666666666665)*(D + 3*F)*kappaP*(2*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2) - 
  2*ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2) + 
  2*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*((m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar)) + 
        2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        ma*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
  m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(ma*(Power(m2,2) - Power(q,2))) - 
(Complex(0,2.6666666666666665)*(D + 3*F)*kappaP*(-(ma*Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  2*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2) + 
  2*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*(-(ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        ma*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*((m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar)) + 
        2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        ma*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
  m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(ma*(Power(m2,2) - Power(q,2))) + 
(16*kappaP*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-(ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (-2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(16*kappaP*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-(ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma - 
(32*Power(D + 3*F,2)*mLambda*(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
  2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Power(q,2),2)) - 
(32*Power(D + 3*F,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   ((m1*ma - 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
     4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   (-2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))))*
Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Power(q,2),2)) + 
(32*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
  (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
     2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Power(q,2),2)) + 
(32*Power(D + 3*F,2)*Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   ((m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  ((m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(9.*(Power(m2,2) - Power(q,2))) + 
(16*kappaP*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(16*kappaP*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma - 
(Complex(0,10.666666666666666)*(D + 3*F)*ma*mLambda*(Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
(Complex(0,10.666666666666666)*(D + 3*F)*(Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (-2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) - 
(4*(D + 3*F)*kappaP*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(Complex(0,-1)*ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  4*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Complex(0,1)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  Complex(0,1)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  4*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  4*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  4*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(3.*ma) + 
(Complex(0,2.6666666666666665)*(D + 3*F)*kappaP*mLambda*(Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*((m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(ma*(Power(m2,2) - Power(q,2))) + 
(Complex(0,2.6666666666666665)*(D + 3*F)*kappaP*mLambda*(Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*(m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(ma*(Power(m2,2) - Power(q,2))) - 
(4*Power(kappaP,2)*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(ma,2) - 
(Complex(0,1.3333333333333333)*(D + 3*F)*kappaP*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(Power(m2,2) - Power(q,2)) + 
(Complex(0,1.3333333333333333)*(D + 3*F)*kappaP*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(ma*(Power(m2,2) - Power(q,2))) - 
(4*Power(kappaP,2)*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
(16*(D + 3*F)*kappaP*Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  ma*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(3.*ma) - 
(Complex(0,2.6666666666666665)*(D + 3*F)*kappaP*mLambda*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-(Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(ma*(Power(m2,2) - Power(q,2))) + 
(Complex(0,2.6666666666666665)*(D + 3*F)*kappaP*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-(Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(ma*(Power(m2,2) - Power(q,2))) - 
(Complex(0,2.6666666666666665)*(D + 3*F)*kappaP*mLambda*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-(Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(ma*(Power(m2,2) - Power(q,2))) + 
(Complex(0,1.3333333333333333)*(D + 3*F)*kappaP*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-(Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(Power(m2,2) - Power(q,2)) + 
(Complex(0,1.3333333333333333)*(D + 3*F)*kappaP*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-(Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
    Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*(-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(ma*(Power(m2,2) - Power(q,2))) - 
(4*Power(kappaP,2)*Power(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   ((-2*m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     2*(m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     (-2*m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) - 
(32*Power(D + 3*F,2)*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(9.*(Power(m2,2) - Power(q,2))) - 
(32*Power(D + 3*F,2)*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-(m1*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     m1*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(9.*(Power(m2,2) - Power(q,2))) + 
(4*Power(kappaP,2)*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  ma*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) - 
(32*Power(D + 3*F,2)*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      ((m1*ma - 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     ((m1*ma - 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
        2*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(9.*(Power(m2,2) - Power(q,2))) - 
(32*Power(D + 3*F,2)*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-(m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      ((m1*ma - 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     ((m1*ma - 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
        2*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(9.*(Power(m2,2) - Power(q,2))) + 
(32*Power(D + 3*F,2)*Power(mLambda,2)*(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      ((m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ((m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(9.*(Power(m2,2) - Power(q,2))) + 
(4*Power(kappaP,2)*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  ma*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) - 
(32*(D + 3*F)*kappaP*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*(-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   ((-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  4*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/(3.*ma) + 
(16*kappaP*Power(mLambda,2)*(2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(16*kappaP*(4*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(16*kappaP*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(16*kappaP*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(32*(D + 3*F)*kappaP*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*(-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/(3.*ma) + 
(32*(D + 3*F)*kappaP*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*(-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/(3.*ma) + 
(32*(D + 3*F)*kappaP*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/(3.*ma) + 
(32*(D + 3*F)*kappaP*mLambda*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/(3.*ma) + 
(32*(D + 3*F)*kappaP*Power(mLambda,2)*(2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   ((-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/(3.*ma) + 
(4*Power(kappaP,2)*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (-4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      ((-(m1*ma) + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        (m1*ma - 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        (-(m1*ma) + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) - 
(4*Power(kappaP,2)*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(-2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        (-(m1*ma) + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      ((m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        (-(m1*ma) + 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
(32*Power(D + 3*F,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      ((m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     (-((m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/(9.*(Power(m2,2) - Power(q,2))) + 
(4*Power(kappaP,2)*mLambda*(-2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
     m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) - 
(4*Power(kappaP,2)*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*m1*ma*Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
             2*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) - 
(16*kappaP*mLambda*(2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(16*kappaP*mLambda*(2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(4*Power(kappaP,2)*Power(mLambda,2)*(Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      ((2*m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        2*(m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        (2*m1*ma + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
(4*Power(kappaP,2)*mLambda*(-2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) - 
(4*Power(kappaP,2)*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*m1*ma*Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
           2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
           2*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
(32*(D + 3*F)*kappaP*mLambda*(-2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   ((-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/(3.*ma) + 
(32*(D + 3*F)*kappaP*mLambda*(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   ((-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/(3.*ma) - 
(16*kappaP*mLambda*(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      ((-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/ma - 
(32*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      (-(m1*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        m1*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     ma*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/(9.*(Power(m2,2) - Power(q,2))) + 
(32*Power(D + 3*F,2)*mLambda*(2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      (-(m1*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        m1*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
     ma*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         (-(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/(9.*(Power(m2,2) - Power(q,2))) + 
(4*Power(kappaP,2)*mLambda*(-2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  4*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
            (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
              2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
(4*Power(kappaP,2)*mLambda*(2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
  4*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     ma*(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
         (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
            (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
              2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
(4*Power(kappaP,2)*(-(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     (4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2) - 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
  2*(4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2)*
      Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      (4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)),2) - 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
      ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
           m1*ma*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
        Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
              Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
            (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
               Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
              2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
               Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
              2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
               Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
           Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
            (-2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
               Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
              Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
               Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))))/Power(ma,2))/2};