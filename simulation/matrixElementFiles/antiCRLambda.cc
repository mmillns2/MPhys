int GF{1};
double numeratorConstants{(Power(ACRLambda,2)*Power(D + 3*F,2)*Power(GF,2)*Power(Vus,2))};

double denominatorConstants{(16*Power(fPi,2)*Power(-Power(mLambda,2) + Pair(Momentum(s,p1,W,theta,thetaStar,phiStar)+
   Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar)+Momentum(s,q,W,theta,thetaStar,phiStar)),2))};

// Matrix element body
double diagonals{(-16*Power(D + 3*F,2)*Power(mLambda,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/9. + 
16*Power(mLambda,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
(16*Power(D + 3*F,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
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
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/9. + 
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
(16*Power(D + 3*F,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/9. + 
(8*Power(D + 3*F,2)*Power(Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)),2)*
(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  (m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Power(q,2),2)) + 
(8*Power(D + 3*F,2)*(Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)),2)*
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
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Power(q,2),2)) - 
(8*Power(D + 3*F,2)*Power(mLambda,2)*(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
   Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(9.*Power(Power(m2,2) - Power(q,2),2)) + 
(Power(kappaP,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (-2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
        8*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
         Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
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
         (4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
            Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
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
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (-4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
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
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     4*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
(Power(kappaP,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
(Power(kappaP,2)*Power(mLambda,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
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
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/Power(ma,2) + 
(Power(kappaP,2)*Power(mLambda,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
(Power(kappaP,2)*(-8*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
(Power(kappaP,2)*(-8*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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

double offDiagonals{(2*(-144*mLambda*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
     ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
16*Power(D + 3*F,2)*mLambda*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
144*mLambda*(m1*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) - 
16*Power(D + 3*F,2)*mLambda*(m1*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (-2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
     ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
48*(D + 3*F)*(Complex(0,-1)*Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
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
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
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
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
48*(D + 3*F)*(Complex(0,1)*Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
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
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
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
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
  Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
    Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
  2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
48*(D + 3*F)*mLambda*(2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
48*(D + 3*F)*mLambda*(2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
   Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
48*(D + 3*F)*Power(mLambda,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
  Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
   (Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
     2*Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) - 
48*(D + 3*F)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
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
16*Power(D + 3*F,2)*(4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
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
144*(4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
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
         Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))) - 
48*(D + 3*F)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
(Complex(0,6)*(D + 3*F)*kappaP*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
    Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) + 
(Complex(0,6)*(D + 3*F)*kappaP*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
    Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) + 
(Complex(0,6)*(D + 3*F)*kappaP*mLambda*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
(Complex(0,6)*(D + 3*F)*kappaP*mLambda*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
    Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(ma*(Power(m2,2) - Power(q,2))) - 
(Complex(0,12)*(D + 3*F)*kappaP*Power(mLambda,2)*(-(m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
        Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
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
    Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(ma*(Power(m2,2) - Power(q,2))) - 
(Complex(0,12)*(D + 3*F)*kappaP*mLambda*(-(m1*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
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
(8*Power(D + 3*F,2)*mLambda*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Power(q,2),2) + 
(8*Power(D + 3*F,2)*mLambda*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
  (2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    (m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
  (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Power(q,2),2) - 
(8*Power(D + 3*F,2)*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
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
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(Power(m2,2) - Power(q,2),2) - 
(Complex(0,24)*(D + 3*F)*ma*mLambda*(Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
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
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) - 
(Complex(0,24)*(D + 3*F)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
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
(Complex(0,48)*(D + 3*F)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
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
(36*kappaP*Power(mLambda,2)*(2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
(24*(D + 3*F)*kappaP*mLambda*(Complex(0,-1)*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma - 
(24*(D + 3*F)*kappaP*Power(mLambda,2)*(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
(36*kappaP*mLambda*(-2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
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
(36*kappaP*mLambda*(2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
(36*kappaP*mLambda*(-2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
(36*kappaP*(2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
(12*(D + 3*F)*kappaP*(4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
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
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma - 
(Complex(0,6)*(D + 3*F)*kappaP*mLambda*(Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
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
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(Complex(0,6)*(D + 3*F)*kappaP*mLambda*(-(Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
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
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma - 
(12*(D + 3*F)*kappaP*mLambda*(Complex(0,-2)*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),
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
(6*(D + 3*F)*kappaP*(Complex(0,2)*m1*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    Complex(0,2)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*ma*Eps(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*ma*Eps(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,2)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,1)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q1,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Complex(0,1)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q3,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
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
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*m1*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Complex(0,2)*ma*Eps(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar),
      Momentum(s,q2,W,theta,thetaStar,phiStar))*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
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
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/ma + 
(8*Power(D + 3*F,2)*mLambda*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
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
(8*Power(D + 3*F,2)*mLambda*(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
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
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/(Power(m2,2) - Power(q,2)) + 
(9*Power(kappaP,2)*mLambda*(2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
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
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
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
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))/Power(ma,2) + 
(9*Power(kappaP,2)*mLambda*(2*m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
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
(36*kappaP*Power(mLambda,2)*(-(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
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
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma - 
(36*kappaP*(4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    4*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
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
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
     (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) - 
    ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma - 
(8*Power(D + 3*F,2)*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
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
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/(Power(m2,2) - Power(q,2)) - 
(8*Power(D + 3*F,2)*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
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
(8*Power(D + 3*F,2)*Power(mLambda,2)*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
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
(8*Power(D + 3*F,2)*mLambda*(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
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
(8*Power(D + 3*F,2)*mLambda*(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
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
(9*Power(kappaP,2)*Power(mLambda,2)*(2*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
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
(9*Power(kappaP,2)*(-4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
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
(36*kappaP*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
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
(36*kappaP*(-2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
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
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma + 
(36*kappaP*(-2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
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
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
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
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     ((Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))) + 
    2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     (-(m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
          Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) - 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (-Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))*
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
       2*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
       m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma + 
(36*kappaP*mLambda*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
(36*kappaP*mLambda*(m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/ma - 
(9*Power(kappaP,2)*mLambda*(Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
     (-4*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       2*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     (m1*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
           (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
             2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))) + 
    m1*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) - 
       Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        (2*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     (2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       2*m1*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
        Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
       ma*(2*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
          Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
           Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
(9*Power(kappaP,2)*mLambda*(Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
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
(9*Power(kappaP,2)*mLambda*(Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
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
(9*Power(kappaP,2)*mLambda*(Power(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
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
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))/Power(ma,2) + 
(24*(D + 3*F)*kappaP*mLambda*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma + 
(24*(D + 3*F)*kappaP*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
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
(18*Power(kappaP,2)*mLambda*(-2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
(36*kappaP*(4*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
(24*(D + 3*F)*kappaP*(-2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
          Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))))/ma - 
(24*(D + 3*F)*kappaP*(-2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
(36*kappaP*(2*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
(36*kappaP*mLambda*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
(36*kappaP*mLambda*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
(36*kappaP*mLambda*(2*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
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
(8*Power(D + 3*F,2)*(m1*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
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
(8*Power(D + 3*F,2)*(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
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
(36*kappaP*(4*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
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
(9*Power(kappaP,2)*(-4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
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
(9*Power(kappaP,2)*(-4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
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
(9*Power(kappaP,2)*(-4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
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
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
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
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
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
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    m1*ma*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p1,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
     Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
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
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
(9*Power(kappaP,2)*(-4*Power(Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)),2)*
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
             Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2) + 
(9*Power(kappaP,2)*(4*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
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
              Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))))/Power(ma,2)))/9};