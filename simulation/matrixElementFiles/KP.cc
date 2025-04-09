int GF{1};
double constants{(Power(AKP,2)*Power(GF,2)*Power(Vus,2))/(16*Power(fPi,2)*Power(Power(m2,2) - 
  Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)),2))};
double matrixElementBody{8*((m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
  Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar)) + 
 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
  Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
 2*m1*ma*Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
  Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
  Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
 2*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s,q,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
    Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))) + 
 m1*ma*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) - 
 Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar))*
  Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)))*
(Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
  (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
    2*Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
 Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
  Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))};