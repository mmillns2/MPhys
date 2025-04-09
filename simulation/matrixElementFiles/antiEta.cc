int GF{1};
double constants{(Power(AEta,2)*Power(D - 3*F,2)*Power(GF,2)*Power(ma,2)*Power(Vus,2))/(4*Power(fPi,2)*Power(-Power(mEta,2)+ 
  Pair(Momentum(s,q,W,theta,thetaStar,phiStar)-Momentum(s,q2,W,theta,thetaStar,phiStar),
  Momentum(s,q,W,theta,thetaStar,phiStar)-Momentum(s,q2,W,theta,thetaStar,phiStar)),2))};
double matrixElementBody{-8*(m1*ma - Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)))*
   (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
      (-4*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     4*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q2,W,theta,thetaStar,phiStar))*
      (2*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
        Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))) + 
     Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      (Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) - 
        2*(2*Pair(Momentum(s,q2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
           Pair(Momentum(s,q3,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)))))};