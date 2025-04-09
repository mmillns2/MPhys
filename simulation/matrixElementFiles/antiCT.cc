int GF{1};
double constants{(Power(ACT,2)*Power(GF,2)*Power(Vus,2))/(4*Power(fPi,2))};
double matrixElementBody{64*(Power(-1 + BCT,2)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))*
   Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q1,W,theta,thetaStar,phiStar)) + 
  (1 + BCT)*((-1 + BCT)*m1*ma*Pair(Momentum(s,p2,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar)) + 
     (1 + BCT)*Pair(Momentum(s,p1,W,theta,thetaStar,phiStar),Momentum(s,p2,W,theta,thetaStar,phiStar))*
      Pair(Momentum(s,q1,W,theta,thetaStar,phiStar),Momentum(s,q3,W,theta,thetaStar,phiStar))))};