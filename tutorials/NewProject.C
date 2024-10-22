#include<iostream>
#include"TF1.h"

auto pi = TMath::Pi();

// function code
double single(double *x, double *par)
{
return pow(sin(pi*par[0]*x[0])/(pi*par[0]*x[0]), 2);
}

double nslit0(double*x, double *par)
{
return pow(sin(pi*par[1]*x[0])/sin(pi*x[0]), 2);
}

double nslit(double *x, double *par)
{
return single(x, par) * nslit0(x, par);
}

// main function
int NewProject()
{
float r, ns;
TCanvas* c = new TCanvas("c", "My Canvas",800 , 600);

// get user input
cout << "slit width /g ? ";
scanf("%f", &r);
cout << "number of slits? ";
scanf("%f", &ns);

// define function and set options
TF1 *Fnslit = new TF1("Fnslit", nslit, -5.00, 5, 2);
Fnslit->SetNpx(500);

// set parameters
Fnslit->SetParameter(0,r);
Fnslit->SetParameter(1,ns);

// draw and print
Fnslit->Draw();
c->Update();
c->Print("slitsPlot.png");


return 0;
}
