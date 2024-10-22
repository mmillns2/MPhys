#include"TF1.h"

void functionPlot()
{
	TCanvas* c = new TCanvas("c", "My Canvas", 800, 600);

	TF1* f1 = new TF1("f1","x",0.,10.);	
	f1->Draw();

	c->Update();
	c->Print("myplot.png");
}
