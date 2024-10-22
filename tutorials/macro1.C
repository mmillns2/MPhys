#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TRandom3.h"

void macro1()
{
	// The values and the errors on the Y axis
	const int n_points=10;
	double x_vals[n_points]={1,2,3,4,5,6,7,8,9,10};
	double y_vals[n_points]={6,12,14,20,22,24,35,45,44,53};
	double y_errs[n_points]= {5,5,4.7,4.5,4.2,5.1,2.9,4.1,4.8,5.43};

	
	// Instance of the graph
	TGraphErrors graph(n_points,x_vals,y_vals,nullptr,y_errs);
    	graph.SetTitle("Measurement XYZ;length [cm];Arb.Units");

	// make the plot estetically better
	graph.SetMarkerStyle(kOpenCircle);
    	graph.SetMarkerColor(kBlue);
    	graph.SetLineColor(kBlue);

	// create canvas
	auto  mycanvas = new TCanvas("c", "my canvas", 2048, 1536);

	// draw graph
	graph.DrawClone("APE");

	// define a linear function
	TF1 f("Linear law","[0]+x*[1]",.5,10.5);

	// Let's make the function line nicer
	f.SetLineColor(kRed);
	f.SetLineStyle(2);

	// Fit it to the graph and draw it
	graph.Fit(&f);
	f.DrawClone("Same");

	// Generate random numbers and input into f
	// Create an instance of TRandom3
	TRandom3 randGen;
	const int n_rand_points = 30; 
	double x_rand_vals[n_rand_points];
	double y_rand_vals[n_rand_points];
	for(int i = 0; i < n_rand_points; i++) 
	{
		double randomInRange = randGen.Uniform(0, 10); // Uniform in [0, 10]
		x_rand_vals[i] = randomInRange;
		y_rand_vals[i] = f.Eval(randomInRange);
	}

	// instance of rand_graph
	TGraphErrors rand_graph(n_rand_points,x_rand_vals,y_rand_vals,nullptr,nullptr);
	rand_graph.SetMarkerStyle(3);
	rand_graph.SetMarkerColor(kGreen + 3);
	rand_graph.SetLineColor(kGreen);
	rand_graph.SetMarkerSize(3.0);
	rand_graph.DrawClone("P SAME");

	// Build and Draw a legend
	TLegend leg(.1,.7,.3,.9,"Test program");
    	leg.SetFillColor(0);
    	graph.SetFillColor(0);
    	leg.AddEntry(&graph,"Exp. Points");
    	leg.AddEntry(&f,"Th. Law");
    	leg.AddEntry(&rand_graph,"Rand. Points");
    	leg.DrawClone("Same");

	// Draw an arrow on the canvas
	TArrow arrow(8,8,6.2,23,0.02,"|>");
    	arrow.SetLineWidth(2);
    	arrow.DrawClone();

	// Add some text to the plot
	TLatex text(8.2,7.5,"#splitline{Maximum}{Deviation}");
    	text.DrawClone();

    	mycanvas->Print("myplot2.png");
}
