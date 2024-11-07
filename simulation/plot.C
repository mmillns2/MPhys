#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TApplication.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>


//  Function to read whitespace-delimited data from a text file
std::vector<std::vector<std::string>> readDataFromTxt(const std::string& filename) {
  std::vector<std::vector<std::string>> data;
  std::ifstream file(filename);

  // Check if the file was opened successfully
  if (!file.is_open()) 
  {                         
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return data;                                     
  }
                                         
  std::string line;
  // Read each line from the file
  while (std::getline(file, line)) 
  {
    std::vector<std::string> row;
    std::istringstream line_stream(line);
    std::string word;

    // Read each word separated by whitespace
    while (line_stream >> word)
    {
      row.push_back(word);
    }

    // Add the row to the data vector if its not empty
    if (!row.empty())
    {
      data.push_back(row);
    }
  }
  file.close();
  return data;
} 


// Function to separate columns into individual vectors
std::vector<std::vector<std::string>> separateColumns(const std::vector<std::vector<std::string>>& data)
{
  std::vector<std::vector<std::string>> columns;

  if (data.empty()) 
  {
    std::cerr << "Error: Data is empty." << std::endl;
    return columns;
  }

  // Determine the number of columns (assumes all rows have the same number of columns)
  size_t num_columns = data[0].size();
  columns.resize(num_columns);

  // Populate each column vector
  for (const auto& row : data) 
  {
    for (size_t i = 0; i < row.size(); ++i)
    { 
      columns[i].push_back(row[i]);
    }
  }
  return columns;
}


// Main function
int plot()
{
  // Import data from txt file (whitespace delimiter)
  std::string filename = "data.txt";
  std::vector<std::vector<std::string>> data = readDataFromTxt(filename);

  // Separate each column into individual vectors
  std::vector<std::vector<std::string>> columns = separateColumns(data);

 //  Create a canvas
  TCanvas *c1 = new TCanvas("c1","Our Plot", 2048, 1536);
 
  c1->SetGrid();
  c1->GetFrame()->SetBorderSize(12);
 
  //  Define the data to be plotted
  const int n_points = columns[0].size();
  Float_t x[n_points]  = columns[0];
  Float_t y[n_points]  = columns[1];
  Float_t ex[n_points] = columns[2];
  Float_t ey[n_points] = columns[3];

  //  Instantiate a graph with errors
  TGraphErrors *graph = new TGraphErrors(n_points, x, y, ex, ey);

  //  Change how the plot looks
  graph->SetTitle("TGraphErrors Example");
  graph->SetMarkerColor(kBlue);
  graph->SetMarkerStyle(kOpenCircle);
 
  //  Plot the graph
  graph->DrawClone("APE");
 
  c1->Update();
  c1->Print("firstPlot.png");
       
  return 0;
}
