

#include "droplet_diffusion.h"
#include "config_file.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>


using namespace std;


struct FunctorC0Out {
  FunctorC0Out(double a, double b, double c, double d, double x0)
    : a(a), b(b), c(c), d(d), x0(x0) {}

  double operator() (double x, double y, double t) {
    return a + b * x + c * tanh( (x - x0) / d ) ; 
  }

  double a, b, c, d, x0;
};

int main()
{

  // read configuration file
  Config params("input.txt");
 

  // set parameter values 
  double D = params.get_parameter<double>("D");
  double Dd1 = params.get_parameter<double>("Dd1");
  double l_gamma = params.get_parameter<double>("l_gamma");
  double c0_in = params.get_parameter<double>("c0_in");
  double c_out = params.get_parameter<double>("c_out");
  double Rco = params.get_parameter<double>("Rco");
  unsigned int number_of_droplets =
      params.get_parameter<unsigned int>("number_of_droplets");
  double Lx = params.get_parameter<double>("Lx");
  double Ly = params.get_parameter<double>("Ly");
  unsigned int Nx =params.get_parameter<unsigned int>("Nx");
  unsigned int Ny =params.get_parameter<unsigned int>("Ny");
  double dt = params.get_parameter<double>("dt");
  unsigned int seed =params.get_parameter<unsigned int>("seed");
  double integration_time =
      params.get_parameter<double>("integration_time");
  double save_every =
      params.get_parameter<double>("save_every");


  double c0_out_a = params.get_parameter<double>("a");
  double c0_out_b = params.get_parameter<double>("b");
  double c0_out_c = params.get_parameter<double>("c");
  double c0_out_d = params.get_parameter<double>("d");
  double c0_out_x0 = params.get_parameter<double>("x0");

  // functor for the equilibrium concentration outide of the drop
  FunctorC0Out f_c0_out(c0_out_a, c0_out_b, c0_out_c,
                        c0_out_d, c0_out_x0);


  // simulation object
  DropletDiffusion<FunctorC0Out>
        droplet_diffusion(D, Dd1, l_gamma, c0_in, c_out,
            f_c0_out, Rco, number_of_droplets, Lx, Ly, dt,
            Nx, Ny, seed);

  cout << "Concentration integration time step: "
       << droplet_diffusion.GetMaxTimeStepConcentration() << endl;
  
  cout << "Number of droplets start: " 
       << droplet_diffusion.GetNumberOfDroplets()  << endl;

  cout << "Total concentration start: " 
       << droplet_diffusion.GetTotalConcentration() << endl;

  unsigned int ti = 0;
  droplet_diffusion.SaveDroplets(
      "data/drops_" + to_string(ti) + ".dat");
  droplet_diffusion.SaveConcentration(
      "data/c_" + to_string(ti) + ".dat");

  while (droplet_diffusion.GetTime() < integration_time) {
    // integrate time
    droplet_diffusion.TimeEvolve(save_every);
    ti++;

    // save data
    droplet_diffusion.SaveDroplets(
        "data/drops_" + to_string(ti) + ".dat");
    droplet_diffusion.SaveConcentration(
        "data/c_" + to_string(ti) + ".dat");

  }

  cout << "Total concentration end: "
       << droplet_diffusion.GetTotalConcentration() << endl;
  cout << "Number of droplets end: " 
       << droplet_diffusion.GetNumberOfDroplets()  << endl;

  // save final droplets and concentration
  droplet_diffusion.SaveDroplets("data/drops.dat");
  droplet_diffusion.SaveConcentration("data/c.dat");


  return 0;
}
