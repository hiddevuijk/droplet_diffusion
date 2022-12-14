#ifndef GUARD_CONCENTRATION_H
#define GUARD_CONCENTRATION_H

#include <vector>
#include <string>
#include <fstream>



class Concentration
{
 public:
  Concentration(double c_init, double D, double Lx, double Ly, unsigned int Nx, unsigned int Ny);

  // Source of dn particles in bin centered ad x,y
  void Source(double x, double y, double dn);
 
  // evolve time by delta_t
  void TimeEvolve(double delta_t); 

  // save concentration field to file named {name}
  void Save(std::string name) const;

  // return the concentration in the bin centered at x,y
  double GetConcentration(double x, double y) const;

  double GetAverageConcentration() const;

  double GetMaxTimeStep() const { return dt_max_;}
 private:

  // diffusion constant
  double D_;
  // box size
  double Lx_, Ly_;
  // number of bins in the x and y directions
  double Nx_, Ny_;

  // bin dimensions
  double dx_, dy_;

  
  std::vector<std::vector<double> > concentration1_;
  std::vector<std::vector<double> > concentration2_;

  // pointer to the current concentration  (either concentration1_ or 
  // concentratoin2_)
  std::vector<std::vector<double> > *current_concentration_ptr_;

  // pointer to concentration array that current_concentration_ptr
  // does not point to
  // used as a temporary container
  std::vector<std::vector<double> > *temp_concentration_ptr_;

  // set to 0.1 time the Von Nuemann stability
  double dt_max_;

  // make a single time step of length dt
  void MakeTimeStep(double dt);
};

Concentration::Concentration(double c_init, double D, double Lx, double Ly,
      unsigned int Nx, unsigned int Ny)
  : D_(D), Lx_(Lx), Ly_(Ly), Nx_(Nx), Ny_(Ny), dx_(Lx/Nx), dy_(Ly/Ny),
    concentration1_(Nx, std::vector<double>(Ny,c_init) ),
    concentration2_(Nx, std::vector<double>(Ny)),
    current_concentration_ptr_(&concentration1_),
    temp_concentration_ptr_(&concentration2_)
{
  if (dx_ < dy_) {
    dt_max_ = 0.1 * dx_ * dx_ / (2 * D_);  
  } else {
    dt_max_ = 0.1 * dy_ * dy_ / (2 * D_);  
  }
}


void Concentration::TimeEvolve(double delta_t)
{
  while (delta_t > dt_max_) {
    MakeTimeStep(dt_max_);
    delta_t -= dt_max_;
  }
  MakeTimeStep(delta_t);
}

void Concentration::MakeTimeStep(double dt)
{

  double Dx = dt * D_ / (dx_ * dx_);
  double Dy = dt * D_ / (dy_ * dy_);
 
  for (unsigned int ix = 0; ix < Nx_; ++ix) {
  for (unsigned int iy = 0; iy < Ny_; ++iy) {
    (*temp_concentration_ptr_)[ix][iy] = (*current_concentration_ptr_)[ix][iy];

    // x derivative
    if (ix == 0) {
      (*temp_concentration_ptr_)[ix][iy] += Dx * (
          (*current_concentration_ptr_)[ix+1][iy]
          - (*current_concentration_ptr_)[ix][iy]);
    }else if (ix == Nx_-1) {
      (*temp_concentration_ptr_)[ix][iy] += Dx * (
          (*current_concentration_ptr_)[ix-1][iy]
          - (*current_concentration_ptr_)[ix][iy]);
    } else {
      (*temp_concentration_ptr_)[ix][iy] += Dx * (
          (*current_concentration_ptr_)[ix+1][iy]
          + (*current_concentration_ptr_)[ix-1][iy]
          - 2 * (*current_concentration_ptr_)[ix][iy]);
    }

    // y derivative
    if (iy == 0) {
      (*temp_concentration_ptr_)[ix][iy] += Dy * (
          (*current_concentration_ptr_)[ix][iy+1]
          - (*current_concentration_ptr_)[ix][iy]);
    } else if (iy == Ny_ - 1) {
      (*temp_concentration_ptr_)[ix][iy] += Dy * (
          + (*current_concentration_ptr_)[ix][iy-1]
          - (*current_concentration_ptr_)[ix][iy]);
    } else {
      (*temp_concentration_ptr_)[ix][iy] += Dy * (
          (*current_concentration_ptr_)[ix][iy+1]
          + (*current_concentration_ptr_)[ix][iy-1]
          - 2 * (*current_concentration_ptr_)[ix][iy]);
    }

  }} 

  // swap pointers
  std::vector<std::vector<double> > *temp_ptr;
  temp_ptr = temp_concentration_ptr_;
  temp_concentration_ptr_ = current_concentration_ptr_;
  current_concentration_ptr_ = temp_ptr;
}

void Concentration::Source(double x, double y, double dn)
{
  // CHECK FOR UNPHYSICAL VALUES
  unsigned int xi = x / dx_; 
  unsigned int yi = y / dy_; 
  (*current_concentration_ptr_)[xi][yi] += dn / (dx_ * dy_);  
}

double Concentration::GetConcentration(double x, double y) const
{
  // CHECK FOR UNPHYSICAL VALUES
  unsigned int xi = x / dx_; 
  unsigned int yi = y / dy_; 
  return (*current_concentration_ptr_)[xi][yi];
}

void Concentration::Save(std::string name) const 
{
  std::ofstream out;
  out.open(name);
  for (unsigned int ix = 0; ix < Nx_; ++ix) {
    for (unsigned int iy = 0; iy < Ny_; ++iy) {
      out << (*current_concentration_ptr_)[ix][iy];
      if (iy < Ny_ -1) out << "\t";
    }
    out << "\n";
  }
}

double Concentration::GetAverageConcentration() const
{

  double c_average = 0.0;
  for (unsigned int ix = 0; ix < Nx_; ++ix) {
  for (unsigned int iy = 0; iy < Ny_; ++iy) {
    c_average += (*current_concentration_ptr_)[ix][iy] * dx_ * dy_;
  }}
  return c_average / (Lx_ * Ly_);
}
#endif
