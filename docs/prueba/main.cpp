#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

struct atom{
  int spin;
  bool spin_b;
  int Sa, Sb;
}hierro ; 

struct system_info{
	int mc_steps;
	double Dq, STemperature, average;
	const int B = 23002;
}sys_info;

double metropolis(atom &, system_info &);
int montecarlo(atom &, system_info &);
int int_randomNumber_generator(const int, const int);
double real_randomNumber_generator(const int, const int);

int main(int argc, char const *argv[])
{
	//remove old files
	if(remove("table.dat") != 0) std::perror("Error deleting file: table.dat");
    else std::puts("The 'table.dat' file was successfully deleted.");

	//put a for loop running over temperature. External loop
  hierro.spin_b = 1;

  sys_info.STemperature = 1; //fixed temperature
  sys_info.average = 0;        //average of spin, (0 or 1)
  sys_info.mc_steps = 100;
  sys_info.Dq = 21000;
  
  std::ofstream outputfileflux ("table.dat", std::ios::out|std::ios::app);
  if (!outputfileflux.is_open()) {
    std::wcerr << "Error abriendo \"table.dat\" para escribir en el." << std::endl;
    abort ();
  }
 
  outputfileflux <<"\n#Temperature:"<<sys_info.STemperature
                 <<"\n#initial Dq:" <<sys_info.Dq 
                 <<"\n#initial spin_b:" << hierro.spin_b
                 <<"\nDq \t spin"   <<std::endl;
  
  for (int i = 21100; i < 24030; i+=1){
  	sys_info.Dq = i;
  	montecarlo(hierro, sys_info);
  	outputfileflux <<i<<"\t"<<sys_info.average<<std::endl;
  	sys_info.average = 0;
  }
  
  outputfileflux <<"\n"<<std::endl;
  sys_info.average = 0;
  for (int i = 24030; i > 21100; i-=1000){
  	sys_info.Dq = i;
  	montecarlo(hierro, sys_info);
  	outputfileflux <<i<<"\t"<<sys_info.average<<std::endl;
  	sys_info.average = 0;
  }

  outputfileflux.close();
  return 0;
}


int montecarlo(atom & hierro, system_info & sys_info)
{
  //the firs case will be a isotherm
  double p_omega = 0, p_calculed = 0, sumation = 0;
  bool tmp_spin;

  for (int i = 0; i < sys_info.mc_steps; ++i){

    tmp_spin = hierro.spin_b;
    double delta_energy = metropolis(hierro, sys_info);

    p_omega = exp(-(delta_energy)/sys_info.STemperature * 0.695);
    p_calculed = real_randomNumber_generator(0, 1);
    tmp_spin = hierro.spin_b;
      
    if (metropolis(hierro, sys_info) <= 0) hierro.spin_b = 1 - hierro.spin_b;

    else if (p_omega >= p_calculed) hierro.spin_b = 1 - hierro.spin_b;
      else hierro.spin_b = tmp_spin;

    sumation += hierro.spin_b;
  }
  sys_info.average = sumation/sys_info.mc_steps;

  return 0;
}

double metropolis(atom & hierro, system_info & sys_info)
{
  double Eex = 0.0, Eb = 0.0;
  
  if (sys_info.Dq < 22612){
    Eb = -4 * sys_info.Dq;
    Eex = -24 * sys_info.Dq + 45224.16;
    
    if(hierro.spin_b == 1)return -(Eex-Eb);
    else if(hierro.spin_b == 0)return -(Eb-Eex);
  }
  else if (sys_info.Dq > 22612){	  
    Eb = -24 * sys_info.Dq + 45224.16;
    Eex = -4 * sys_info.Dq;

    if(hierro.spin_b == 1)return -(Eex-Eb);
    else if(hierro.spin_b == 0)return -(Eb-Eex);
  }
  return 0;
}

int int_randomNumber_generator(const int min/*!min*/,
                               const int max/*!max*/) /*!to generate a real number between <const min:const max>*/
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis_int(min, max);
  return dis_int(gen);
}

double real_randomNumber_generator(const int min/*!min*/,
                                   const int max/*!max*/)/*! to generate a real number between <const min:const max> */
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis_real(min, max);
  return dis_real(gen);
}
