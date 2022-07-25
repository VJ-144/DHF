///////////////////////////////////////////////////////////////////////////////////
//    Parameters.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>

//////////////////
class Parameters
{
 public:
  static std::map<std::string,std::string> string_par;
  static std::map<std::string,double> double_par;
  static std::map<std::string,int> int_par;
  static std::map<std::string,std::vector<std::string>> vec_par;

  Parameters(){};
  Parameters(int, char**);
  void ParseCommandLineArgs(int, char**);
  void PrintOptions();
  std::string s(std::string);
  double d(std::string);
  int i(std::string);
  std::vector<std::string> v(std::string);
};

std::map<std::string,std::string> Parameters::string_par = {
  {"orbitals", "rel-s2"},        // Orbitals constructing the model space. nonrel-s(num)-p(num)- ... is also possible
  {"atom","He"},        // Target atom 
  {"radial_function_type","LSpinor"}, 
  {"eeintegral_mesh_type","Legendre"}, 
  {"filename_coulomb",""}, 
  {"filename_summary",""},
};


std::map<std::string,double> Parameters::double_par = {
  {"zeta_inv",2},        // Underlying parameter for basis
  {"eeintegral_rmax",100}, 
};

std::map<std::string,int> Parameters::int_par = {
  {"NMesh",2000}
};

std::map<std::string,std::vector<std::string>> Parameters::vec_par = {
};


// The constructor
Parameters::Parameters(int argc, char** argv)
{
  ParseCommandLineArgs(argc, argv);
}

void Parameters::ParseCommandLineArgs(int argc, char** argv)
{
  std::cout << "====================  Parameters (defaults set in Parameters.hh) ===================" << std::endl;
  for (int iarg=1; iarg<argc; ++iarg)
  {
    std::string arg = argv[iarg];
    size_t pos = arg.find("=");
    std::string var = arg.substr(0,pos);
    std::string val = arg.substr(pos+1);
    if (string_par.find(var) != string_par.end() )
    {
      if (val.size() > 0)
        std::istringstream(val) >> string_par[var];
      std::cout << var << " => " << string_par[var] << std::endl;
    }
    else if (double_par.find(var) != double_par.end() )
    {
      if (val.size() > 0)
        std::istringstream(val) >> double_par[var];
      std::cout << var << " => " << double_par[var] << std::endl;
    }
    else if (int_par.find(var) != int_par.end() )
    {
      if (val.size() > 0)
        std::istringstream(val) >> int_par[var];
      std::cout << var << " => " << int_par[var] << std::endl;
    }
    else if (vec_par.find(var) != vec_par.end() )
    {
      if (val.size() > 0)
      {
        std::istringstream ss(val);
        std::string tmp;
        while( getline(ss,tmp,',')) vec_par[var].push_back(tmp);
      }
      std::cout << var << " => ";
      for (auto x : vec_par[var]) std::cout << x << ",";
      std::cout << std::endl;
    }
    else
    {
      std::cout << "Unkown parameter: " << var << " => " << val << std::endl;
    }

  }
  std::cout << "====================================================================================" << std::endl;
}

std::string Parameters::s(std::string key)
{
  return string_par[key];
}
double Parameters::d(std::string key)
{
  return double_par[key];
}
int Parameters::i(std::string key)
{
  return int_par[key];
}
std::vector<std::string> Parameters::v(std::string key)
{
  return vec_par[key];
}

void Parameters::PrintOptions()
{
  std::cout << "Input parameters and default values: " << std::endl;
  for (auto& strpar : string_par)
  {
    std::cout << "\t" << std::left << std::setw(30) << strpar.first << ":  " << strpar.second << std::endl;
  }
  for (auto& doublepar : double_par)
  {
    std::cout <<  "\t" << std::left << std::setw(30) <<doublepar.first << ":  " << doublepar.second << std::endl;
  }
  for (auto& intpar : int_par)
  {
    std::cout <<  "\t" << std::left << std::setw(30) <<intpar.first << ":  " << intpar.second << std::endl;
  }
  for (auto& vecpar : vec_par)
  {
    std::cout <<  "\t" << std::left << std::setw(30) << vecpar.first << ":  ";
    for (auto& op : vecpar.second) std::cout << op << ",";
    std::cout << std::endl;
  }

}
