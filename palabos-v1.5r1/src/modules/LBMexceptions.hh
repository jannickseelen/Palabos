#ifndef LBMEXCEPTIONS_HH
#define LBMEXCEPTIONS_HH

#include <exception>

class MachException: public std::exception{
  virtual const char* what() const throw()
  {
    return "Compressibility Error too Large. Lower your u0lb in parameters.xml section LBM";
  }
}machEx;

class LocalMachException: public std::exception{
  virtual const char* what() const throw()
  {
    return "Local Compressibility Error too Large. Do one of the following...\n"
	"1. Decrease u0lb \n"
	"2. Decrease the max Grid refinement Level \n"
	"in parameters.xml \n";
  }
}localMachEx;

class SuperException: public std::exception{
  virtual const char* what() const throw()
  {
    return "The model does not support supersonic flows. Do one of the following...\n"
	"1. Increase u0lb \n"
	"2. Decrease the reference Resolution \n"
	"3. Decrease the max Grid refinement Level \n"
	"in parameters.xml \n";
  }
}superEx;

class MarginException: public std::exception{
  virtual const char* what() const throw()
  {
    return "Margin cannot be smaller then borderWidth, please correct the xml file";
  }
}marginEx;

class ResolutionException: public std::exception{
  virtual const char* what() const throw()
  {
    return "The resolution cannot be zero, please correct the xml file";
  }
}resolEx;
#endif // define LBMEXCEPTIONS_HH
