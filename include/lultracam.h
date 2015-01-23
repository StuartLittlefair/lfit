/*
 *  lultracam.h
 *  lfit_xcode
 *
 *  Created by Stuart Littlefair on 10/02/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <cmath>
#include <map>
#include <string>
#include <vector>
#include "trm/subs.h"
#include "trm/position.h"
#include "trm/telescope.h"

namespace LFIT{
    
    //holds information needed to properly process each individual logfile
    struct FileInfo {
        std::string fname; // filename
        std::string redfilt;     // which redfilt. 
        int targap,compap; // aperture number for target & comparison
        double comprMag,compgMag,compuMag; // magnitude of comparison star
		virtual ~FileInfo(){
			//std::cout << "deleting FileInfo" << std::endl;
		}
    };
    // in and out operators to allow FileInfo structures to be put in
    // a Buffer object
    std::istream& operator>>(std::istream& s, const FileInfo& el);
    std::ostream& operator<<(std::ostream& os, const FileInfo& el);
    
    // holds the information from the input file
    struct Input {
		std::string format;
        Subs::Telescope site;
        Subs::Position objPos;
        Subs::Buffer1D<FileInfo> files;
        double xmin, xmax; // phase range of interest
        double hjd0, per; // ephemeris of star
		virtual ~Input(){
			//std::cout << "deleting input" << std::endl;
		}
    };    

    //! Holds all the data for a single point of a light curve
    class Datum {
        
	public:
		
		Datum() : time(0.0), expose(0.0), flux(0.0), ferr(0.0), bad(false){}
		virtual ~Datum(){}
        //! The time (originally in MJDs, later converted to orbital phase)
        double time;        
        //! The exposure length in the same units as the time
        double expose;       
        //! The flux
        double flux;        
        //! The uncertainty on the flux in the same units
        double ferr;        
        //! logical mask to indicate goodness or otherwise
        bool bad;
    };
    // Dummy ASCII input operator to allow Datum's to be put in Buffer1D's
    std::istream& operator>>(std::istream& s, const Datum& p);
    // Dummy ASCII output operator to allow use of Buffer1D
    std::ostream& operator<<(std::ostream& s, const Datum& p);

    // enum object allows us to test for the appropriate CCD of a lightcurve
    // useful for proper extinction corrections etc...
    enum CCD {
        RED=1,
        GRN=2,
        BLU=3
    };
    
    // Holds a light curve
    class Data : public Subs::Buffer1D<Datum>{
        
    public:
        
        //! Default constructor
      Data(){}
      
      //! Constructor with pre-defined size
      Data(int n){
	this->resize(n);
      }
      
      virtual ~Data(){
	//std::cout << "Destroying data objects " << std::endl;
      }
      
      Data& operator=(const Data& obj){
	if (this != &obj){
	  this->clear();
	  for(size_t i=0; i<obj.size(); i++){
	    this->push_back(obj[i]);
	  }
	  this->chip = obj.chip;
	}
	return *this;
      }
              
      // returns min x value
      double get_xmin(){
	return xmin;
      }
      // returns max x value
      double get_xmax(){
	return xmax;
      }
      // returns min y value
      double get_ymin(){
	return ymin;
      }
      // returns max y value
      double get_ymax(){
	return ymax;
      }
      
      // sets min x value
      void set_xmin(const double& in){
	xmin=in;
      }
      // set max x value
      void set_xmax(const double& in){
	xmax=in;
      }
      // sets min y value
      void set_ymin(const double& in){
	ymin=in;
      }
      // sets max y value
      void set_ymax(const double& in){
	ymax=in;
      }
      void resetLimits(){
	double x1=1.0e30,x2=-1.0e30,y1=1.0e30,y2=-1.0e30;
	for (size_t i=0; i < this->size(); i++ ){
	  if(this->buff[i].time < x1) x1 = this->buff[i].time;
	  if(this->buff[i].time > x2) x2 = this->buff[i].time;
	  if(this->buff[i].flux < y1) y1 = this->buff[i].flux;
	  if(this->buff[i].flux > y2) y2 = this->buff[i].flux;
	}
	this->set_xmin(x1);
	this->set_xmax(x2);
	this->set_ymin(y1);
	this->set_ymax(y2);
      }

      CCD chip;
    private:
      double xmin,xmax,ymin,ymax;
    };
    // Dummy ASCII input operator to allow Data objectss to be put in Buffer1D's
    std::istream& operator>>(std::istream& s, const Data& p);
    // Dummy ASCII output operator to allow use of Buffer1D
    std::ostream& operator<<(std::ostream& s, const Data& p);

    // set the limit properties of a Data object
    void setlims(Data& data);
    
    // define less than for data types which compares them by time (for sorting)
    bool operator<(Datum x, Datum y);
    
    // workhorse routine and the main entrance into these functions
    // call this if you want to read in ultracam data from logfiles
    // uses a file called input.dat to determine how to process data
    void ucam_read(LFIT::Data& red, LFIT::Data& grn, LFIT::Data& blu);
    
    // deletes unwanted phases form lightcurves
    void trim(Data& data, const Input& input);
    
    // sorts by phase
    void psort(Data& data);
    
    // splits strings into two by first occurence of delim
    void split_on(const std::string& str, const std::string& delim,
                  std::string& str1, std::string& str2);

    // does the grunt work of airmass correction, dividing by std.
    void process_ucam_data(const LFIT::Input& input, const LFIT::Data& targ, 
                           const LFIT::Data& comp, const int& fileNum,
                           LFIT::Data& out);

    // reads in a lightcurve from specified ccd and aperture from logfile
    bool loadultracam(LFIT::Data& lcurve, int nccd, int naper, std::string file); 

    // reads input.dat
    void read_infile(LFIT::Input& input);

	// reads an lcurve format file
	bool read_lcurve(LFIT::Data& lcurve, std::string file);
	
	void process_ascii_data(const LFIT::Input& input, const LFIT::Data& targ, LFIT::Data& out);

}
