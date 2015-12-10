// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2012 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk. 
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia.h"
#include "LHAFortran.h"

using namespace Pythia8; 

extern "C" {
  extern struct {
    double vetoscaletp;
    double vetoscaletm;
    double vetoscalewp;
    double vetoscalewm;
    int canveto; // 1: can veto, #1 cannot
  } resonancevetos_;
  extern struct {
    int py8tune;
    int nohad;
  } cpy8tune_;
}

class MyUserHooks : public UserHooks {

public:

  MyUserHooks() {cout << "Setting up Hook";}

  // Destructor deletes anti-kT jet finder.
  ~MyUserHooks() {;}

  // Allow process cross section to be modified..

  virtual bool canSetResonanceScale() {
    cout<< "called cansetresonancescale" <<endl;
    cout<< "can veto resonances? 1=true, 0=false: "<<resonancevetos_.canveto<<endl;
    if(resonancevetos_.canveto == 1) return true;
    else return false; }

  virtual double scaleResonance( const int iRes, const Event& event) {
//     cout << "in scaleResonance: resonance "; cout << event[iRes].id();
//     cout << " " <<resonancevetos_.vetoscaletp << " " <<resonancevetos_.vetoscaletm;
//     cout << "\n";
    switch(event[iRes].id())
      {
      case 6  : return  resonancevetos_.vetoscaletp;
      case -6 : return  resonancevetos_.vetoscaletm;
      case 24 : return  resonancevetos_.vetoscalewp;
      case -24: return  resonancevetos_.vetoscalewm;
      }
  }    

private:


};




Pythia pythia;
LHAupFortran LHAinstance;
MyUserHooks* MyHook=new MyUserHooks();



extern "C" {
  // F77 interface to pythia8
  void pythia_option0_(char *string) {
    pythia.readString(string);
  }

  void pythia_init_() {
  // Generator. Process selection. LHC initialization. Histogram.
    //  pythia.readString("Beams:eCM = 8000.");    
    //  pythia.readString("HardQCD:all = on");    
    //  pythia.readString("PhaseSpace:pTHatMin = 20.");

    // shower vetoing
    pythia.readString("SpaceShower:pTmaxMatch = 1");
    pythia.readString("TimeShower:pTmaxMatch = 1");

    pythia.readString("SpaceShower:QEDshowerByQ = off"); // From quarks.        
    pythia.readString("SpaceShower:QEDshowerByL = off"); // From Leptons.       
    pythia.readString("TimeShower:QEDshowerByQ = off"); // From quarks.         
    pythia.readString("TimeShower:QEDshowerByL = off"); // From Leptons.

    // tune
    if(cpy8tune_.py8tune == 14) {
      pythia.readString("Tune:pp = 14"); // Monash2013 tune
      cout << "pythia8F77: setting pythia tune 14";
      //      pythia.readString("Tune:pp = 5"); // Default
    }

    pythia.readString("PartonLevel:MPI = off");
    if(cpy8tune_.nohad == 1) {
      pythia.readString("HadronLevel:All = off");
    }

    pythia.readString("111:mayDecay = off"); 
    pythia.readString("521:mayDecay = off"); 
    pythia.readString("-521:mayDecay = off"); 
    pythia.readString("511:mayDecay = off"); 
    pythia.readString("-511:mayDecay = off"); 
    pythia.readString("531:mayDecay = off"); 
    pythia.readString("-531:mayDecay = off"); 
    pythia.readString("5222:mayDecay = off"); 
    pythia.readString("5112:mayDecay = off"); 
    pythia.readString("5232:mayDecay = off"); 
    pythia.readString("-5132:mayDecay = off"); 
    pythia.readString("5132:mayDecay = off"); 
    pythia.readString("541:mayDecay = off"); 
    pythia.readString("-541:mayDecay = off"); 
    pythia.readString("553:mayDecay = off"); 
    pythia.readString("-5112:mayDecay = off"); 
    pythia.readString("-5222:mayDecay = off"); 
    pythia.readString("-5122:mayDecay = off"); 
    pythia.readString("5332:mayDecay = off"); 
    pythia.readString("-5232:mayDecay = off"); 
    pythia.readString("-5332:mayDecay = off"); 
    pythia.readString("5122:mayDecay = off"); 


    pythia.readString("HadronLevel:All = off");

    
    pythia.setUserHooksPtr(MyHook); 
    pythia.init(&LHAinstance);
  }

  void pythia_next_(int & iret){
  // Begin event loop. Generate event. Skip if error. List first one.
    iret = pythia.next();
  }

  void pythia_to_hepevt_(const int &nmxhep, int & nhep, int * isthep,
			 int * idhep,
			 int  (*jmohep)[2], int (*jdahep)[2],
			 double (*phep)[5], double (*vhep)[4]) {
    nhep = pythia.event.size();
    if(nhep>nmxhep) {cout << "too many particles!" ; exit(-1); }
    for (int i = 0; i < pythia.event.size(); ++i) {
      *(isthep+i) = pythia.event.statusHepMC(i);
      *(idhep+i) = pythia.event[i].id();
      // All pointers should be increased by 1, since here we follow
      // the c/c++ convention of indeces from 0 to n-1, but i fortran
      // they are from 1 to n.
      (*(jmohep+i))[0] = pythia.event[i].mother1() + 1;
      (*(jmohep+i))[1] = pythia.event[i].mother2() + 1;
      (*(jdahep+i))[0] = pythia.event[i].daughter1() + 1;
      (*(jdahep+i))[1] = pythia.event[i].daughter2() + 1;
      (*(phep+i))[0] = pythia.event[i].px();
      (*(phep+i))[1] = pythia.event[i].py();
      (*(phep+i))[2] = pythia.event[i].pz();
      (*(phep+i))[3] = pythia.event[i].e();
      (*(phep+i))[4] = pythia.event[i].m();
    }
    // override mother of very first event, set to 0
    *(jmohep)[0] = 0 ;
    *(jmohep)[1] = 0 ;
  }

}


// User-written routine that does the event generation and fills hepeup.
// We do nothing here; we fill the common block on the fortran side
// before we start pythia and before we call each event.

bool LHAupFortran::fillHepRup() {
  return true;
}

bool  LHAupFortran::fillHepEup(){
  return true;
}

