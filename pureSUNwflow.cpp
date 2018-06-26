#include "fermiqcd.h"              

#include "ploop3.cpp"
#include "wilsonflow_rk.cpp"
#include "readmilcascii.cpp"
#include "field_strength.cpp"

// for KH's lattice output format
#include <ctime>


int main(int argc, char** argv) {

  int verbose=false;
  int ndim=4,nc=3;
  int L[10]={4,4,4,4,4,4,4,4,4,4};
  mdp_real beta=6.0;
  int steps=0;
  char input[1024]="";
  char output[1024]="";
  int mode=0;
  int warms=0;
  int trajecs=1;
  int meas=1;
  register int i, j;
  double plaq;
  double staple_plaq=0;
  int mu, parity;
  long seed=-666;
  mdp_complex aveP;
  // saveformat determines the suffix of the saved
  // lattices.
  // saveformat=0: JH format-- .n (trajectory number)
  // saveformat=1: KH format-- .YYYYMMDDHHMMSS (really??)
  int savestep=-666;
  time_t rawtime;
  struct tm *timeinfo;
  char timestamp[80];
  char latfile[1024];
  char savestyle[3];
  int saveformat=0; // JH style
  float wf_t=0.05, wf_T=0, wf_Tmax=0, wf_savestep=0;
  int wf_expn=6;
  int UPDATEOFF=0;
  double Espace, Etime, topoQ;


  // //////////////////////////////
  // Parsing command line arguments
  // //////////////////////////////

  for(int i=1; i<argc; i++) {
     //     cout << i <<" "<< argv[i] << endl;
    if(strncmp(argv[i],"-verbose",8)==0)      verbose=true;
    else if(strncmp(argv[i],"-beta",5)==0)    sscanf(argv[i+1],"%f",&beta);    
    else if(strncmp(argv[i],"-cold",5)==0)    mode=0;
    else if(strncmp(argv[i],"-hot",4)==0)     mode=1;
    else if(strncmp(argv[i],"-savestyle",10)==0) sscanf(argv[i+1],"%s",&savestyle);  
    else if(strncmp(argv[i],"-input",6)==0)   {mode=2; sscanf(argv[i+1],"%s",input); }
    else if(strncmp(argv[i],"-readmilcascii",14)==0) {mode=3; sscanf(argv[i+1],"%s",input); }
    else if(strncmp(argv[i],"-seed",5)==0)    sscanf(argv[i+1],"%ld",&seed);
    else if(strncmp(argv[i],"-output",7)==0)  sscanf(argv[i+1],"%s",output);
    else if(strncmp(argv[i],"-warms",6)==0)   sscanf(argv[i+1],"%i",&warms);  
    else if(strncmp(argv[i],"-trajecs",8)==0) sscanf(argv[i+1],"%i",&trajecs);  
    else if(strncmp(argv[i],"-updateOFF",10)==0) UPDATEOFF=1;
    else if(strncmp(argv[i],"-meas",5)==0)    sscanf(argv[i+1],"%i",&meas);  
    else if(strncmp(argv[i],"-nc",3)==0)      sscanf(argv[i+1],"%i",&nc);
    else if(strncmp(argv[i],"-steps",6)==0)   sscanf(argv[i+1],"%i",&steps);
    else if(strncmp(argv[i],"-L",2)==0)       ndim=sscanf(argv[i+1],"%ix%ix%ix%ix%ix%ix%ix%ix%ix%i",
							  L,L+1,L+2,L+3,L+4,L+5,L+6,L+7,L+8,L+9);
    else if(strncmp(argv[i],"-savestep",9)==0)   sscanf(argv[i+1],"%i",&savestep);
    else if(strncmp(argv[i],"-wf_t",5)==0)      sscanf(argv[i+1],"%f",&wf_t);
    else if(strncmp(argv[i],"-wf_Tmax",8)==0)      sscanf(argv[i+1],"%f",&wf_Tmax);
    else if(strncmp(argv[i],"-wf_savestep",12)==0)   sscanf(argv[i+1],"%f",&wf_savestep);
    // help
    else if(strncmp(argv[i],"--help",5)==0) {
       mdp << "\nusage: [defaults]\n\n";
       mdp << "\n===Simulation Parameters:\n";
       mdp << "(works for any nc >=2, and up to 10 dimensions)\n";
       mdp << "\t-L NxNxNxN [4x4x4x4]\t \"lattice dimensions\"\n";
       mdp << "\t-beta [6.0]\t\t \"beta coupling\"\n";
       mdp << "\t-warms [0]\t\t \"number of warm-up (equilibration) sweeps\n";
       mdp << "\t-trajecs [1]\t\t \"number of \'production\' sweeps\"\n";
       mdp << "\t-updateOFF [0]\t\t \"No updates in trajectory loop\"\n";
       mdp << "\t-meas [1]\t\t \"trajecs between measurements\"\n";
       mdp << "\t-nc [3]\t\t\t \"number of colors: N of SU(N)\"\n";
       mdp << "\n===Initialization:\n";
       mdp << "\t-hot / -cold [cold]\t \"initial links\"\n";
       mdp << "\t-seed [time]\t\t \"random number seed\"\n";
       mdp << "\t-input fqcdfilename\t \"read starting FQCD file\"\n";
       mdp << "\t-readmilcascii milcfilename\t \"read starting MILC ascii file\"\n";
       mdp << "\n===Measurements and Saving gauge field U(x,mu):\n";
       mdp << "currently:   MEAS beta plaq P.re P.im\n";
       mdp << "\t-output filename_stub [fqcdlattice]\t \"lattice save filename stub (prefix)\"\n";
       mdp << "\t-savestep [none]\t \"trajecs between lattice saves\"\n";
       mdp << "\t    to save just the last lattice, set -savestep N = -trajecs N\n";
       mdp << "\t-savestyle [JH]\t  \"JH or KH suffix to saved lattice files\"\n";
       mdp << "\n===Wilson Flow Parameters, Measurements, and Saving flowed field W(x,mu):\n";
       mdp << "At meas steps in the Markov chain generating U, the Wilson Flow module\n";
       mdp << "begins, if wf_Tmax > 0. Wilson flow\'ed lattices are saved at multiples\n";
       mdp << "of wf_savestep (i.e. T = wf_savestep, 2*wf_savestep, 3*..., to files named\n";
       mdp << "wf+(Tnn)+(tnn)+filename_stub+... according to the above output format.\n";
       mdp << "\t-wf_t [0.05] \"time step for Wilson Flow: exp(itQ)*W\"\n";
       mdp << "\t-wf_Tmax [0] \"maximum accumulated Wilson Flow time\"\n";
       mdp << "\t-wf_expn [6] \"number of terms in Taylor expansion of exp(itQ)\"\n";
       mdp << "\t-wf_savestep [none=0]\t \"Wilson flow'ed lattices saved every wf_savestep in T\"\n";
       exit(1);
    }
    //    else error("Error in command line options");
  }


  mdp.open_wormholes(argc,argv);    // open communications

  
  // //////////////////////////////
  // Output parameters
  // //////////////////////////////



  // Lattice Fields and Parameters /////////////////////////
  
  if(!verbose)  mdp.print=false;  // eventualy print off
  mdp_lattice   lattice(ndim,L); // declare lattice
  mdp_site      x(lattice);      // declare site variable
  gauge_field	U(lattice,nc);   // declare SU(3) field
  coefficients	gauge;		 // declare coefficients

  // tmp storage matrix fields for Polyakov loops, etc.
  mdp_matrix_field  T1(lattice,U.nc,U.nc);
  mdp_matrix_field  T2(lattice,U.nc,U.nc);
  mdp_matrix_field  T3(lattice,U.nc,U.nc);

  // for the Staples (all dirs)
  mdp_nmatrix_field  S(lattice,U.ndim,U.nc,U.nc);
  mdp_nmatrix_field  AccumRK(lattice,U.ndim,U.nc,U.nc);

  // for the Wilson Flow
  gauge_field	W(lattice,nc);


  gauge["beta"]=beta;            // set beta



  // Initialize Lattice ////////////////////////////////////////

  if(mode==2) {
    mdp_field_file_header header;
    if(is_file(input)) header=get_info(input);
    else error("Unable to access input gauge configuration\n");
    ndim=header.ndim;
    nc=(int) sqrt((double) header.bytes_per_site/(ndim*sizeof(mdp_complex)));
    for(i=0; i<ndim; i++) L[i]=header.box[i];
  }

  switch(mode) {  
    case 0:set_cold(U); break;
    case 1:set_hot(U); break;       
    case 2:U.load(input); break;
    case 3: read_ascii_lat(input, U); break;
  }


  // Get random seed
  if(seed == -666) {
     // default: seed related to time started
     seed = (long)(time(NULL)-1495056447);
     lattice.initialize_random(seed);
     if(ME==0) { cout << "Initialized lattice random number generator with seed = "
		      << seed << endl; }
  } else {
     lattice.initialize_random(seed);
     if(ME==0) { cout << "Initialized lattice random number generator with seed = "
		      << seed << endl; }
  }

  // lattice save fromat
  if(strncmp(savestyle,"JH",2)==0) { saveformat = 0; } 
  else if(strncmp(savestyle,"KH",2)==0) { saveformat = 1; }
	


  //// report parameters as we know them
  if(ME==0) {
     cout <<"Parameters"<<endl;
     cout <<"----------------------------"<<endl;
     cout <<"L "; for(i=0; i<ndim; i++) cout << L[i] <<" ";
     cout << endl;
     cout <<"nc "<< nc <<endl;
     cout <<"beta "<< beta <<endl;
     cout <<"warms "<< warms <<endl;
     cout <<"trajecs "<< trajecs <<endl;
     if(UPDATEOFF==1) cout << "Updates are OFF" << endl;
     cout <<"meas "<< meas <<endl;
     cout <<"seed "<< seed <<endl;
     cout <<"input/output "<< input <<"/"<< output <<endl;
     if(savestep>0) {
	cout <<"saving lattices every "<< savestep <<" trajectories"<<endl;
     } else {
	cout <<"not saving intermediate lattices"<<endl;
     }
     cout <<"startmode "<< mode <<endl;
     switch(mode) {
       case 0: cout << "Creating a cold configuration\n"; break;
       case 1: cout << "Creating a hot configuration\n"; break;
       case 2: cout << "Reading MDP input file = " << input << '\n'; break;
       case 3: cout << "Reading MILC ascii input file = " << input << '\n'; break;
     }
     cout << "Performing HEATBATH updates" << endl;
     cout << "Wilson Flow:" << endl;
     if(wf_Tmax == 0) {cout << "wf_Tmax = 0" << endl;}
     else{
	cout << "time step: wf_t = "<< wf_t << endl;
	cout << "max flow time: wf_Tmax = "<< wf_Tmax << endl;
	cout << "Taylor expansion order: wf_expn = "<< wf_expn << endl;
	cout << "WF savestep: "<< wf_savestep << endl;
     }
     cout <<"----------------------------"<<endl;
  }

  /////////////////////////////////////////////////////////////
  // Start Simulation
  /////////////////////////////////////////////////////////////

  // Inital Measurement ///////////////////////////////////
  plaq = average_plaquette(U);
  aveP = ave_polyakov_loop(U, U.ndim-1,T1,T2,T3);
  if(ME==0) {
     cout <<"INIT "<< beta <<" "<< plaq 
	  <<" "<< aveP.real() <<" "<< aveP.imag() << endl;
  }




  // Thermalize ///////////////////////////////////////////
  if((UPDATEOFF==1)&&(warms>0)) {
     if(ME==0) {
	cout << "UPDATEOFF=1 (updates are OFF) -and- warms>0" << endl;
     }
     exit(1);
  }
  for(i=0; i<warms; i++) {
     WilsonGaugeAction::heatbath(U,gauge,1); //heatbath 
  }


  plaq = average_plaquette(U);
  aveP = ave_polyakov_loop(U, U.ndim-1, T1, T2, T3);
  if(ME==0) {
     cout << warms << " Heatbath steps completed" << "\n";
     cout <<"MEAS "<< beta <<" "<< plaq 
	  <<" "<< aveP.real() <<" "<< aveP.imag() << endl;
  }



  // Do the measurement trajectories ////////////////////
  ///////////////////////////////////////////////////////
  for(i=1; i<=trajecs; i++) {
     if(UPDATEOFF!=1) WilsonGaugeAction::heatbath(U,gauge,1); //heatbath 
     //     if(ME==0) cout << "update" << endl;
     
     // Measurements ////////////////////////////////////
     if(i%meas==0) {
	plaq = average_plaquette(U);
	aveP = ave_polyakov_loop(U, U.ndim-1,T1,T2,T3);
	if(ME==0) {
	   cout <<"MEAS "<< beta <<" "<< plaq 
		<<" "<< aveP.real() <<" "<< aveP.imag() << endl;
	}
	// Save U(x,mu) lattice every 'savestep' trajectories
	if((savestep>0)&&(i%savestep==0)) {
	   if(strlen(output)==0) { sprintf(output,"fqcdlattice"); }
	   if(saveformat==0) {
	      sprintf(latfile,"%s.%d",output,i);
	   } else if(saveformat==1) {
	      time(&rawtime);
	      timeinfo = localtime(&rawtime);
	      strftime (timestamp,80,".%G%m%d%H%M%S",timeinfo);
	      strcpy(latfile,output);
	      strcat(latfile,timestamp);
	   }
	   if(ME==0) cout <<"Saving lattice at trajectory "<< i <<" in file: "
			  << latfile <<endl;
	   U.save(latfile);      // save file
	}

	

	// Wilson Flow module //////////////////////////
	if(wf_Tmax > 0) {
	   if(ME==0) {
	      cout<<"Starting Wilson Flow Module"<<endl;
	      cout<<"WFFMT beta wf_T plaq Es Et Q T^2*(Es+Et)"
		  << endl;
	   }
	   forallsites(x) {
	      for(mu=0; mu<W.ndim; mu++) W(x,mu) = U(x,mu);
	   }
	}
	while((wf_Tmax > 0) && (wf_T <= wf_Tmax)) {
	   W.update();
	   wilsonFlow_RK(W, AccumRK, S, wf_t);
	   wf_T += wf_t;

	   // WFLOW Measurements /////////////////
	   plaq = average_plaquette(W);
	   em_field_strength(W, &Espace, &Etime, &topoQ);


	   if(ME==0) { cout <<"WFLOW "<< beta <<" "<< wf_T <<" "<< plaq 
			    <<" "<< Espace <<" "<< Etime <<" "<< topoQ
			    <<" "<< wf_T*wf_T*32*(1-plaq) <<endl; }
	   
	   // Wilson Flow SAVES //////////////////
	   // Save W(x,mm) lattice every 'wf_savestep' in wf_T
	   if((wf_savestep>0) && (fabs(fmod(wf_T, wf_savestep))<1e-6)) {
	      if(strlen(output)==0) { sprintf(output,"fqcdlattice"); }
	      if(saveformat==0) {
		 sprintf(latfile,"wfT%.2ft%1.2f%s.%d",wf_T, wf_t, output,i);
	      } else if(saveformat==1) {
		 time(&rawtime);
		 timeinfo = localtime(&rawtime);
		 strftime (timestamp,80,".%G%m%d%H%M%S",timeinfo);
		 sprintf(latfile,"wfT%.2ft%1.2f%s",wf_T, wf_t, output);
		 strcat(latfile,timestamp);
	      }
	      if(ME==0) cout <<"Saving WF lattice at time T= "<< wf_T <<" in file: "
			     << latfile <<endl;
	      W.save(latfile);      // save file
	   }
	} 
	wf_T = 0; // end Wilson Flow Module 

     } // end MEASurement
  } // end trajecs


  // Output lattice //////////////////////////////////////
  if(string(output)!="") U.save(output);      // save file



  // Shutdown gracefully /////////////////////////////////
  mdp.close_wormholes();
  return 0;
}
