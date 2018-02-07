// Utility to read in a MILC save_ascii file

int read_ascii_lat(char *fname, gauge_field &U) {
   char line[100];
   FILE *fp;
   int x,y,z,t,i,a,b,dir, Nt,Nx,Ny,Nz;
   mdp_site siteX(U.lattice());
   mdp_matrix M(U.nc,U.nc);
   float Ur, Ui;

   Nt = U.lattice().size(0);
   Nx = U.lattice().size(1);
   Ny = U.lattice().size(2);
   Nz = U.lattice().size(3);

   i = 0;
   fp = fopen(fname, "r");

   cout << "reading ASCII MILC lattice from " << fname << endl;

   // discard header
   //   cout << "discarding header lines\n";
   while (i<3) {
      if( fgets (line, 100 , fp)!=NULL ) i++; 
      //      cout << "header line " << i << endl;
   }

   //   ### Define Nt, Nx, etc.

   for(t=0; t<Nt; t++)
      for(z=0; z<Nz; z++)
	 for(y=0; y<Ny; y++)
	    for(x=0; x<Nx; x++) {
	       //      ### need to set a site for U(x)(a,b)
	       siteX.set(t,x,y,z);
	       cout <<  "set siteX: " <<t<<" "<<x<<" "<<y<<" "<<z<<endl;
	       for(dir=0; dir<4; dir++) {
		  cout << "dir: "<< dir << endl;
		  // get x= y= z= t= line
		  //		     fscanf(fp,"%s\n", line);
		  if (fgets(line, 100, fp) != NULL) {
		     cout <<"x= line: "<< line << endl;
		     //		     break;
		  }
                  //read a link
		  for(a=0; a<3; a++) for(b=0; b<3; b++) {
		     fscanf(fp,"%e %e\n", &Ur, &Ui);
		     M(a,b) = mdp_complex(Ur, Ui);
		     //		     if(i==0) {
		     //			printf("coords %d %d %d %d dir %d Ur: %f Ui: %f\n", x,y,z,t,dir,Ur, Ui);
		     //			L[t][z][y][x][dir] = Ur + I*Ui;
		  }
		  cout << M << endl;
		  U(siteX, dir) = M;
		  cout << "U: siteX=" << siteX <<"\n"<< U(siteX, dir);
	       }
   }

   fclose(fp);

}
		  
