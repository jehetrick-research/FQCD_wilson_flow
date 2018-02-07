/////////////////////////
// Polyakov routines



// T1, T2, and T3 are temp matrix fields. Can use the staple field S(x) for example.

mdp_complex ave_polyakov_loop(gauge_field& U, int dir, 
			      mdp_matrix_field &T1, mdp_matrix_field &T2, mdp_matrix_field &T3) {

   int t, nt; 

   nt = U.lattice().size(dir);
   //   cout << "nt = "<< nt << endl;
   //   cout << "dir = "<< dir <<endl;

   mdp_site x(U.lattice());
   mdp_complex z=0;

   forallsites(x) { T1(x) = U(x,dir); }
   T1.update();

   for(t=1; t<nt; t++){
      forallsites(x) if(x(dir)==0) {
	 if(t==1){
	    T2(x) = U(x,dir) * T1(x+dir);
	 }
	 else {
	    T2(x) = T2(x) * T1(x+dir);
	 }
	 /**
	 cout <<"t = "<< t <<endl;
	 cout <<"x = "<< x <<endl;
	 cout <<"T2=U/T2 * T1(x+dir)\n";
	 cout <<"U(x):\n"; cout << U(x,dir);
	 cout <<"T1:\n"; cout << T1(x) <<endl;
	 cout <<"T1+dir:\n"; cout << T1(x+dir) <<endl;
	 cout <<"T2:\n"; cout << T2(x) <<endl;
	 **/
      }
      forallsites(x) { 
	 T3(x) = T1(x+dir);
      }
      forallsites(x) { 
	 T1(x) = T3(x);
	 //	 cout <<"x="<< x <<"T1 now = "<< endl;
	 //	 cout << T1(x);
      }

      T1.update();
      T2.update();
   }

   forallsites(x) if(x(dir)==0) z += trace(T2(x));
   mdp.add(z);
   z /= U.lattice().nvol_gl/nt;
   //   cout <<"vol/nt = " << U.lattice().nvol_gl/nt << endl;

   return z;
}



// T is a tmp SU(N) field
// Pls is where the result will be put



mdp_complex ploop_less_slice(gauge_field &U, mdp_matrix_field &Pls, int dir,
			     mdp_matrix_field &T1, mdp_matrix_field &T2, int parity) {
   mdp_site x(U.lattice());
   int opp_parity;
   int nt=U.lattice().size(dir);
   int t;
   opp_parity = parity==0 ? 1 : 0; 


   forallsites(x) T1(x) = U(x,dir);
   T1.update();

   forallsitesofparity(x,parity) {
      for(t=1; t<nt; t++) {
         if(t==1){
            Pls(x) = U(x+dir,dir);
         }
         else {
            Pls(x) = Pls(x) * T1(x+dir);
         }
      }

      forallsites(x) T2(x) = T1(x+dir);
      forallsites(x) T1(x) = T2(x);
      T1.update();
      //      T2.update();
      Pls.update();
   }


   // report ave Ploop from parity * tslice Ploops
   mdp_complex z = 0;
   forallsitesofparity(x, parity)
      z += trace(U(x,dir) * Pls(x));

   mdp.add(z);
   z /= U.lattice().nvol_gl/2;

   return z;

}


mdp_complex ploop_less_slice_t(gauge_field &U, mdp_matrix_field &Pls, int dir, int tslice, 
			     mdp_matrix_field &T1, mdp_matrix_field &T2, int parity) {
   mdp_site x(U.lattice());
   int opp_parity;
   int nt=U.lattice().size(dir);
   int t;
   opp_parity = parity==0 ? 1 : 0; 


   forallsites(x) T1(x) = U(x,dir);
   T1.update();

   for(t=1; t<nt; t++) {
      //      forallsitesofparity(x,parity) if(x(dir)==tslice) {
      forallsites(x) if(x(dir)==tslice) {
      
         if(t==1){
            Pls(x) = U(x+dir,dir);
         }
         else {
            Pls(x) = Pls(x) * T1(x+dir);
         }
	 /**
	 cout <<"MAKE:PLS ts="<<tslice<<" loop-t="<< t <<" x = "<< x <<" "<< 
	    real(Pls(x)[0,0]) <<endl;
	 cout <<"U(x)="<<real(U(x,dir)[0,0])<<" U(x+dir)="<<real(U(x+dir,dir)[0,0])
	      <<" T1(x)="<<real(T1(x)[0,0])<<" T1(x+dir)="<<real(T1(x+dir)[0,0])
	      <<" PLS(x)= "<< real(Pls(x)[0,0]) <<endl;
	 **/
      }
      //      cout <<endl;

      //      forallsitesofparity(x,opp_parity) T2(x) = T1(x+dir);
      forallsites(x) T2(x) = T1(x+dir);
      T2.update();
      //      forallsitesofparity(x,parity) T1(x) = T2(x);
      forallsites(x) T1(x) = T2(x);
      T1.update();
      Pls.update();
   }

   // report ave Ploop from parity * tslice Ploops
   mdp_complex z = 0;
   forallsitesofparity(x, parity) if(x(dir)==tslice)
      z += trace(U(x,dir) * Pls(x));

   mdp.add(z);
   z /= U.lattice().nvol_gl/2/nt;

   return z;

}

