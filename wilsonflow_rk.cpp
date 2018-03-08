/////////////////////////
// Routines to do Wilson Flow on a gauge field.
// Takes a field W1 and flows it to W2, with timestep t
// using W2 = (1 + itQ + 1/2*(itQ)^2 + 1/6*(itQ)^3 +...) W1


// external storage required:
//
//     mdp_matrix_field S[4]    :to store staples
//
// 7 Feb 2018: added need for S[4] to do all staple dirs at once. 



int load_staples_alldirs(gauge_field &W, mdp_nmatrix_field &S) {

   site x(W.lattice());
   int mu,dir;
   mdp_matrix staple(W.nc, W.nc);

   for(dir=0; dir<W.ndim; dir++) {
      forallsites(x) {
	 staple = 0;
	 for(mu=0; mu<W.ndim; mu++) {
	    if(mu != dir) {
	       staple = staple +
		  // upper staple contribution
		  W(x+dir, mu) * hermitian(W(x+mu, dir)) * hermitian(W(x, mu)) +
		  // lower staple contribution
		  hermitian(W(x+dir-mu, mu)) * hermitian(W(x-mu, dir)) * W(x-mu, mu);
	    }
	 }
	 S(x,dir) = staple;
      }
      //      S.update();
   }
   return 0;
}



// Stout Smearing, a la MILC
// Tmp is an "accumulator field for the AH projection
//
void stoutSmear(gauge_field &W, mdp_nmatrix_field &Tmp, mdp_nmatrix_field &S, float c1, float c2) {

   site x(W.lattice());
   int N=W.nc;
   int mu;
   int k;
   int expn=8;
   mdp_matrix Omega, Q, Qn, expQ;

   load_staples_alldirs(W, S);

   for(mu=0; mu<W.ndim; mu++) {
      forallsites(x) {
	 //	 cout << "S["<<x<<","<<mu<<"]\n" << S(x,mu) <<endl;

	 Omega = W(x,mu) * S(x,mu); // MILC staple is Herm.Conj. of ours.

	 //	 cout<<"Omega: "<< x <<":"<< mu <<endl;
	 //	 cout<< Omega <<endl;


	 //	 Q = I*(hermitian(Omega) - Omega)/2 - I*trace(hermitian(Omega) - Omega)/2/N;

	 Q = (Omega - hermitian(Omega))/2 - trace(Omega - hermitian(Omega))/2/N;
	 //	 cout<<"Q: "<< x <<":"<< mu <<endl;
	 //	 cout<< Q <<endl;

	 // Tmp(x,mu) = MILC Acum field
	 Tmp(x,mu) += c1*Q;
	 //	 cout<<"Tmp: "<< x <<":"<< mu <<endl;
	 //	 cout<< Tmp(x,mu) <<endl;

	 // build exp(c2*Tmp)
	 /**/
	 expQ = mdp_identity(N);
	 Qn = mdp_identity(N);
	 for(k=1; k<=expn; k++) {
	    Qn = c2*Tmp(x,mu)*Qn/k;
	    expQ += Qn;
	 }
	 /**/
	 //	 expQ = exp(c2*Tmp(x,mu));
	 //	 cout<<"expQ: "<< x <<":"<< mu <<endl;
	 //	 cout<< expQ <<endl;
	 
	 W(x,mu) = expQ*W(x,mu); 
      }
   }
   //   W.update();
}

	    
// Wilson flows W with timestep t  
void wilsonFlow_RK(gauge_field &W, mdp_nmatrix_field &Tmp,  mdp_nmatrix_field &S, float t) {

   site x(W.lattice());
   mdp_matrix Zero=mdp_zero(W.nc);
   int mu;

   for(mu=0;mu<W.ndim;mu++) forallsites(x) { Tmp(x,mu) = Zero; }
   stoutSmear(W, Tmp, S, 17./36*t, -9./17);
   stoutSmear(W, Tmp, S, -8./9*t, 1);
   stoutSmear(W, Tmp, S,  3./4*t, -1);

}
