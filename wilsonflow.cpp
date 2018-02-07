/////////////////////////
// Routines to do Wilson Flow on a gauge field.
// Takes a field W1 and flows it to W2, with timestep t
// using W2 = (1 + itQ + 1/2*(itQ)^2 + 1/6*(itQ)^3 +...) W1


// external storage required:
//
//     mdp_matrix_field S    :to store staples




int load_staples(gauge_field &W1, mdp_matrix_field &S, int dir) {

   site x(W1.lattice());
   int mu;
   mdp_matrix staple(W1.nc, W1.nc);

   forallsites(x) {
      staple = 0;
      for(mu=0; mu<W1.ndim; mu++) {
	 if(mu != dir) {
           staple = staple +
	      // upper staple contribution
	      W1(x+dir, mu) * hermitian(W1(x+mu, dir)) * hermitian(W1(x, mu)) +
	      // lower staple contribution
              hermitian(W1(x+dir-mu, mu)) * hermitian(W1(x-mu, dir)) * W1(x-mu, mu);
	 }
      }
      S(x) = staple;
   }
   S.update();

   return 0;
}



// Wilson flows W1 -> W2, with timestep t, using exp() Taylor series of order n  
void wilsonFlow(gauge_field &W1,gauge_field &W2, mdp_matrix_field &S, double t, int n) {

   site x(W1.lattice());
   int N=W1.nc;
   int mu;
   int parity;
   int k;
   mdp_matrix Omega, Q, Qn, expQ;

   for(parity=0; parity<2; parity++) {
      for(mu=0; mu<W1.ndim; mu++) {
	 load_staples(W1, S, mu);
	 forallsitesofparity(x, parity) {
	    Omega = S(x) * hermitian(W1(x,mu));
	    Q = I*(hermitian(Omega) - Omega)/2 - I*(hermitian(Omega) - Omega)/2/N;
	    // build expQ
	    expQ = mdp_identity(N);
	    Qn = mdp_identity(N);
	    for(k=1; k<=n; k++) {
	       Qn = I*t*Q*Qn/k;
	       expQ += Qn;
	    }
	    W2(x,mu) = expQ*W1(x,mu); 
	 }
      }
   }
}

	    
