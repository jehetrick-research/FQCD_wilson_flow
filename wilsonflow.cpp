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
      S.update();
   }
   return 0;
}



// Wilson flows W1 -> W2, with timestep t, using exp() Taylor series of order n  
void wilsonFlow(gauge_field &W1,gauge_field &W2, mdp_nmatrix_field &S, double t, int n) {

   site x(W1.lattice());
   int N=W1.nc;
   int mu;
   int parity;
   int k;
   mdp_matrix Omega, Q, Qn, expQ;

   load_staples_alldirs(W1, S);
   for(mu=0; mu<W1.ndim; mu++) {
      forallsites(x) {
	 Omega = W1(x,mu) * hermitian(S(x,mu));
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
   W2.update();
}

	    
