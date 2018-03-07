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
void stoutSmear(gauge_field &W1, mdp_matrix_field &Tmp, mdp_nmatrix_field &S, float c1, float c2) {

   site x(W1.lattice());
   int N=W1.nc;
   int mu;
   int k;
   int expn=8;
   mdp_matrix Omega, Q, Qn, expQ;

   load_staples_alldirs(W1, S);

   for(mu=0; mu<W1.ndim; mu++) {
      forallsites(x) {
	 //	 cout << "S["<<x<<","<<mu<<"]\n" << S(x,mu) <<endl;

	 /* Update the accumulation matrix A += c1*proj(U*S) */
	 /*
	 mult_su3_na( U, &(s->staple[dir]), &tempS1 );
	 anti_hermitian_traceless_proj( &tempS1, &tempA1 );
	 scalar_mult_add_ah( Acur, &tempA1, c1, Acur );
	 */
	 /* Update the links U = exp(c2*A)*U */
	 /*
	 scalar_mult_ah( Acur, c2, &tempA1 );
	 exp_anti_hermitian( &tempA1, &tempS1, exp_order );
	 mult_su3_nn( &tempS1, U, &tempS2 );
	 su3mat_copy( &tempS2, U );
	 */

	 Omega = W1(x,mu) * S(x,mu); // MILC staple is Herm.Conj. of ours.
	 //	 Q = I*(hermitian(Omega) - Omega)/2 - I*trace(hermitian(Omega) - Omega)/2/N;
	 Q = (hermitian(Omega) - Omega)/2 - trace(hermitian(Omega) - Omega)/2/N;
	 Tmp(x) += c1 * Q;

	 // build expQ
	 expQ = mdp_identity(N);
	 Qn = mdp_identity(N);
	 for(k=1; k<=expn; k++) {
	    Qn = I*c2*Q*Qn/k;
	    expQ += Qn;
	 }
	 //	 expQ = exp(I*c2*Q);
	 W1(x,mu) = expQ*W1(x,mu); 
      }
   }
   //   W1.update();
}

	    
// Wilson flows W1 with timestep t  
void wilsonFlow_RK(gauge_field &W1, mdp_matrix_field &Tmp,  mdp_nmatrix_field &S, float t) {

   site x(W1.lattice());
   mdp_matrix Zeero=mdp_zero(W1.nc);

   forallsites(x) { Tmp(x) = Zeero; }
   stoutSmear(W1, Tmp, S, 17./36*t, -9./17);
   stoutSmear(W1, Tmp, S, -8./9*t, 1);
   stoutSmear(W1, Tmp, S,  3./4*t, -1);

}
