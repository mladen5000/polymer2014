
#include <stdio.h>


void gaussj(MatDoub_IO &a, MatDoub_IO &b)
	//Input matrix is a
	//b is input containing the m rhs vectors
	//on output, a is replaced by the inverse
	//b is replaced by the solutions set
{
	Int i,icol,irow,j,k,l,ll,n=a.nrows(),m=b.ncols();
	Doub big,dum,pivinv;
	VecInt indxc(n), indxr(n), ipiv(n); //integer arrays used for bookkeepingduring pivoting

	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) { //The main loop
		big = 0.0;
		for (j=0;j<n;j++) { 
			//Begin search for pivot element
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (abs(a[j][k]) >= big) {
							big = abs(a[j][k]);
							irow = j;
							icol = k;
							}
					}
			}
		++(ipiv[icol]);
		/* By this point you have found a pivot element */

		if ( irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
			for (l=0;l<m;l++) SWAP(a[irow][l],a[icol][l]);
		}

		indxr[i] = irow;
		indxc[i] = icol;
		if ( a[icol][icol] == 0.0) throw("gaussj: singular matrix");
		pivinv = 1.0/a[icol][icol];
		a[icol][icol] = 1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
		}

		for (l=n-1;l>=0;l--){
			if (indxr[l] != indxc[l])
				for (k=0;k<n;k++)
					SWAP(a[k][indexr[l]], a[k][indexc[l]]);
		}
	}

		
			

}


int jacobi(double** a, int n, double d[], double** v)
{
	int nrot=0;
	double tresh,theta,tau,t,sm,s,h,g,c,b[n],z[n];

	//Initialize v to an identity matrix.
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
			v[i][j] = 0.0;
		v[i][i] = 1.0;
	}



