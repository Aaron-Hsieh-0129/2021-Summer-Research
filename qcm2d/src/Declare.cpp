#include "Declare.h"

Array::Array() {
	for (int k = 0; k <= nz-1; k++) {
		for (int i = 0; i <= nx-1; i++){
			pip[i][k] = pi[i][k] = pim[i][k] = 0.;
			thp[i][k] = th[i][k] = thm[i][k] = 0.;
			qvp[i][k] = qv[i][k] = qvm[i][k] = 0.;
			qcp[i][k] = qc[i][k] = qcm[i][k] = 0.;
			qrp[i][k] = qr[i][k] = qrm[i][k] = 0.;
			up[i][k] = u[i][k] = um[i][k] = 0.;
			wp[i][k] = w[i][k] = wm[i][k] = 0.;
		}
		tb[k] = 0.; 
		tb_zeta[k] = 0.;
		rhou[k] = 0.;
		rhow[k] = 0.;
		pib[k] = 0.;
		qvb[k] = 0.;
		qvsb[k] = 0.;
		tvb[k] = 0.;
		pb[k] = 0.;
		ubm[k] = 0.;
	}
}

void Array::BoundaryProcess(double tmp[][nz]) {
	for (int k = 1; k <= nz-2; k++) {
		tmp[0][k] = tmp[nx-2][k];
		tmp[nx-1][k] = tmp[1][k];
	}
	for (int i = 0; i <= nx-1; i++) {
		tmp[i][0] = tmp[i][1];
		tmp[i][nz-1] = tmp[i][nz-2];
	}
}

void Array::BoundaryProcessZETA(double tmp[][nz]) {
	for (int k = 1; k <= nz-2; k++) {
		tmp[0][k] = tmp[nx-2][k];
		tmp[nx-1][k] = tmp[1][k];
	}
	for (int i = 0; i <= nx-1; i++) {
		tmp[i][0] = tmp[i][1];
		tmp[i][nz-1] = 0.;
	}
}

void Array::BoundaryProcessDouble(double tmp[][nz]) {
	for (int k = 1; k <= nz-2; k++) {
		tmp[0][k] = tmp[nx-2][k];
		tmp[nx-1][k] = tmp[1][k];
	}
	for (int i = 0; i <= nx-1; i++) {
		tmp[i][0] = tmp[i][nz-2];
		tmp[i][nz-1] = tmp[i][1];
	}
}
