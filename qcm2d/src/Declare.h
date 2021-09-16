#include <iostream>
#include <cmath>
#include "Const.h"

using namespace std;

class Array {
	public:
			// constructor
			Array();

			double tb[nz], tb_zeta[nz], rhou[nz], rhow[nz], pib[nz], qvb[nz], qvsb[nz], tvb[nz], pb[nz];
			double up[nx][nz], u[nx][nz], um[nx][nz];
			double wp[nx][nz], w[nx][nz], wm[nx][nz];
			double pip[nx][nz], pi[nx][nz], pim[nx][nz];
			double thp[nx][nz], th[nx][nz], thm[nx][nz];
			double qvp[nx][nz], qv[nx][nz], qvm[nx][nz];
			double qcp[nx][nz], qc[nx][nz], qcm[nx][nz];
			double qrp[nx][nz], qr[nx][nz], qrm[nx][nz];
			double addflx[nx];
			double ubm[nz];

			static void BoundaryProcess(double tmp[][nz]);
			static void BoundaryProcessZETA(double tmp[][nz]);
			static void BoundaryProcessDouble(double tmp[][nz]);
};