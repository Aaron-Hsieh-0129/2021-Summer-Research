#include "Iteration.h"

void Iteration::pu_pt(Array& myArray) {

	// calculate ubm
	for (int k = 1; k <= nz-2; k++) {
		double ubar = 0.;
		for (int i = 0; i <= nx-1; i++) {
			ubar += myArray.um[i][k];
		}
		myArray.ubm[k] = ubar / nx;
	}
	myArray.ubm[0] = myArray.ubm[1];
	myArray.ubm[nz-1] = myArray.ubm[nz-2];

	double puzeta_px = 0., pwzeta_pz = 0., g_tb_pth_px = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			myArray.up[i][k] = myArray.um[i][k] - 0.25 * d2t * rdx * (pow(myArray.u[i+1][k] + myArray.u[i][k], 2) - pow(myArray.u[i][k] + myArray.u[i-1][k], 2))
												- 0.25 * d2t * rdz * (myArray.rhow[k+1] * (myArray.w[i][k+1] + myArray.w[i-1][k+1]) * (myArray.u[i][k+1] + myArray.u[i][k]) 
																	  - myArray.rhow[k] * (myArray.w[i][k] + myArray.w[i-1][k]) * (myArray.u[i][k] + myArray.u[i][k-1])) / myArray.rhou[k]
												- d2t * rdx * C_p * myArray.tvb[k] * (myArray.pi[i][k] - myArray.pi[i-1][k]);

			// Advection test
			// #if defined(ADVECTIONU)
			// 	g_tb_pth_px = 0.;
			// #elif defined(ADVECTIONW)
			// 	g_tb_pth_px = 0.;
			// #else
			// 	g_tb_pth_px = g / myArray.tb_zeta[k] * (0.5*(myArray.th[i][k] + myArray.th[i][k-1]) - 0.5*(myArray.th[i-1][k] + myArray.th[i-1][k-1])) * rdx;
			// #endif

			// Add diffusion
			#ifdef DIFFUSION
				myArray.up[i][k] += d2t * Kx * rdx2 * (myArray.um[i+1][k] - 2. * myArray.um[i][k] + myArray.um[i-1][k]) + 
									d2t * Kz * rdz2 * (myArray.um[i][k+1] - 2. * myArray.um[i][k] + myArray.um[i][k-1]) - 
									d2t * Kz * rdz2 * (myArray.ubm[k+1] - 2. * myArray.ubm[k] + myArray.ubm[k-1]);
			#endif
		}
	}
	#if defined(ADVECTIONU)
		myArray.BoundaryProcessDouble(myArray.up);
	#elif defined(ADVECTIONW)
		myArray.BoundaryProcessDouble(myArray.up);
	#else
		myArray.BoundaryProcess(myArray.up);
	#endif

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.u[i][k] += TIMETS * (myArray.up[i][k] - 2 * myArray.u[i][k] + myArray.um[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::pw_pt(Array& myArray) {
	double puzeta_px = 0., pwzeta_pz = 0., g_tb_pth_px = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			#if defined(NoBouyance)
				myArray.wp[i][k] = myArray.wm[i][k] - 0.25 * d2t * rdx * ((myArray.u[i+1][k] + myArray.u[i+1][k-1]) * (myArray.w[i+1][k] + myArray.w[i][k]) - (myArray.u[i][k] + myArray.u[i][k-1]) * (myArray.w[i][k] + myArray.w[i-1][k]))
										- 0.25 * d2t * rdz * (myArray.rhou[k] * pow(myArray.w[i][k+1] + myArray.w[i][k], 2) - myArray.rhou[k-1] * pow(myArray.w[i][k] + myArray.w[i][k-1], 2)) / myArray.rhow[k]
										- 0.5 * d2t * rdz * C_p * (myArray.tvb[k] + myArray.tvb[k-1]) * (myArray.pi[i][k] - myArray.pi[i][k-1]);
			#else
				myArray.wp[i][k] = myArray.wm[i][k] - 0.25 * d2t * rdx * ((myArray.u[i+1][k] + myArray.u[i+1][k-1]) * (myArray.w[i+1][k] + myArray.w[i][k]) - (myArray.u[i][k] + myArray.u[i][k-1]) * (myArray.w[i][k] + myArray.w[i-1][k]))
										- 0.25 * d2t * rdz * (myArray.rhou[k] * pow(myArray.w[i][k+1] + myArray.w[i][k], 2) - myArray.rhou[k-1] * pow(myArray.w[i][k] + myArray.w[i][k-1], 2)) / myArray.rhow[k]
										- 0.5 * d2t * rdz * C_p * (myArray.tvb[k] + myArray.tvb[k-1]) * (myArray.pi[i][k] - myArray.pi[i][k-1])
										+ d2t * g * 0.5 * (myArray.th[i][k] / myArray.tb[k] + myArray.th[i][k-1] / myArray.tb[k-1]);
			#endif
			// Advection test
			// #if defined(ADVECTIONU)
			// 	g_tb_pth_px = 0.;
			// #elif defined(ADVECTIONW)
			// 	g_tb_pth_px = 0.;
			// #else
			// 	g_tb_pth_px = g / myArray.tb_zeta[k] * (0.5*(myArray.th[i][k] + myArray.th[i][k-1]) - 0.5*(myArray.th[i-1][k] + myArray.th[i-1][k-1])) * rdx;
			// #endif

			// Add water 
			#if defined(WATER)
				myArray.wp[i][k] += d2t * (0.61 * myArray.qv[i][k] - myArray.qc[i][k] - myArray.qr[i][k]);
			#endif

			// Add diffusion
			#ifdef DIFFUSION
				myArray.wp[i][k] += d2t * Kx * rdx2 * (myArray.wm[i+1][k] - 2. * myArray.wm[i][k] + myArray.wm[i-1][k]) + 
									d2t * Kz * rdz2 * (myArray.wm[i][k+1] - 2. * myArray.wm[i][k] + myArray.wm[i][k-1]);
			#endif
		}
	}
	#if defined(ADVECTIONU)
		myArray.BoundaryProcessDouble(myArray.wp);
	#elif defined(ADVECTIONW)
		myArray.BoundaryProcessDouble(myArray.wp);
	#else
		myArray.BoundaryProcessZETA(myArray.wp);
	#endif

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.w[i][k] += TIMETS * (myArray.wp[i][k] - 2 * myArray.w[i][k] + myArray.wm[i][k]);
			}
		}
	#endif
	return;
}


void Iteration::pth_pt(Array& myArray) {
	double puth_px = 0., prhowth_pz_rho = 0., wptb_pz = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puth_px = (myArray.u[i+1][k] * 0.5*(myArray.th[i+1][k] + myArray.th[i][k]) - myArray.u[i][k] * 0.5*(myArray.th[i][k] + myArray.th[i-1][k])) * rdx;
			prhowth_pz_rho = (myArray.rhow[k+1] * myArray.w[i][k+1] * 0.5*(myArray.th[i][k+1] + myArray.th[i][k]) - 
							  myArray.rhow[k] * myArray.w[i][k] * 0.5*(myArray.th[i][k] + myArray.th[i][k-1])) * rdz / myArray.rhou[k];
			wptb_pz = 0.5*(myArray.w[i][k+1] + myArray.w[i][k]) * (0.5*(myArray.tb[k+1] + myArray.tb[k]) - 0.5*(myArray.tb[k] + myArray.tb[k-1])) * rdz;

			myArray.thp[i][k] = myArray.thm[i][k] + d2t * (-puth_px - prhowth_pz_rho - wptb_pz);

			#ifdef DIFFUSION
				myArray.thp[i][k] += d2t * Kx * rdx2 * (myArray.thm[i+1][k] - 2. * myArray.thm[i][k] + myArray.thm[i-1][k]) + 
									 d2t * Kz * rdz2 * (myArray.thm[i][k+1] - 2. * myArray.thm[i][k] + myArray.thm[i][k-1]);
			#endif
		}
	}
	#if defined(ADVECTIONU)
		myArray.BoundaryProcessDouble(myArray.thp);
	#elif defined(ADVECTIONW)
		myArray.BoundaryProcessDouble(myArray.thp);
	#else
		myArray.BoundaryProcess(myArray.thp);
	#endif

	#if defined(HEATFLUX)
		for (int i = 1; i <= nx-2; i++) {
			if (i <= nx / 2 - 1) heatflux(myArray, i, 1, 1);
		}
	#endif

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.th[i][k] += TIMETS * (myArray.thp[i][k] - 2 * myArray.th[i][k] + myArray.thm[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::ppi_pt(Array& myArray) {
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			myArray.pip[i][k] = myArray.pim[i][k] - pow(cs, 2) / (myArray.rhou[k] * C_p * pow(myArray.tvb[k], 2)) * 
								(d2t * rdx * myArray.rhou[k] * myArray.tvb[k] * (myArray.u[i+1][k] - myArray.u[i][k])
								 + 0.5 * d2t * rdz * (myArray.rhow[k+1] * myArray.w[i][k+1] * (myArray.tb[k+1] + myArray.tb[k]) 
								 - myArray.rhow[k] * myArray.w[i][k] * (myArray.tb[k] + myArray.tb[k-1])));
			
			#ifdef DIFFUSION
				myArray.pip[i][k] += d2t * Kx * rdx2 * (myArray.pim[i+1][k] - 2. * myArray.pim[i][k] + myArray.pim[i-1][k]) + 
									 d2t * Kz * rdz2 * (myArray.pim[i][k+1] - 2. * myArray.pim[i][k] + myArray.pim[i][k-1]);
			#endif
		}
	}
	myArray.BoundaryProcess(myArray.pip);

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.pi[i][k] += TIMETS * (myArray.pip[i][k] - 2 * myArray.pi[i][k] + myArray.pim[i][k]);
			}
		}
	#endif


	return;
}


void Iteration::pqv_pt(Array & myArray) {
	double puqv_px = 0., prhowqv_pz_rho = 0., wpqvb_pz = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puqv_px = (myArray.u[i+1][k] * 0.5*(myArray.qv[i+1][k] + myArray.qv[i][k]) - myArray.u[i][k] * 0.5*(myArray.qv[i][k] + myArray.qv[i-1][k])) * rdx;
			prhowqv_pz_rho = (myArray.rhow[k+1] * myArray.w[i][k+1] * 0.5*(myArray.qv[i][k+1] + myArray.qv[i][k]) - 
							  myArray.rhow[k] * myArray.w[i][k] * 0.5*(myArray.qv[i][k] + myArray.qv[i][k-1])) * rdz / myArray.rhou[k];
			wpqvb_pz = 0.5*(myArray.w[i][k+1] + myArray.w[i][k]) * (0.5*(myArray.qvb[k+1] + myArray.qvb[k]) - 0.5*(myArray.qvb[k] + myArray.qvb[k-1])) * rdz;
			myArray.qvp[i][k] = myArray.qvm[i][k] + d2t * (-puqv_px - prhowqv_pz_rho - wpqvb_pz);

			// diffusion
			#ifdef DIFFUSION
				myArray.qvp[i][k] += d2t * Kx * rdx2 * (myArray.qvm[i+1][k] - 2. * myArray.qvm[i][k] + myArray.qvm[i-1][k]) + 
									 d2t * Kz * rdz2 * (myArray.qvm[i][k+1] - 2. * myArray.qvm[i][k] + myArray.qvm[i][k-1]);
			#endif

			// negative qv process: (source)
			if (myArray.qvp[i][k] + myArray.qvb[k] < 0.) myArray.qvp[i][k] = -myArray.qvb[k];
		}
	}
	myArray.BoundaryProcess(myArray.qvp);

	// saturation process: if cloudless.
	#ifdef CLOUDLESS
		for (int i = 1; i <= nx-2; i++) {
			for (int k = 1; k <= nz-2; k++) {
				double pc = 380. / (pow(myArray.pib[k], C_p / Rd) * P0);   // coefficient
				double pth = myArray.thp[i][k] + myArray.tb[k];
				double qvs = pc * exp(17.27 * (myArray.pib[k] * pth - 273.) / (myArray.pib[k] * pth - 36.));
				if (myArray.qvp[i][k] + myArray.qvb[k] > qvs) {
					myArray.qvp[i][k] = qvs - myArray.qvb[k];
				}
			}
		}
		myArray.BoundaryProcess(myArray.qvp);
	#endif


	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.qv[i][k] += TIMETS * (myArray.qvp[i][k] - 2 * myArray.qv[i][k] + myArray.qvm[i][k]);
			}
		}
	#endif

	return;
}

void Iteration::pqc_pt(Array & myArray) {
	double puqc_px = 0., prhowqc_pz_rho = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puqc_px = (myArray.u[i+1][k] * 0.5*(myArray.qc[i+1][k] + myArray.qc[i][k]) - myArray.u[i][k] * 0.5*(myArray.qc[i][k] + myArray.qc[i-1][k])) * rdx;
			prhowqc_pz_rho = (myArray.rhow[k+1] * myArray.w[i][k+1] * 0.5*(myArray.qc[i][k+1] + myArray.qc[i][k]) - 
							  myArray.rhow[k] * myArray.w[i][k] * 0.5*(myArray.qc[i][k] + myArray.qc[i][k-1])) * rdz / myArray.rhou[k];

			myArray.qcp[i][k] = myArray.qcm[i][k] + d2t * (-puqc_px - prhowqc_pz_rho);
			// negative qc process
			if (myArray.qcp[i][k] < 0.) myArray.qcp[i][k] = 0.;

			// saturation process: sink and source (qv <--> qc)
			condensation(myArray, i, k);

			#ifdef DIFFUSION
				myArray.qcp[i][k] += d2t * Kx * rdx2 * (myArray.qcm[i+1][k] - 2. * myArray.qcm[i][k] + myArray.qcm[i-1][k]) + 
									 d2t * Kz * rdz2 * (myArray.qcm[i][k+1] - 2. * myArray.qcm[i][k] + myArray.qcm[i][k-1]);
			#endif
		}
	}
	myArray.BoundaryProcess(myArray.qcp);

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.qc[i][k] += TIMETS * (myArray.qcp[i][k] - 2 * myArray.qc[i][k] + myArray.qcm[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::pqr_pt(Array & myArray) {
	double puqr_px = 0., prhowVTqr_pz_rho = 0.;
	for (int i = 1; i <= nx-2; i++) {
		for (int k = 1; k <= nz-2; k++) {
			puqr_px = (myArray.u[i+1][k] * 0.5*(myArray.qr[i+1][k] + myArray.qr[i][k]) - myArray.u[i][k] * 0.5*(myArray.qr[i][k] + myArray.qr[i-1][k])) * rdx;
			// TODO: VT
			double VT = 6.;
			prhowVTqr_pz_rho = (myArray.rhow[k+1] * (myArray.w[i][k+1] - VT) * 0.5*(myArray.qr[i][k+1] + myArray.qr[i][k]) - 
							    myArray.rhow[k] * (myArray.w[i][k] - VT) * 0.5*(myArray.qr[i][k] + myArray.qr[i][k-1])) * rdz / myArray.rhou[k];
			myArray.qrp[i][k] = myArray.qrm[i][k] + d2t * (-puqr_px - prhowVTqr_pz_rho);
			// negative qr process
			if (myArray.qrp[i][k] < 0.) myArray.qrp[i][k] = 0.;
			
			#ifdef DIFFUSION
				myArray.qrp[i][k] += d2t * Kx * rdx2 * (myArray.qrm[i+1][k] - 2. * myArray.qrm[i][k] + myArray.qrm[i-1][k]) + 
									 d2t * Kz * rdz2 * (myArray.qrm[i][k+1] - 2. * myArray.qrm[i][k] + myArray.qrm[i][k-1]);
			#endif

			autoconversion(myArray, i, k);
			accretion(myArray, i, k);
			evaporation(myArray, i, k);
		}
	}
	myArray.BoundaryProcess(myArray.qrp);

	// Time filter
	#ifdef TIMEFILTER
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.qr[i][k] += TIMETS * (myArray.qrp[i][k] - 2 * myArray.qr[i][k] + myArray.qrm[i][k]);
			}
		}
	#endif
	return;
}

void Iteration::condensation(Array & myArray, int i, int k) {
	double pc = 380. / (pow(myArray.pib[k], C_p / Rd) * P0); 	// coefficient
	double pth = myArray.thp[i][k] + myArray.tb[k];
	double qvs = pc * exp(17.27 * (myArray.pib[k] * pth - 273.) / (myArray.pib[k] * pth - 36.));
	double phi = qvs * (17.27 * 237. * Lv) / (C_p * pow(pth * myArray.pib[k] - 36., 2));

	double C = (myArray.qvp[i][k] + myArray.qvb[k] - qvs) / (1 + phi); 

	// C should less than qc
	if (fabs(C) > myArray.qcp[i][k] && C < 0) C = -myArray.qcp[i][k];
	
	myArray.qvp[i][k] = myArray.qvp[i][k] - C;
	myArray.qcp[i][k] = myArray.qcp[i][k] + C;
	myArray.thp[i][k] = myArray.thp[i][k] + Lv / (C_p * myArray.pib[k]) * C;
}

// autoconversion of qc to qr
void Iteration::autoconversion(Array & myArray, int i, int k) {
	double autort = 0.001, autotr = 0.001; // autocon rate [1/sec], autocon threshold [kg/kg]
	double qcplus = max(0., myArray.qcp[i][k]);
	double ar = autort * (qcplus - autotr);
	ar = max(0., ar);
	double arcrdt = min(ar * d2t, qcplus);
	myArray.qcp[i][k] = myArray.qcp[i][k] - arcrdt;
	myArray.qrp[i][k] = myArray.qrp[i][k] + arcrdt;
	return;
}

// accretion of qc by qr
void Iteration::accretion(Array & myArray, int i, int k) {
	double accrrt = 2.2; // accretion rate [1/sec]
	double qcplus = max(0., myArray.qcp[i][k]);
	double qrplus = max(0., myArray.qrp[i][k]);

	double cr = myArray.rhou[k] * accrrt * qcplus * pow(qrplus, 0.875);
	double arcrdt = min(cr * d2t, qcplus);

	myArray.qcp[i][k] = myArray.qcp[i][k] - arcrdt;
	myArray.qrp[i][k] = myArray.qrp[i][k] + arcrdt;
	return;
}

// evaporation of rain water
void Iteration::evaporation(Array & myArray, int i, int k) {
	double qrplus = max(0., myArray.qrp[i][k]);
	double qvplus = max(0., myArray.qvp[i][k] + myArray.qvb[k]);

	double pc = 380. / (pow(myArray.pib[k], C_p / Rd) * P0);	 // coefficient
	double pth = myArray.thp[i][k] + myArray.tb[k];
	double qvs = pc * exp(17.27 * (myArray.pib[k] * pth - 273.) / (myArray.pib[k] * pth - 36.));	// Tetens equation

	double coef = 1.6 + 30.39 * pow((myArray.rhou[k] * qrplus), 0.2046);	// ventilation coef.
	double deficit = max((1. - qvplus / qvs), 0.);							// saturation dificit (RH < 100%)

	double er = coef * deficit * (pow(myArray.rhou[k] * qrplus, 0.525)) / 
				((2.03e4 + 9.584e6 / (myArray.pb[k] * qvs)) * myArray.rhou[k]);
	double erdt = min(qrplus, max(0., er * d2t));

	myArray.qrp[i][k] = myArray.qrp[i][k] - erdt;
	myArray.qvp[i][k] = myArray.qvp[i][k] + erdt;
	myArray.thp[i][k] = myArray.thp[i][k] - Lv * erdt / (C_p * myArray.pib[k]);
	return;
}

void Iteration::heatflux(Array & myArray, int i, int k, int ishflux) {
	if (ishflux == 1 && k == 1) {
		double cdh = 7e-3;
		double tground = 303.;
		double tdif = max(tground - (myArray.thm[i][k] + myArray.tb[k]), 0.);
		double avgu = 0.5 * abs(myArray.u[i+1][k] + myArray.u[i][k]);
		avgu = max(avgu, 2.);
		double wnetc = 2. * sqrt(tdif);      // "convective velocity" adjustment
		double vel = sqrt(pow(avgu, 2) + pow(wnetc, 2));
		myArray.thp[i][k] = myArray.thp[i][k] + d2t * cdh * vel * myArray.addflx[i] * tdif * rdz;
	}
}

void Iteration::LeapFrog(Array & myArray) {
	int n = 0;
	double timenow = 0.;
	double temp = TIMEEND / dt;
	int nmax = (int) temp;
	while (n < nmax) {
		cout << n << endl;
		// output
		if (n % WRITEFILESTEP == 0 || n == 1) {
			Output::output_pi(n, myArray);
			Output::output_th(n, myArray);
			Output::output_u(n, myArray);
			Output::output_w(n, myArray);
			Output::output_qv(n, myArray);
			Output::output_qc(n, myArray);
			Output::output_qr(n, myArray);
		}
		n++;
		timenow = n * dt;

		// calculate
		pu_pt(myArray);
		pw_pt(myArray);
		pth_pt(myArray);
		ppi_pt(myArray);

		pqv_pt(myArray);
		pqc_pt(myArray);
		pqr_pt(myArray);

		// next step
		for (int i = 0; i <= nx-1; i++) {
			for (int k = 0; k <= nz-1; k++) {
				myArray.um[i][k] = myArray.u[i][k];
				myArray.u[i][k] = myArray.up[i][k];

				myArray.wm[i][k] = myArray.w[i][k];
				myArray.w[i][k] = myArray.wp[i][k];

				myArray.thm[i][k] = myArray.th[i][k];
				myArray.th[i][k] = myArray.thp[i][k];

				myArray.pim[i][k] = myArray.pi[i][k];
				myArray.pi[i][k] = myArray.pip[i][k];

				myArray.qvm[i][k] = myArray.qv[i][k];
				myArray.qv[i][k] = myArray.qvp[i][k];

				myArray.qcm[i][k] = myArray.qc[i][k];
				myArray.qc[i][k] = myArray.qcp[i][k];

				myArray.qrm[i][k] = myArray.qr[i][k];
				myArray.qr[i][k] = myArray.qrp[i][k];
			}
		}

	}
	return;

}


