#include "Outputfile.h"

void Output::printInit(Array & myArray) {
	double z;
	cout << "z          tb        rhou       rhow       qvb	 RH      pib" << endl;
	for (int k = 0; k <= nz-1;k++){
		z = (double) (k - 0.5) * dz ;
		cout << fixed << setprecision(2) << z << "    " << myArray.tb[k] << "    " << myArray.rhou[k] << "     " 
		<< myArray.rhow[k] << "   	 " << myArray.qvb[k] * 1000 << "    " << myArray.qvb[k] / myArray.qvsb[k] << "    "
		<< myArray.pib[k] << endl;
	}
	fstream initout;
	string initName = "../outputs/init.txt";
	initout.open(initName, ios::out);
	for (int k = 0; k <= nz-1; k++) {
		z = (double) (k - 0.5) * dz ;
		initout << z << "    " << myArray.tb[k] << "    " << myArray.rhou[k] << "     " 
		<< myArray.rhow[k] << "   	 " << myArray.qvb[k] << "    " << myArray.qvsb[k] << "    " << myArray.qvb[k] / myArray.qvsb[k] << "    "
		<< myArray.pib[k] << endl;
	}
	return;
}

void Output::output_pi(int n, Array & myArray) {
	fstream foutzeta;
	string zetaName = "../outputs/pi/pi_" + to_string(n) + ".txt";
	foutzeta.open(zetaName, ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutzeta << myArray.pi[i][k] << " ";
		}
	}
}

void Output::output_th(int n, Array & myArray) {
	fstream foutth;
	string thName = "../outputs/th/th_" + to_string(n) + ".txt";
	foutth.open(thName, ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutth << myArray.th[i][k] << " ";
		}
	}
}

void Output::output_u(int n, Array & myArray) {
	fstream foutu;
	string uName = "../outputs/u/u_" + to_string(n) + ".txt";
	foutu.open(uName, ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutu << myArray.u[i][k] << " ";
		}
	}
}

void Output::output_w(int n, Array & myArray) {
	fstream foutw;
	string wName = "../outputs/w/w_" + to_string(n) + ".txt";
	foutw.open(wName, ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutw << myArray.w[i][k] << " ";
		}
	}
}

void Output::output_qv(int n, Array & myArray) {
	fstream foutqv;
	string qvName = "../outputs/qv/qv_" + to_string(n) + ".txt";
	foutqv.open(qvName, ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutqv << myArray.qv[i][k] + myArray.qvb[k] << " ";
		}
	}
}

void Output::output_qc(int n, Array & myArray) {
	fstream foutqc;
	string qcName = "../outputs/qc/qc_" + to_string(n) + ".txt";
	foutqc.open(qcName, ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutqc << myArray.qc[i][k] << " ";
		}
	}
}

void Output::output_qr(int n, Array & myArray) {
	fstream foutqr;
	string qrName = "../outputs/qr/qr_" + to_string(n) + ".txt";
	foutqr.open(qrName, ios::out);
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nx; i++) {
			foutqr << myArray.qr[i][k] << " ";
		}
	}
}
