#include <fstream>
#include <iomanip>
#include "Init.hpp"
#include "/Users/Aaron/miniconda3/include/python3.8/Python.h"
#ifndef WITHOUT_NUMPY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "/Users/Aaron/miniconda3/lib/python3.8/site-packages/numpy/core/include/numpy/arrayobject.h"
#endif

class Output {
	public:
			static void printInit(vvmArray &);
			static void output_zeta(int, vvmArray &);
			static void output_th(int, vvmArray &);
			static void output_u(int, vvmArray &);
			static void output_w(int, vvmArray &);
			static void output_qv(int, vvmArray &);
			static void output_qc(int, vvmArray &);
			static void output_qr(int, vvmArray &);
};


class Plot {
    public:
        static int *initNumpy(void);
        static void initPython(void);
        static void plotInit(vvmArray &);
        static void plot_zeta(int, vvmArray &);
        static void plot_th2d(int, vvmArray &);
        static void plot_u(int, vvmArray &);
        static void plot_w(int, vvmArray &);
        static void plot_qv(int, vvmArray &);
        static void plot_qc(int, vvmArray &);
        static void plot_qr(int, vvmArray &);
        static void plot_qc_qr(int, vvmArray &);
        static void plot_qv_qc(int, vvmArray &);
        static void plot_qc_qr_th_u_w(int, vvmArray &);
        static void plot_qr_th_u_w(int, vvmArray &);
};
