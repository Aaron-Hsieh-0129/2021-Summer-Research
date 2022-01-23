#include <fstream>
#include <iomanip>
#include "Init.hpp"

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