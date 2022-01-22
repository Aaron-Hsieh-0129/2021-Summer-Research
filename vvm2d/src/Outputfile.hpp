#include <fstream>
#include <iomanip>
#include "Init.h"

class Output {
	public:
			static void printInit(Array &);
			static void output_zeta(int, Array &);
			static void output_th(int, Array &);
			static void output_u(int, Array &);
			static void output_w(int, Array &);
			static void output_qv(int, Array &);
			static void output_qc(int, Array &);
			static void output_qr(int, Array &);
};