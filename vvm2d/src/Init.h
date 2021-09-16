#include <iostream>
#include <random>
#include "Const.h"
#include "Declare.h"

using namespace std;

class Init {
	public:
			static void Init1d(Array &);
			static void Init2d(Array &);
			
	private:
			static double GetTB(int);
			static double GetTHRAD(int, int);
			static double GetTH(int, int);
			static double GetQVB(int);
};