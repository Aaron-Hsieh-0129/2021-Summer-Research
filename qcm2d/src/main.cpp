#include <iostream>
#include "Iteration.h"

Array myArray;

int main(void) {
	Init::Init1d(myArray);
	Init::Init2d(myArray);
	Output::printInit(myArray);
	Iteration::LeapFrog(myArray);
	return 0;
}
