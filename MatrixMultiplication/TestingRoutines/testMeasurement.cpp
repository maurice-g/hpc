#include "../HelperClasses/Measurement.hpp"
#include <iostream>
#include <fstream>

int main() {
	Measurement m("just a test",0,0,15.3);
	std::cout << m;
	
	std::ofstream myfile;
	myfile.open("ex.txt");
	myfile << m;
	myfile.close();

return 0;
}
