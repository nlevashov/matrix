#include <iostream>
#include "matrix.h"
using namespace std;

int main () {
	matrix<double> a;
	cin >> a;
	cout << a[0][1] << endl;
	matrix<int> b = a;
	cout << b.inv();
	matrix<float> c;
	c = b.inv();
	cout << c;
	cout << c.transp() << c.trace() << endl << c.det() << endl;
	return 0;
}
