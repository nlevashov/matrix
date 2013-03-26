#include <iostream>
#include "matrix.h"
using namespace std;

int main () {
	matrix<double> a;
	cin >> a;
	cout << a.inv();

	return 0;
}
