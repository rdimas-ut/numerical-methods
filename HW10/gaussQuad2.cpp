#include <iostream>
#include <cmath>

using namespace std;

double gaussQuad2(double f(double), double a, double b){
	double result = 0.0;
	double gauss1 [2] = {sqrt(1.0/3.0), -sqrt(1.0/3.0)};
	double ba = (b - a);
	double ab = (b + a);

	for (int i =0; i < 2; i++){
		result += f(((ba*gauss1[i])+ab)/2.0)*(ba/2.0);
	}

	return result;
}

double test(double x){
	return log(x);
}

int main(){
	double a, b;
	cout << "Please enter a, b: ";
	cin >> a >> b;
	cout << gaussQuad2(test, a, b) << endl;
	return 0;
}
