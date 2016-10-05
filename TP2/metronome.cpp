#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "consts.h"

using namespace std;

// Forward declare functions

void solve(int tau, bool dampened);
float getNext(vector<float> theta, float step, float a, float r);
float getNextWithDampening(vector<float> theta, float step, float a, float r);
float derivative(vector<float> theta, float step);
string recordResult(vector<float> theta);

// Application entry - solves undampened case, then dampened, both for tau = 1000s
int main()
{
	solve(1000, false);
	solve(1000, true);
	return 0;
}

/*
    Solves equation of motion.
    parameter: dampened (bool) - is there dampening?
    parameter: tau (int) - the time are we modeling to
    parameter: _fileNumber (int) - the problem in the TP
*/
void solve(int tau, bool dampened)
{
	ofstream file;
	
	string fileName = "soln" + (dampened ? "dampened" : "undampened") + ".txt";

	file.open(fileName);

	float step = .01;
	float r = step / (2.0 * tau);
	float a = 0.0001;

	
	vector<float> theta = { 0, theta_t0, theta_tmoinsd_t0, 0, 0 };
	
	string write;

	for (int t = 0; t <= tau; t++)
	{
		theta[4] = t;

		theta[0] = dampened ? getNextWithDampening(theta, step, a, r) : 
				      getNext(theta, step, a, r);
		
		theta[3] = derivative(theta, step);

		file << recordResult(theta);

		float temp = theta[1];
		theta[1] = theta[0];
		theta[2] = temp;
		a += dampened ? .0000009 : 0.0;
	}

	file.close();
}

string recordResult(vector<float> theta)
{
	string write = to_string(theta[1]) +
		"\t" +
		to_string(theta[3]) +
		"\t" +
		to_string(theta[4]) +
		"\n";

	return write;
}

float derivative(vector<float> theta, float step)
{
	return (theta[0] - theta[2]) / (2.0 * step);
}

float getNext(vector<float> theta, float step, float a, float r)
{
	return (2.0 * theta[1] - theta[2] + pow(step, 2) * ((g * sin(theta[1]) / l) -
		(C / m1) * theta[1] + a * sin(2 * PI * NU * theta[3])));
}

float getNextWithDampening(vector<float> theta, float step, float a, float r)
{
	return (2.0 * theta[1] + (r - 1)* theta[2] + pow(step, 2) * ((g * sin(theta[1]) / l) -
		(C / m2) * theta[1] + a * sin(2 * PI * NU * theta[3]))) / (r + 1);
}
