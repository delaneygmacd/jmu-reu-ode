#include "psmGeneralFunctions.h"
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

bool debug = 0;

int PSM::getMaxDegree(){
	return maxDegree;
}

void PSM::setMaxDegree(int maxDegree){
	this->maxDegree = maxDegree;
} 

template <class T>
T PSM::nthCoefficientProduct(const vector<T> &x, const vector<T> &y, int n){
	/* This method computes the nth coefficient of x*y, runs in O(n)
	 */
	T sum = 0;
	for(int i  = 0; i<=n; i++){
		sum+=x[i]*y[n-i];
	}
	return sum;
}

template <class T>
T PSM::eval(const vector<T> &coeff, T x){
	/* Evaluates a Polynomial given its coefficients
	 * Inputs:
	 * coeff: coefficients of the polynomial to be evaluated in increasing order
	 * x    : number to be evaluated
	 * Output:
	 * p(x) 
	 */
	
	//Horners Algorithm
	T val = 0;
	for(int i = coeff.size()-1; i>=0; i--){
		val = val*x+coeff[i];
	}
	return val;
}

template <class T>
vector<T> PSM::evaluateAll(const vector<vector<T>> &coeff, T x){
	/* Evaluates a set of polynomials at x
	 */
	int n = coeff.size();
	vector<T> sol(n);
	for(int i = 0; i<n; i++){
		sol[i] = PSM::eval(coeff[i],x);
	}
	return sol;
}

template <class T>
vector<T> PSM::calculateDerivative(const vector<T> &f){
	/*This method calculates the derivative of a polynomial*/
	vector<T> deriv(f.size()-1,0);
	for(unsigned int i = 0; i< deriv.size(); i++){
		deriv[i] = f[i+1]*(i+1);
	}
	return deriv;
}



//nn means degree
template <class T>
vector<vector<T>> PSM::computeCoefficients(const vector<T> &parameters, const vector<T> &initialConditions, int nn){
	/* Inputs
	 * parameters: parameters of the equation
	 * initialConditions: initial conditions
	 * nn : order of the polynomial approximation

	 * Output:
	 * coefficients of the nth polynomial approximation in increasing order
	 * 
	 * Runs in O(n^2), memory O(n)
	 */
     //Initial Conditions
	int k = initialConditions.size();
	vector<vector<T>> coefficients(k,vector<T> (nn,0));
	
	//Add initial conditions
	for(int i = 0; i<k; i++){
		coefficients[i][0] = initialConditions[i];
	}
	/*Calculate solution*/
	for(int n = 0; n<nn-1; n++){
		//This function has to be defined for every problem
		updateNthCoefficient(parameters,coefficients,n);
	}
	return coefficients;
}


template <class T>
PSM::Solution PSM::findSolution(const vector<T> &parameters,
        const vector<T> &initialConditions, double step, double end, int n, bool forward){
	/*This function is the general PSM  solver
	 * Input:
	 * parameters: vector of parameters of the ODE
	 * initialConditions: vector of initial conditions, it's a vector for systems of ODEs 
	 * step: size of the step
	 * end : end of the interval
	 * forward: direction in which we should move
	 */
	int numberOfEquations = initialConditions.size();
	vector<vector<T>> solution(numberOfEquations); //We might want to include any auxiliary variables in here as well
	vector<double> steps;
	steps.push_back(0);
	//Initialize initial conditions
	vector<T> currentInitialConditions(numberOfEquations);
	for(int i = 0; i<numberOfEquations; i++){
		solution[i].push_back(initialConditions[i]);
		currentInitialConditions[i] = initialConditions[i];
	}
	//Start Stepping
	for(int i = 1; step*i<= end; i++){
		vector<vector<T>> coeff = PSM::computeCoefficients(parameters,currentInitialConditions,n);
		if(forward){
			steps.push_back(step*i);
			currentInitialConditions = PSM::evaluateAll(coeff,step);
		}
		else{
			steps.push_back(-step*i);
			currentInitialConditions = PSM::evaluateAll(coeff,-step);
		}
		for(int j = 0; j<numberOfEquations; j++){
			solution[j].push_back(currentInitialConditions[j]);
		}	
	}
	Solution sol {steps, solution};
	return sol;
	//templated struct
} 


//Adaptive Stuff



template <class T>
T PSM::findRadiusOfConvergence(const vector<T> &a,bool EvenRadius){
	// This method returns 1/r, where r is an approximation to the radius of convergence
	// returns -1 if n<8
	int n = a.size();
	if(n<8){
		T lastNonZeroCoefficient = max(fabs(a[n-1]),fabs(a[n-2]));
		T radius = pow(lastNonZeroCoefficient,1.0/(n-1)); //This is an estimation to the radius of convergence
		
		return radius;
	}
	if((n%2 == 0 && EvenRadius) || (n%2 != 0 && !EvenRadius)){
		n-=2;
	}
	else{
		n--;
	}
	T xn= pow(fabs(a[n]),1.0/n);
	T xn1 = pow(fabs(a[n-2]),1.0/(n-2));
	T xn2 = pow(fabs(a[n-4]),1.0/(n-4));
	T bn = xn-pow(xn-xn1,2)/((xn-xn1)-(xn1-xn2));
	return bn;
}


template <class T> 
T PSM::approximateRadiusOfConvergence(const vector<T> &coeff){
	return max(findRadiusOfConvergence(coeff,1),findRadiusOfConvergence(coeff,0));	

}  

template <class T>
T PSM::approximateRadiusOfConvergence(const vector<vector<T>> &coeff){
	int n = coeff.size();
	T invR = approximateRadiusOfConvergence(coeff[0]);
	for(int i = 1; i<n; i++){
		invR = max(invR,approximateRadiusOfConvergence(coeff[i]));
	}
	return 1.0/invR;
}  

/*As for now I'll just have this as a global variable
I should change it later on
*/


//The following are 3 adaptive methods


template <class T>
PSM::Solution PSM::findSolutionAdaptive(const vector<T> &parameters,
        const vector<T> &initialConditions, double end, bool forward,double eps){
			
	int numberOfEquations = initialConditions.size();
	vector<vector<T>> solution(numberOfEquations); //We might want to include any auxiliary variables in here as well
	vector<double> steps;
	steps.push_back(0);
	//Initialize initial conditions
	vector<T> currentInitialConditions(numberOfEquations);
	for(int i = 0; i<numberOfEquations; i++){
		solution[i].push_back(initialConditions[i]);
		currentInitialConditions[i] = initialConditions[i];
	}
	
	//Start Stepping
	double cur = 0;
	double step = 0;
	while(cur< end){
		v-ector<vector<double>> coeff = PSM::computeCoefficients(parameters,currentInitialConditions,maxDegree);
		step = pow(eps,1.0/(coeff[0].size()-1))*approximateRadiusOfConvergence(coeff); //This is an estimation to the radius of convergence
		cur+= step;
		if(forward){
			steps.push_back(cur);
			currentInitialConditions = PSM::evaluateAll(coeff,step);
		}
		else{
			steps.push_back(-cur);
			currentInitialConditions = PSM::evaluateAll(coeff,-step);
		}
		for(int j = 0; j<numberOfEquations; j++){
			solution[j].push_back(currentInitialConditions[j]);
		}	
	}
	Solution sol {steps, solution};
	return sol;				
}

template <class T>
PSM::Solution PSM::findAdaptiveSolutionTruncation(const vector<T> &parameters,
        const vector<T> &initialConditions, double end, bool forward,double eps){
	
	int numberOfEquations = initialConditions.size();
	vector<vector<T>> solution(numberOfEquations); //We might want to include any auxiliary variables in here as well
	vector<double> steps;
	steps.push_back(0);
	//Initialize initial conditions
	vector<T> currentInitialConditions(numberOfEquations);
	for(int i = 0; i<numberOfEquations; i++){
		solution[i].push_back(initialConditions[i]);
		currentInitialConditions[i] = initialConditions[i];
	}
	
	//Start Stepping
	double cur = 0;
	double step = 0;
	while(cur< end){
		vector<vector<T>> coeff = PSM::computeCoefficients(parameters,currentInitialConditions,maxDegree);
		T lastNonZeroCoefficient = max(fabs(coeff[0][maxDegree-1]),fabs(coeff[0][maxDegree-2]));
		step = pow(fabs(eps/(2*lastNonZeroCoefficient)),1.0/(maxDegree-1)); //This is an estimation to the radius of convergence
		cur+= step;
		if(forward){
			steps.push_back(cur);
			currentInitialConditions = PSM::evaluateAll(coeff,step);
		}
		else{
			steps.push_back(-cur);
			currentInitialConditions = PSM::evaluateAll(coeff,-step);
		}
		for(int j = 0; j<numberOfEquations; j++){
			solution[j].push_back(currentInitialConditions[j]);
		}	
	}
	Solution sol {steps, solution};
	return sol;		
	
}


template <class T>
PSM::Solution PSM::findAdaptiveSolutionJorbaAndZou(const vector<T> &parameters,
        const vector<T> &initialConditions, double end, bool forward,double eps){
	
	int numberOfEquations = initialConditions.size();
	vector<vector<T>> solution(numberOfEquations); //We might want to include any auxiliary variables in here as well
	vector<double> steps;
	steps.push_back(0);
	//Initialize initial conditions
	vector<T> currentInitialConditions(numberOfEquations);
	for(int i = 0; i<numberOfEquations; i++){
		solution[i].push_back(initialConditions[i]);
		currentInitialConditions[i] = initialConditions[i];
	}
	
	//Start Stepping
	double cur = 0;
	double step = 0;
	while(cur< end){
		double A = 1;
		int degree = ceil(-log(eps/A)/2.0 -1);
		vector<vector<T>> coeff = PSM::computeCoefficients(parameters,currentInitialConditions,degree);
		step = approximateRadiusOfConvergence(coeff)/exp(2.0);//This is an estimation to the radius of convergence
		cur+= step;
		if(forward){
			steps.push_back(cur);
			currentInitialConditions = PSM::evaluateAll(coeff,step);
		}
		else{
			steps.push_back(-cur);
			currentInitialConditions = PSM::evaluateAll(coeff,-step);
		}
		for(int j = 0; j<numberOfEquations; j++){
			solution[j].push_back(currentInitialConditions[j]);
		}	
	}
	Solution sol {steps, solution};
	return sol;		
	
}

