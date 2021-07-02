#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;

ofstream out;

double C, M;

//Cauchy Product
template <class T>
T nthCoefficientProduct(vector<T> &x, vector<T> &y, int n){
	/* This method computes the nth coefficient of x*y, runs in O(n)
	 */
	T sum = 0;
	for(int i  = 0; i<=n; i++){
		sum+=x[i]*y[n-i];
	}
	return sum;
}

template <class T>
void updateNthCoefficient(vector<T> &y, vector<T> &u, vector<T> &w, int n){
    y[n+1] = u[n]/(n+1);
    u[n+1] = (nthCoefficientProduct(u,w,n))/(n+1);
    w[n+1] = -nthCoefficientProduct(u,u,n)/(n+1);
}

//nn means degree
template <class T>
vector<T> coefficients(T x0, int nn){
	/* Inputs
	 * nn : order of the polynomial approximation
	 * x0: initial conditions 
	 * Output:
	 * coefficients of the nth polynomial approximation in increasing order
	 * 
	 * Runs in O(n^2), memory O(n)
	 */
	 vector<T> y(nn,0);
	 vector<T> u(nn,0);
	 vector<T> w(nn,0);
     
     //Initial Conditions
	 y[0] = x0;  
	 u[0] = sin(x0);
     w[0] = cos(x0);

	 for(int n = 0; n<nn-1; n++){
        updateNthCoefficient(y,u,w,n);
	 }
	 return y;
}

template <class T>
T eval(vector<T> coeff, T x){
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
vector<T> derivative(vector<T> f){
	vector<T> deriv(f.size()-1,0);
	for(int i = 0; i< deriv.size(); i++){
		deriv[i] = f[i+1]*(i+1);
	}
	return deriv;
}

/* template <class T>
vector<T> solve(T x0,T step, T end, int n, bool forward){
    vector<T> solution;
    solution.push_back(x0);
    int step_count = 0;
    T newx0 = x0;
    for(int i = 1; step*i<= end; i++){
        //if use_double ....
        vector<T> coeff = coefficients(newx0,n);
        if(forward){
            newx0 = eval(coeff,step);
        }
        else{
            newx0 = eval(coeff,-step);
        }
        solution.push_back(newx0);
        step_count++;
    }
    return solution;
} */

//newx0 = eval(coefficients(newx0,n), (forward ? step ; -step));
//newx0 = eval(coefficients((float)newx0,n), (forward ? step ; -step));

template <class T>
vector<T> findDerivative(T x0,T step, T end, int n, bool forward){
    vector<T> deriv; //Keeps track of the derivative 
    deriv.push_back(sin(x0));
    T newx0 = x0;
    
    for(int i = 1; step*i<= end; i++){
        vector<T> coeff = coefficients(newx0,n);
        if(forward){
            newx0 = eval(coeff, step);
        }
        else{
            newx0 = eval(coeff, -step);
        }
        T dv = eval(derivative(coeff),newx0);
        deriv.push_back(dv);
    }
    return deriv;
}

template <class T>
T acot(T x) {
  return (2*atan(1)) - atan(x);
}

template <class T>
T cot(T x){
    return 1/tan(x);
}

template <class T>
T evaluateFunction(T t,T x0){
    T res = 2*acot(exp(-t)*cot(x0/2));
    return res;
}

int main(int argc, const char* argv[]){

    if (argc != 6) {
        cout << "Usage: " << argv[0] << " <x0> <tn> <dt> <n> <type>" << endl;
        return EXIT_FAILURE;
    }

    double x0 = stod(argv[1]);
    double step = stod(argv[2]);
    double end = stod(argv[3]);
    int n = stod(argv[4]);
    bool use_double = stod(argv[5]);
    bool forward;

    if (use_double){
        double newx0 = x0;
        int step_count = 0;
        vector<double> solution;
        solution.push_back(newx0);
        for(int i = 1; step*i<= end; i++){
        //if use_double ....
            vector<double> coeff = coefficients(newx0,n);
        if(forward){
            newx0 = eval(coeff,step);
        }
        else{
            newx0 = eval(coeff,-step);
        }            solution.push_back(newx0);
            step_count++;
        }      
    }
    else{
        float newx0 = x0;
        int step_count = 0;
        vector<float> solution;
        solution.push_back(newx0);
        for(int i = 1; step*i<= end; i++){
        //if use_double ....
            vector<float> coeff = coefficients(newx0,n);
        if(forward){
            newx0 = eval(coeff,step);
        }
        else{
            newx0 = eval(coeff,-step);
        }            solution.push_back(newx0);
            step_count++;
        }
    }

    }
    //vector<double> dSolF = findDerivative(x0,step,end,n,1);
    //vector<double> solB = solve(x0, step,end,n, 0);
    //vector<double> dsolB = findDerivative(x0,step,end,n,0);

    //vector<float> solF = solve(x0, step,end,n, 1);
    //vector<float> dSolF = findDerivative(x0,step,end,n,1);
    //vector<float> solB = solve(x0, step,end,n, 0);
    //vector<float> dsolB = findDerivative(x0,step,end,n,0);


    double Terror = 0;
    double error = 0;   
    double difference = 0;
    //out.open("out.dat");
    /* for(int i = solB.size()-1; i>0; i--){
        C = M = 1;
        if(fabs(solB[i-1])>1){
            C = M = (solB[i-1]);
        }
        Terror = (C*pow(M*step,n+1))/(1-fabs(M*step));
        difference = fabs(sin(solB[i])-dsolB[i]); //|sin(y)-y'|
        error = abs(evaluateFunction(-i*step,x0)-solB[i]);
        cout<<setw(15)<<-i*step<<setw(15)<<solB[i]<<setw(15)<<error<<setw(15)<<difference<<setw(15)<<error<<"\n";       
    }
    cout<<setw(15)<<0<<setw(15)<<solF[0]<<setw(15)<<0<<"\n";       
    for(unsigned int i = 1; i<solF.size();i++){
        C = M = 1;
        if(fabs(solF[i-1])>1){
            C = M = (solF[i-1]);
        }
        Terror = (C*pow(M*step,n+1))/(1-fabs(M*step));
        difference = fabs(sin(solF[i])-dSolF[i]); //|sin(y)-y'|
        error = abs(evaluateFunction(i*step,x0)-solF[i]);
        cout<<setw(15)<<i*step<<setw(15)<<solF[i]<<setw(15)<<Terror<<setw(15)<<difference<<setw(15)<<error<<"\n";
    } */
   // out.close();
    return 0;
}