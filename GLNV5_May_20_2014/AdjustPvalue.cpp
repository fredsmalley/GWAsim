// AdjustPvalue.cpp 
//
// Curtis Luce, Joe Song
// Created: 2007
// Last modified: September 27, 2008
// Last modified: October 4, 2011.  File name changed from "adjusted_pvalue.cpp"
// Last modified: May 19, 2014. fix factorial function to calculate huge factorials.

#include <iostream>
#include <sstream>
using std::cout;
using std::endl;

#include <cmath>
#include <stdio.h>
#include <math.h>

#include "GLN.h"
#include "TransitionTable.h"
//#include "FactorialByGMP.h"//Hua added, May 19 2014

#if defined _WIN32 || defined _WIN64
#define isinf(x) (!_finite(x))
#endif

#include <boost/math/special_functions/factorials.hpp> //Hua added, Nov 10 2014

using namespace boost::math;

////
//Hua modified, Mar 9, 2015
bool compare_abs(const int & a, const int & b){
    return abs(a)<abs(b);
}

double factorial_vector(vector<int> v){//Hua Add, May 2 2014
    //Hua added, Nov 10 2014, calculate factorial by boost
    sort(v.begin(), v.end(), compare_abs);
    double value=1.0;
    for (int i=0; i<(int)v.size(); i++) {
        if(v[i]>0){
            value = value * factorial<double>(abs(v[i]));
        }else if(v[i]<0){
            value = value / factorial<double>(abs(v[i]));
        }else{
            
        }
    }
    return value;
    ////
}
////

/*
int nfac[10]={1,1,2,6,24,120,720,};

double factorial (int n)
{
    double answer = 1.0;

    if(n <= NFACSIZE)
    switch (n) {
        case 0:
        case 1:
            answer = 1.0;
            break;
        case 2:
            answer = 2.0;
            break;
        case 3:
            answer = 6.0;
            break;
            
        default:
            for (int i = 1; i <= n; i++)
                answer *= i;
            break;
    }
    
	return answer;
}
*/

double factorial (double n)
{
	if(n==0)
		return 1;
    
	double answer = 1;

	for (int i = 1; i <= n; i++) {
		answer *= i;
        if(std::isinf(answer)) {
            break;
        }
    }
	return answer;
}

void adjusted_pvalue(GeneralizedLogicalNetwork & net, vector<TransitionTable> & tts, int numParents)
{
	double adjust = 0.0;
	int fact1;
	int fact2;
	int fact3;
	double adjustedfinal = 0.0;

	fact1 = (int) factorial(net.size());
	//cout << "max offset :  " << net.max_Markov_order << endl;
	//cout << "numParents:  " << numParents << endl;
	fact2 = (int) factorial(net.getMaxMarkovOrder() * -1);
	fact3 = (int) factorial((net.size()) -(net.getMaxMarkovOrder()*-1));

	//cout << fact1 << ":" << fact2 <<":" <<fact3 << endl;

	//adjust = (float)(factorial(net.numNodes)/(factorial(net.max_Markov_order*-1)*factorial(net.numNodes-(net.max_Markov_order*-1))));
	adjust= (fact1 * 1.0)/((fact2*1.0) * (fact3 *1.0));
	//cout << " adjust:  " << adjust << endl;
	
	for (int i = 1; i < numParents; i++)
		adjustedfinal += adjust * pow(net.getMaxMarkovOrder()*-1.0, i*1.0);
	//adjust = 100.0;

	//cout << "num nodes: " << net.numNodes << endl;
	for (size_t i = 0; i < net.size(); i++){
		if(net.getNodeType(i) == 'i'){		
			//cout << "changing p_value" << endl;
			//cout << "old pvalue :  " << net.m_nodes[i].getPvalue() << endl;
			// net.getNodes()[i].setPvalue(net.getNodes()[i].getPvalue() * adjustedfinal);
			
			// net.setNodePvalue((int) i, net.getNodes()[i].getPvalue() * adjustedfinal);
			
			tts[i].setAdjustedpValue(tts[i].getpValue() * adjustedfinal);

			//cout << "new pvalue :  " << net.m_nodes[i].getPvalue() << endl;
		}
	}

	//cout << "adjust pvalue by:  " << adjustedfinal << endl;
}

void print_adjusted(const GeneralizedLogicalNetwork & net, const vector<TransitionTable> & tts, double alpha)
{
	cout << "Adjusted p-values:" << endl;

	cout << "node.id\t" << "parents\t" << "adjusted.p.value" << endl;

	for (size_t i = 0; i < net.size(); i++) {

		cout << i+1 << "\t";
		for (size_t j = 0; j < tts[i].getParents().size(); j++) {

			if(net.getNodeType(i) == 'i'){

				if(tts[i].getAdjustedpValue() <= alpha) {
					// cout << net.getNodes()[i].getParent(j)+1 << ",";
					//cout << tts[i].getParents()[j]+1 << ','; 
					cout << tts[i].getParents()[j] << ",";
				}
			}
		}

		cout << "\t";
		cout <<  tts[i].getAdjustedpValue() << endl;
	}

}

