// binom.cpp -- compute binomial distribution and conduct statistical 
//              testing regarding binomial distribution
// Joe Song
// Created: June 16, 2006
// Modified: March 25, 2014. "register" type removed because it is deprecated.

#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>


using namespace std;

inline long factorial(int k) 
{
	long f = 1;
	for(int i=1; i<=k; i++) {
		f *= i;
	}
	return f;
}

inline double ln_factorial(int k) 
{
    double f = 0;
	for(int i=1; i<=k; i++) {
		f += log((double)i);
	}
	return f;
}

inline long product(int a, int b) 
{
    long f=1;
	int lower = min(a, b);
	int upper = max(a, b);

	for(int i=lower; i<=upper; i++) {
		f *= i;
	}

	return f;
}

inline double ln_product(int a, int b) 
{
    double f=0;
	int lower = min(a, b);
	int upper = max(a, b);

	for(int i=lower; i<=upper; i++) {
		f += log((double)i);
	}

	return f;
}


inline long choose(int n, int x) 
{
	if(x == 0 || x == n) {
		return 1;
	}

	if(x > n || x < 0) {

		return 0;

	} else {
		if(x > n-x) {

			return product(n, x+1) / factorial(n-x);
			
		} else {

			return product(n, n-x+1) / factorial(x);

		}
	}
}

inline double ln_choose(int n, int x) 
{
	if(x == 0 || x == n) {
		return log(1.0);
	}

	if(x > n || x < 0) {

		return 0;

	} else {
		if(x > n-x) {

			return ln_product(n, x+1) - ln_factorial(n-x);
			
		} else {

			return ln_product(n, n-x+1) - ln_factorial(x);

		}
	}
}

double dbinom(int x, int n, double p) 
{
	// compute the binomial distribution
	return choose(n, x) * pow(p, x) * pow(1-p, n-x);
}

double ln_dbinom(int x, int n, double p) 
{
	// compute the binomial distribution
	return ln_choose(n, x) + x * log(p) +  (n-x) * log(1-p);
}

double pbinom(int x, int n, double p) 
{
	// compute the accumulative probability over [0, x].
	double F=0;
	for(int i=0; i<=x; i++) {
		F += dbinom(i, n, p);
	}
	return F;
}

double pbinom_via_log(int x, int n, double p) 
{
	// compute the accumulative probability over [0, x].
	double F=0;
	for(int i=0; i<=x; i++) {
		F += exp(ln_dbinom(i, n, p));
	}
	return F;
}

double BinomialTestOneSample(int x, int n, double psuccess, const char alternative[]) 
{
	double pvalue = 1;
	double p = psuccess;

	if(! strcmp(alternative, "less")) { 
		
		pvalue = pbinom_via_log(x, n, p);

	} else if(! strcmp(alternative, "greater") ) {	
	
		pvalue = 1 - pbinom_via_log(x - 1, n, p);	
		
	} else if(! strcmp(alternative, "two.sided") ) {		

		if (p == 0) {

			pvalue = (x == 0) ? 1 : 0;

		} else if (p == 1) { 

			pvalue = (x == n) ? 1 : 0;

		} else {

            double relErr = 1 + 1e-07;
            double d = exp(ln_dbinom(x, n, p));
            double m = n * p;
			int y=0;

			if (x == m) {

                pvalue = 1;

			} else if (x < m) {
				
				for(int i=(int)ceil(m); i<=n; i++) {
					if(exp(ln_dbinom(i, n, p)) <= d * relErr) {
						y ++;
					}
					// y <- sum(dbinom(i, n, p) <= d * relErr)
				}

                pvalue = pbinom_via_log(x, n, p) + 1 - pbinom_via_log(n - y, n, p);
            }
            else {
                // i <- seq(from = 0, to = floor(m))
                // y <- sum(dbinom(i, n, p) <= d * relErr)
                // pvalue = pbinom(y - 1, n, p) + pbinom(x - 1, n, p, lower = FALSE)
				for(int i = 0; i <= floor(m); i++) {
					if(exp(ln_dbinom(i, n, p)) <= d * relErr) {
						y ++;
					}
				}
                pvalue = pbinom_via_log(y - 1, n, p) + 1 - pbinom_via_log(x - 1, n, p);
            }
        }
    }

	return pvalue;
}

void test_binom()
{
	cout << "4!=" << factorial(4) << " (24 expected)." << endl;

	cout << "0!=" << factorial(0) << " (1 expected)." << endl;

	cout << "5 choose 3 = " << choose(5, 3) << " (10 expected)." << endl;

	cout << "5 choose 0 = " << choose(5, 0) << " (1 expected)." << endl;

	cout << "dbinom(2, 3, 0.7)= " << dbinom(2, 3, 0.7) << " (0.441 expected)." << endl;

	cout << "dbinom(0, 3, 0.7)= " << dbinom(0, 3, 0.7) << " (0.027 expected)." << endl;

	cout << "dbinom(3, 3, 0.7)= " << dbinom(3, 3, 0.7) << " (0.343 expected)." << endl;

	cout << "pbinom(2, 3, 0.7)= " << pbinom(2, 3, 0.7) << " (0.657 expected)." << endl;

    cout << "pbinom(3, 3, 0.7)= " << pbinom(3, 3, 0.7) << " (1 expected)." << endl;

	cout << "binom.test(2,3,0.7) (less) = " << BinomialTestOneSample(2, 3, 0.7, "less") << " (0.657 expected)." << endl;

	cout << "binom.test(2,3,0.7) (greater) = " << BinomialTestOneSample(2, 3, 0.7, "greater") << " (0.784 expected)." << endl;

	cout << "binom.test(7,10,0.8) (two.sided) = " << BinomialTestOneSample(7, 10, 0.8, "two.sided") << " (0.4296 expected)." << endl;

	cout << "binom.test(0,3,0.7) (less) = " << BinomialTestOneSample(0, 3, 0.7, "less") << " (0.027 expected)." << endl;

	cout << "binom.test(0,3,0.7) (greater) = " << BinomialTestOneSample(0, 3, 0.7, "greater") << " (1 expected)." << endl;

	cout << "binom.test(0,3,0.7) (two.sided) = " << BinomialTestOneSample(0, 3, 0.7, "two.sided") << " (0.027 expected)." << endl;

	cout << "binom.test(3,3,0.7) (less) = " << BinomialTestOneSample(3, 3, 0.7, "less") << " (1 expected)." << endl;

	cout << "binom.test(3,3,0.7) (greater) = " << BinomialTestOneSample(3, 3, 0.7, "greater") << " (0.343 expected)." << endl;

	cout << "binom.test(3,3,0.7) (two.sided) = " << BinomialTestOneSample(3, 3, 0.7, "two.sided") << " (0.559 expected)." << endl;

}
