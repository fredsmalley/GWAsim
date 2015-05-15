//
//  CompareChisq.cpp
//  gln
//
//  Created by Joe Song on 9/12/11.
//  Copyright 2011 New Mexico State University. All rights reserved.
//
// Modified:
//  MS 2/14/2013. Add absolute value for all chisq comparison to handle negative
//    chisq.

#include <cmath>
#include <iostream>

using namespace std;

#include "StatDistributions.h"

// double qchisq(int df, double x);

bool compare_chisq(double p_1, int numofp_1, double chisquare_1, 
				   double p_2, int numofp_2, double chisquare_2)
{
	return	
	( p_1 < p_2 )
	||	( ( p_1 == p_2 ) 
		 &&  ( numofp_1 < numofp_2 ) )
	||	( ( p_1 == p_2 ) 
		 && ( numofp_1 == numofp_2 )
		 && ( abs(chisquare_1) > abs(chisquare_2) ) );
}

bool compare_chisq_imp(double p1, int df1, double chisq1, 
                       double p2, int df2, double chisq2)
{
	bool improved;
	
	if(df1 > df2) {
		
		if (abs(chisq1) <= abs(chisq2)) {
			
			improved = false;
			
		} else {
			double chisq_imprv = abs(chisq1) - abs(chisq2);
			int df_imprv = df1 - df2;
			double p_imprv = ChisqPvalue(chisq_imprv, df_imprv); // qchisq(df_imprv, chisq_imprv);
			
			improved = compare_chisq(p_imprv, df_imprv, chisq_imprv, p2, df2, chisq2);
		}
		
	} else if(df1 == df2) {
		
		improved = ( abs(chisq1) > abs(chisq2) );
		
	} else {  // df1 < df2
		
		if ( abs(chisq1) >= abs(chisq2) ) {
			
			improved = true;
			
		} else {
			double chisq_imprv = abs(chisq2) - abs(chisq1);
			int df_imprv = df2 - df1;
			double p_imprv = ChisqPvalue(chisq_imprv, df_imprv); // qchisq(df_imprv, chisq_imprv);
			
			improved = compare_chisq(p1, df1, chisq1, p_imprv, df_imprv, chisq_imprv);
		}		
		
	}
	
	return improved;
    
}

double compute_F_imp(int df1, double chisq1, int df2, double chisq2)
// return F improvement
{
    double improvement;

    if (df1 == 0 && df2 == 0) {

        improvement = 0;
    
    } else if(df2 == 0) {
        
        improvement = abs(chisq1) / df1;
        
    } else if(df1 == df2) { 
        
        improvement = ( abs(chisq1) - abs(chisq2) ) / df1;
        
    } else if(df1 > df2) {
    
        // ">" but not ">=" implies a preference for the smaller df2 
        improvement = ( (abs(chisq1) - abs(chisq2) ) / (df1 - df2) ) - ( abs(chisq2) / df2 );

    } else if(df1 == 0) {
        
        improvement = - ( abs(chisq2) - abs(chisq1) ) / (df2 - df1);
    
    } else { // df1 < df2 and df1 != 0
        
        // "<=" but not "<" implies a preference for the smaller df1
        improvement = ( abs(chisq1) / df1 ) - ( (abs(chisq2) - abs(chisq1) ) / (df2 - df1) );

    }
    
    return improvement;
}

bool compare_F_imp(int df1, double chisq1, int df2, double chisq2)
// return true if & only if the (chisq1, df1) is better than (chisq2, df2)
{
    bool improved;
    
    if (df2 == 0) {
        
        improved = true;

    } else if(df1 == df2) { 
        
        improved = abs(chisq1) > abs(chisq2);
        
    } else if(df1 > df2) {
        
        // ">" but not ">=" implies a preference for the smaller df2 
        improved = ( ( abs(chisq1) - abs(chisq2) ) / (df1 - df2) ) > ( abs(chisq2) / df2 );
        
    } else if(df1 == 0) {
    
        improved = false;
        
    } else { // df1 < df2 and df1 != 0
        
        // "<=" but not "<" implies a preference for the smaller df1
        improved = ( (abs(chisq2) - abs(chisq1) ) / (df2 - df1) ) <= ( abs(chisq1) / df1 );
        
    }
    
    return improved;
}


bool compare_F_imp_over_chisq(double p1, int df1, double chisq1, double p2, int df2, double chisq2)
{
	bool improved;
	
	if(df1 > df2) {
		
		if ( abs(chisq1) <= abs(chisq2) ) {
			
			improved = false;
			
		} else {
            
            if(df2 > 0 && abs(chisq2) > 0) {

                double p_imprv = FPvalue( ( abs(chisq1) - abs(chisq2) )
                                          / (df1-df2) / (abs(chisq2)/df2),
                                         df1-df2, df2);
                
                /*
                Fdist F(df1 - df2, df2);
                double p_imprv = F.cdf((chisq1 - chisq2)/(df1-df2)/(chisq2/df2));
                */
                
                improved = p_imprv < p2;
                
            } else {
                
                improved = true;
                
            }
		}
	} else if(df1 == df2) {
		
		improved = ( abs(chisq1) > abs(chisq2) );
		
	} else {  // df1 < df2
		
		if ( abs(chisq1) >= abs(chisq2) ) {
			
			improved = true;
			
		} else {
            
            if(df1 > 0 && abs(chisq1) > 0) {

                double p_imprv = FPvalue( (abs(chisq2) - abs(chisq1))
                                         / (df2-df1) / ( abs(chisq1) /df1),
                                         df2 - df1, df1);
                
                /*
                Fdist F(df2 - df1, df1);
                double p_imprv = F.cdf((chisq2 - chisq1)/(df2-df1)/(chisq1/df1));
                */
                
                improved = p_imprv < p2;
                
            } else {
                
                improved = false;
                
            }           
            
		}		
		
	}
	
	return improved;
    
}
