//
//  EnvGLNRecParentSets.cpp
//  gln
//
//  Created by Joe Song on 10/4/11.
//  Copyright 2011 New Mexico State University. All rights reserved.
//

#include <iostream>
#include <algorithm>

#include "EnvGLNRec.h"

bool compare_chisq(double p_1, int numofp_1, double chisquare_1, 
				   double p_2, int numofp_2, double chisquare_2);

bool compare_chisq_imp(double p_1, int df_1, double chisq_1, 
                       double p_2, int df_2, double chisq_2);

bool compare_F_imp_over_chisq(double p_1, int df_1, double chisq_1, 
                              double p_2, int df_2, double chisq_2);

bool compare_F_imp(int df_1, double chisq_1, int df_2, double chisq_2);


bool EnvGLNRec::compareParentSets(const TransitionTable & tt1, const TransitionTable & tt2) 
{
	// return true if tt1 is "better" than tt2; return false otherwise
    
    bool improved = false;
    
    switch (m_gc.m_methodRecParCmp) {
            
        case BY_PVAL: // compare p-values directly
            improved = compare_chisq(tt1.getpValue(), 
                                     tt1.getDf(), //m_parents.size(), 
                                     tt1.getChisq(), 
                                     tt2.getpValue(), 
                                     tt2.getDf(), // m_parents.size(), 
                                     tt2.getChisq() 
                                     );
            break;
            
        case BY_PVAL_IMP: // compare p-values of absolute improvement
            improved = compare_chisq_imp(tt1.getpValue(), 
                                         tt1.getDf(), //m_parents.size(), 
                                         tt1.getChisq(), 
                                         tt2.getpValue(), 
                                         tt2.getDf(), // m_parents.size(), 
                                         tt2.getChisq()
                                         );
            break;
            
        case BY_PVAL_F_OVER_CHISQ: // compare p-values of relative improvement
            improved = compare_F_imp_over_chisq(tt1.getpValue(), 
                                                tt1.getDf(), //m_parents.size(), 
                                                tt1.getChisq(), 
                                                tt2.getpValue(), 
                                                tt2.getDf(), // m_parents.size(), 
                                                tt2.getChisq()
                                                );
            break;
            
        case BY_IMP_PER_DF: // compare relative improvement per degrees of freedom
            
            improved = compare_F_imp(tt1.getDf(), tt1.getChisq(),
                                     tt2.getDf(), tt2.getChisq() );
            
            break;

		case BY_IMP_AND_PVAL: // compare improvement only on subsets, otherwise using pval
				//check whether the new parent set contains currect best tt2
				if(includes(tt1.getParents().begin(), tt1.getParents().end(),
								tt2.getParents().begin(),tt2.getParents().end()))
					{

            			improved = compare_F_imp(tt1.getDf(), tt1.getChisq(),
												tt2.getDf(), tt2.getChisq() );

 				}else{

						improved = compare_chisq(tt1.getpValue(), 
														tt1.getDf(), //m_parents.size(), 
														tt1.getChisq(), 
														tt2.getpValue(), 
														tt2.getDf(), // m_parents.size(), 
														tt2.getChisq() 
														);		

				}           

			break;

        default:
            cerr << "ERROR: wrong comparison method!" << endl;
            break;
    }
    
	/*
     return	
     ( tt1.m_pValue < tt2.m_pValue )
     ||	( ( tt1.m_pValue == tt2.m_pValue ) 
     &&  ( tt1.m_parents.size() < tt2.m_parents.size() ) )
     ||	( ( tt1.m_pValue == tt2.m_pValue ) 
     && ( tt1.m_parents.size() == tt2.m_parents.size() )
     && ( tt1.m_chisq > tt2.m_chisq ) );
	 */
    
    return improved;
}
