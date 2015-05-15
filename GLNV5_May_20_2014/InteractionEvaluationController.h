//
//  InteractionEvaluationController.h
//  GLN
//
//  Created by Tyler Hunt on 7/31/12.
//
//  Singleton
//  The Object of this class will control the systems actions based on the pval mode
//  given in the arguments
//

#ifndef GLN_InteractionEvaluationController_h
#define GLN_InteractionEvaluationController_h
#define DPVALMODE 3 // default pvalue mode for the system

#include "TransitionTable.h"
#include "GLNConfig.h"
#include "ConditionalAssociation.h"

class InteractionEvaluationController{
public:
    // the only way to get a reference to the controller
    static InteractionEvaluationController& getController();
    
    // set the P-value mode (changes bool condition also)
    void setPvalMode(int mode);
    
    // apply the Statistic based on the pValue mode
    void applyStatistic(TransitionTable &table);
    
    // return the pvalue given a statistic and degrees of freedome (currently only supporst chisquare/gtest)
    double getPval(double statisicResult, int df);
    
    //Set options in this class based on a GLNConfig
    void setGLNConfig(GLNConfig gc){m_gc = gc; setPvalMode(m_gc.m_pValueMode); };
    
    // Apply statistics to a table(make a decicion about conditional/non conditional fill first)
    //TransitionTable evaluate(int child, const vector<int>  &parents, const vector<int> &delays, TrajectoryCollection trajCol);
    
    // Same as above only takes imput as ready made table (results are then stored in said table)
    void evaluate(TransitionTable &table, const TrajectoryCollection &trajCol);
    
protected:
    int pvalMode;
    GLNConfig m_gc;
    bool conditional;
    static InteractionEvaluationController m_controller;
    vector<int> findPossibleParents(const TrajectoryCollection &trajCol);

private:
    InteractionEvaluationController(){ setPvalMode(DPVALMODE); };  // Private to enforce singleton
    InteractionEvaluationController(InteractionEvaluationController const&){}; //singleton cannot be copied 

};

#endif
