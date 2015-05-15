//
//  InteractionEvaluationController.cpp
//  GLN
//
//  Created by Tyler Hunt on 7/31/12.
//
//

#include "InteractionEvaluationController.h"
#include "GLNGlobals.h"
#include "StatDistributions.h"

InteractionEvaluationController InteractionEvaluationController::m_controller = InteractionEvaluationController();

void InteractionEvaluationController::setPvalMode(int mode){
    if((mode >= 1 && mode <= 10) || mode==15){
        conditional = false;
        pvalMode = mode;
    }
    else if(mode >= 101 && mode <= 107){
        pvalMode = mode - 100;
        conditional = true;
    }
    else
        GLNExit(EXIT_FAILURE, "ERROR: unknown p-value mode! (InteractionEvaluationController");
}


void InteractionEvaluationController::applyStatistic(TransitionTable &table){
    switch(pvalMode){
        case 1:
        case 2:
        case 3:
        case 5:
        case 6:
        case 7:
		case 8:
        case 9:
        case 10: // case 10 is added by Joe Song 2/28/2015
        case 15: // case 15 is added by Hua 3/16/2015
            table.applyChisquareTest(pvalMode);
            break;
        case 4:
            table.applyChisqPermTest(pvalMode,m_gc.m_permStatTable,m_gc.m_max_Markov_order,m_gc.m_min_Markov_order);
            break;
        default:
            GLNExit(EXIT_FAILURE, "ERROR: unknown p-value mode! (InteractionEvaluationController");
    }
}


double InteractionEvaluationController::getPval(double statisticResult, int df){
    switch(pvalMode){
        case 1:
        case 2:
            return 0; // need to add
            break;
        case 3:
        case 5:
        case 6:
        case 10: // case 10 is added by Joe Song 2/28/2015
            return ChisqPvalue(statisticResult, df);
            break;
        case 7:
		case 8:
        case 9:
        case 15: // case 15 is added by Hua 3/16/2015
        case 4:
            break; // need to add
        default:
           GLNExit(EXIT_FAILURE, "ERROR: unknown p-value mode! (InteractionEvaluationController"); 
    }
    return 0;
}

InteractionEvaluationController& InteractionEvaluationController::getController(){
    return m_controller;
}

/*TransitionTable InteractionEvaluationController::evaluate(int child, const vector<int>  &parents, const vector<int> &delays, TrajectoryCollection trajCol){
    if(conditional){
        ConditionalAssociation conditional(parents, child, findPossibleParents(trajCol), delays, trajCol);
        return conditional.evaluate(trajCol);
    }
    TransitionTable table(child, parents, delays, trajCol);
    table.fill(trajCol);
    applyStatistic(table);
    return table;
}*/

void InteractionEvaluationController::evaluate(TransitionTable &table, const TrajectoryCollection &trajCol){
    if(conditional){
        // ConditionalAssociation conditional(table, findPossibleParents(trajCol));
        ConditionalAssociation::evaluate(table, findPossibleParents(trajCol), trajCol);
		  //table.clear();
        table.fill(trajCol);
        return;
    }
    table.fill(trajCol);
    applyStatistic(table);
}

vector<int> InteractionEvaluationController::findPossibleParents(const TrajectoryCollection & trajCol){
    vector<int> possibleParents;
    vector<bool> beParent = trajCol.getBeParent();
    //int possibleParentsSize = trajCol.getIntNodeNames().size(); 
    for(size_t i = 0; i < beParent.size(); i++){
        if(beParent[i]){
            possibleParents.push_back((int)i);
        
        }
    }
    return possibleParents;
}
