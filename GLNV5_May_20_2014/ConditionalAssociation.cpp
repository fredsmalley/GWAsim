//
//  ConditionalAssociation.cpp
//  GLN
//
//  Created by Tyler Hunt on 8/1/12.
//  Copyright (c) 2012 Tyler Hunt. All rights reserved.
//

#include <algorithm>
#include "ConditionalAssociation.h"

void ConditionalAssociation::evaluate(TransitionTable & table,
                                      const vector<int> possibleParents,
                                      const TrajectoryCollection& trajCol)
{
    double statRes = 0;
    double pVal = 0;
    int df = 0;
    vector<int> parents = table.getParents();

    // iterate through the possible parents
    for(size_t i = 0; i < possibleParents.size(); i++)
        for(int lev = 0; lev < trajCol.getBases()[i]; lev++){
            // A condition cannot also be the child
            if(possibleParents[i] == table.getChild()) continue; // need to modify to handle delays
            if(std::find(parents.begin(), parents.end(), possibleParents[i]) != parents.end()) {
                continue; // The condition is not a parent
            }
             // Populate the table
            populate(possibleParents[i], lev, table, trajCol);
            InteractionEvaluationController::getController().applyStatistic(table);
            statRes += table.getChisq();
            df += table.getDf();
        }
    pVal = InteractionEvaluationController::getController().getPval(statRes, df);
    table.setChisq(statRes);
    table.setpValue(pVal);
    table.setDf(df);  
}


/*ConditionalAssociation::ConditionalAssociation(const vector<int> &parents, int child, const vector<int> &possibleParents, const vector<int> &delays, const TrajectoryCollection & trajCol){
    m_table = new TransitionTable(child, parents, delays, trajCol);
    m_possibleParents = possibleParents;
    vector<int> m_possibleParentsDelays = *new vector<int>(possibleParents.size(), 0); // Will need changing to include actual delays for conditions
}


ConditionalAssociation::ConditionalAssociation(TransitionTable &table, const vector<int> &possibleParents){
    m_table = &table;
    m_possibleParents = possibleParents;
    //vector<int> m_possibleParentsDelays = *new vector<int>;
}
*/


void ConditionalAssociation::populate(int condition, int level, TransitionTable & table, const TrajectoryCollection & trajCol){
    table.reset();
    table.condFill(trajCol, condition, 0, level);
}
