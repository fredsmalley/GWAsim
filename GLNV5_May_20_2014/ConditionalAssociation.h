//
//  ConditionalAssociation.h
//  GLN
//
//  Created by Tyler Hunt on 8/1/12.
//  Copyright (c) 2012 Tyler Hunt. All rights reserved.
//

#ifndef GLN_ConditionalAssociation_h
#define GLN_ConditionalAssociation_h

#include <vector>

#include "TransitionTable.h"
#include "InteractionEvaluationController.h"

class ConditionalAssociation{

public:
    // This constructor creates the table for you
    //ConditionalAssociation(const vector<int> &parents, int child, const vector<int> & possibleParents,  const vector<int> &delays,
    //                       const TrajectoryCollection &trajCol);
    
    // Constructor for when table is constructed beforehand
    //ConditionalAssociation(TransitionTable& table, const vector<int> &possibleParents);
    
    // Fill the table and evaluate it for all conditions results are stored in m_table which is also returned
    static void evaluate(TransitionTable & table, const vector<int> possibleParents, const TrajectoryCollection& trajCol);
    
protected:
    //vector<int> m_possibleParents; //< A vector containing the indexes of all the nodes of the network that can be parents
    //vector<int> m_possibleParentDelays; //< A network containing the delays of all parents (currently not utilized)
    //TransitionTable* m_table; //< Used to perform calculations and store results
    static void populate(int condition, int level, TransitionTable & table, const TrajectoryCollection & trajCol);

private:
    ConditionalAssociation(){}; // Default construction impossible
    // void clearTable(){ m_table->clear(); };
};

#endif
