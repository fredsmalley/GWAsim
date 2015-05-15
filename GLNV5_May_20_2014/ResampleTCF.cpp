//
//  ResampleTCF.cpp
//  ResampleGLNV5
//
//  Created by Yang Zhang on 4/23/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <vector>
#include "TrajectoryCollection.h"

#include "GLNGlobals.h"

/*
vector<TrajectoryCollection>
resample(const vector<TrajectoryCollection> & trajCols,
         ResampleType strategy)
{
    vector<TrajectoryCollection> trjColsBS = trajCols;
    
    for(size_t trjIndex = 0; trjIndex < trajCols.size(); trjIndex++)
    {
        const vector<vector<vector<int> > > & trajs = trajCols[trjIndex].getTrajTables();
        
        vector<vector<vector<int> > > trajsFinal(trajs.size());

        vector<vector<int> > rows;  //reduce 3-d TCF into 2-d, for bootstrapping purpose
        
        for(size_t i=0; i<trajs.size(); i++)
        {
            for(size_t j=0; j<trajs[i].size(); j++)
            {
                rows.push_back(trajs[i][j]);
            }
        }
        
        int rownum = rows.size();
        
        for(size_t i=0; i<trajsFinal.size(); i++)
        {
            switch (strategy) {
                    
                case BOOTSTRAP_PLUS_ONE:
                    trajsFinal[i].resize(trajs[i].size()+1);
                    break;
                    
                default:
                    trajsFinal[i].resize(trajs[i].size());
                    break;
            }

            for(size_t j=0; j<trajsFinal[i].size();j++)
            {
                double randomnumber = randomNumberGenerator(0,rownum-1);
                trajsFinal[i][j] = rows[randomnumber];
            }
        }
        
        trjColsBS[trjIndex].setTrajTables(trajsFinal);
    }
    
    return trjColsBS;
}*/


//no need to create rows, we just need index to resample
vector<TrajectoryCollection>
resample(const vector<TrajectoryCollection> & trajCols,
         ResampleType strategy)
{
    vector<TrajectoryCollection> trjColsBS = trajCols;
    
    //each row contains three number: condition number, trajectory number and time number
    vector<vector<int> > rowindices;
    vector<int> onerow(3);//store condition number, trajectory number and time number

    
    for(size_t condIndex = 0; condIndex < trajCols.size(); condIndex++)
    {
        const vector<vector<vector<int> > > & trajs = trajCols[condIndex].getTrajTables();
        
        for(size_t trjIndex=0; trjIndex<trajs.size();trjIndex++)
        {
            for(size_t timeIndex=0; timeIndex<trajs[trjIndex].size(); timeIndex++)
            {
                onerow[0] = condIndex;
                onerow[1] = trjIndex;
                onerow[2] = timeIndex;
                rowindices.push_back(onerow);
            }
        }
    }    
    

    for(size_t trjIndex = 0; trjIndex < trajCols.size(); trjIndex++)
    {
        const vector<vector<vector<int> > > & trajs = trajCols[trjIndex].getTrajTables();
                
        int rownum = rowindices.size();
        
        for(size_t i=0; i<trajs.size(); i++)
        {
            switch (strategy) {
                    
                case BOOTSTRAP_PLUS_ONE:
                    trjColsBS[trjIndex].
                    setTrajectory(i, vector< vector<int> >(trajs[i].size() + 1,
                                            vector<int>(trajs[i][0].size()) ) );
                    break;
                    
                default:
                    break;
            }
            
            for(size_t j=0; j<trjColsBS[trjIndex].getTrajTables()[i].size();j++)
            {
                double randomnumber = randomNumberGenerator(0,rownum-1);
                int tempCondIndx = rowindices[randomnumber][0];
                int tempTRJIndx = rowindices[randomnumber][1];
                int tempTimeIndx = rowindices[randomnumber][2];
                
                trjColsBS[trjIndex].setState(i, j,
                                             trajCols[tempCondIndx].getState(tempTRJIndx, tempTimeIndx) );

            }
        }

        /*
        switch (strategy) {
            case BOOTSTRAP_PLUS_ONE:
            {
                vector<int> trajLengths(trajs.size());
                for( size_t i=0; i<trajs.size(); i++ ) {
                    trajLengths[i] = trajs[i].size() + 1;
                }
                trjColsBS[trjIndex].setTrajectoryLength(trajLengths);
                break;
            }
            default:
                break;
        }
         */
    }

    return trjColsBS;
}