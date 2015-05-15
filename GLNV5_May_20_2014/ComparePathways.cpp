//
//  ComparePathways.cpp -- compare pathways by permutation (originally written by Yang Zhang)
//  
//
//  Extracted from main.cpp
//  Created by Joe Song on 9/30/11.
//  Copyright 2011 New Mexico State University. All rights reserved.
//
//  Obsolete, now using resampling based approximation rather than permutation.

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <sstream>

using namespace std;

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "GLN.h"
#include "GLNConfig.h"
#include "PathwayStats.h"
#include "EnvGLNCmp.h"
//#include "Topology.h"

EnvGLNCmp CompareGLNs(const GLNConfig & gc);

void parseResultFile(const string & resultFileLocation, vector<double> & Chisq, vector<int> & Df, vector<double> & Pvalue);

vector<vector<vector<int> > > concatenateTrjs(vector<vector<vector<int> > > & firstTrj, vector<vector<vector<int> > > & secondTrj);
vector<vector<vector<int> > > mergeTrjs(vector<vector<vector<int> > > & firstTrj, vector<vector<vector<int> > > & secondTrj);
vector<vector<vector<int> > > splitTrjsByRow(vector<vector<vector<int> > > & trjTable, int begin, int end);
vector<vector<vector<int> > > splitTrjsByCol(vector<vector<vector<int> > > & trjTable, int begin, int end);

void ComparePathways(GLNConfig & gc, int permuteStrategy, const vector<string> & input_files) 
{
    // srand( (unsigned)time(NULL) );
	
#ifdef ENABLE_MPI
    
    int rank =0;	
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    vector<TrajectoryCollection> m_trajCols;
    vector<vector<double> > Chisquares;  //chi_d, chi_c, chi_t, chi_p_z, chi_c_z
    vector<vector<int> > Degrees;  // df_d, df_c, df_t, df_p_z, df_c_z
    vector<vector<double> > Pvalues; // p_d, p_c, p_t, p_p_z, p_c_z
    
    string twzPEAResultFile;
    
    if(gc.m_pathways.size()==1)
    {
       twzPEAResultFile = gc.m_pathways[0].getFile();
    }else if(gc.m_candidateTopologys.size()==1)
    {
        twzPEAResultFile = gc.m_candidateTopologys[0].getFile();
    }
    twzPEAResultFile = twzPEAResultFile.substr(0,twzPEAResultFile.find(".txt"));
    twzPEAResultFile = twzPEAResultFile + "result.txt";
    
    for(int i=0; i<gc.m_numPermutations+1; i++)
    {
        if(i!=0) // do not permute the original trjs
        {
            if(rank == 0)
            {
                //scan trjs				
                m_trajCols.clear();
                m_trajCols.resize( gc.m_listTrajFiles.size() );
                for(size_t i=0; i<gc.m_listTrajFiles.size(); i++) {
                    m_trajCols[i].scan(gc.m_listTrajFiles[i]);
                }
                
                if(permuteStrategy==0)  // permute on time(row) across conditions 
                {
                    vector<vector<vector<int> > > concatenatedTrj = m_trajCols[0].getTrajTables() ;
                    for(size_t j=1; j<m_trajCols.size(); j++) {
                        vector<vector<vector<int> > > trjtable = m_trajCols[j].getTrajTables();
                        concatenatedTrj = concatenateTrjs(concatenatedTrj,trjtable);
                    }
                    TrajectoryCollection concatenateTrjCol = m_trajCols[0];
                    concatenateTrjCol.setTrajTables(concatenatedTrj);
                    concatenateTrjCol.permuteTrjtableOnTime();
                    for(size_t j=0; j<m_trajCols.size(); j++) {
                        vector<vector<vector<int> > > tempTrj = concatenateTrjCol.getTrajTables();
                        int trjLength = tempTrj.size()/m_trajCols.size();
                        m_trajCols[j].setTrajTables(splitTrjsByRow(tempTrj,trjLength*j,trjLength*(j+1)));
                        string fileName = "temp";
                        stringstream trjIndex;
                        trjIndex << j;
                        fileName = fileName + trjIndex.str() + ".trj";
                        m_trajCols[j].save(fileName.c_str());
                        gc.m_listTrajFiles[j] = fileName;	
                    }
                }else if(permuteStrategy == 1) // permute on genes(column)
                {
                    vector<vector<vector<int> > > concatenatedTrj = m_trajCols[0].getTrajTables() ;
                    for(size_t j=1; j<m_trajCols.size(); j++) {
                        vector<vector<vector<int> > > trjtable = m_trajCols[j].getTrajTables();
                        concatenatedTrj = mergeTrjs(concatenatedTrj,trjtable);
                    }
                    TrajectoryCollection concatenateTrjCol = m_trajCols[0];
                    concatenateTrjCol.setTrajTables(concatenatedTrj);
                    concatenateTrjCol.permuteTrjtableSwapingColumn();
                    for(size_t j=0; j<m_trajCols.size(); j++) {
                        vector<vector<vector<int> > > tempTrj = concatenateTrjCol.getTrajTables();
                        int trjWidth = tempTrj[0][0].size()/m_trajCols.size();
                        m_trajCols[j].setTrajTables(splitTrjsByCol(tempTrj,trjWidth*j,trjWidth*(j+1)));
                        string fileName = "temp";
                        stringstream trjIndex;
                        trjIndex << j;
                        fileName = fileName + trjIndex.str() + ".trj";
                        m_trajCols[j].save(fileName.c_str());
                        gc.m_listTrajFiles[j] = fileName;	
                    }
                }
                else if(permuteStrategy == 2) // permute on condition
                {
                    vector<vector<vector<int> > > concatenatedTrj = m_trajCols[0].getTrajTables() ;
                    for(size_t j=1; j<m_trajCols.size(); j++) {
                        vector<vector<vector<int> > > trjtable = m_trajCols[j].getTrajTables();
                        concatenatedTrj = mergeTrjs(concatenatedTrj,trjtable);
                    }
                    //currently supports 2 conditions only, need to change later on
                    for (int k = 0; k < m_trajCols[0].nNodes(); ++k) {
                        int m = rand() % 2;//m_trajCols.size();
                        //do the actual permute job 
                        if(m==1)
                        {
                            for(size_t trjIndex = 0; trjIndex < concatenatedTrj.size(); trjIndex++)
                            {
                                for(size_t rowIndex = 0; rowIndex < concatenatedTrj[trjIndex].size(); rowIndex++)
                                {
                                    int temp = concatenatedTrj[trjIndex][rowIndex][k];
                                    concatenatedTrj[trjIndex][rowIndex][k] = concatenatedTrj[trjIndex][rowIndex][k+m_trajCols[0].nNodes()];
                                    concatenatedTrj[trjIndex][rowIndex][k+m_trajCols[0].nNodes()] = temp;
                                }
                            }
                        }
                    }
                    
                    for(size_t j=0; j<m_trajCols.size(); j++) {
                        vector<vector<vector<int> > > tempTrj = concatenatedTrj;
                        int trjWidth = tempTrj[0][0].size()/m_trajCols.size();
                        m_trajCols[j].setTrajTables(splitTrjsByCol(tempTrj,trjWidth*j,trjWidth*(j+1)));
                        string fileName = "temp";
                        stringstream trjIndex;
                        trjIndex << j;
                        fileName = fileName + trjIndex.str() + ".trj";
                        m_trajCols[j].save(fileName.c_str());
                        gc.m_listTrajFiles[j] = fileName;	
                    }
                }
            }
            
            gc.m_listTrajFiles[0] = "temp0.trj";
            gc.m_listTrajFiles[1] = "temp1.trj";	
            
            //cout << gc.m_listTrajFiles[0] << endl;
        }//if(i!=0)
        MPI_Barrier (MPI_COMM_WORLD);
        
        CompareGLNs(gc);
        
        MPI_Barrier (MPI_COMM_WORLD);
        
        if(rank == 0)
        {
            vector<double> Chisq; //chi_d, chi_c, chi_t, chi_p_z, chi_c_z
            vector<int> Df; //df_d, df_c, df_t, df_p_z, df_c_z
            vector<double> Pvalue; //p_d,p_c,p_t,p_p_z,p_c_z
            Chisq.resize(5);
            Df.resize(5);
            Pvalue.resize(5);
            
            parseResultFile(twzPEAResultFile, Chisq, Df, Pvalue);
            
            Chisquares.push_back(Chisq);
            Degrees.push_back(Df);
            Pvalues.push_back(Pvalue);
        }
        gc.m_listTrajFiles = input_files;
    }
    
    if(rank == 0)
    {
        ofstream out(twzPEAResultFile.c_str());
        out<<"";
        out.close();
        
        double rank_p_d = 0.0; 
        double rank_p_c = 0.0; 
        double rank_p_t = 0.0; 
        double rank_p_p_z = 0.0; 
        double rank_p_c_z = 0.0; 
        for ( size_t i=1; i<Pvalues.size(); i++)
        {
            if(Pvalues[i][0]<Pvalues[0][0])
            {
                rank_p_d += 1;
            }
            if(Pvalues[i][1]<Pvalues[0][1])
            {
                rank_p_c += 1;
            }
            if(Pvalues[i][2]<Pvalues[0][2])
            {
                rank_p_t += 1;
            }
            if(Pvalues[i][3]<Pvalues[0][3])
            {
                rank_p_p_z += 1;
            }
            if(Pvalues[i][4]<Pvalues[0][4])
            {
                rank_p_c_z += 1;
            }
        }
        
        MultiplePathwayStats mps;
        
        mps.m_het = Chisq(Chisquares[0][0], Degrees[0][0], rank_p_d/gc.m_numPermutations);
        mps.m_hom = Chisq(Chisquares[0][1], Degrees[0][1], rank_p_c/gc.m_numPermutations);
        mps.m_tot = Chisq(Chisquares[0][2], Degrees[0][2], rank_p_t/gc.m_numPermutations);
        mps.m_wzParents = Chisq(Chisquares[0][3], Degrees[0][3], rank_p_p_z/gc.m_numPermutations);
        mps.m_wzChildren = Chisq(Chisquares[0][4], Degrees[0][4], rank_p_c_z/gc.m_numPermutations);
        mps.save(twzPEAResultFile);
    }
    
#else	
    
    vector<TrajectoryCollection> m_trajCols;
    vector<vector<double> > Chisquares;  //chi_d, chi_c, chi_t, chi_p_z, chi_c_z
    vector<vector<int> > Degrees;  // df_d, df_c, df_t, df_p_z, df_c_z
    vector<vector<double> > Pvalues; // p_d, p_c, p_t, p_p_z, p_c_z
    
	string base;
    if(gc.m_pathways.size()==1)
	{
		base = gc.m_pathways[0].getFile().substr(0, gc.m_pathways[0].getFile().find(".txt"));
	}else if(gc.m_candidateTopologys.size()==1){
		base = gc.m_candidateTopologys[0].getFile().substr(0, gc.m_candidateTopologys[0].getFile().find(".txt"));
	}
    
    string twzPEAResultFile = base + "-PWStats.txt";
     
    for(int i=0; i<gc.m_numPermutations+1; i++)
    {
        if(i!=0) // do not permute the original trjs
        {
            //scan trjs				
            m_trajCols.clear();
            m_trajCols.resize( gc.m_listTrajFiles.size() );
            
            for(size_t i=0; i<gc.m_listTrajFiles.size(); i++) {
                m_trajCols[i].scan(gc.m_listTrajFiles[i]);
            }
            
            if(permuteStrategy==0)  // permute on time(row) across conditions 
            {
                vector<vector<vector<int> > > concatenatedTrj = m_trajCols[0].getTrajTables() ;
                for(size_t j=1; j<m_trajCols.size(); j++) {
                    vector<vector<vector<int> > > trjtable = m_trajCols[j].getTrajTables();
                    concatenatedTrj = concatenateTrjs(concatenatedTrj,trjtable);
                }
                TrajectoryCollection concatenateTrjCol = m_trajCols[0];
                concatenateTrjCol.setTrajTables(concatenatedTrj);
                concatenateTrjCol.permuteTrjtableOnTime();
                for(size_t j=0; j<m_trajCols.size(); j++) {
                    vector<vector<vector<int> > > tempTrj = concatenateTrjCol.getTrajTables();
                    int trjLength = tempTrj.size()/m_trajCols.size();
                    m_trajCols[j].setTrajTables(splitTrjsByRow(tempTrj,trjLength*j,trjLength*(j+1)));
                    string fileName = "temp";
                    stringstream trjIndex;
                    trjIndex << j;
                    fileName = fileName + trjIndex.str() + ".trj";
                    m_trajCols[j].save(fileName.c_str());
                    gc.m_listTrajFiles[j] = fileName;	
                }
            }else if(permuteStrategy == 1) // permute on genes(column)
            {
                vector<vector<vector<int> > > concatenatedTrj = m_trajCols[0].getTrajTables() ;
                for(size_t j=1; j<m_trajCols.size(); j++) {
                    vector<vector<vector<int> > > trjtable = m_trajCols[j].getTrajTables();
                    concatenatedTrj = mergeTrjs(concatenatedTrj,trjtable);
                }
                TrajectoryCollection concatenateTrjCol = m_trajCols[0];
                concatenateTrjCol.setTrajTables(concatenatedTrj);
                concatenateTrjCol.permuteTrjtableSwapingColumn();
                for(size_t j=0; j<m_trajCols.size(); j++) {
                    vector<vector<vector<int> > > tempTrj = concatenateTrjCol.getTrajTables();
                    int trjWidth = tempTrj[0][0].size()/m_trajCols.size();
                    m_trajCols[j].setTrajTables(splitTrjsByCol(tempTrj,trjWidth*j,trjWidth*(j+1)));
                    string fileName = "temp";
                    stringstream trjIndex;
                    trjIndex << j;
                    fileName = fileName + trjIndex.str() + ".trj";
                    m_trajCols[j].save(fileName.c_str());
                    gc.m_listTrajFiles[j] = fileName;	
                }
            }
            else if(permuteStrategy == 2) // permute on condition
            {
                vector<vector<vector<int> > > concatenatedTrj = m_trajCols[0].getTrajTables() ;
                for(size_t j=1; j<m_trajCols.size(); j++) {
                    vector<vector<vector<int> > > trjtable = m_trajCols[j].getTrajTables();
                    concatenatedTrj = mergeTrjs(concatenatedTrj,trjtable);
                }
                //currently supports 2 conditions only, need to change later on
                for (int k = 0; k < m_trajCols[0].nNodes(); ++k) {
                    int m = rand() % 2;//m_trajCols.size();
                    //do the actual permute job 
                    if(m==1)
                    {
                        for(size_t trjIndex = 0; trjIndex < concatenatedTrj.size(); trjIndex++)
                        {
                            for(size_t rowIndex = 0; rowIndex < concatenatedTrj[trjIndex].size(); rowIndex++)
                            {
                                int temp = concatenatedTrj[trjIndex][rowIndex][k];
                                concatenatedTrj[trjIndex][rowIndex][k] = concatenatedTrj[trjIndex][rowIndex][k+m_trajCols[0].nNodes()];
                                concatenatedTrj[trjIndex][rowIndex][k+m_trajCols[0].nNodes()] = temp;
                            }
                        }
                    }
                }
                
                for(size_t j=0; j<m_trajCols.size(); j++) {
                    vector<vector<vector<int> > > tempTrj = concatenatedTrj;
                    int trjWidth = tempTrj[0][0].size()/m_trajCols.size();
                    m_trajCols[j].setTrajTables(splitTrjsByCol(tempTrj,trjWidth*j,trjWidth*(j+1)));
                    string fileName = "temp";
                    stringstream trjIndex;
                    trjIndex << j;
                    fileName = fileName + trjIndex.str() + ".trj";
                    m_trajCols[j].save(fileName.c_str());
                    gc.m_listTrajFiles[j] = fileName;	
                }
            }
            
            gc.m_listTrajFiles[0] = "temp0.trj";
            gc.m_listTrajFiles[1] = "temp1.trj";	
            
            //cout << gc.m_listTrajFiles[0] << endl;
        }//if(i!=0)
        
        CompareGLNs(gc);
        
        vector<double> Chisq; //chi_d, chi_c, chi_t, chi_p_z, chi_c_z
        vector<int> Df; //df_d, df_c, df_t, df_p_z, df_c_z
        vector<double> Pvalue; //p_d,p_c,p_t,p_p_z,p_c_z
        Chisq.resize(5);
        Df.resize(5);
        Pvalue.resize(5);
        
        parseResultFile(twzPEAResultFile, Chisq, Df, Pvalue);
        
        Chisquares.push_back(Chisq);
        Degrees.push_back(Df);
        Pvalues.push_back(Pvalue);
        
        gc.m_listTrajFiles = input_files;
    }
    
    ofstream out(twzPEAResultFile.c_str());
    out<<"";
    out.close();
    
    double rank_p_d = 0.0; 
    double rank_p_c = 0.0; 
    double rank_p_t = 0.0; 
    double rank_p_p_z = 0.0; 
    double rank_p_c_z = 0.0; 
    for ( size_t i=1; i<Pvalues.size(); i++)
    {
        if(Pvalues[i][0]<Pvalues[0][0])
        {
            rank_p_d += 1;
        }
        if(Pvalues[i][1]<Pvalues[0][1])
        {
            rank_p_c += 1;
        }
        if(Pvalues[i][2]<Pvalues[0][2])
        {
            rank_p_t += 1;
        }
        if(Pvalues[i][3]<Pvalues[0][3])
        {
            rank_p_p_z += 1;
        }
        if(Pvalues[i][4]<Pvalues[0][4])
        {
            rank_p_c_z += 1;
        }
    }
    
    MultiplePathwayStats mps;
    
    mps.m_het = Chisq(Chisquares[0][0], Degrees[0][0], rank_p_d/gc.m_numPermutations);
    mps.m_hom = Chisq(Chisquares[0][1], Degrees[0][1], rank_p_c/gc.m_numPermutations);
    mps.m_tot = Chisq(Chisquares[0][2], Degrees[0][2], rank_p_t/gc.m_numPermutations);
    mps.m_wzParents = Chisq(Chisquares[0][3], Degrees[0][3], rank_p_p_z/gc.m_numPermutations);
    mps.m_wzChildren = Chisq(Chisquares[0][4], Degrees[0][4], rank_p_c_z/gc.m_numPermutations);
    mps.save(twzPEAResultFile);
    
#endif	
    
}

