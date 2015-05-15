// JobGLN.cpp
// 
// Curtis Luce, Joe Song
// Created: June 2007
// Modified: September 5, 2008. Joe Song. Incorporated minimum number of parents
// Last modified: September 27, 2008.  Joe Song. Name changed from JobLognet.cpp
#ifdef ENABLE_MPI			
#include <mpi.h>
#endif

#include <new>
#include <iostream>
#include <fstream>
#include <sstream>
using std::ios;
using std::bad_alloc;
using std::endl;
using std::cerr;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ofstream;

#include <cstdlib>

#include "JobGLN.h"
#include "GLNGlobals.h"

void removeDuplicates(std::vector<int>& vec);

int JobGLN::pack(unsigned char* &jobBuffer) const 
{
	int elementsSize = 1 + 1 + getSize() + getSize();

	int totalSize = sizeof(int) * elementsSize;
	
	try {

		jobBuffer = new unsigned char [totalSize];

	} catch (bad_alloc & ba) {

		cerr << "ERROR: Failed to allocate a buffer for job!" << endl;
		cerr << "\t" << ba.what() << endl;
	}

	unsigned char * p = jobBuffer;

	*(int *) p = m_nodeId;
	p += sizeof(int);
	
	*(int *) p = m_size;
	p += sizeof (int);
	
	for (int i = 0; i < m_size; i++){
		*(int *) p = m_beginParentList[i];
		p += sizeof(int);
	}

	for (int i = 0; i < m_size; i++){
		*(int *) p = m_endParentList[i];
		p += sizeof(int);
	}
	
	return totalSize;
}

void JobGLN::unpack(unsigned char * jobBuffer, int length)
{
	unsigned char * p = 0;
	p = jobBuffer;
	m_nodeId = *(int *) p;
	p += sizeof (int);
	m_size = *(int *) p;
	p += sizeof (int);

	m_beginParentList.resize(m_size);
	m_endParentList.resize(m_size);

	for (int i = 0; i < m_size; i++){
		m_beginParentList[i] = *(int *) p;
		p += sizeof(int);
	}
	
	for (int i = 0; i < m_size; i++){
		m_endParentList[i] = *(int *) p;
		p += sizeof(int);
	}
}

void JobGLN::first(const JobEnvironment &env)
{
	int nNodes = env.getGLNSize();

	//int nParents;
    int numOfParentsAllowed;

    //With the given 'min_parents_allowed' being the size initially,
    //'beginParents' records the starting set (parent IDs)
    //and 'endParents' records the last possbile parent set (parent IDs).

    //For example, in a network with 5 nodes and minimum number of parents being 3
    // 'beginparents' will be "1, 2, 3" and
    // 'endParents" will be "3, 4, 5".
    // JobGLN::preocess() will iterate through combinations "1,2,3", "1,2,4", "1,2,5", ... "2,4,5", "3,4,5"
    // with allowed delays to find out the combinatioin of parents with largest p-value.
	vector<int> beginParents;
	vector<int> endParents;

	int child = env.m_gc.m_startAtNode - 1; // index from 0


    //If all nodes are external nodes, i.e. nodes that can not have parents, return immediately
    /*************************************************************/
	while ( child < nNodes && env.getNodeType(child) == 'e' ) {
		++child;
	}
	if ( child == nNodes ) {
		terminal();
		return;
	}
	/***********************************************************/



	// nParents = 1; 
	//nParents = env.m_gc.m_min_parents_allowed; // Joe Song, September 5, 2008
    numOfParentsAllowed = env.m_gc.m_min_parents_allowed;
	
    
    //beginParents.resize(nParents);
	//endParents.resize(nParents);
    beginParents.resize(numOfParentsAllowed);
	endParents.resize(numOfParentsAllowed);

	// Generate initial set of parents
	//for(int i=0; i < nParents; i++){
    for(int i=0; i<numOfParentsAllowed; ++i){
		beginParents[i] = i + 1;  //Parents ID starts from 1
	}//end for

	// endParents[0] = penv -> m_gln.size();
	// endParents[0] = env.getGLNSize();
	//for(int i=0; i < nParents; i++){ // Joe Song November 30, 2008
    int indexDiff = env.getGLNSize() - numOfParentsAllowed ;
	for(int i=0; i<numOfParentsAllowed; ++i){	
        endParents[i] = indexDiff + i + 1;
	} //end for

	initialize(child+1, beginParents, endParents, numOfParentsAllowed);
	
}

void JobGLN::next(const JobEnvironment &env)
{
    //For one node, if already test the maximum number of parents allowed, then proceed to next node.
    //Otherwise, increase the number of parents and JobGLN::process() will iterate through all the combinations.

    //For example, if maximum number of parents allowed is 4.
    //And JobGLN::process() has already tested all parent combinations of size 3.
    //JobGLN::next() will generate 'beginParentlist' to be "1,2,3,4" and
    // 'endParentList" to be "2,3,4,5" (assuming network size is 5).


	if( isTerminal() )  return;

	//int child = getNodeId() - 1;
    int child = m_nodeId - 1;

	//size_t nParents = getSize();
    // nParents: number of parents allowed for this iternation
    size_t nParents = m_size;

	int nNodes = env.getGLNSize();

	//vector<int> beginParentList = getBeginParents();
	//vector<int> endParentList = getEndParents();
    vector<int> beginParentList = m_beginParentList;
    vector<int> endParentList = m_endParentList;

    if((nParents == env.getMaxNumDiffParents()) || (child >= nNodes)) { // Joe Song Nov 18, 2011
        
	// if((nParents == env.m_gc.m_max_parents_allowed) || (child >= nNodes)) {

		++child; //Move to next child node

		if( child >= nNodes ) {
			terminal();
			return;
		}

		while ( child < nNodes && env.getNodeType(child) == 'e' ) {
			++child;
		}
		if ( child == nNodes ) {
			terminal();
			return;
		}

		// nParents = 1;

        //Reset the number of parents to the minimum and increase it up tp the upper limit to find the set of parents with largerst p-value.
		nParents = env.m_gc.m_min_parents_allowed; // Joe Song. September 5, 2005.

	} else { //For the same child, increase the possible number of parents, and calculate all combinations p-value again.
		++nParents;
	}

	if(! isTerminal()) {

		beginParentList.resize(nParents);
		endParentList.resize(nParents);

		/**************************************
		****Generate initial set of parents****
		***************************************/
		for(size_t i=0; i < nParents; i++){
			beginParentList[i] = i+1;
		}//end for

		if(nParents == 1) {
			endParentList[0] = nNodes;
		} else {
			for(int i=nParents-1, d=0; i>=0; --i,++d){
				endParentList[i] = nNodes - d;
			}//end for
		}//end else
	}

	initialize(child+1, beginParentList, endParentList, nParents);
}

/*bool JobGLN::process(JobEnvironment &env) const
{
	bool improved = false;
	if(env.m_gc.m_pathway.getPathway().size()<=0)
	{
		improved = env.enumerateParents(getNodeId()-1, getBeginParents(), 
		getEndParents(), env.getGLNSize(), env.m_gc);
		
		//added by YangZhang 2/24/2009 to record all results with p-value less than threshhold
		if(env.m_gc.m_recordResultFile.length()!=0)
		{
			ResutlsToFile(env);
		}
 	}
	else
	{
		vector<int> parentlist = (env.m_gc.m_pathway.getPathway())[getNodeId()-1];	
		improved = env.enumerateParents(getNodeId()-1, parentlist,parentlist, env.getGLNSize(), env.m_gc);
	}
	return improved;
}*/

bool JobGLN::process(JobEnvironment &env) const
{
    bool improved = false;
    
    //int child = getNodeId()-1;  //Which node, i.e. child, to work on.
    int child = m_nodeId - 1;  // 'child' is 0-based

    //int nParents = getBeginParents().size();
    // nParents: number of parents allowed for this iternation
    int nParents = m_beginParentList.size();

    //Record IDs of nodes
    vector<int> candidates;
    vector<int> begin;
    vector<int> end;
    
    int nNodes = 0;  // Number of nodes to consider in the enumeration
    
    if( env.m_gc.m_pathways.size() ==1) { // only provide one pathway file, shared by both condition
        
        begin = end = (env.m_gc.m_pathways[0].get())[child][0];
    
        nNodes = env.getGLNSize();
        
    } else if( env.m_gc.m_candidateTopologys.size() == 1) {
        
        //begin = getBeginParents();
        begin = m_beginParentList;

        candidates = (env.m_gc.m_candidateTopologys[0].get())[child][0];
        
        nNodes = candidates.size();
        
        if( nNodes >= nParents) {
            
            end.resize(nParents);
            
            for(int i=0; i < nParents; i++){
                end[i] = (candidates.size() - nParents) + i + 1;
            }
        }
        
    } else if (env.m_gc.m_pathways.size() > 1 || env.m_gc.m_candidateTopologys.size() > 1 ) {
        //individual topology for each condition
        //treat them as if candidate topology file although we should use exactly the topology provided
        begin = m_beginParentList;
        
        //env.m_gc.m_max_parents_allowed = 2*(env.m_gc.m_max_parents_allowed);
        
        if(env.m_gc.m_pathways.size() > 1)
        {
            for(size_t i=0; i<env.m_gc.m_pathways.size(); i++ )
            {
                vector<int> tempParents = (env.m_gc.m_pathways[i].get())[child][0];
                if(tempParents.size()>0)
                {
                    for(size_t pid=0; pid<tempParents.size(); pid++)
                    {
                        candidates.push_back(tempParents[pid]);
                    }
                }
            }
            
        }else if(env.m_gc.m_candidateTopologys.size()>1)
        {
            for(size_t i=0; i<env.m_gc.m_candidateTopologys.size(); i++ )
            {
                vector<int> tempParents = (env.m_gc.m_candidateTopologys[i].get())[child][0];
                if(tempParents.size()>0)
                {
                    for(size_t pid=0; pid<tempParents.size(); pid++)
                    {
                        candidates.push_back(tempParents[pid]);
                    }
                }
            }
        }
        
        removeDuplicates(candidates);

        nNodes = candidates.size();
        
        if( nNodes >= nParents) {
            
            end.resize(nParents);
            
            for(int i=0; i < nParents; i++){
                end[i] = (candidates.size() - nParents) + i + 1;
            }
        }
    } else {

        //begin = getBeginParents();
        begin = m_beginParentList;

        if(env.m_gc.m_HP[child]) {
            for(size_t i=0; i<env.getGLNSize(); ++i) {
                if(env.m_gc.m_BP[i]) {
                    candidates.push_back(i+1);
                }
            }

            nNodes = candidates.size();
            if(nNodes >= nParents) {
                
                end.resize(nParents);
                int indexDiff = nNodes - nParents;

                for(int i=0; i<nParents; ++i) {
                    end[i] = indexDiff + i + 1;
                }
            }
        }
    }//end of else
    

    // 'mustIncludeParents' records all the parents must be include from child 'child'
    // Note that 'child' is 0-based and parents ID in 'mustIncludeParents' are 1-based.
    vector<int> mustIncludeParents(0);
    if( env.m_gc.m_mustIncludeParents.size() > 0)
        mustIncludeParents = env.m_gc.m_mustIncludeParents.get()[child][0];



    if(nNodes >= nParents) {
        improved = env.enumerateParentSubsets(child, begin, end, nNodes, env.m_gc, candidates, mustIncludeParents);
    }
    


    if(env.m_gc.m_recordResultFile.length()!=0) {
        //added by YangZhang 2/24/2009 to record all results with p-value less than threshhold
        ResutlsToFile(env);
    }



	return improved;
}


void JobGLN::ResutlsToFile(JobEnvironment &env) const
{
#ifdef ENABLE_MPI

    //m_resultrecord.push_back(_result);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //if(rank == 0)
    //{
    stringstream out;  //child node id
    //convert int to string
    out << rank;
    
	string file = env.m_gc.m_recordResultFile + out.str() + "SummarizedResults.txt";

#else  // #ifndef ENABLE_MPI

    string file = env.m_gc.m_recordResultFile + "SummarizedResults.txt";

#endif

    ofstream ofs(file.c_str(), ios::app);

    if(! ofs.is_open() ) {
        cerr << "ERROR: openning file \"" << env.m_gc.m_recordResultFile << "\"" << endl;
        GLNExit(EXIT_FAILURE);
    }

    for(size_t i=0; i<env.m_recordresult[getNodeId()-1].size(); i++)
    {
        ofs << env.m_recordresult[getNodeId()-1][i] << endl;
    }

	ofs.close();

    //added by Yang Zhang 11.09.2011 remove previous results to avoid duplicate results
    env.m_recordresult[getNodeId()-1].clear();
    
// #endif
    
    //**************************************************************
    //** Record contingency tables found in m_recordresults
    //** only the counts are recorded and placed in to separate
    //** files for easy of extraction.  NOTE: that the number in
    //** the filename correlates to the order in which a result is
    //** listed in m_recordresult (CMHJ 4/18/09)
    /*FILE * continTableFile;
     
     for(size_t i=0;i<env.m_recordcounts[getNodeId()-1].size();i++)
     {
     stringstream continTableFilename;
     continTableFilename << env.m_gc.m_recordResultFile;
     continTableFilename << "ContinTable_NodeID";
     continTableFilename << getNodeId()-1;
     continTableFilename << "_";
     
     stringstream lookingforname;
     size_t increamentor = 0;
     bool done = false;
     do{
     lookingforname.str("");
     lookingforname << continTableFilename.str();
     
     lookingforname << i + increamentor;
     lookingforname << ".txt";
     continTableFile = fopen(lookingforname.str().c_str(), "r");
     if(continTableFile)
     {
     fclose(continTableFile);
     increamentor++;
     }
     else
     {
     continTableFile = fopen(lookingforname.str().c_str(), "w");
     done = true;
     }
     }while(!done);
     
     if(continTableFile==NULL)
     {
     cerr << "ERROR: openning file \"" << lookingforname.str().c_str() << "\"" << endl;
     GLNExit(EXIT_FAILURE);
     }
     fprintf(recordResultFile,"# %s\n",env.m_recordresult[getNodeId()-1][i].c_str());
     fprintf(continTableFile,"%s\n",env.m_recordcounts[getNodeId()-1][i].c_str());
     fclose(continTableFile);
     }*/
    //**************************************************************
}


bool JobGLN::processToBeDeleted(JobEnvironment &env) const
{
	bool improved = false;
	if(env.m_gc.m_acceptableTopologies == SAME) // 1)
	{
		/* modify enumerating parent strategy by Yang Zhang 7.25.2010
         if(env.m_gc.m_pathway.getPathway().size()<=0)
         {
         improved = env.enumerateParents(getNodeId()-1, getBeginParents(), 
         getEndParents(), env.getGLNSize(), env.m_gc);
         
         //added by YangZhang 2/24/2009 to record all results with p-value less than threshhold
         if(env.m_gc.m_recordResultFile.length()!=0)
         {
         ResutlsToFile(env);
         }
         }
         else
         {
         vector<int> parentlist = (env.m_gc.m_pathway.getPathway())[getNodeId()-1];	
         improved = env.enumerateParents(getNodeId()-1, parentlist,parentlist, env.getGLNSize(), env.m_gc);
         }
         */
		if(env.m_gc.m_pathways.size() == 1)
		{
			vector<int> parentlist = (env.m_gc.m_pathways[0].get())[getNodeId()-1][0];	
			improved = env.
            // enumerateParentSubsets
            enumerateParents
            (getNodeId()-1, parentlist,parentlist, env.getGLNSize(), env.m_gc, vector<int>());
            
		}
		else if(env.m_gc.m_candidateTopologys.size()==1)
		{
			vector<int> parentlist = (env.m_gc.m_candidateTopologys[0].get())[getNodeId()-1][0];
			int NumbOfParent = getBeginParents().size();
			if( (int) parentlist.size() >= NumbOfParent)
			{
				vector<int> endParentsInParentlist;
				endParentsInParentlist.resize(NumbOfParent);
				for(int i=0; i < NumbOfParent; i++){
					endParentsInParentlist[i] = (parentlist.size() - NumbOfParent) + i + 1;
				}	
				improved = env.
                //enumerateParentSubsets
                enumerateParents
                (getNodeId()-1, getBeginParents(), endParentsInParentlist, parentlist.size(), env.m_gc, parentlist);
			}
		}
		else 
		{
			vector<int> parentlist;
            
			if(env.m_gc.m_HP[getNodeId()-1])
			{
				for(size_t i=0; i<env.getGLNSize(); i++)
				{
					if(env.m_gc.m_BP[i])
					{
						parentlist.push_back(i+1);
					}
				}
                
				int NumbOfParent = getBeginParents().size();
				if((int) parentlist.size()>=NumbOfParent)
				{
					vector<int> endParentsInParentlist;
					endParentsInParentlist.resize(NumbOfParent);
					for(int i=0; i < NumbOfParent; i++){
						endParentsInParentlist[i] = (parentlist.size() - NumbOfParent) + i + 1;
					}	
					improved = env.
                    // enumerateParentSubsets
                    enumerateParents
                    (getNodeId()-1, getBeginParents(),endParentsInParentlist, parentlist.size(), env.m_gc, parentlist);
				}
			}
			/*improved = env.enumerateParents(getNodeId()-1, getBeginParents(), 
             getEndParents(), env.getGLNSize(), env.m_gc, vector<int>());*/
            
		}
		//added by YangZhang 2/24/2009 to record all results with p-value less than threshhold
		if(env.m_gc.m_recordResultFile.length()!=0)
		{
			ResutlsToFile(env);
		}
	} else if(env.m_gc.m_acceptableTopologies == DIFFERENT) // 2)
	{
		/*modify enumerating parent strategy by Yang Zhang 7.25.2010
         if(env.m_gc.m_pathway.getPathway().size()<=0)
         {
         improved = env.enumerateParents2(getNodeId()-1, getBeginParents(), 
         getEndParents(), env.getGLNSize(), env.m_gc);
         
         //added by YangZhang 2/24/2009 to record all results with p-value less than threshhold
         if(env.m_gc.m_recordResultFile.length()!=0)
         {
         ResutlsToFile(env);
         }
         }
         else
         {
         vector<int> parentlist = (env.m_gc.m_pathway.getPathway())[getNodeId()-1];	
         improved = env.enumerateParents2(getNodeId()-1, parentlist,parentlist, env.getGLNSize(), env.m_gc);
         }
         */
		if(env.m_gc.m_pathways.size() == 1)
		{
			vector<int> parentlist = (env.m_gc.m_pathways[0].get())[getNodeId()-1][0];	
			improved = env.
            //enumerateParentSubsets
            enumerateParents2
            (getNodeId()-1, parentlist,parentlist, env.getGLNSize(), env.m_gc, vector<int>());
		}
		else if(env.m_gc.m_candidateTopologys.size()==1)
		{
			vector<int> parentlist = (env.m_gc.m_candidateTopologys[0].get())[getNodeId()-1][0];	
			//vector<int> parentlistMask;
            
			int NumbOfParent = getBeginParents().size();
			if((int) parentlist.size() >= NumbOfParent)       
			{
				vector<int> endParentsInParentlist;	
				endParentsInParentlist.resize(NumbOfParent);
				for(int i=0; i < NumbOfParent; i++){
					endParentsInParentlist[i] = (parentlist.size() - NumbOfParent) + i + 1;
				}
				//improved = env.enumerateParents2(getNodeId()-1, getBeginParents(), 
                //getEndParents(), env.getGLNSize(), env.m_gc, vector<int>());
				improved = env.
                // enumerateParentSubsets
                enumerateParents2
                (getNodeId()-1, getBeginParents(),endParentsInParentlist, parentlist.size(), env.m_gc, parentlist);
			}
		}
		else
		{
			vector<int> parentlist;
            
			if(env.m_gc.m_HP[getNodeId()-1])
			{
				for(size_t i=0; i<env.getGLNSize(); i++)
				{
					if(env.m_gc.m_BP[i])
					{
						parentlist.push_back(i+1);
					}
				}
                
				int NumbOfParent = getBeginParents().size();
				if((int) parentlist.size()>=NumbOfParent)
				{
					vector<int> endParentsInParentlist;
					endParentsInParentlist.resize(NumbOfParent);
					for(int i=0; i < NumbOfParent; i++){
						endParentsInParentlist[i] = (parentlist.size() - NumbOfParent) + i + 1;
					}	
					improved = env.
                    // enumerateParentSubsets
                    enumerateParents2
                    (getNodeId()-1, getBeginParents(),endParentsInParentlist, parentlist.size(), env.m_gc, parentlist);
				}
			}
            
			/*improved = env.enumerateParents2(getNodeId()-1, getBeginParents(), 
             getEndParents(), env.getGLNSize(), env.m_gc, vector<int>());*/
		}
        
		//added by YangZhang 2/24/2009 to record all results with p-value less than threshhold
		if(env.m_gc.m_recordResultFile.length()!=0)
		{
			ResutlsToFile(env);
		}
        
	}
	return improved;
}

/* 
void JobGLN::generateNextJob(JobEnvironment &env)
{
	JobEnvironmentGLN *penv = (JobEnvironmentGLN *) &env;

	/ *
	if(m_beginParentList) {
		delete [] m_beginParentList;
		m_beginParentList = 0;
	}

	if(m_endParentList) {
		delete [] m_endParentList;
		m_endParentList = 0;
	}
	* /

	// ****************************cout << "generating job" << endl;
	if((penv -> m_nParents == penv -> m_estPar.m_max_parents_allowed) || (penv -> m_child >= penv -> m_gln.size())) {

		penv -> m_child++;
		m_nodeId++;
		size = penv -> m_nParents;

		// **********************cout << "m_child: " << penv -> m_child << endl;
		// ************************cout << "nodes" << penv -> m_gln.size() << endl;
		// ***************************cout << "m_nParents" << penv -> m_nParents << endl;

		if(penv -> m_child >= penv -> m_gln.size()){

			penv -> indexParents.resize(penv -> m_nParents);

			/ *
			try{
				m_endParentList = new int [penv -> m_nParents];
			}
			catch(bad_alloc xa){
				cout << "Allocation error for array 'm_beginParentList' in function 'generateNextJob'" << endl; // Joe Song
			}//end catch
			* /
			m_beginParentList.resize(penv -> m_nParents);
			m_endParentList.resize(penv -> m_nParents);

			for(int i = 0; i < penv -> m_nParents; i++){
				m_beginParentList[i] = -1;
				m_endParentList[i] = -1;
				penv -> indexParents [i] = -1;
			}
			size = penv -> m_nParents;
			return;
		}
		//int flag2 = 0;
		while (true){
			if(penv -> m_child < penv -> m_gln.size()){
				if(penv -> m_gln.getNodes()[penv -> m_child].getType() == 'e'){
					penv -> m_child++;
					m_nodeId++;
				} else {
					break;
				}
			} else {
				penv -> indexParents.resize(penv -> m_nParents);
				penv -> m_nParents = penv -> m_estPar.m_max_parents_allowed;
				
				/ *
				try{
					m_beginParentList = new int [penv -> m_nParents];
				}
				catch(bad_alloc xa){
					cout << "Allocation error for array 'm_beginParentList' in function 'generateNextJob'" << endl; // Joe Song
				}//end catch
				* /
				m_beginParentList.resize(penv -> m_nParents);

				/ *
				try{
					m_endParentList = new int [penv -> m_nParents];
				}
				catch(bad_alloc xa){
					cout << "Allocation error for array 'm_beginParentList' in function 'generateNextJob'" << endl; // Joe Song
				}//end catch
				* /
				m_endParentList.resize(penv -> m_nParents);

				for(int i = 0; i < size; i++){
					m_beginParentList[i] = -1;
					m_endParentList[i] = -1;
					penv -> indexParents [i] = -1;
				}
				size = penv -> m_nParents;
				return;
				//flag2 =1;
				//break;
			}
		}
		//if(flag2 == 0 ){
		size = penv -> m_nParents;  // ?? this -> size = 1;   

		// penv -> m_nParents = 1;
		penv -> m_nParents = penv -> m_estPar.m_min_parents_allowed; // Joe Song. September 5, 2005.

		//}
	} else {
		penv -> m_nParents++;
		size = penv -> m_nParents;
	}

	// **********************cout << "m_child in generate next job    " << penv -> m_child << endl;
	// **********************cout << "Testing node type: " << penv -> m_gln.getNodes()[penv -> m_child].getType()  << endl;
	if(isTerminal(env) == false){
		penv -> position = penv -> m_nParents-1;			//Adjust position to end of array
		penv -> first_max = penv -> m_gln.size() - (penv -> m_nParents-1);    	//maximum number the array[0] can be
		penv -> trans_table=0;			//reset pointer to zero
		penv -> delete_array = true;		// flag for if we should delete the array at the very moment
		penv -> tt_rows = 1;			// truth table rows, reset to 1 to figure out the real tt_rows
		penv -> tt_cols =  penv -> m_trajCol.base(penv -> m_child);  	//number of truth table columns is based on the current child's base

		penv -> indexParents.resize(penv -> m_nParents);

		/ *
		m_beginParentList = 0;
		m_endParentList = 0;

		try{	
			m_beginParentList = new int[penv -> m_nParents];
		}//end try
		catch(bad_alloc xa){
			cout << "Allocation error for array 'm_beginParentList' in function 'generateNextJob'" << endl; // Joe Song
		}//end catch

		try{	
			m_endParentList = new int[penv -> m_nParents];
		}//end try
		catch(bad_alloc xa){
			cout << "Allocation error for array 'm_endParentList' in function 'generateNextJob'" << endl; // Joe Song
		}//end catch
		* /

		m_beginParentList.resize(penv -> m_nParents);
		m_endParentList.resize(penv -> m_nParents);

		// **************************************
		// ****Generate initial set of parents****
		// ***************************************
		for(int i=0; i < penv -> m_nParents; i++){
			penv -> indexParents[i]= i+1;
			m_beginParentList[i] = i+1;
		}//end for

		if(penv -> m_nParents == 1)
			m_endParentList[0] = penv -> m_gln.size();
		else{
			int d = 0;
			//this -> m_endParentList[penv -> m_nParents-1] = penv -> m_gln.size();
			for(int i = penv -> m_nParents-1; i>=0; i--){
				m_endParentList[i] = penv -> m_gln.size()-d;
				d++;
			}//end for
		}//end else

		penv -> tt_cols =  penv -> m_trajCol.base( penv -> m_child );	//Readjust tt_cols
		penv -> tt_rows=1;

		for(int i=0; i < penv -> m_nParents; i++){
			penv -> tt_rows = penv -> m_trajCol.base( penv -> indexParents[i]-1 ) * penv -> tt_rows;
		}//end for

		//penv -> trans_table = createTable(penv -> tt_rows, penv -> tt_cols);	//create a new truth table
		penv -> delete_array = true;
		size = penv -> m_nParents;
		// **********************cout << "m_child in generate next job  " << penv -> m_nParents << endl;
	}
}
*/

/*
JobResult* JobGLN::processJob(JobEnvironment &env)
{
	JobEnvironmentGLN *penv = (JobEnvironmentGLN*) &env;

	JobResultRec *result = new JobResultRec;
	bool setResult = false;
	int flag = 0;
	do{
	//setResult = true;
	penv -> tt_cols =  penv -> m_trajCol.base( penv -> m_child );	//Readjust tt_cols
	penv -> tt_rows=1;

	for(int i=0; i < penv -> m_nParents; i++){
		penv -> tt_rows = penv -> m_trajCol.base( penv -> indexParents[i]-1 ) * penv -> tt_rows;
	}//end for

	//penv -> trans_table = createTable(penv -> tt_rows, penv -> tt_cols);	//create a new truth table
	penv -> delete_array = true;	// Joe Song 08/04/06 Critical
	// ********************cout << "-----------------m_nParents" << penv -> m_nParents << endl;
	if(penv -> m_nParents > penv -> m_gln.size()){
		for(int i = 0; i < penv -> m_nParents; i++){
			if (this -> m_beginParentList[i] == this -> m_endParentList[i])
				flag = 1;
			else{
				flag = 0;
				break;
			}
		}
	}
	if (flag == 0){
	// ****************cout << "process Node id:  " << this -> m_nodeId << endl;
	// ******************cout << "process begin list:  ";
	// ************for(int i = 0; i < penv -> m_nParents; i++){
	// ***************	cout << this -> m_beginParentList[i] << " ";
	// ***************}
	// **************cout << endl;
	// *****************cout << "process end list:  ";
	// ************for(int i = 0; i < penv -> m_nParents; i++){
	// ************	cout  << this -> m_endParentList[i] << " ";
	// ************}
	// ******************cout << endl;
	
	}

	int j = penv -> m_nParents;

	if(! penv -> m_estPar.allowSelfCycle()) {
		for(j=0; j < penv -> m_nParents ; j++) {
			if(penv -> m_child == penv -> indexParents[j]-1) {
				break;
			}
		}
	}

	if( j == penv -> m_nParents ) {
		
		offset(penv->m_estPar, penv->m_trajCol, penv->trans_table, penv->indexParents, 
			penv->indexParents[penv->m_nParents-1]-1, penv->m_child, penv->m_nParents, 
			penv->tt_rows,  penv->tt_cols, penv->indexParents,  penv->p_val, penv->chisq, 
			penv->p_values, penv->chisqs, penv->delete_array, penv->listsOfParents, 
			penv->best_trans_tables, penv->num_parents, penv->sizeTruthTables, 
			penv->degree_of_freedoms, penv->listsOfOffsets, setResult);
	}

	penv -> position = penv -> m_nParents - 1;

	while(penv -> position >= 0){
		//Check to see if current position has reached it's max
		if(penv -> indexParents[penv->position] < (penv->m_gln.size() -(penv->m_nParents-penv->position-1)) ){
			break;						
		}//end if

		penv -> position--;	
		}//end while

		if(penv -> position >= 0){
		penv -> indexParents[penv -> position]++;
		this -> m_beginParentList [penv -> position]++;
			if(penv -> position < penv -> m_nParents-1) {
				//fill in the rest of the indexParents
				for(int i = penv -> position+1; i < penv->m_nParents; i++){
					penv->indexParents[i] = penv->indexParents[i-1]+1;
					this -> m_beginParentList[i] = this-> m_beginParentList[i-1]+1;
				}//end for
			}//end if
		}
		//result -> setResult(penv);
	}while(penv -> position >=0 || penv -> indexParents[0] != penv -> first_max);
	//cout << "11111111111111111111111111111111111111111set Result bool is:  "  << setResult << endl;
	if(setResult == true){
		result -> setResult(penv);
		return (JobResult*) result;
	}
	else{
		if(result)
			delete result;
		result = 0;
		return (JobResult*) result;
	}
	//return (JobResult*) result;
}
*/

/*
void JobGLN::transfer(JobEnvironment &env)
{
	JobEnvironmentGLN *penv = (JobEnvironmentGLN *) &env;
	penv -> m_reconstruct.transfer(* this );
}
*/
