//
//  Topology.cpp
//  gln
//
//  Created by Joe Song on 10/6/11.
//  Copyright 2011 New Mexico State University. All rights reserved.
//
//  Modified:
//    Feb 26, 2012. MS fixed a bug reporting total number of lines in 
//      Topology::save().
//    March 11, 2013. MS Moved Pathway::randomGeneratePathway() to
//      Topology::randomGenerate()

#include <iostream>
using std::cerr;
using std::ios;
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <cstdlib>
#include <algorithm>

using std::find;
using std::sort;
using std::unique;

#include "Topology.h"

vector<int> Topology::getAllNodes() const
// added by Yang Zhang 7.23.2010
// Return indices of all nodes that show up on the topology 
//   regardless of having or not having parents.
//   All IDs start with 1, to be consistent throughout the project
{
	vector <int> AllNodes;
	vector <int> NodesOnPathway;
    
	AllNodes.resize(m_topology.size(),0);
    
	for(size_t i=0; i<m_topology.size(); i++ )
	{
		if(m_topology[i][0].size() != 0)
		{
			AllNodes[i] = 1;
			for(size_t j=0; j<m_topology[i][0].size(); j++)
			{
				AllNodes[ m_topology[i][0][j] - 1] = 1;
			}
		}
	}
	
	for (size_t index =0; index <AllNodes.size(); index++)
	{
		if(AllNodes[index]==1)
		{
			NodesOnPathway.push_back(index+1);
		}
	}
    
	return NodesOnPathway;
}

vector<int> Topology::getAllChildNodes() const
// Return indices of all child nodes that have parents on the topology
{
	vector <int> allChildNodes;
    
	for(size_t i=0; i<m_topology.size(); i++ )
	{
		if(m_topology[i][0].size() != 0)
		{
			allChildNodes.push_back(i + 1);
		}
	}
	
	return allChildNodes;
}

void Topology::scan(const string & file)
{
    if( file.empty() ) return;
            
	ifstream ifs(file.c_str(), ios::in); 	//open file to load data
    
	if(!ifs.is_open()) { //check to see if file is open
		cerr << "ERROR: cannot open file \"" << file << "\"!" << endl;
		exit(EXIT_FAILURE);
	}
    
    size_t nlines = 0;
	string line; 
    
	//added by YangZhang 09/05/09 to be clarified with pathway parents file
	getline(ifs, line);
    
	//Get the number of lines in the topology file, 
    //  It is not necessarily the size of the Network
	if(getline(ifs, line)) {
		nlines = atoi(line.c_str());
	}
    
    if(nlines > 0 && size() < nlines) {
        cout << "Warning: Topology file many contain more nodes than the trajectory file!" << endl;
        // exit(EXIT_FAILURE);
    }
    
	vector<string> strAdjacencyList(nlines);
    
    size_t numNonEmptyLines = 0;
    
	while( ! ifs.eof() ) { 	//Copy each line of non-empty text to vector strAdjacencyList.

        if(numNonEmptyLines < nlines) {
            
            getline(ifs, strAdjacencyList[ numNonEmptyLines ]);

            numNonEmptyLines = (strAdjacencyList[ numNonEmptyLines ].size() > 0) ? numNonEmptyLines+1 : numNonEmptyLines;
            
        } else {
            
            getline(ifs, line);
            
            if(line.size() > 0) {
                cerr << "ERROR: The number of non-empty lines in file \"" << file 
                    << "\" is more than the total number of lines reported in the file!" << endl;
                exit(EXIT_FAILURE);
            }
        }
	}
    
    if(numNonEmptyLines != nlines) {

		cerr << "ERROR: The number of non-empty lines in file \"" << file 
            << "\" is different from the total number of lines reported in the file!" << endl;
		exit(EXIT_FAILURE);
        
    }
    
	ifs.close(); 	// Close the candidate file
    
    size_t nMatched = 0;
    
	for (size_t i=0; i < strAdjacencyList.size(); i++ ) {
        
        int m = strAdjacencyList[i].find('*');
        
		string nameChild = strAdjacencyList[i].substr(0, m);
        
        size_t idChild = find(m_nodeNames.begin(), m_nodeNames.end(), nameChild) - m_nodeNames.begin();
        
        if(idChild >= m_nodeNames.size()) {
            cout << "Warning: Variable " << nameChild 
            << " appeared in candidate file \"" << file 
            << "\" but is not a variable in the TCF file!" << endl;
            continue;
            
		} else {
            
            nMatched ++;
            
        }

        m ++;

        string strParents = strAdjacencyList[i].substr(m);
        
        /*
        string strParentscopy = strParents;
        
		int nParents = 0; //number of parents

		while( strParentscopy.find(',',0)!=string::npos )
		{
			strParentscopy = strParentscopy.substr( strParentscopy.find(',', 0) + 1 );
			nParents ++;
		}
        */
        

        int nParents = 0; //number of parents

		while( strAdjacencyList[i].find(',', m) != string::npos ) {
            m = strAdjacencyList[i].find(',', m) + 1;
			nParents ++;
		}
        
		vector<int> idParents;
		vector<int> delayParents;
        
		for(int j=0; j < nParents; j++) {
            //modified by Yang Zhang on 5.20.2013
            //now one topology file only corresponds to one condition
			string strParentAndDelay = strParents.substr(0, strParents.find(',', 0) ); //like GeneA -1,
			string strParentName = strParentAndDelay.substr(0, strParentAndDelay.find('\t', 0)); // GeneA
			            
			size_t pid = find(m_nodeNames.begin(), m_nodeNames.end(), strParentName) - m_nodeNames.begin();
            
            if( pid < m_nodeNames.size() ) { // if parent id is found
                
				idParents.push_back(pid + 1); // record parent id after adjusting to base 1
                
				if(strParentAndDelay == strParentName) { //no delay is provided
                    
					delayParents.push_back(9); //using 9 to represent no delays requirements
                    
				} else {
                    
						int delay = atoi(strParentAndDelay.substr( strParentAndDelay.find('\t', 0) + 1 ).c_str());
						delayParents.push_back(delay);
                }		
			}
			strParents = strParents.substr(strParents.find(',', 0) + 1);
		}
        
        //idChild: 0-based
		m_topology[idChild][0] = idParents;
		m_topology[idChild][1] = delayParents;
        
	}
         
	for(size_t r = 0; r < m_topology.size(); r ++) {  
        // Check for duplicated entries.  Should check for delays as well.  A duplicate should be
        //   defined as a node having the same parent id and delay with another node.
        
		vector<int> idParents = m_topology[r][0];
		
		sort(idParents.begin(), idParents.end());
        
		// idParents.erase(unique(idParents.begin(), idParents.end()), idParents.end());
        
        unique( idParents.begin(), idParents.end() );
        
		if( idParents.size() != m_topology[r][0].size() ) {
            
			cerr << "ERROR: found duplicated entries for " << m_nodeNames[r] 
                << ". Check the pathway file " << file << "!" << endl;

			exit(EXIT_FAILURE);
		}
	}
    
    if(nMatched < strAdjacencyList.size()) {
        cout << "Warning: " << nMatched << " out of " << strAdjacencyList.size() 
            << " nodes in the topology file matched the trajectory file." << endl;
    }
}

void Topology::save(const string & file) const
{  
    ofstream ofs(file.c_str(), ios::out);
	if(!ofs.is_open()) 
	{
		cerr << "ERROR: cannot open file \"" << file << "\" for saving!" << endl;
		exit(EXIT_FAILURE);
	}
    
	ofs << "Pathway/Candidate File" << endl;
    
    // Feb 26, 2012. MS added the following 7 lines
    size_t nNodesWithParents=0;
	for (size_t i=0; i < m_topology.size(); i++) {
		if(m_topology[i][0].size() > 0) {
            nNodesWithParents ++; 
        }
    }
	ofs << nNodesWithParents << endl;
    
	for (size_t i=0; i < m_topology.size(); i++)	{
        
		if(m_topology[i][0].size() > 0) {
			ofs << m_nodeNames[i] << "*";  
            for(size_t j=0; j < m_topology[i][0].size(); j++ )
            {
                ofs << m_nodeNames[ m_topology[i][0][j]-1 ] << "\t";
                for(int k=0; k < m_nConditions; k++ )
                {
                    if(m_topology[i][1].size()>0)
                    {
                        if(k!=m_nConditions-1)
                        {
                            ofs << m_topology[i][1][m_nConditions*j+k] << "\t";
                        }else
                        {
                            ofs << m_topology[i][1][m_nConditions*j+k] << ",";
                        }
                    }else
                    {
                        if(k!=m_nConditions-1)
                        {
                            ofs << "9" << "\t";
                        }else
                        {
                            ofs << "9" << ",";
                        }
                    }
                }
            }
            ofs << endl;
        }
        
    }
	
	ofs.close();  // Close file
}

//added by Yang Zhang 6/14/2009
void Topology::randomGenerate(vector < vector<int> > ParentsVector,
                              int NumberOfNodes,const char * fileTopology)
{
    // srand( (unsigned)time(NULL) );
    
    if (NumberOfNodes >= (int)ParentsVector.size())
    {
        cout<<"Pathway size can not be bigger than network size!"<<endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        ofstream out(fileTopology);
        out<<ParentsVector.size()<<endl;
        
        vector<int> indexOfNodes;
        //indexOfNodes.resize(NumberOfNodes);
        for(int i=0;i<NumberOfNodes;i++)
        {
            int randomnumber = (int)(((double)(rand())/RAND_MAX)*(ParentsVector.size()-1));
            for(size_t j=0; j<indexOfNodes.size(); j++)
            {
                if(randomnumber==indexOfNodes[j])
                {
                    randomnumber = (int)(((double)(rand())/RAND_MAX)*(ParentsVector.size()-1));
                }
            }
            indexOfNodes.push_back(randomnumber);
        }
        
        for(size_t i=0;i<indexOfNodes.size();i++)
        {
            out<<indexOfNodes[i]+1<<"*";
            for(size_t j=0;j<ParentsVector[i].size();j++)
            {
                out<<ParentsVector[i][j]<<",";
            }
            out<<endl;
        }
        out.close();
    }
}
