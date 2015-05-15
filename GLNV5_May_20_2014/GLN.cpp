// GLN.cpp -- The GeneralizedLogicalNetwork class
//
// Eric Lance, Chris Lewis, Curtis Luce, Joe Song
// Created: 2006
// Modified: September 5, 2008
// Last modified: November 30, 2008

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <fstream>
using std::ofstream;

#include <sstream>
using std::stringstream;

#include <string>
using std::string;

#include <utility>
using std::pair;

#include <cstdlib>

#include "GLNGlobals.h"
#include "GLN.h"
#include "PermStatTable.h"

const vector<string> GeneralizedLogicalNetwork::getNodeNames(char type) const
{
    vector<string> names;
    
    for(size_t j=0; j<size(); j++) {
        if(m_nodes[j].getType() == type)
        {
            names.push_back(m_nodes[j].getName());
        }
    }
    return names;
}

const vector<string> GeneralizedLogicalNetwork::getNodeNames(const vector<int> & ids) const
{
    vector<string> names(ids.size());
    
    for(size_t j=0; j<ids.size(); j++) {
        
        names[j] = m_nodes[j].getName();
        
    }
    return names;    
}


int GeneralizedLogicalNetwork::detectCycle(const vector< vector<int> > & trajectory) const
{
	for(int t=m_current_time-1;t >= 0; t--) {
		size_t i;
		for(i=0; i < size(); i++)
			if (trajectory[t][i] != trajectory[m_current_time][i])
				break;
		if (i == size())
			return (m_current_time - t); //found
	}
	return -1; //not found
}

void GeneralizedLogicalNetwork::initialize(int networksize) 
{
	m_current_time=0;
	m_ready = false;			//No nodes have been inserted yet so default is false
	m_nodes.resize(networksize);
	m_gtts.resize(networksize);
}//end initialize

void GeneralizedLogicalNetwork::initialize(int nNodes, const vector<int> & extNodes,
	const vector<string> & extNodeNames, const vector<string> & intNodeNames)
{
	initialize(nNodes);

	for(int i=0; i < nNodes; i++){  // Joe Song 8/5/06
		m_nodes[i].setType('i'); 
	}

	for(int m=0; m < (int)extNodes.size(); m++){  // Joe Song 8/5/06
		m_nodes[extNodes[m]-1].setType('e'); 
	}

	for(int j=0, i=0, m=0; j<nNodes; j++) {
		switch(m_nodes[j].getType()) {
			
			case 'i':
				if(intNodeNames.size() > 0) {
					// char node_name[50]; // Joe Song 8/8/2008
					// sprintf(node_name, "%s", intNodeNames[i].c_str()); // Joe Song 8/8/2008
					// m_nodes[j].setName(node_name); // Joe Song 8/8/2008
					m_nodes[j].setName(intNodeNames[i].c_str());  // Joe Song 8/8/2008
					i ++;
				}
				else{
					//this is for the default case when no node name is
					//specified, the default name will be node%d, j will fill
					//in %d which is current node,  sprintf makes a string node_name

					stringstream numberStream;
					numberStream << j+1;
					string nodeName = "Node" + numberStream.str();
					
					// char node_name[50];
					// sprintf (node_name, "node%d", j+1);

					m_nodes[j].setName(nodeName.c_str()); //set node name

				}
				break;
			case 'e':
				if(extNodeNames.size() > 0) {
					//char node_name[50]; // Joe Song 8/8/2008
					//sprintf (node_name, "%s", extNodeNames[m].c_str()); // Joe Song 8/8/2008
					//m_nodes[j].setName(node_name); // Joe Song 8/8/2008
					m_nodes[j].setName(extNodeNames[m].c_str()); // Joe Song 8/8/2008
					m ++;
				}
				else{
					//this is for the default case when no node name is
					//specified, the default name will be node%d, j will fill
					//in %d which is current node,  sprintf makes a string node_name
					// char node_name[50];
					// sprintf (node_name, "node%d", j+1);

					stringstream numberStream;
					numberStream << j+1;
					string nodeName = "Node" + numberStream.str();

					m_nodes[j].setName(nodeName.c_str());  //set node name

				}
				break;
			default:
				cerr << "ERROR: Unknown node type " << m_nodes[j].getType() << "!" << endl;
				GLNExit(EXIT_FAILURE);
		}
	}
}

//added by YangZhang 4/24/2009
void GeneralizedLogicalNetwork::initialize(int nNodes, const vector<int> & extNodes,const vector<string> & extNodeNames, 
	const vector<string> & intNodeNames,const vector<bool> & beParent, const vector<bool> & haveParent)
{
	initialize(nNodes);

	for(int i=0; i < nNodes; i++){  // Joe Song 8/5/06
		m_nodes[i].setType('i'); 
		m_nodes[i].setBP(beParent[i]);
		m_nodes[i].setHP(haveParent[i]);
	}

	for(int m=0; m < (int)extNodes.size(); m++){  // Joe Song 8/5/06
		m_nodes[extNodes[m]-1].setType('e'); 
	}

	for(int j=0, i=0, m=0; j<nNodes; j++) {
		switch(m_nodes[j].getType()) {
			
			case 'i':
				if(intNodeNames.size() > 0) {
					// char node_name[50]; // Joe Song 8/8/2008
					// sprintf(node_name, "%s", intNodeNames[i].c_str()); // Joe Song 8/8/2008
					// m_nodes[j].setName(node_name); // Joe Song 8/8/2008
					m_nodes[j].setName(intNodeNames[i].c_str());  // Joe Song 8/8/2008
					i ++;
				}
				else{
					//this is for the default case when no node name is
					//specified, the default name will be node%d, j will fill
					//in %d which is current node,  sprintf makes a string node_name

					stringstream numberStream;
					numberStream << j+1;
					string nodeName = "Node" + numberStream.str();
					
					// char node_name[50];
					// sprintf (node_name, "node%d", j+1);

					m_nodes[j].setName(nodeName.c_str()); //set node name

				}
				break;
			case 'e':
				if(extNodeNames.size() > 0) {
					//char node_name[50]; // Joe Song 8/8/2008
					//sprintf (node_name, "%s", extNodeNames[m].c_str()); // Joe Song 8/8/2008
					//m_nodes[j].setName(node_name); // Joe Song 8/8/2008
					m_nodes[j].setName(extNodeNames[m].c_str()); // Joe Song 8/8/2008
					m ++;
				}
				else{
					//this is for the default case when no node name is
					//specified, the default name will be node%d, j will fill
					//in %d which is current node,  sprintf makes a string node_name
					// char node_name[50];
					// sprintf (node_name, "node%d", j+1);

					stringstream numberStream;
					numberStream << j+1;
					string nodeName = "Node" + numberStream.str();

					m_nodes[j].setName(nodeName.c_str());  //set node name

				}
				break;
			default:
				cerr << "ERROR: Unknown node type " << m_nodes[j].getType() << "!" << endl;
				GLNExit(EXIT_FAILURE);
		}
	}
}

//added by yangzhang 11.15.2008,return number of internal nodes to decide how to generate random initial values
int GeneralizedLogicalNetwork::getNumofInternalNode() const
{
	int numberofinternalnode = 0;

	for(size_t i=0; i<m_nodes.size(); i++)
	{
		if(m_nodes[i].getType()=='i')
		{
			numberofinternalnode++;
		}
	}
	return numberofinternalnode;
}

/*
void GeneralizedLogicalNetwork::initializeNode(int nodenum, unsigned short initialValue, 
									  unsigned short radix, int num_parents,
									  vector<int> parentindexes,Predictor pred,vector<int> parentoffsets, 
									  double p_val, int node_id, string node_name,unsigned char node_type)
{
	m_ready = true;					//network has at least one node
	vector<Node *> parents;
	vector<int> offsets;
	vector<Node>::iterator net_iter = m_nodes.begin();
	for(int i=0;i< num_parents; i++) {
		parents.push_back(&(*(net_iter+parentindexes[i]))); //does it copy?
		offsets.push_back(((parentoffsets[i])));
	}

	m_nodes[nodenum].initialize(nodenum, initialValue,radix,num_parents,parents,pred,node_name,node_type, offsets, node_id,p_val);
}//end initializeNode
*/

void GeneralizedLogicalNetwork::initializeNode(int nodeID, unsigned short initialValue, 
											   unsigned short radix, int node_id, 
											   string node_name, unsigned char node_type)
{
	m_ready = true;					//network has at least one node
	m_nodes[nodeID].initialize(initialValue, radix, node_name, node_type, node_id);
}//end initializeNode

/*
//called when no parent offsets are provided initialize to zero
void GeneralizedLogicalNetwork::initializeNode(int nodenum, unsigned short initialValue, 
									  unsigned short radix, int num_parents,
									  vector<int> parentindexes,Predictor pred, double p_val, string node_name,unsigned char node_type)
{
	vector<int> parentoffsets;
	//initialize all parent offsets to 0
	for (int i = 0; i< num_parents; i++)
		parentoffsets.push_back(0);

	initializeNode (nodenum, initialValue, radix, num_parents, parentindexes, pred, parentoffsets, p_val, -1, node_name, node_type);

}//end initializeNode
*/

//called when no parent offsets are provided initialize to zero
void GeneralizedLogicalNetwork::initializeNode(int nodeID, unsigned short initialValue, 
									  unsigned short radix, string node_name,unsigned char node_type)
{
	initializeNode(nodeID, initialValue, radix, -1, node_name, node_type);
}//end initializeNode

int GeneralizedLogicalNetwork::getMaxMarkovOrder() const 
{
	return m_max_Markov_order;
}

void GeneralizedLogicalNetwork::setMaxMarkovOrder(int order) 
{
	m_max_Markov_order = order;
}

vector<int> GeneralizedLogicalNetwork::getBases() const
{
	vector<int> bases( size() );

	for(size_t i=0; i < size(); i++)
	{
		bases[i] = m_nodes[i].getBase();
	}

	return bases;
}


/*
added by yang zhang 4.3.2009
in order to compare two truth tables,children must be in same order
*/
//bool GeneralizedLogicalNetwork::compareTruthTables(const vector<GeneralizedTruthTable> & m_gtts1,const vector<GeneralizedTruthTable> & m_gtts2)
//{
//	bool compResult = true;
//	if(m_gtts1.size()!=m_gtts2.size())
//	{
//		compResult = false;
//	}
//	else
//	{
//		for(int i=0;i<m_gtts1.size();i++)
//		{
//			if(m_gtts1[i].getChild()!=m_gtts2[i].getChild())
//			{
//				compResult = false;
//			}
//			else
//			{
//				if(m_gtts1[i].getParents().size()!=m_gtts2[i].getParents().size())
//				{
//					compResult = false;
//				}
//				else
//				{
//					for(int j=0;j<m_gtts1[i].getParents().size();j++)
//					{
//						if ( (m_gtts1[i].getParents())[j] != (m_gtts2[i].getParents())[j] )
//						{
//							compResult = false;
//						}
//						else
//						{
//							for(int k=0;k<m_gtts1[i].getTruthValues().size();k++)
//							{
//								if( (m_gtts1[i].getTruthValues())[k] != (m_gtts2[i].getTruthValues())[k] )
//								{
//									compResult = false;
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	return compResult;
//}

bool GeneralizedLogicalNetwork::compareTruthTables(const GeneralizedLogicalNetwork & gln)
{
	vector<GeneralizedTruthTable> & m_gtts1 = m_gtts;
	const vector<GeneralizedTruthTable> & m_gtts2 = gln.getGTTs();
	bool compResult = true;
	if(m_gtts1.size()!=m_gtts2.size())
	{
		compResult = false;
	}
	else
	{
		for(int i=0;i<(int)m_gtts1.size();i++)
		{
			if(m_gtts1[i].getChild()!=m_gtts2[i].getChild())
			{
				compResult = false;
			}
			else
			{
				if(m_gtts1[i].getParents().size()!=m_gtts2[i].getParents().size())
				{
					compResult = false;
				}
				else
				{
					for(int j=0;j<(int)m_gtts1[i].getParents().size();j++)
					{
						if ( (m_gtts1[i].getParents())[j] != (m_gtts2[i].getParents())[j] )
						{
							compResult = false;
						}
						else
						{
							for(int k=0;k<(int)m_gtts1[i].getTruthValues().size();k++)
							{
								if( (m_gtts1[i].getTruthValues())[k] != (m_gtts2[i].getTruthValues())[k] )
								{
									compResult = false;
								}
							}
						}
					}
				}
			}
		}
	}
	return compResult;
}

/*
void GeneralizedLogicalNetwork::compareAllTruthTables(const GeneralizedLogicalNetwork & gln, const string & fileIntxCmp)
{
	vector<char> interactionComparison(gln.size());

	const vector<GeneralizedTruthTable> & m_gtts1 = m_gtts;

	const vector<GeneralizedTruthTable> & m_gtts2 = gln.getGTTs();

	bool compResult = true;

	if(m_gtts1.size()!=m_gtts2.size())
	{
		compResult = false;
	}
	else
	{
		for(int i=0;i<(int)m_gtts1.size();i++)
		{
			interactionComparison[i] = 'C';

			if(m_gtts1[i].getChild()!=m_gtts2[i].getChild())
			{
				interactionComparison[i] = 'D';
			}
			else
			{
				if(m_gtts1[i].getParents().size()!=m_gtts2[i].getParents().size())
				{
					interactionComparison[i] = 'D';
				}
				else
				{
					for(int j=0;j<(int)m_gtts1[i].getParents().size();j++)
					{
						if ( (m_gtts1[i].getParents())[j] != (m_gtts2[i].getParents())[j] )
						{
							interactionComparison[i] = 'D';
						}
						else
						{
							for(int k=0;k<(int)m_gtts1[i].getTruthValues().size();k++)
							{
								if( (m_gtts1[i].getTruthValues())[k] != (m_gtts2[i].getTruthValues())[k] )
								{
									interactionComparison[i] = 'D';
									break;
								}
							}
						}
					}
				}
			}
		}
	}

	ofstream ofs(fileIntxCmp.c_str());

	for(int i=0; i<(int)interactionComparison.size(); i++) {
		ofs << interactionComparison[i] << endl;
	}

	ofs.close();

}

*/

void GeneralizedLogicalNetwork::compareAllTruthTables(const GeneralizedLogicalNetwork & gln, const string & fileIntxCmp)
{
	const vector<GeneralizedTruthTable> & m_gtts1 = m_gtts;

	const vector<GeneralizedTruthTable> & m_gtts2 = gln.getGTTs();

	if(m_gtts1.size() == m_gtts2.size()) {

		vector<char> interactionComparison( gln.size() );
		vector<double> distances( gln.size() );
		vector<bool> sameParents( gln.size() );

		for(int i=0;i<(int)m_gtts1.size();i++)	{

			interactionComparison[i] = compareGenTruthTables(m_gtts1[i], m_gtts2[i]) ? 'C' : 'D';
			distances[i] = distanceTransProbTables(m_gtts1[i], m_gtts2[i]);
			sameParents[i] = sameChildParentsDelays(m_gtts1[i], m_gtts2[i]);

		}

        ofstream ofs(fileIntxCmp.c_str());

		ofs << "SameGTT\tSameParents\tTransProbTableDistance\tParents1\tParents2" << endl;

    	for(int i=0; i<(int)interactionComparison.size(); i++) {
	    	ofs << interactionComparison[i] << '\t' << sameParents[i] << '\t' << distances[i] << '\t';

			for(size_t j=0; j < m_gtts1[i].getParents().size(); j++) {

				if(j > 0) ofs << ',';

				ofs << m_gtts1[i].getParents()[j];
			}

			ofs << '\t'; 

			for(size_t j=0; j < m_gtts2[i].getParents().size(); j++) {

				if(j > 0) ofs << ',';

				ofs << m_gtts2[i].getParents()[j];
			}

			ofs << endl;
    	}

    	ofs.close();

	} else {

		cerr << "WARNING: Two GLN files are not the same size. They are not compared." << endl;

	}

}


//Added by Haizhou Wang, March 23, 2012
void GeneralizedLogicalNetwork::eliminateFictitiousParents(const Topology & mustIncludeParents )
{
	for(size_t i=0; i<m_gtts.size(); i++)
    {
		if( m_nodes[i].getType() == 'i') //Only check internal node, since external node does NOT have truth table
        {
            if( mustIncludeParents.size() != 0)
            {
                m_gtts[i].getRidOfFictitious( mustIncludeParents.get()[i][0] );
            }
            else
            {
                m_gtts[i].getRidOfFictitious();
            }
        }
    }
}


//Added by Haizhou Wang, Oct 25, 2012
bool GeneralizedLogicalNetwork::operator==( const GeneralizedLogicalNetwork & gln ) const
{
    if( m_nodes != gln.m_nodes ||
        m_min_Markov_order != gln.m_min_Markov_order ||
        m_max_Markov_order != gln.m_max_Markov_order ||
        m_current_time != gln.m_current_time ||
        m_version != gln.m_version ||
        m_gtts != gln.m_gtts ||
        m_ready != gln.m_ready
        )
        return false;


    return true;
}
//End of Haizhou's code


//Added by Haizhou Wang, March 18, 2013
vector< pair<int,int> > GeneralizedLogicalNetwork::getEdgeList() const
{
    vector< pair<int,int> > result;
    for( const auto & node: m_gtts ){
        for( const auto & parent: node.getParents() ){
            //edge pair: <parent, child>
            result.push_back( pair<int, int>( parent, node.getChild()+1) );  //edge<1-based, 1-based>  
        }
    }

    return result;
}


void GeneralizedLogicalNetwork::saveEdgeList(const string & outputFileName) const
{
    ofstream ou(outputFileName);
    ou<<"Parent\tChild"<<endl;

    for( const auto & node: m_gtts ){
        for( const auto & parent: node.getParents() ){
            //edge pair: <parent, child>
            //result.push_back( pair<int, int>( parent, node.getChild()+1) );  //edge<1-based, 1-based>  
            ou<<m_nodes[parent-1].getName()<<"\t"<<m_nodes[node.getChild()].getName()<<endl;
        }
    }

    ou.close();
}
