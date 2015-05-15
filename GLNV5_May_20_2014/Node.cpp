// Node.cpp
//
// Joe Song
// Created: May, 2006

#include <iostream>
#include <vector>
#include <string>
#include "Node.h"

using namespace std;

void Node::initialize(unsigned short initialValue, unsigned short radix, 
					  string node_name, unsigned char node_type, int node_id)
{
	// index_of_node=index;
	m_base = radix;

	m_value = initialValue;

	//check to see if default value
	if(node_name.compare("!!") != 0) // if(node_name != "!!")
		m_name = node_name;
	//check to see if default value
	if(node_type != 'n')
		m_type = node_type;

	//check to see if default value
	if(node_id != -1)
		m_nodeId = node_id;

	
	// if(p_val != -1.0)
	//	p_value = p_val; 

	/*
	NumParents=num_parents;
	pred_func = pred;
	parentNodes=parents;
	parent_time_offset = offsets;
	*/

}//end

/*
void Node::initialize(int index,unsigned short initialValue, unsigned short radix, int num_parents, vector<Node *> parents, Predictor pred, string node_name, unsigned char node_type, vector<int> offsets, int node_id, double p_val){
	index_of_node=index;
	base=radix;
    	NumParents=num_parents;
   	value = initialValue;
    	pred_func = pred;
   	parentNodes=parents;
	//check to see if default value
	if(node_name.compare("!!") != 0) // if(node_name != "!!")
		m_name = node_name;
	//check to see if default value
	if(node_type != 'n')
		m_type = node_type;
	parent_time_offset = offsets;
	//check to see if default value
	if(node_id != -1)
		nodeId = node_id;
	if(p_val != -1.0)
		p_value = p_val; 
}//end
*/

// int Node::getThisIndex() const {
//	return index_of_node;
//}//end getThisIndex

void Node::setValue(unsigned short val){
    m_value=val;
}//end setValue

unsigned short Node::getValue() const {
	return m_value;
}//end getValue

unsigned short Node::getBase() const {
	return m_base;
}//end getBase

/*
int Node::getNumParents() const {
	return NumParents;
}//end getNumParents

int Node::getParent(int number) const {
//This only fetches a parent given by the index, should be used together with a for loop
	if(number >=  NumParents){
		cout << "Error in function 'Node::getParent' : 'number' >= 'NumParents'" << endl;
		return -1;
	}//end if

	return parentNodes[number]->getThisIndex();
}//end getParents

vector<int> Node::getOffset() const {
	//return the vector, offsets are in order of parents
	return parent_time_offset;
}//end getOffset

vector<unsigned short> Node::getBaseOfParents() const {
	return pred_func.getParentBases();
}//end getBaseOfParent

vector<unsigned short> Node::getTruthTable() const {
	return pred_func.getTruthTable();
}//end getBaseOfParent


unsigned short Node::predict(vector< vector<unsigned short> >& trajectory){
    vector<Node *>::iterator iter = parentNodes.begin();
    vector<unsigned short> vals;
    //only first order offsets use this, case for when a initial trajectory file not provided
    if(trajectory.size() == 0){
	for(int i=0;i < NumParents;i++) {
        	vals.push_back((*iter)->value);
        	iter++;
    	}//end for
   }//end if
   
   //first to k-th order use this, case for initial trajectory file
   else{
    	for(int i=0;i < NumParents;i++){
		vals.push_back(trajectory[trajectory.size()+parent_time_offset[i]] [((*iter)->getId())-1]);
		//cout << trajectory[trajectory.size()+parent_time_offset[i]][((*iter)->getId())-1] << " " << ((*iter)->getId())-1 << " " << trajectory.size()<< " " <<parent_time_offset[i] << endl;
		iter++;
	}//end for
   }//end else
    
    return pred_func.evaluate(vals);

}//end predict
*/

/*  Removed. Joe Song. 10/6/2011.  This is a very general function and does not seem to belong to here. 
//added by Yang Zhang 10/28/2009 to convert genenames to id for GLN program to use
int Node::returnIndex(const string & node, const vector<string> & nodeNames)
{
	for(size_t i=0; i < nodeNames.size(); i++)
	{
		if(node==nodeNames[i])
		{
			return (int) i;
		}
	}
	return -1;
}
*/



//Added by Haizhou Wang, Oct 25, 2012
bool Node::operator==(const Node & node) const
{
    //NOTE: m_value maybe different
    if( m_type != node.m_type ||
        m_name != node.m_name ||
        m_nodeId != node.m_nodeId ||
        m_base != node.m_base ||
        m_bp != node.m_bp ||
        m_hp != node.m_hp
        )
        return false;

    return true;
}
//End of Haizhou's code