// Node.h
//
// Eric Lance, Chris Lewis, Joe Song
// Created: May 2006
// Last modified: September 27, 2008. Joe Song

#pragma once

#include <list>
using std::list;

#include <vector>
using std::vector;

#include <string>
using std::string;

// #include "predictor.h"

class Node {
//private:
public:
	unsigned char m_type;	// 'e': external, 'i': internal, Joe Song 08/05/06
	string m_name; // Joe Song 08/05/06
	int m_nodeId;
	unsigned short m_value;
	unsigned short m_base;			//radix


	//added by Yang Zhang 4-24-2009
	bool m_bp; //variable indicate whether or not a node can be parent
	bool m_hp; //variable indicate whether or not a node can have parent
		
	// double p_value;

	/*
	int index_of_node;			//Holds the current index of this node

	vector<Node *> parentNodes;		//Linked list of parent nodes
	vector<int> parent_vector;
	vector<int> parent_time_offset;		//holds parents time offsets (t-1, t-2, t-n)
	vector<int> temp_parents;		//Temporary vector of parent nodes
	Predictor pred_func;          		//Unique predictor function for node
	int NumParents;

	*/

	//Input output functions.
	//operator>>();
	//operator<<();

public:

	void initialize(unsigned short initialValue, unsigned short radix, 
		string node_name, unsigned char node_type, int node_id);

	int getThisIndex() const;
	unsigned short getValue() const;		 //Returns "value"
	unsigned short getBase() const;
	void setValue(unsigned short val);   //Set the value incase we suddendly need to change it manually

	void setType(char t) { m_type = t; }  // Joe Song 08/05/06
	char getType() const { return m_type; }	// Joe Song 08/05/06
	void setName(const char * name) { m_name = name; }  // Joe Song 08/05/06
	const string & getName() const { return m_name; }	// Joe Song 08/05/06
	vector<int> getOffset() const;
	int getId() const {return m_nodeId;}

	void setBP(bool bp) { m_bp = bp; }	//Yang Zhang 4/24/2009
	void setHP(bool hp) { m_hp = hp; }  //Yang Zhang 4/24/2009

	bool getBP() const { return m_bp; } //Yang Zhang 4/24/2009
	bool getHP() const { return m_hp; } //Yang Zhang 4/24/2009
	
	// static int returnIndex(const string & node, const vector<string> & nodeNames);
	/*

	void initialize(int index, unsigned short initialValue,
		unsigned short radix,
		int num_parents,
		vector<Node *> parents,
		Predictor pred,string node_name = "!!", unsigned char node_type = 'n',
		vector<int> offsets = vector <int> (0), int node_id = 0, double p_val = -1.0
		);

	double getPvalue() const { return p_value; }
	void setPvalue(double val) { p_value = val; }

	int getNumParents() const;
	int getParent(int number) const;
	vector<unsigned short> getBaseOfParents() const;
	vector<unsigned short> getTruthTable() const;
	unsigned short predict(vector< vector<unsigned short> >& trajectory);			
	*/


    //Added by Haizhou Wang, Oct 25, 2012
    bool operator==( const Node & node ) const;
    //End of Haizhou's code
};     

