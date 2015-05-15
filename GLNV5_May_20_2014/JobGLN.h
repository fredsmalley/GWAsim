// JobGLN.h
//
// Curtis Luce, Joe Song
// Created 2007
// Last modified: September 7, 2008

#pragma once

#include "JobEnvironment.h"
#include "Job.h"

//added by Yang Zhang, to obtain BP(can be parent) or HP(can have parent) meta information
#include "TrajectoryCollection.h"

#include <iostream>
using std::cout;
using std::endl;

#include <vector>
using std::vector;

class JobGLN: public Job {

private:

	int m_nodeId;
	vector<int> m_beginParentList;
	vector<int> m_endParentList;
	int m_size;

public:

	JobGLN():m_nodeId(0),m_beginParentList(), m_endParentList(), m_size(0){}

	void initialize(int id, const vector<int> & beginParents, const vector<int> & endParents, int size) {
		m_nodeId = id;
		m_beginParentList = beginParents;
		m_endParentList = endParents;
		m_size = size;
	}

	virtual void first(const JobEnvironment &env);
	virtual void next(const JobEnvironment &env);

	void terminal() { m_nodeId = -1; m_size = 0; m_beginParentList.clear(); m_endParentList.clear(); }

	// void transfer(JobEnvironment &env);

	// JobResult* process(JobEnvironment &env);
	virtual bool process(JobEnvironment &env) const;

    virtual bool processToBeDeleted(JobEnvironment &env) const;

	// bool isTerminal(JobEnvironment &env);
	virtual bool isTerminal() const { return -1 == m_nodeId; }

	virtual int pack(unsigned char* &jobBuffer) const;
	virtual void unpack(unsigned char * jobBuffer, int length);

	void setBeginList(int index, int value) { m_beginParentList[index] = value; }
	int getBeginList(int index) const {	return m_beginParentList[index]; }
	
	int getEndList (int index) const { return m_endParentList[index]; }

	const vector<int> & getBeginParents() const { return m_beginParentList; }
	const vector<int> & getEndParents() const { return m_endParentList; }

	void setNodeId(int id) { m_nodeId = id; }
	int getNodeId () const { return m_nodeId; }

	int getSize() const	{ return m_size; }

	virtual void print() const
	{
		cout << "Node=" << m_nodeId
		 << ", Parents(begin)=";
		for (int i = 0; i < m_size; i++)  
			cout << m_beginParentList [i] << ":";
		cout << ", Parents(end)=";
		for (int i = 0; i < m_size; i++)  
			cout << m_endParentList [i] << ":";
		cout << '\r';
	}

	virtual void brief() const
	{
		// cout << "N" << getNodeId() << '-' << getSize() << ',';
	}

	//Managed by Yang and CMHJ:  Prints a description file of all signification contingency table
	//and their counts.
	void ResutlsToFile(JobEnvironment &env) const;
};

