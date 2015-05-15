//
//  EnvGLNRecIO.cpp
//  gln
//
//  Created by Joe Song on 10/4/11.
//  Copyright 2011 New Mexico State University. All rights reserved.
//
// Modified:
//   July 13, 2014. MS. Added saveTopology() member function


#include <iostream>
#include <sstream>
using namespace std;

#include "EnvGLNRec.h"

void print_adjusted(const GeneralizedLogicalNetwork & net, const vector<TransitionTable> & tts, double alpha); 

void EnvGLNRec::print(void) const
{
	cout << endl << endl;
    
	cout << "id\t";
	cout << "name\t";
	cout << "type\t";
	cout << "num.parents\t";
	cout << "parents\t";
	cout << "offsets\t";
	cout << "p.value\t";
	cout << "chisq\t";
	cout << "df" << endl;
    
	for(size_t i=0; i < m_gln.size(); i++){
        
		cout << i+1 << "\t";
		cout << m_gln.getNodeName(i) << "\t";
		cout << m_gln.getNodeType(i) << "\t";
        
		m_transTables[i].heading();
        
	}
    
	cout << "Overall p-value of the reconstructed generalized logical network = " 
        << m_p_val << endl << endl;
    
	print_adjusted(m_gln, m_transTables, m_gc.m_alpha);	
    
    cout << "Significant transition tables:" << endl;

	for(size_t i=0; i < m_gln.size(); i++){
        
        // cout << i+1 << '(' << m_gln.getNodeName(i) << ')' <<": ";
        
		if(m_transTables[i].nrow()==0){
			cout << "None." << endl;
		} else {
			m_transTables[i].print();
		}
	}
    
	m_gln.printGTTs();
}

void EnvGLNRec::recordResults(int child, const TransitionTable & tt)
{
    //m_resultrecord.push_back(_result);
    stringstream out1;  //child node id
    //convert int to string
    out1 << child + 1;
    string _result = out1.str();
    
    _result = _result + " " + m_gln.getNodeName(child);
    _result = _result + " " + m_gln.getNodeType(child);
    
    stringstream out2; //parent size
    out2 << tt.getParents().size();
    _result = _result + " " + out2.str() + " ";
    
    for(size_t j=0; j < tt.getParents().size(); j++) {
        //convert int to string
        stringstream out3; //parent node id
        out3 << tt.getParents()[j];
        _result = _result + out3.str() + ":" + m_gln.getNodeName(tt.getParents()[j]-1) + ",";
        //_result = _result + out3.str() +  ",";
    }
    
    stringstream out4;
    out4<< tt.getpValue();
    _result = _result + " " + out4.str();
    
    stringstream out5;
    out5<< tt.getChisq();
    _result = _result + " " + out5.str();
    
    stringstream out6;
    out6<<tt.getDf();
    _result = _result + " " + out6.str();
    
    m_recordresult[child].push_back(_result.c_str());
    m_recordcounts[child].push_back(tt.countsToString());
}

void EnvGLNRec::saveTopology(const string & fileTopology) const
{ // Save all interactions with p-value less than the given alpha level
    
	if(fileTopology.empty()) return;
    
    Topology topology( m_trajCol.getIntNodeNames(), 1);
    
	for(int i=0; i < (int) getGLNSize(); i++) {
        if (m_transTables[i].getpValue() <= m_gc.m_alpha) {
            topology.set(i, m_transTables[i].getParents());
        }
	}
	topology.save( fileTopology );
}


/*
//m_resultrecord.push_back(_result);
stringstream out1;  //child node id
//convert int to string
out1 << child + 1;
string _result = out1.str();

_result = _result + " " + m_gln.getNodes()[child].getName();
_result = _result + " " + m_gln.getNodes()[child].getType();

stringstream out2; //parent size
out2 << tt.getParents().size();
_result = _result + " " + out2.str() + " ";

for(size_t j=0; j < tt.getParents().size(); j++) {
    //convert int to string
    stringstream out3; //parent node id
    out3 << tt.getParents()[j];
    _result = _result + out3.str() + ":" + m_gln.getNodes()[(tt.getParents()[j]-1)].getName() + ",";
    //_result = _result + out3.str() +  ",";
}

stringstream out4;
out4<< tt.getpValue();
_result = _result + " " + out4.str();

stringstream out5;
out5<< tt.getChisq();
_result = _result + " " + out5.str();

stringstream out6;
out6<<tt.getDf();
_result = _result + " " + out6.str();

m_recordresult[child].push_back(_result.c_str());
m_recordcounts[child].push_back(tt.countsToString());
*/
