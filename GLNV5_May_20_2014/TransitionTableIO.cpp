//
//  TransitionTableIO.cpp
//  gln
//
//  Created by Joe Song on 5/21/13.
//
//  Updated:
//    Feb 3, 2014. MS. Added two general member functions scan(ifstream)
//      and save(ofstream)

#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include "TransitionTable.h"

bool TransitionTable::scan4by2(const string & file)
{
    bool success = true;

    // parents
    vector<int> parents(1);
    parents[0]=1;//, parents[1]=2;
    
    vector<int> parentBases(1);
    parentBases[0]=0;//, parentBases[1]=2;
    
    vector<string> nameParents(1);
    nameParents[0]="P1";//, nameParents[1]="P2";
    
    vector<int> delays(1, 0);
    
    // child
    int child = 0, base = 0;
    string nameChild="Child";

    Transition::initialize(child, base, parents, parentBases, delays, nameChild, nameParents);

    allocate();

    m_transitionTable.clear();
    
    ifstream ifs(file);
    
    if (ifs.is_open()) {
        
        string line;
        while(getline(ifs, line)){
            istringstream iss(line);
            int n;
            std::vector<int> v;
            while (iss >> n)
            {
                v.push_back(n);
            }
            m_transitionTable.push_back(v);
        }
        
    } else {
        success = false;
    }
    
    parentBases[0] = m_transitionTable.size();
    base = m_transitionTable[0].size();
    
    return success;
}

void TransitionTable::saveTable(const string & file) const
{
    ofstream ofs(file);
    
    if (ofs.is_open()) {
        
        for (size_t i=0; i<m_transitionTable.size(); i++) {
            for (size_t j=0; j<m_transitionTable[i].size(); j++) {
                if(j != 0) ofs << "\t";
                ofs << m_transitionTable[i][j];
            }
            ofs << endl;
        }
    } else {
        cerr << "ERROR: saving transition table" << endl;
    }
}

template <class T>
void read_and_split_line(ifstream & ifs, vector<T> & v)
{
    string line;

    getline(ifs, line);
    istringstream iss(line);
    T element;

    while (iss >> element) {
        v.push_back(element);
    }
}

bool TransitionTable::scan(ifstream & ifs)
{
    if( ! ifs.is_open() ) return false;
    
    ifs >> m_child;
    ifs >> m_nameChild;
    ifs >> m_base;
    
    string line;
    getline(ifs, line);
    
    read_and_split_line(ifs, m_parents);
    read_and_split_line(ifs, m_parentBases);
    read_and_split_line(ifs, m_nameParents);
    read_and_split_line(ifs, m_delays);
    
    allocate();

    for (size_t i=0; i<m_transitionTable.size(); i++) {
        for (size_t j=0; j<m_transitionTable[i].size(); j++) {
            ifs >> m_transitionTable[i][j];
        }
    }
    
    return true;
}

void TransitionTable::save(ofstream & ofs) const
{
    if (ofs.is_open()) {
 
        ofs << m_child << '\t' << m_nameChild << '\t' << m_base << endl;
        
        for (size_t m=0; m<m_parents.size(); m++) {
            ofs << m_parents[m];
            if (m < m_parents.size()-1) {
                ofs << '\t';
            } else {
                ofs << endl;
            }
        }

        for (size_t m=0; m<m_parentBases.size(); m++) {
            ofs << m_parentBases[m];
            if (m < m_parentBases.size()-1) {
                ofs << '\t';
            } else {
                ofs << endl;
            }
        }
        
        for (size_t m=0; m<m_nameParents.size(); m++) {
            ofs << m_nameParents[m];
            if (m < m_nameParents.size()-1) {
                ofs << '\t';
            } else {
                ofs << endl;
            }
        }

        for (size_t m=0; m<m_delays.size(); m++) {
            ofs << m_delays[m];
            if (m < m_delays.size()-1) {
                ofs << '\t';
            } else {
                ofs << endl;
            }
        }

        for (size_t i=0; i<m_transitionTable.size(); i++) {
            for (size_t j=0; j<m_transitionTable[i].size(); j++) {
                if(j != 0) ofs << "\t";
                ofs << m_transitionTable[i][j];
            }
            ofs << endl;
        }
    } else {
        cerr << "ERROR: saving transition table" << endl;
    }
}

void TransitionTable::saveDOT(ofstream & ofs) const
{
    if (ofs.is_open()) {
        for (size_t m=0; m<m_nameParents.size(); m++) {
            ofs << "\"" << m_nameParents[m] << "\"" << "->" << "\""
            << m_nameChild << "\""<< ';' << endl;
        }
    } else {
        cerr << "ERROR: saving transition table" << endl;
    }
}
