/*
 * TrjBuilder.h
 *
 *  Created on: Apr 10, 2015
 *      Author: nanderso
 */

#ifndef TABLEBUILDER_H_
#define TABLEBUILDER_H_

#include "FileStats.h"
#include <fstream>

namespace std {

class TrjBuilder {
public:
    TrjBuilder(string inFileName, string outFileName);
	virtual ~TrjBuilder(void);
    
    void run(vector<string> &children);
    
private:
    string inputFileName;
    string outputFileName;
    FileStats fs;
    vector<string> readBuf;
    vector<long> childLocations;
    
    void seekBody(ifstream &ifs);
    
    void readLine(ifstream &ifs);
    
    void writeLine(ofstream &ofs);

    void writeHeader(ofstream &ofs);
    
    string getTableValue(string buffer);
};

} /* namespace std */

#endif /*  TABLEBUILDER_H_ */