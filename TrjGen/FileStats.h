/*
 * FileStats.h
 *
 *  Created on: Mar 31, 2015
 *      Author: nanderso
 */

#define myDebug 0

#ifndef FILESTATS_H_
#define FILESTATS_H_

#include <vector>
#include <sstream>

namespace std {

struct FileStats {
	string name;
	string fileName;
	long geneLocCount;
   vector<string> columns;
	vector<long> fileKeys;
	vector<string> participants;

	string toString(){
		stringstream output;
		output << "Total data points: " << fileKeys.size() << endl;
        
        if (myDebug){
            int i = 1;
            for (long key : fileKeys)
                output << i++ << "-" << key << " ";

            output << endl;
        }
        
        output << "Gene range columns: " << geneLocCount << endl;
        
        output << "Participant count: " << participants.size() << endl;
        
        if (myDebug){
            int i = 1;
            for (string name : participants)
                output << i++ << "-" << name << " ";
            
            output << endl;
        }
        
		return output.str();
	}
};
} /* namespace std */
#endif /* FILESTATS_H_ */
