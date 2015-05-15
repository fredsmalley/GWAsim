/*
 * FileStatsBuilder.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: nanderso
 */

#include "FileStatsBuilder.h"
#include <iostream>
#include <sstream>
#include <limits>

namespace std{

FileStatsBuilder::FileStatsBuilder(string file) : fileName(file) { }

FileStatsBuilder::~FileStatsBuilder() { }

bool FileStatsBuilder::gotoNextComment(ifstream &ifs){
    char ch;
    
    do {
        ifs >> ch;
    } while (ch != '#' );
    
    if (ch == '#' && ifs.peek() == '#'){
        ifs >> ch;
        return true;
    }
    
    return false;
}

void FileStatsBuilder::generateFileStats(FileStats& fs) {
    ifstream ifs(fileName);
    char ch;
    string buffer;
    stringstream ss1;
    stringstream ss2;
    long key;
    long loc1;
    long loc2;
    
    // Look through all header comments
    while (gotoNextComment(ifs)){
        // Read comment line
        while (ifs.peek() != '\n'){
            ifs >> buffer;
            if (buffer.compare("Positions:") == 0){
               ifs >> buffer;
               ss1.clear();
               ss1 << buffer;
               ss1 >> fs.geneLocCount;
            }
        }
    }
    
    // Loop until end of line getting location
    do{
        ifs >> buffer;
        fs.columns.emplace_back(buffer);
    } while (ifs.peek() != '\n');
    
    // Loop through making file keys
    for (string s : fs.columns){
        ss1.clear();
        ss2.clear();
        
        ss1 << s.substr(0, s.find("_"));
        ss2 << s.substr(s.find("_") + 1, s.length());
        
        ss1 >> loc1;
        ss2 >> loc2;
        key = loc2 - loc1;
        
        fs.fileKeys.emplace_back(key);
    }
    
    // Loop to count participants
    do{
        ifs >> buffer;
        
        // Break to prevent the double count at end of file
        if (ifs.peek() == EOF || buffer.compare("") == 0)
            break;
        
        if (isalpha(buffer.at(0))){
            fs.participants.emplace_back(buffer);
        }
        // Loop until end of line or end of file
        while (ifs.peek() != '\n' && ifs.peek() != EOF) {
            ifs >> ch;
        }
    } while (ifs.peek() != EOF);
	
	cout << fs.toString();

    ifs.close();
    return;
}
} /* namespace std */