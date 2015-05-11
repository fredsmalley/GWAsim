/*
 * TableBuilder.cpp
 *
 *  Created on: Apr 10, 2015
 *      Author: nanderso
 */

#include "TrjBuilder.h"
#include "FileStatsBuilder.h"

#include <iostream>
#include <algorithm> // replace
#include <string>

#define debug 0
namespace std {

TrjBuilder::TrjBuilder(string inFile, string outFile) {
    FileStatsBuilder fsb(inFile);
    fsb.generateFileStats(fs);
    inputFileName = inFile;
    outputFileName = outFile;
}

TrjBuilder::~TrjBuilder() { }

void TrjBuilder::run(vector<string> &children){
    ifstream ifs(inputFileName);
    ofstream ofs(outputFileName);
    readBuf.resize(fs.columns.size());
    
    // Find locations of child nodes.
    for (string str : children){
        for (long i = 0; i < fs.columns.size(); ++i){
            if (fs.columns.at(i).compare(str) == 0)
                childLocations.emplace_back(i);
        }
    }
    
    // move input file to start of table
    seekBody(ifs);
    
    writeHeader(ofs);
    
    // Compute tables.
    for (long i = 0; i < fs.participants.size(); ++i){
        readLine(ifs);
        writeLine(ofs);
    }
    
    ifs.close();
    ofs.close();
}
    
void TrjBuilder::seekBody(ifstream &ifs){
    char ch;
    bool inComment = true;
    
    // Get to start of the body
    while (inComment){
        do {
            ifs >> ch;
        } while (ch != '#' );
        
        if (ch == '#' && ifs.peek() == '#'){
            ifs >> ch;
            
            while (ifs.peek() != '\n'){
                ifs >> ch;
            }
        } else {
            inComment = false;
        }  
    }
    
    // Skip the column titles
    while (ifs.peek() != '\n'){
        ifs >> ch;
    }
}
  
void TrjBuilder::readLine(ifstream &ifs) {
    string buffer;
    // Skip participent name
    ifs >> buffer;
    
    // Buffer the gene data string
    for (long i = 0; i < fs.columns.size(); ++i){
        ifs >> buffer;
        readBuf[i] = buffer;
    }
}

void TrjBuilder::writeLine(ofstream &ofs) {
    // Write row count to file
    ofs << endl << endl << "1" << endl;
    
    // Write child nodes first
    if (readBuf.at(fs.columns.size()-1).at(0) == '~'){
		ofs << "0 ";
        cout << readBuf.at(fs.columns.size()-1) << endl;
    } else {
		ofs << "1 ";
        cout << readBuf.at(fs.columns.size()-1) << endl;
    }
    
   for (long j = 0; j < fs.columns.size() - 3; ++j){
      ofs << getTableValue(readBuf.at(j));
   }
    
   if (readBuf.at(fs.columns.size()-3).at(0) == '~')
		ofs << "0 ";
	else
		ofs << "1 ";

	if (readBuf.at(fs.columns.size()-2).at(0) == 'W' || readBuf.at(fs.columns.size()-2).at(0) == 'C')
		ofs << "0";
	else if (readBuf.at(fs.columns.size()-2).at(0) == '~')
		ofs << "-1";
	else
		ofs << "1";
}
    
void TrjBuilder::writeHeader(ofstream &ofs) {
   ofs << "TRAJECTORY_VER2.1" << endl;
   ofs << fs.participants.size() << " " << fs.columns.size() << " " << "0" << endl;
   ofs << "2 ";
	// For homozygous, heterozygous or reference.
	for (int i = 0; i < fs.geneLocCount; ++i)
		ofs << "3 ";

	ofs << "2 " << "2" << endl << endl; // For south west parameter and race
	
	ofs << fs.columns.at(fs.columns.size()-1) << " ";
	// List gene locations
	for (unsigned int i = 0; i < fs.columns.size()-1; ++i)
		ofs << fs.columns.at(i) << " ";
	ofs.seekp((long)ofs.tellp() - 1); // Take off last space
	
	ofs << endl << "No ";
	for (unsigned int i = 0; i < fs.columns.size() - 1; ++i)
		ofs << "Yes ";
	ofs.seekp((long)ofs.tellp() - 1); // Take off last space

	ofs << endl << "Yes ";
	for (unsigned int i = 0; i < fs.columns.size() - 1; ++i)
		ofs << "No ";
	ofs.seekp((long)ofs.tellp() - 1); // Take off last space

}

string TrjBuilder::getTableValue(string buffer){
   //cout << buffer << endl;
   string strArray[20];
   replace(buffer.begin(), buffer.end(), ':', ' ');
   replace(buffer.begin(), buffer.end(), '.', ' ');
   replace(buffer.begin(), buffer.end(), ',', ' ');

   stringstream ss(buffer);
   string outstr;

   for (int i = 0; i < 5; ++i)
	   ss >> strArray[i];

   if (debug){
	   string test = strArray[1] + strArray[2] + strArray[3] + strArray[4];
	   if (test.length() != 0 && test.length() != 4)
		   cout << " " << "Gene data not standard! = " << strArray[0] << " " << test << " ";
   }

   switch (strArray[0].at(0)){
   case 'n':
   case 'j':	outstr = "-1 ";
			   break;

   case '~':	outstr = "0 ";
			   break;

   case 's':	if (strArray[1].compare(strArray[2]) != 0 && strArray[3].compare(strArray[4]) == 0)
				   outstr = "1 "; // Heterozygous
			   else if (strArray[1].compare(strArray[2]) == 0 && strArray[3].compare(strArray[4]) != 0)
				   outstr = "1 "; // Heterozygous
			   else if (strArray[1].compare(strArray[2]) != 0 && strArray[3].compare(strArray[4]) != 0)
				   outstr = "2 "; // Homozygous
			   else
				   outstr = "0 ";

			   if (debug)
				   cout << buffer << " ";
			   break;

   case 'd':	if (strArray[1].compare(strArray[2]) != 0 && strArray[3].compare(strArray[4]) == 0)
				   outstr = "1 "; // Heterozygous
			   else if (strArray[1].compare(strArray[2]) == 0 && strArray[3].compare(strArray[4]) != 0)
				   outstr = "1 "; // Heterozygous
			   else if (strArray[1].compare(strArray[2]) != 0 && strArray[3].compare(strArray[4]) != 0)
				   outstr = "2 "; // Homozygous
			   else
				   outstr = "0 ";

			   if (debug)
				   cout << buffer << " ";
			   break;

   case 'r':	outstr = "0 ";
			   break;

   case 'S':	outstr = "1 ";
			   break;

   case 'i':	outstr = "1 ";
			   break;

   default:   cout << endl << buffer << endl;
   }
   return outstr;
}
}
