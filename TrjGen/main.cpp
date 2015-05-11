/*
 * main.cpp
 *
 *  Created on: Mar 24, 2015
 *      Author: nanderso
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <limits>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <ctype.h>
#include "TrjBuilder.h"

#define debug 0

#define inFileName "TRPM8largeTable.tsv"
// /Users/nathananderson/Documents/fredSmalley/GeneVar/files/TRPM8largeTable.tsv"
#define outFileName "TRPM8largeTable.trj"
///Users/nathananderson/Documents/TrjGen/TRPM8largeTable.trj"

using namespace std;

int main(int argc, char** argv){
	string buffer;
	string strArray[5];
	vector<string> strBuf;
	vector<int> filekey;

	//ifstream ifs (inFileName, ifstream::in);
	//ofstream ofs (outFileName, ofstream::out);
    
    vector<string> children;
    
    TrjBuilder tb(inFileName, outFileName);
    
    tb.run(children);
    
	//FileStatsBuilder fsb(inFileName);

	//FileStats fs;
	//fsb.generateFileStats(fs);

	//ofs << "TRAJECTORY_VER2.1" << endl;

    //ifs.close();
    //ofs.close();
    
    return 0;
    /*ifs.close();
	ifs.open(inFileName, ifstream::in);

	ifs.ignore(numeric_limits<streamsize>::max(), '\t'); // clear out row name

	while (ifs.peek() != EOF && ifs.peek() != '\n'){
		ifs >> buffer;
	    strBuf.emplace_back(buffer);//.substr(buffer.find("_") + 1, buffer.length())); // trim column name to just the trailing number and store in vector.
	}

	cout << endl << "number of datapoint " << strBuf.size() << endl;
	int j = 1;
	for (string s : strBuf)
		cout << j++ << "-"<< s << " ";

	//ofs << fs.participants.size() << " " << fs.fileKeys.size() << " " << "0" << endl;

	ofs << "2 ";
	// For homozygous, heterozygous or reference.
	for (int i = 0; i < strBuf.size() - 3; ++i)
		ofs << "3 ";

	ofs << "2 " << "2" << endl << endl; // For south west parameter and race

	ofs << strBuf.at(strBuf.size()-1) << " ";
	// List gene locations
	for (unsigned int i = 0; i < strBuf.size()-1; ++i)
		ofs << strBuf.at(i) << " ";
	ofs.seekp((long)ofs.tellp() - 1); // Take off last space

	ofs << endl << "No ";
	for (unsigned int i = 0; i < strBuf.size() - 1; ++i)
		ofs << "Yes ";
	ofs.seekp((long)ofs.tellp() - 1); // Take off last space

	ofs << endl << "Yes ";
	for (unsigned int i = 0; i < strBuf.size() - 1; ++i)
		ofs << "No ";
	ofs.seekp((long)ofs.tellp() - 1); // Take off last space

	ofs  << endl << endl << "1" << endl;

	string out;
	if (debug){
		ifs >> buffer;
		cout << endl << buffer << endl;
		ofs << endl << endl << buffer;
	} else {
		ifs >> out;//ifs.ignore(numeric_limits<streamsize>::max(), '\t'); // clear out row name
	}

	unsigned int column = -1;
	unsigned int row = 0;
	stringstream outss;

	while (ifs.peek() != EOF){
		if (ifs.peek() == '\n'){
			if (debug){
				ifs >> buffer;
				if (ifs.peek() != EOF){
					cout << endl << buffer << endl;
					outss << endl << endl << buffer;
				}
			} else {
				ifs >> out;//ifs.ignore(numeric_limits<streamsize>::max(), '\t'); // clear out row name
			}

				outss.seekp((long)outss.tellp() - 1); // Take off last space
				if (ifs.peek() != EOF)
					outss << endl << endl << "1" << endl;

				while (!outss.eof()){
					//outss >> out;
					switch(outss.peek()){
					case ' ': ofs << " "; break;
					case '\n': ofs << "\n"; break;
					case '-': ofs << "-"; break;
					default: if (isdigit(outss.peek())){ ofs << (char)outss.peek();} break;
					}
					outss.ignore();
				}
				outss.clear();
			//}
			column = -1;
		} else {
			if (column++ == strBuf.size() - 4){

				ifs >> buffer;

				if (buffer.at(0) == '~')
					outss << "0 ";
				else
					outss << "1 ";

				ifs >> buffer;
				if (buffer.at(0) == 'W' || buffer.at(0) == 'C')
					outss << "0 ";
				else if (buffer.at(0) == '~')
					outss << "-1 ";
				else
					outss << "1 ";

				ifs >> buffer;

				if (buffer.at(0) == '~')
					ofs << "0 ";
				else
					ofs << "1 ";

				outss << endl;

				column = -1;
			} else {
				ifs >> buffer;
				replace(buffer.begin(), buffer.end(), ':', ' ');
				replace(buffer.begin(), buffer.end(), '.', ' ');
				replace(buffer.begin(), buffer.end(), ',', ' ');

				stringstream ss(buffer);

				for (int i = 0; i < 5; ++i)
					ss >> strArray[i];

				if (debug){
					string test = strArray[1] + strArray[2] + strArray[3] + strArray[4];
					if (test.length() != 0 && test.length() != 4)
						cout << " " << "Gene data not standard! = " << strArray[0] << " " << test << " ";
				}

				switch (strArray[0].at(0)){
				case 'n':
				case 'j':	outss << "-1 ";
							break;

				case '~':	outss << "0 ";
							break;

				case 's':	if (strcmp(strArray[1].c_str(), strArray[2].c_str()) != 0 && strcmp(strArray[3].c_str(), strArray[4].c_str()) == 0)
								outss << "1 "; // Heterozygous
							else if (strcmp(strArray[1].c_str(), strArray[2].c_str()) == 0 && strcmp(strArray[3].c_str(), strArray[4].c_str()) != 0)
								outss << "1 "; // Heterozygous
							else if (strcmp(strArray[1].c_str(), strArray[2].c_str()) != 0 && strcmp(strArray[3].c_str(), strArray[4].c_str()) != 0)
								outss << "2 "; // Homozygous
							else
								outss << "0 ";

							if (debug)
								cout << buffer << " ";
							break;

				case 'd':	if (strcmp(strArray[1].c_str(), strArray[2].c_str()) != 0 && strcmp(strArray[3].c_str(), strArray[4].c_str()) == 0)
								outss << "1 "; // Heterozygous
							else if (strcmp(strArray[1].c_str(), strArray[2].c_str()) == 0 && strcmp(strArray[3].c_str(), strArray[4].c_str()) != 0)
								outss << "1 "; // Heterozygous
							else if (strcmp(strArray[1].c_str(), strArray[2].c_str()) != 0 && strcmp(strArray[3].c_str(), strArray[4].c_str()) != 0)
								outss << "2 "; // Homozygous
							else
								outss << "0 ";

							if (debug)
								cout << buffer << " ";
							break;

				case 'r':	outss << "0 ";
							break;

				case 'S':	outss << "1 ";
							break;

				case 'i':	outss << "1 ";
							break;

				default:   cout << endl << buffer << endl;
				}
			}
		//cout << buffer << " ";
		}
	}

	ifs.close();
	ofs.close();
	if (debug)
		cout << strBuf.size();
*/
	return 0;
}



