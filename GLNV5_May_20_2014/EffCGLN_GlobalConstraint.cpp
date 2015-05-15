#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "EffCGLN_GlobalConstraint.h"

EffCGLN_GlobalConstraint::EffCGLN_GlobalConstraint(void){}
EffCGLN_GlobalConstraint::~EffCGLN_GlobalConstraint(void){}


EffCGLN_GlobalConstraint::EffCGLN_GlobalConstraint( const string fileName )
{
	string token;
	ifstream in(fileName.c_str());

	if( ! in.is_open())
	{
		cerr<<"File open failed\n"<<endl;
		exit(-1);
	}


	if (in)
	{
		in >> token;
	}

	//Parse Global constraint file line by line
	while(in)
	{
		if (!isSectionName(token))
		{
			cerr<<"ERROR: file format is invalid";
			exit(-1);
		}

		string secName( token.substr(1, token.size() - 2) );

		if (secName == "Node")
		{
			while(in >> token)
			{
				if (isSectionName(token))
				{
					break;
				}
				else
				{
					m_nodes.push_back(token);
				}
			}
		}
		if (secName == "Radix")
		{
			while(in >> token)
			{
				if (isSectionName(token))
				{
					break;
				}
				else
				{
					m_radixs.push_back( stoi(token) );
				}
			}
		}
		else if (secName == "DynamicPattern")
		{
            EffCGLN_DynamicPattern tmp;
            tmp.m_DynamicPattern.push_back(vector<int>(m_nodes.size(), -1));

			while(in >> token)
			{
				if (isSectionName(token))
				{
					break;
				}
                else
                {
					while(true)
					{
						if (token == "DONE")
						{
							m_DynamicPatternList.push_back(tmp);
                            tmp.m_DynamicPattern.clear();
                            tmp.m_DynamicPattern.push_back(vector<int>(m_nodes.size(), -1));
                            break;
						}
                        else if(token == "*")
                        {
                            tmp.m_DynamicPattern.push_back(vector<int>(m_nodes.size(), -1));
                            in>>token;
                        }

						size_t nodeID = static_cast<size_t>(findNodeID(token));
                        assert( nodeID != static_cast<size_t>(-1) );

						in >> token;

						tmp.m_DynamicPattern.back()[nodeID] = stoi(token);
						in >> token;
					}
				}
			}
		}
	}//end of while(), reading in global constraint file
	
	in.close();
}


void EffCGLN_GlobalConstraint::generateGlobalConstraintTestFile()
{
	ofstream out("./GlobalConstraint.mod");

	out  << "[Node]" << endl;
	out  << "A	B	C" << endl;
	out  << "\n" << endl;
	out  << "[Radix]" << endl;
	out  << "2	2	2" << endl;
	out  << "\n" << endl;
	out  << "[DynamicPattern]" << endl;
	out  << "A	0" << endl;
	out  << "B	0" << endl;
	out  << "C	0" << endl;
	out  << "DONE" << endl;
	out  << "\n" << endl;
	out  << "A	1" << endl;
	out  << "B	1" << endl;
	out  << "C	1" << endl;
	out  << "DONE" << endl;
	out  << "\n" << endl;
	out  << "A	0" << endl;
	out  << "B	0" << endl;
	out  << "C	1" << endl;
	out  << "DONE" << endl;
	
	out.close();
}




void EffCGLN_GlobalConstraint::showDynamicPattern() const
{
	cout<<"Global constraint/Dynamic patten is: "<<endl;
	for(size_t x=0; x<m_DynamicPatternList.size(); x++)
	{
		for(size_t y=0; y<m_DynamicPatternList[x].m_DynamicPattern.size();y++)
		{
			for(size_t z=0; z<m_DynamicPatternList[x].m_DynamicPattern[y].size(); z++)
                if(m_DynamicPatternList[x].m_DynamicPattern[y][z] == -1)
				    cout<<"*";
			    else
				    cout<<m_DynamicPatternList[x].m_DynamicPattern[y][z];
            cout<<endl;
		}
		cout<<endl;
	}
	cout<<endl;
}


void EffCGLN_GlobalConstraint::saveGlobalConstraint(const string output_file_name) const
{
    ofstream ou( output_file_name.c_str() );

    ou<<"[Node]\n";
    for( const auto & name: m_nodes)
        ou<<name<<"\t";
    ou<<endl<<endl;

    ou<<"[Radix]\n";
    for( const auto & radix: m_radixs)
        ou<<radix<<"\t";
    ou<<endl<<endl;

    ou<<"[DynamicPattern]\n";
    for(size_t row=0; row<m_DynamicPatternList[0].m_DynamicPattern.size()-1; ++row)
    {
        for(size_t col=0; col<m_DynamicPatternList[0].m_DynamicPattern[row].size(); col++)
        {
            ou<<m_nodes[col]<<"\t"<<m_DynamicPatternList[0].m_DynamicPattern[row][col]<<endl;
        }
        ou<<"*\n";
    }

    //The last dynamic pattern is followed by "DONE" instread of "*"
    for(size_t col=0; col<m_DynamicPatternList[0].m_DynamicPattern.back().size(); col++)
    {
        ou<<m_nodes[col]<<"\t"<<m_DynamicPatternList[0].m_DynamicPattern.back()[col]<<endl;
    }
    ou<<"DONE";

    ou.close();
}



void EffCGLN_GlobalConstraint::saveAsTraColFile( const string & outputFileName) const
{
    int nodeNUM = m_DynamicPatternList.front().m_DynamicPattern.front().size();

    ofstream ou( outputFileName );

    ou<<"TRAJECTORY_VER2.1"<<endl;
    ou<<"1\t"<<nodeNUM<<"\t0"<<endl;
    for(int i=0; i<nodeNUM; ++i)
        ou<<"2\t";
    ou<<endl;
    
    for(int i=0; i<nodeNUM; ++i)
        ou<<"Node"+to_string(i+1)<<"\t";
    ou<<endl<<endl;
    
    for(int i=0; i<nodeNUM; ++i)
        ou<<"Yes\t";
    ou<<endl;

    for(int i=0; i<nodeNUM; ++i)
        ou<<"Yes\t";
    ou<<endl;

    ou<<m_DynamicPatternList.front().m_DynamicPattern.size()<<endl;
    for( const auto & row: m_DynamicPatternList.front().m_DynamicPattern)
    {
        for( const auto & bit: row)
            ou<<bit<<"\t";
        ou<<endl;
    }
    

    ou.close();
}

