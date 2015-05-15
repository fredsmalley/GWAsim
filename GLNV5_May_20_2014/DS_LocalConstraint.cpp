#include <cassert>
#include <algorithm>
using std::sort;

#include <sstream>
using std::istringstream;

#include <fstream>
using std::ofstream;

#include "GLNGlobals.h"

#include "DS_LocalConstraint.h"


DS_LocalConstraint::DS_LocalConstraint(void){}
DS_LocalConstraint::~DS_LocalConstraint(void){}

DS_LocalConstraint::DS_LocalConstraint( const string fileName )
{
    scan( fileName );
}

void DS_LocalConstraint::scan( const string & fileName )
{
	string token;
	ifstream in(fileName.c_str());

	if (in)
	{
		in >> token;
	}

	while(in)
	{
		if (!isSectionName(token))
		{
			GLNExit(EXIT_FAILURE, "ERROR: file format is invalid");
		}

		string secName(token.substr(1, token.size() - 2));

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
					// m_nodes.push_back(Node());
					// m_nodes.back().setName(token.c_str());
				}
			}

			m_compulsiveParents.resize(m_nodes.size());

		}
		else if (secName == "Radix")
		{
			while(in >> token)
			{
				if (isSectionName(token))
				{
					break;
				}
				else
				{
					istringstream tmp(token);
					int intV;
					tmp >> intV;
					m_radixs.push_back(intV);
				}
			}
		}
		else if (secName == "InitialValue")
		{
			while(in >> token)
			{
				if (isSectionName(token))
				{
					break;
				}
				else
				{
					istringstream tmp(token);
					int intV;
					tmp >> intV;
					m_initValue.push_back(intV);
				}
			}

			m_compulsiveParents.resize(m_nodes.size());

		}
        else if( secName == "ExcludedParent" )
        {
            while(in >> token)
            {
                if( isSectionName( token ) )
                {
                    break;
                }
                else
                {//Store name of nodes who can not be parent
                    m_ExcludedParent.push_back( token );
                }
            }
        }
        else if( secName == "ExcludedChild" )
        {
            while(in >> token)
            {
                if( isSectionName( token ) )
                {
                    break;
                }
                else
                {//Store name of nodes who can not have parent
                    m_ExcludedChild.push_back( token );
                }
            }
        }
		else if (secName == "CompulsiveInteraction")
		{
			while(in >> token)
			{
				if (isSectionName(token))
				{
					break;
				}
				else
				{
					//nodes' ID start from 1 but program count from 0
        		                // findNodeID() will return 0-base index
					int parentID = findNodeID(token);

					in >> token;
					size_t childID = static_cast<size_t>(findNodeID(token));

					in >> token;
					//int relation = 0;
                    double relation = 0;
					if (token == "+")
					{
						relation = 1.0;
					}
					else if (token == "-")
					{
						relation = -1.0;
					}
                    else  //relation is now double, Oct 28, 2012
                    {
                        relation = stod( token );
                    }

					m_compulsiveParents[childID].m_nodeList.push_back(parentID + 1); //parent's index 1-based.
					m_compulsiveParents[childID].m_polarityRelations.push_back(relation);
				}
			}
		}
	

	}

	in.close();

	if (!m_initValue.empty())
	{
		if(m_nodes.size() != m_initValue.size())
		{
			GLNExit(EXIT_FAILURE, "ERROR: The length of [Node] and [InitialValue] should be equal.");
		}
	}
}

valarray<int> DS_LocalConstraint::convert(const vector<int> & vec)
{
        valarray<int> ret(vec.size());

        for (size_t i = 0; i < vec.size(); i++)
        {
                ret[i] = vec[i];
        }

        return ret;
}



float DS_LocalConstraint::correlation( const vector<int> & x, const vector<int> & y )
{
	return correlation(convert(x), convert(y));
}

float DS_LocalConstraint::correlation(const valarray<int> & x, const valarray<int> & y)
{
	return signCorrelation(x, y) / sqrt(expectedValue(x * x) - pow(expectedValue(x), 2)) * sqrt(expectedValue(y * y) - pow(expectedValue(y), 2));
}


float DS_LocalConstraint::signCorrelation(const valarray<int> & x, const valarray<int> & y)
{
	return (expectedValue(x * y) - expectedValue(x) * expectedValue(y));
}

float DS_LocalConstraint::signCorrelation( const vector<int> & x, const vector<int> & y )
{
	return signCorrelation(convert(x), convert(y));
}


float DS_LocalConstraint::expectedValue(const valarray<int> & v)
{
	return static_cast<float>(v.sum())/static_cast<float>(v.size());
}



void DS_LocalConstraint::generateLocalConstraintTestFile()
{
	ofstream out("./LocalConstraint.mod");

	out  << "[Node]" << endl;
	out  << "A	B	T" << endl;
	out  << "\n" << endl;
	out  << "[Radix]" << endl;
	out  << "2	2	2" << endl;
	out  << "\n" << endl;
//	out  << "[InitialValue]" << endl;
//	out  << "1	1	1" << endl;
	out  << "\n" << endl;
	out  << "[CompulsiveInteraction]" << endl;
	out  << "T	A	-" << endl;
	out  << "T	B	+" << endl;
	out  << "\n" << endl;
	out  << "[MaxParentsForEnumeration]" << endl;
	out  << "2" << endl;
	out  << "\n" << endl;
	out  << "[Counter]" << endl;
	out  << "COUNTER" << endl;
	out  << "\n" << endl;
	out  << "[RecognizeState]" << endl;
	out  << "T	1" << endl;
	out  << "\n" << endl;
	out  << "[Feedback]" << endl;
	out  << "T	1" << endl;

	out.close();
}

void DS_LocalConstraint::generateLocalConstraintTestGLNFile()
{
	ofstream out("./LocalConstraintTestGLN.txt");

	out << "GENERALIZED_LOGICAL_NETWORK_VER1" << endl;
	out << "3" << endl;
	out << "1\nA\n2\n1\ni\n2" << endl;
	out << "2 3" << endl;
	out << "-1 -1" << endl;
	out << "2 2" << endl;
	out << "1 1 1 0" << endl;
	out << "*" << endl;
	out << "2\nB\n2\n1\ni\n1" << endl;
	out << "3" << endl;
	out << "-1" << endl;
	out << "2" << endl;
	out << "0 1" << endl;
	out << "*" << endl;
	out << "3\nT\n2\n1\ni\n1" << endl;
	out << "1" << endl;
	out << "-1" << endl;
	out << "2" << endl;
	out << "0 1" << endl;
	out << "*" << endl;
	out << "done" << endl;

	out.close();
}




//Added by Haizhou Wang, Oct 15, 2012
void DS_LocalConstraint::saveLocalConstraint(const string output_file_name) const
{
    ofstream ou(output_file_name.c_str());
    ou<<"[Node]\n";
    for( const auto & name: m_nodes)
        ou<<name<<"\t";
    ou<<endl<<endl;

    ou<<"[Radix]\n";
    for( const auto & radix: m_radixs)
        ou<<radix<<"\t";
    ou<<endl<<endl;

    ou<<"[CompulsiveInteraction]\n";
    for(size_t child=0; child<m_compulsiveParents.size(); ++child)
    {
        assert( m_compulsiveParents[child].m_nodeList.size() == m_compulsiveParents[child].m_polarityRelations.size() );
        for(size_t parent=0; parent<m_compulsiveParents[child].m_nodeList.size(); ++parent)
        {
            ou<<m_nodes[m_compulsiveParents[child].m_nodeList[parent]-1]<<"\t"<<m_nodes[child]<<"\t";
            
            ou<<m_compulsiveParents[child].m_polarityRelations[parent]<<endl;
            /*
            if( m_compulsiveParents[child].m_polarityRelations[parent] == 1 )
                ou<<"+"<<endl;
            else if( m_compulsiveParents[child].m_polarityRelations[parent] == -1 )
                ou<<"-"<<endl;
            else
            {
                cerr<<"Error! When output relationship between parent and child, met with unknown relationship\n"<<endl;
                exit(-1);
            }
            */
        }
    }
    ou<<endl;

    ou.close();
}
//End of Haizhou's Code











void DS_LocalConstraint::saveConstraints(const string output_file_name) const
{
    ofstream ou(output_file_name.c_str());

    for(size_t child=0; child<m_compulsiveParents.size(); ++child)
    {
        assert( m_compulsiveParents[child].m_nodeList.size() == m_compulsiveParents[child].m_polarityRelations.size() );
        for(size_t parent=0; parent<m_compulsiveParents[child].m_nodeList.size(); ++parent)
        {
            ou<<m_nodes[m_compulsiveParents[child].m_nodeList[parent]-1]<<"\t"<<m_nodes[child]<<"\t";
            
            ou<<m_compulsiveParents[child].m_polarityRelations[parent]<<endl;
            /*
            if( m_compulsiveParents[child].m_polarityRelations[parent] == 1 )
                ou<<"+"<<endl;
            else if( m_compulsiveParents[child].m_polarityRelations[parent] == -1 )
                ou<<"-"<<endl;
            else
            {
                cerr<<"Error! When output relationship between parent and child, met with unknown relationship\n"<<endl;
                exit(-1);
            }
            */
        }
    }
    ou<<endl;

    ou.close();
}
//End of Haizhou's Code



























//Added by Haizhou Wang, Jan 23, 2013
int DS_LocalConstraint::getMaxNumParents() const
{
    int result = -1;

    for(const auto & node: m_compulsiveParents)
    {
        if( (int)node.m_polarityRelations.size() > result )
            result = node.m_polarityRelations.size();
    }

    return result;
}



void DS_LocalConstraint::saveAsTopologyFile( const string & out_put_file_name) const
{
    ofstream ou( out_put_file_name);

    ou<<"Local Constraint File\n";

    int numOfLines = 0;
    for(const auto & i: m_compulsiveParents)
        if( i.m_polarityRelations.size() != 0 )
            ++numOfLines;

    ou<<numOfLines<<endl;

    for(size_t i=0; i<m_nodes.size(); ++i)
    {
        if( m_compulsiveParents[i].m_nodeList.size() != 0 )
        {
            ou<<m_nodes[i]<<"*";
            for(const auto & parent: m_compulsiveParents[i].m_nodeList)
            {
                ou<<m_nodes[parent-1]<<",";  //parent is 1-based
            }
            ou<<endl;
        }
    }

    ou.close();
}



void DS_LocalConstraint::saveLocalConstraintPair(const string output) const
{
     ofstream ou( output );

    for(size_t i=0; i<m_nodes.size(); ++i)
    {
        if( m_compulsiveParents[i].m_nodeList.size() != 0 )
        {
            for(const auto & parent: m_compulsiveParents[i].m_nodeList)
            {
                ou<<m_nodes[i]<<" "<<m_nodes[parent-1]<<endl;  //'parent' is 1-based
            }
        }
    }

    ou.close();
}

