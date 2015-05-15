#pragma once

#include <vector>
using std::vector;

#include <string>
using std::string;

class DS_ConstraintFile
{
public:
	DS_ConstraintFile(void);
	~DS_ConstraintFile(void);

	bool isSectionName(const string & token);
	int findNodeID(const string & nodeName);

public:
	vector<string> m_nodes;
	vector<int> m_radixs;
};

