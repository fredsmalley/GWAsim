#include "DS_ConstraintFile.h"


DS_ConstraintFile::DS_ConstraintFile(void)
{
}


DS_ConstraintFile::~DS_ConstraintFile(void)
{
}

bool DS_ConstraintFile::isSectionName(const string & token)
{
	if ( (token[0] != '[') || (token[token.size()-1] != ']'))
	{
		return false;
	}

	return true;
}


int DS_ConstraintFile::findNodeID(const string & nodeName)
{
	for (size_t i = 0; i < m_nodes.size(); i++)
	{
		// if (m_nodes[i].getName() == nodeName)
		if (m_nodes[i] == nodeName)
		{
			return (int)i;
		}
	}

	return -1;
}