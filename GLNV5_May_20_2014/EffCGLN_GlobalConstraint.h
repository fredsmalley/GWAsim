#pragma once

#include <string>
#include <vector>
using std::vector;

#include "DS_ConstraintFile.h"
#include "EffCGLN_DynamicPattern.h"

class EffCGLN_GlobalConstraint : public DS_ConstraintFile
{
public:
	EffCGLN_GlobalConstraint(void);
	~EffCGLN_GlobalConstraint(void);

	EffCGLN_GlobalConstraint(const string fileName);

	static void generateGlobalConstraintTestFile();
	void showDynamicPattern() const;
    void saveGlobalConstraint(const string output_file_name) const;

public:
	vector< EffCGLN_DynamicPattern > m_DynamicPatternList;

    void saveAsTraColFile( const string & outputFileName) const;
};
