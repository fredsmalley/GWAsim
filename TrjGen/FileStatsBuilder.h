/*
 * FileStatsBuilder.h
 *
 *  Created on: Mar 31, 2015
 *      Author: nanderso
 */

#ifndef FILESTATSBUILDER_H_
#define FILESTATSBUILDER_H_

#include <fstream>
#include "FileStats.h"

namespace std {

class FileStatsBuilder {
public:
    FileStatsBuilder(string file);
    virtual ~FileStatsBuilder(void);
    
    void generateFileStats(FileStats& fs);

private:
    string fileName;
    
    bool gotoNextComment(ifstream &ifs);

};

} /* namespace std */

#endif /* FILESTATSBUILDER_H_ */
