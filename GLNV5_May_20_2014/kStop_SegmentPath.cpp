#include "kStop_SegmentPath.h"


#include <iostream>
using std::cout;
using std::endl;


#include <cassert>

using std::vector;

//Given segment ID, check to see if all vertices (except the first and last one) are marked with this ID
//If not, it indicates that such a vertex is also used by other segment.
bool SegmentPath::pure(const vector<kStop_Node> & G,  int segID) const
{
    if( m_SegmentPath.size() == 2 ) //Source state is connected directly to target state, no need to check purity
        return true;

    assert( m_SegmentPath.size() > 2 );
    for(unsigned t=1; t<m_SegmentPath.size()-1; ++t) //Note that it doesn't check the first and last node
    {
        if( G[ m_SegmentPath[t] ].m_mark != segID )
            return false;
    }
    
    return true;
}


void SegmentPath::showSegmentpath() const
{
    for(unsigned t=0; t<m_SegmentPath.size(); t++)
        cout<<m_SegmentPath[t]<<" "<<endl;

    cout<<endl;
}
