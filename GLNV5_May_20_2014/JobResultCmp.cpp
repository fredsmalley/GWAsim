// JobResultCmp.cpp -- The job result class for comparison based modeling
// 
// Joe Song
// Created: November 24, 2008


#include <iostream>
using std::endl;
using std::cerr;
using std::bad_alloc;

#include "EnvGLNCmp.h"
#include "JobGLN.h"
#include "JobResultCmp.h"

bool compare_chisq(double p_1, int numofp_1, double chisquare_1, 
				   double p_2, int numofp_2, double chisquare_2);

void JobResultCmp::extract(const Job & job, const JobEnvironment & env) 
{
	const EnvGLNCmp & envGLN = (const EnvGLNCmp &) env;
	int child = ((JobGLN &) job).getNodeId() - 1;

	m_diffTransTable = envGLN.getDiffTransTable( child );
	m_pooledTransTable = envGLN.getPooledTransTable( child );

	//added by Yang Zhang to implement parent selecting strategy
	m_TransTables = envGLN.getTransTables( child );
	m_chisq_t = envGLN.getChisq_t( child );
	m_df_t = envGLN.getDf_t( child );
	m_pchisq_t = envGLN.getpChisq_t( child );
	m_type = envGLN.getType( child );

	m_chisq_z = envGLN.getChisq_z( child );
	m_df_z = envGLN.getDf_z( child );
	m_pchisq_z = envGLN.getpChisq_z( child );

}

//void JobResultCmp::merge(JobEnvironment &env) const 
//{
//	EnvGLNCmp & envGLN = (EnvGLNCmp &) env;
//
//	int child = m_diffTransTable.getChild();
//
//	if( compare(m_diffTransTable, envGLN.getDiffTransTable( child ) ) )
//	{
//		envGLN.setDiffTransTable(child, m_diffTransTable);
//	}
//
//	child = m_pooledTransTable.getChild();
//	if( compare(m_pooledTransTable, envGLN.getPooledTransTable( child ) ) )
//	{
//		envGLN.setPooledTransTable(child, m_pooledTransTable);
//	}
//}


void JobResultCmp::merge(JobEnvironment &env) const 
{
	// bool improved = false;

	EnvGLNCmp & envGLN = (EnvGLNCmp &) env;

    int child = m_pooledTransTable.getChild();

    TransitionTable ttTot;
    ttTot.setChisq(m_chisq_t);
    ttTot.setDf(m_df_t);
    ttTot.setpValue(m_pchisq_t); 

    envGLN.select(child, m_TransTables, ttTot, m_pooledTransTable, m_diffTransTable);
    
    
    // bool improved = envGLN.compareParentSets(child, m_type, m_TransTables, ttTot,  m_pooledTransTable, m_diffTransTable);
  
/*    
   This part below needs to be consistent with the one used in EnvGLNCmp::parentSelect
 
	if(envGLN.getType(child)==NULL_INTX && m_type!= NULL_INTX) // if(envGLN.getType(child)==0&&m_type!=0)
	{
		improved = true;
        
	} 
    
	if(envGLN.getType(child)== REL_DIFF) // 1)
	{
		//if(type==1&&diff.getpValue()<=m_pooledTransTables[child].getpValue())
		if(m_type== REL_DIFF // 1 
           // && compare(m_diffTransTable, envGLN.getPooledTransTable( child )))
           && compare_chisq(m_diffTransTable.getpValue(), m_diffTransTable.getDf(), m_diffTransTable.getChisq(), 
                            envGLN.getPooledTransTable( child ).getpValue(), 
                            envGLN.getPooledTransTable( child ).getDf(),
                            envGLN.getPooledTransTable( child ).getChisq() ))
		{
			improved = true;
		}
		if(m_type==ABS_DIFF ||m_type==CONSERVED) // if(m_type==2||m_type==3)
		{
			improved = true;
		}
	}
	if(envGLN.getType(child)== ABS_DIFF) // 2)
	{
		if(m_type== ABS_DIFF // 2 
			&& compare_chisq(m_pchisq_t, 
							 (int) m_pooledTransTable.getParents().size(),
							 m_chisq_t,
							 envGLN.getpChisq_t( child ), 
							 (int) envGLN.getPooledTransTable( child ).getParents().size(),
							 envGLN.getChisq_t( child )))
		{
			improved = true;
		}
		if(m_type== CONSERVED // 3
			&& compare_chisq(m_pooledTransTable.getpValue(), 
							 (int) m_pooledTransTable.getParents().size(),
							 m_pooledTransTable.getChisq(),
							 envGLN.getpChisq_t( child ), 
							 (int) envGLN.getPooledTransTable( child ).getParents().size(),
							 envGLN.getChisq_t( child )))
		{
			improved = true;
		}
	}
	//so far best parents show a homogeneity interaction
	if(envGLN.getType(child)== CONSERVED) // 3)
	{
		if(m_type== CONSERVED // 3
			&& compare(m_pooledTransTable, envGLN.getPooledTransTable( child )))
		{
			improved = true;
		}
		if(m_type== ABS_DIFF // 2 
			&& compare_chisq(m_pchisq_t, (int) m_pooledTransTable.getParents().size(),m_chisq_t,
				envGLN.getPooledTransTable( child ).getpValue(), 
				(int) envGLN.getPooledTransTable( child ).getParents().size(),
				envGLN.getPooledTransTable( child ).getChisq()))
		{
			improved = true;
		}
	}
*/
    /*
	if(improved) {
        
		envGLN.setPooledTransTable( child,m_pooledTransTable );
		envGLN.setDiffTransTable( child, m_diffTransTable );
		envGLN.setTransTables( child , m_TransTables );
		
		envGLN.setChisq_t( child , m_chisq_t);
		envGLN.setDf_t( child , m_df_t );
		envGLN.setpChisq_t( child ,m_pchisq_t);

		envGLN.setType( child , m_type);

		envGLN.setChisq_z( child , m_chisq_z);
		envGLN.setDf_z( child , m_df_z );
		envGLN.setpChisq_z( child ,m_pchisq_z);
	}
    */
    
}

int JobResultCmp::pack(unsigned char* & buffer) const
{

	size_t ttbytes = 0;
	for(size_t i=0;i<m_TransTables.size();i++)
	{
		ttbytes += m_TransTables[i].numBytes();
	}

	size_t size = m_diffTransTable.numBytes() + m_pooledTransTable.numBytes() + sizeof(int)
				+ ttbytes + sizeof(double) + sizeof(int) + sizeof(double) + sizeof(int) 
				+ sizeof(double) + sizeof(int) + sizeof(double);

	try {

		buffer = new unsigned char [size];

	} catch (bad_alloc & ba) {

		cerr << "ERROR: Failed to allocate a buffer for job result!" << endl;
		cerr << "\t" << ba.what() << endl;
	}

	size_t actualSize = m_diffTransTable.pack(buffer, size);

	actualSize += m_pooledTransTable.pack(buffer + actualSize, size - actualSize);

	unsigned char * p = buffer;

	p = buffer + actualSize;
	*(int *) p = (int) m_TransTables.size();
	actualSize += sizeof(int);

	for(size_t i=0;i<m_TransTables.size();i++)
	{
		actualSize += m_TransTables[i].pack(buffer + actualSize, size - actualSize);
	}
	
	p = buffer + actualSize;

	*(double *) p = m_chisq_t;
	actualSize += sizeof(double);
	p += sizeof(double);

	*(int *) p = m_df_t;
	actualSize += sizeof(int);
	p += sizeof(int);

	*(double *) p = m_pchisq_t;
	actualSize += sizeof(double);
	p += sizeof(double);

	*(int *) p = m_type;
	actualSize += sizeof(int);
	p += sizeof(int);

	*(double *) p = m_chisq_z;
	actualSize += sizeof(double);
	p += sizeof(double);

	*(int *) p = m_df_z;
	actualSize += sizeof(int);
	p += sizeof(int);

	*(double *) p = m_pchisq_z;
	actualSize += sizeof(double);
	p += sizeof(double);

	return (int) actualSize;
}

void JobResultCmp::unpack(unsigned char * buffer, int length)
{
	unsigned char * p = buffer;	

	size_t actualSize = m_diffTransTable.unpack(buffer, length);

	actualSize += m_pooledTransTable.unpack(buffer + actualSize, length - actualSize);

	p = buffer + actualSize;
	int num_TransTable = *(int*) p;
	actualSize += sizeof(int);

	m_TransTables.resize(num_TransTable);

	for(size_t i=0;i<m_TransTables.size();i++)
	{
		actualSize += m_TransTables[i].unpack(buffer + actualSize, length - actualSize);
	}

	p = buffer + actualSize;
	m_chisq_t = *(double*) p;
	actualSize += sizeof(double);

	p = buffer + actualSize;
	m_df_t = *(int*) p;
	actualSize += sizeof(int);

	p = buffer + actualSize;
	m_pchisq_t = *(double*) p;
	actualSize += sizeof(double);

	p = buffer + actualSize;
	m_type = *(int*) p;
	actualSize += sizeof(int);

	p = buffer + actualSize;
	m_chisq_z = *(double*) p;
	actualSize += sizeof(double);

	p = buffer + actualSize;
	m_df_z = *(int*) p;
	actualSize += sizeof(int);

	p = buffer + actualSize;
	m_pchisq_z = *(double*) p;
	actualSize += sizeof(double);

}
