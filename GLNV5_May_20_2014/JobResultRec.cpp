// JobResultRec.cpp
// 
// Curtis Luce, Joe Song
// Created: June 2007
// Last modified: September 27, 2008.  Joe Song. Name changed from JobResultLognet.cpp

#include <iostream>
using std::cout;
using std::endl;
using std::bad_alloc;
using std::cerr;

#include "EnvGLNRec.h"
#include "JobGLN.h"
#include "JobResultRec.h"

void JobResultRec::extract(const Job & job, const JobEnvironment &env) 
{
	const EnvGLNRec & envGLN = (const EnvGLNRec &) env;
	int child = ((JobGLN &) job).getNodeId() - 1;
	* (TransitionTable *) this =  envGLN.getTransTable( child );
}

void JobResultRec::merge(JobEnvironment &env) const
{
	EnvGLNRec & envGLN = (EnvGLNRec &) env;

	int child = getChild();

	if( envGLN.compareParentSets((TransitionTable) *this, envGLN.getTransTable( child ) ) )
	{
		envGLN.setTransTable(child, (TransitionTable) *this);
	}
}	

int JobResultRec::pack(unsigned char* & buffer) const
{
	size_t size = ((TransitionTable *) this) -> numBytes();

	try {

		buffer = new unsigned char [size];

	} catch (bad_alloc & ba) {

		cerr << "ERROR: Failed to allocate a buffer for job result!" << endl;
		cerr << "\t" << ba.what() << endl;
	}

	size_t actualSize = ((TransitionTable *) this) -> pack(buffer, size);

	return (int) actualSize;
}

void JobResultRec::unpack(unsigned char * buffer, int length)
{
	// size_t actualSize = ((TransitionTable *) this) -> unpack(buffer, length);
    ((TransitionTable *) this) -> unpack(buffer, length);
}

/*

int JobResultRec::pack(unsigned char * & jobResultBuffer) const
{
	//cout << "in pack function" << endl;
	// EnvGLNRec *penv = (EnvGLNRec *) &env;

	int num_parents = (int) m_parents.size();
	int ttrows= (int) m_transitionTable.size();

	//child + base + bestparents + parentbases + bestoffsets + numparents + row + transtable + df
	int intelementsSize = 1 + 1 + num_parents*3 + 1 + 1 + (m_base * ttrows) + 1;
	
	int doubleelementsSize = 1 + 1;//chisq + pvalue
	int totalSize = (sizeof(int) * intelementsSize) + (sizeof(double) * doubleelementsSize);
	
	try {

		jobResultBuffer = new unsigned char [totalSize];

	} catch (bad_alloc & ba) {

		cerr << "ERROR: Failed to allocate a buffer for job result!" << endl;
		cerr << "\t" << ba.what() << endl;
	}

	unsigned char* p = jobResultBuffer;

	*(int*) p = m_child;
	p = p + sizeof(int);

	*(int*) p = m_base;
	p = p + sizeof(int);


	*(int*) p = num_parents;
	p = p + sizeof(int);

	*(int*) p = ttrows;

	for (int i = 0; i < ttrows; i++){
		for (int j = 0; j < m_base; j++){
			p = p + sizeof(int);
			*(int*) p = m_transitionTable[i][j];
		}
	}

	for (int i = 0; i < num_parents; i++){
		p = p + sizeof(int);
		*(int*) p = m_parents[i];
	}

	for (int i = 0; i < num_parents; i++){
		p = p + sizeof(int);
		*(int*) p = m_parentBases[i];
	}


	for (int i = 0; i < num_parents; i++){
		p = p + sizeof(int);
		*(int*) p = m_delays[i];
	}
	p = p + sizeof(int);

	*(double*) p = m_chisq;
	p += sizeof(double);

	*(int *) p = m_df;
	p += sizeof(int);

	*(double*) p = m_pValue;

	return totalSize;
	
}

*/

/*
void JobResultRec::unpack(unsigned char * jobResultBuffer, int length)
{
	// EnvGLNRec *penv = (EnvGLNRec *) &env;

	unsigned char* p = 0;
	p = jobResultBuffer;

	m_child = *(int*) p;
	p = p+sizeof(int);

	m_base = *(int*) p;
	p = p+sizeof(int);

	int num_parents = *(int*) p;
	p = p+sizeof(int);

	int ttrows = *(int*) p;
	
	// trans_table = createTable(size_transtable, penv->m_trajCol.base(nodeId-1) );
	m_transitionTable.resize(ttrows);
	for(int i=0; i<ttrows; i++) {
		m_transitionTable[i].resize( m_base );
	}

	for (int i = 0; i < ttrows; i++){
		for (int j = 0; j < m_base; j++){
			p = p + sizeof(int);
			m_transitionTable[i][j] =*(int*) p;
		}
	}

	m_parents.resize(num_parents);
	
	for (int i = 0; i < num_parents; i++){
		p = p + sizeof(int);
		m_parents[i] = *(int*) p;
	}

	m_parentBases.resize(num_parents);
	
	for (int i = 0; i < num_parents; i++){
		p = p + sizeof(int);
		m_parentBases[i] = *(int*) p;
	}

	m_delays.resize(num_parents);

	for (int i = 0; i < num_parents; i++){
		p = p + sizeof(int);
		m_delays[i] = *(int*) p;
	}
	p = p + sizeof(int);

	m_chisq = *(double*) p;
	p += sizeof(double);

	m_df = *(int *) p;
	p += sizeof(int);

	m_pValue = *(double*) p;

}
*/

//void JobResultRec::mergeJobResults(JobEnvironment &env)
//{
//	EnvGLNRec *penv = (EnvGLNRec*) &env;

//	penv -> m_reconstruct.merge(* this);

	/*
	// if(p_val <= penv->alpha && (p_val < penv->p_values[nodeId-1] 
	//		|| ( (p_val == penv->p_values[nodeId-1]) && chisq > penv->chisqs[nodeId-1] ) )) {

	if( isResultBetter( penv-> m_estPar.m_alpha, m_p_val, m_chisq, getNumParents(), 
		penv -> p_values[m_nodeId-1], penv -> chisqs[m_nodeId-1], penv->num_parents[m_nodeId-1]) ) {

			//delete_array = keepBestValues(penv->listsOfParents, penv->trans_table, penv->best_trans_tables,this->tt_rows, this->num_parents, this->sizeTruthTables, this->indexParents,penv-> p_values, this->p_val, this-> m_nParents, this->m_child, penv->listsOfOffsets, this->offset);

			//penv -> num_parents[this -> m_child] = this ->num_parents[this -> m_child]; 	

			penv -> chisqs[m_nodeId-1] = m_chisq;	// Joe Song
			penv -> degree_of_freedoms[m_nodeId-1] = (size_transtable-1) * (penv->m_trajCol.base( m_nodeId-1 ) - 1); // Joe Song. September 23, 2006; September 5, 2008. Changed size_transtable to (size_transtable -1)

			penv->listsOfOffsets[nodeId-1].resize( getNumParents() );

			for (int i = 0; i < getNumParents(); i++) {
				penv-> listsOfOffsets[m_nodeId-1][i] = m_best_offsets[i];
			}

			penv->listsOfParents[m_nodeId-1].resize( getNumParents() );

			for(int i=0; i < getNumParents(); i++){
				penv->listsOfParents[m_nodeId-1][i] = m_best_parents[i];
			}//end for

			//Update the best truth table
			if(penv->best_trans_tables[m_nodeId-1]!=0){

				// need to know the size of previous truth table to properly delete
				deleteArray(penv-> best_trans_tables[m_nodeId-1], penv->sizeTruthTables[m_nodeId-1]);
				penv->best_trans_tables[m_nodeId-1] = 0;
			} //end if

			penv->best_trans_tables[m_nodeId-1] = createTable(size_transtable, penv->m_trajCol.base( nodeId-1 ) );
			for (int i = 0; i < size_transtable; i++){
				for (int j = 0; j < penv->m_trajCol.base( nodeId-1 ); j++){
					penv->best_trans_tables[nodeId-1][i][j] = m_trans_table[i][j];
					//cout << "transtable at : " << i << " : " << j << "  " <<trans_table[i][j] <<"\t";
				}
				//cout << endl;
			}
			//cout << endl;
			//}

			penv->sizeTruthTables[m_nodeId-1] = getTransTableSize();	//update the number of rows in the current "best truth table"
			penv->num_parents[nodeId-1] = getNumParents(); 	//update the number of parents
			penv->p_values[nodeId-1] = m_p_val; 	//update lowest p val
			//delete_array = true;


	}//end if	
	
	/ *
	if(trans_table!=0){
		deleteArray( trans_table, size_transtable);
		trans_table = 0;
	}//end if

	* /
}
*/

/*
void JobResultRec::unpack(unsigned char * jobResultBuffer, int length, JobEnvironment & env)
{
	
	// cout << "JobResultRec::unpack() {" << endl;

	EnvGLNRec *penv = (EnvGLNRec *) &env;
	unsigned char* p = 0;
	p = jobResultBuffer;

	m_nodeId = *(int*) p;
	p = p+sizeof(int);

	int num_parents = *(int*) p;
	p = p+sizeof(int);

	int size_transtable = *(int*) p;
	
	// trans_table = createTable(size_transtable, penv->m_trajCol.base(nodeId-1) );
	m_trans_table.resize(size_transtable);
	for(int i=0; i<size_transtable; i++) {
		m_trans_table[i].resize( penv->m_trajCol.base(m_nodeId-1) );
	}

	for (int i = 0; i < size_transtable; i++){
		for (int j = 0; j < penv->m_trajCol.base( m_nodeId-1 ); j++){
			p = p + sizeof(int);
			m_trans_table[i][j] =*(int*) p;
		}
	}

	m_best_parents.resize(num_parents);

	for (int i = 0; i < num_parents; i++){
		p = p + sizeof(int);
		m_best_parents[i] = *(int*) p;
	}

	m_best_offsets.resize(num_parents);

	for (int i = 0; i < num_parents; i++){
		p = p + sizeof(int);
		m_best_offsets[i] = *(int*) p;
	}
	p = p + sizeof(int);

	m_chisq = *(double*) p;;
	
	p = p + sizeof(double);
	
	m_p_val = *(double*) p;

	// print(); 

	// cout << "}" << endl;

}
*/

// void JobResultRec::unpack(unsigned char * jobResultBuffer, int length, const JobEnvironment & env)

/*	
int JobResultRec::pack(unsigned char* &jobResultBuffer, JobEnvironment & env)
{
	//cout << "in pack function" << endl;
	EnvGLNRec *penv = (EnvGLNRec *) &env;

	int num_parents = (int) getNumParents();
	int size_transtable = (int) getTransTableSize();

	//bestparents + bestoffsets + numparents + sizetable + transtabl+nodeId
	int intelementsSize = num_parents + num_parents + 1 + 1 +1+ (penv -> m_trajCol.base( m_nodeId-1 ) * size_transtable);
	
	int doubleelementsSize = 2;
	int totalSize = (sizeof(int) * intelementsSize) + (sizeof(double) * doubleelementsSize);
	
	try {

		jobResultBuffer = new unsigned char [totalSize];

	} catch (bad_alloc & ba) {

		cerr << "ERROR: Failed to allocate a buffer for job result!" << endl;
		cerr << "\t" << ba.what() << endl;
	}

	unsigned char* p = jobResultBuffer;

	*(int*) p = m_nodeId;
	p = p + sizeof(int);

	*(int*) p = num_parents;
	p = p + sizeof(int);

	*(int*) p = size_transtable;

	for (int i = 0; i < size_transtable; i++){
		for (int j = 0; j < penv-> m_trajCol.base( m_nodeId-1 ); j++){
			p = p + sizeof(int);
			*(int*) p = m_trans_table[i][j];
		}
	}

	for (int i = 0; i < num_parents; i++){
		p = p + sizeof(int);
		*(int*) p = m_best_parents[i];
	}

	for (int i = 0; i < num_parents; i++){
		p = p + sizeof(int);
		*(int*) p = m_best_offsets[i];
	}
	p = p + sizeof(int);

	*(double*) p = m_chisq;
	
	p = p + sizeof(double);
	
	*(double*) p = m_p_val;

	return totalSize;
	
}
*/

// int JobResultRec::pack(unsigned char* &jobResultBuffer, JobEnvironment & env)

//void JobResultRec::setResult(EnvGLNRec *env)
/*

void JobResultRec::setResult(JobEnvironment & environment)
{
	EnvGLNRec * env = (EnvGLNRec *) & environment;

	* this = * (JobResultRec *) ( env -> m_reconstruct.extract() );

-------------------------
	nodeId = env-> m_child+1;

	chisq = env->chisqs[env->m_child];

	if(env->listsOfParents[env->m_child].size() > 0){
		for (int i = 0; i < env->m_nParents; i++) {
			best_parents[i] = env->listsOfParents[env->m_child][i];
		}
	}

	if(env->listsOfOffsets[env->m_child].size() > 0){
		for (int i = 0; i < env->m_nParents; i++) {
			best_offsets[i] = env->listsOfOffsets[env->m_child][i];
		}
	}

	//Update the best truth table
	if(trans_table !=0 ){

		// need to know the size of previous truth table to properly delete 
		deleteArray(trans_table, size_transtable);
		trans_table = 0;
	}//end if

	//trans_table = env->best_trans_tables[env->m_child];		//Store pointer to the current best truth table

	size_transtable = env->sizeTruthTables[env->m_child];	//update the number of rows in the current "best truth table"
	//int col = env->m_nParents;
	if(env->best_trans_tables[env->m_child]!=0){
		//int col = env->m_nParents;
		//if(env->m_nParents == 1){
		//trans_table = createTable(size_transtable,env->m_nParents+1);
		//for (int i = 0; i < size_transtable; i++){
		//	for (int j = 0; j < env->  m_nParents+1; j++){  // m_trajCol.base_lookup[env->m_child]
		//		trans_table[i][j] = env->best_trans_tables[env->m_child][i][j];
		//	}
		//}
		//}


		//else{
		trans_table = createTable(size_transtable,env->m_trajCol.base( env->m_child ) );
		for (int i = 0; i < size_transtable; i++){
			for (int j = 0; j < env->m_trajCol.base( env->m_child ); j++){
				trans_table[i][j] = env->best_trans_tables[env->m_child][i][j];
			}
		}
	}
	//}


	// if(){
	//      deleteArray( trans_table, size_transtable);
	//    trans_table = 0;
	//}//end if


	num_parents = env->m_nParents; 	//update the number of parents
	p_val = env->p_values[env->m_child]; 	//update lowest p val

}
	*/

