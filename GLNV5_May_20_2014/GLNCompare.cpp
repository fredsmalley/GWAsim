// GLNCompare.cpp - generalized logical network comparison
//
// Joe Song
// Created: November 28, 2008

#include <iostream>
using std::cout;

#include "EnvGLNCmp.h"
#include "JobGLN.h"
#include "JobResultCmp.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#ifndef ENABLE_MPI

void SingleProcessorRunJobs(JobEnvironment & envman, Job & job, JobResult & managerresult, 
							JobEnvironment ** workerEnvs, Job ** workerJobs, JobResult ** workerResults,
							int nWorkers);

void CompareGLNs(EnvGLNCmp & envman)
{
	JobGLN job;
	JobResultCmp managerresult;

	int numWorkers = 3; // Number of workers

	EnvGLNCmp * workerEnvs = new EnvGLNCmp [numWorkers];
	JobGLN * workerJobs = new JobGLN [numWorkers];
	JobResultCmp * workerResults = new JobResultCmp [numWorkers];

	EnvGLNCmp ** pworkerEnvs = new EnvGLNCmp * [numWorkers];
	JobGLN ** pworkerJobs = new JobGLN * [numWorkers];
	JobResultCmp ** pworkerResults = new JobResultCmp * [numWorkers];


	for(int i=0; i<numWorkers; i++) {
		workerEnvs[i] = envman; // .initialize(gc);
		pworkerEnvs[i] = workerEnvs + i;
		pworkerJobs[i] = workerJobs + i;
		pworkerResults[i] = workerResults + i;
	}

	SingleProcessorRunJobs((JobEnvironment &) envman, (Job &) job,
                           (JobResult &) managerresult,
                           (JobEnvironment  **) pworkerEnvs,
                           (Job **) pworkerJobs,
                           (JobResult **) pworkerResults, numWorkers);

	// envman.finalize();

	if(pworkerResults) delete [] pworkerResults;
	if(pworkerJobs) delete [] pworkerJobs;
	if(pworkerEnvs) delete [] pworkerEnvs;

	if(workerEnvs) delete [] workerEnvs;
	if(workerJobs) delete [] workerJobs;
	if(workerResults) delete [] workerResults;

}

/*
void CompareGLNs(const GLNConfig & gc)
{
	EnvGLNCmp envman;
	JobGLN job;
	JobResultCmp managerresult;
    
	int numWorkers = 3; // Number of workers
    
	EnvGLNCmp * workerEnvs = new EnvGLNCmp [numWorkers];
	JobGLN * workerJobs = new JobGLN [numWorkers];
	JobResultCmp * workerResults = new JobResultCmp [numWorkers];
    
	EnvGLNCmp ** pworkerEnvs = new EnvGLNCmp * [numWorkers];
	JobGLN ** pworkerJobs = new JobGLN * [numWorkers];
	JobResultCmp ** pworkerResults = new JobResultCmp * [numWorkers];
    
	envman.initialize(gc);
    
	for(int i=0; i<numWorkers; i++) {
		workerEnvs[i].initialize(gc);
		pworkerEnvs[i] = workerEnvs + i;
		pworkerJobs[i] = workerJobs + i;
		pworkerResults[i] = workerResults + i;
	}
    
	SingleProcessorRunJobs((JobEnvironment &) envman, (Job &) job, (JobResult &) managerresult,
                           (JobEnvironment  **) pworkerEnvs, (Job **) pworkerJobs, (JobResult **) pworkerResults, numWorkers);
    
	envman.finalize();
    
	if(pworkerResults) delete [] pworkerResults;
	if(pworkerJobs) delete [] pworkerJobs;
	if(pworkerEnvs) delete [] pworkerEnvs;
    
	if(workerEnvs) delete [] workerEnvs;
	if(workerJobs) delete [] workerJobs;
	if(workerResults) delete [] workerResults;
    
}
*/

#else

int MPIRunJobs(JobEnvironment &env, Job & job, JobResult & result);

/*
void CompareGLNs(const GLNConfig & gc)
{
	EnvGLNCmp env;
	JobGLN job;
	JobResultCmp result;

	env.initialize(gc);

	//call MPI file to run parallel processing
	MPIRunJobs((JobEnvironment &) env, (Job &) job, (JobResult &) result);
}
*/

void CompareGLNs(EnvGLNCmp & env)
{
	JobGLN job;
	JobResultCmp result;
    
	//call MPI file to run parallel processing
	MPIRunJobs((JobEnvironment &) env, (Job &) job, (JobResult &) result);
}

#endif


EnvGLNCmp CompareGLNs(const GLNConfig & gc)
{
	EnvGLNCmp envman;
	envman.initialize(gc);
    CompareGLNs(envman);
	envman.finalize();
    return envman;
}


EnvGLNCmp BuildTopologyThenCompareGLNs(const GLNConfig & gc)
{
    GeneralizedLogicalNetwork ReconstructGLN(const GLNConfig & gc);
    
    vector< GeneralizedLogicalNetwork > glns(gc.m_listTrajFiles.size());

    GLNConfig gck(gc);

    // change gck for reconstruction of the k-th network    
    // parent selection method must be specified for reconstruction 
    gck.m_methodRecParCmp = BY_PVAL;

    gck.m_fileGLN.clear(); // no output file
    gck.m_fileDOT.clear();
    gck.m_adjustedpValFile.clear();
    gck.m_adjustedDOTFile.clear();
    
    gck.m_alpha = 1.5;  // So that we keep the maximal topology of the network
    gck.m_listTrajFiles.resize(1);
    
    if(gck.m_pathways.size()>1) {
        gck.m_pathways.resize(1);
    } else if(gck.m_candidateTopologys.size()>1) {
        gck.m_candidateTopologys.resize(1);
    }
    
    // Reconstruct each topology
    for ( size_t k=0; k<glns.size(); k++ ) {
        // only the k-th TCF
        gck.m_listTrajFiles[0] = gc.m_listTrajFiles[k];
        if(gc.m_pathways.size()==1) {
            gck.m_pathways[0] = gc.m_pathways[0];
        } else if(gc.m_candidateTopologys.size()==1) {
            gck.m_candidateTopologys[0] = gc.m_candidateTopologys[0];
        } else if(gc.m_pathways.size()>1) {
            gck.m_pathways[0] = gc.m_pathways[k];
        } else if(gc.m_candidateTopologys.size()>1) {
            gck.m_candidateTopologys[0] = gc.m_candidateTopologys[k];
        }
        glns[k] = ReconstructGLN(gck);
    }

    EnvGLNCmp envman;

#ifdef ENABLE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
#endif
        
    // Apply comparative analysis using the given topology

    // Do the following only when rank is 0

    envman.initialize(gc);

    vector< vector<int> > parentSets( glns.size() );
    vector< vector<int> > delaySets( glns.size() );
    
    for (size_t child=0; child < envman.getGLNSize(); child++) {
        
        for ( size_t k=0; k<glns.size(); k++ ) {
            parentSets[k] = glns[k].getGTT(child).getParents();
            delaySets[k] = glns[k].getGTT(child).getDelays();
        }
        
        envman.processOneEnumeration(child, parentSets, delaySets);
    }

    envman.finalize();

#ifdef ENABLE_MPI
  }
#endif
    return envman;
}
