// GLNReconstruct.cpp - generalized logical network reconstruction
//
// Joe Song
// Created: November 28, 2008 -- extracted from GLN_EstSP.cpp and GLN_EstMPI.cpp

#include <iostream>
using std::cout;

#include "EnvGLNRec.h"
#include "JobGLN.h"
#include "JobResultRec.h"

#ifndef ENABLE_MPI

void SingleProcessorRunJobs(JobEnvironment & envman, Job & job, JobResult & managerresult, 
							JobEnvironment ** workerEnvs, Job ** workerJobs, JobResult ** workerResults,
							int nWorkers);

GeneralizedLogicalNetwork ReconstructGLN(const GLNConfig & gc)
{
	EnvGLNRec envman;
	JobGLN job;
	JobResultRec managerresult;

	int numWorkers = 3; // Number of workers

	EnvGLNRec * workerEnvs = new EnvGLNRec [numWorkers];
	JobGLN * workerJobs = new JobGLN [numWorkers];
	JobResultRec * workerResults = new JobResultRec [numWorkers];

	EnvGLNRec ** pworkerEnvs = new EnvGLNRec * [numWorkers];
	JobGLN ** pworkerJobs = new JobGLN * [numWorkers];
	JobResultRec ** pworkerResults = new JobResultRec * [numWorkers];

	envman.initialize(gc);

	for(int i=0; i<numWorkers; ++i) {
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

	return envman.getGLN();
}

#else

int MPIRunJobs(JobEnvironment &env, Job & job, JobResult & result);

GeneralizedLogicalNetwork ReconstructGLN(const GLNConfig & gc)
{
	EnvGLNRec env;
	JobGLN job;
	JobResultRec result;

	env.initialize(gc);

	//call MPI file to run parallel processing
	MPIRunJobs((JobEnvironment &) env, (Job &) job, (JobResult &) result);

	return env.getGLN();

}

#endif

