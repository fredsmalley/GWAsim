// RunSPJobs.cpp -- Single processor job and environment interactions 
//
// Joe Song
// Created: November 28, 2008.  Extracted from GLN_EstSP.cpp
// Last modified: October 1, 2011. Name changed from "SPRunJobs.cpp"

#ifndef ENABLE_MPI

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <new>
using std::bad_alloc;

#include "JobEnvironment.h"
#include "Job.h"
#include "JobResult.h"

void SingleProcessorRunJobs(JobEnvironment & envman, Job & managerJob, JobResult & managerresult, 
							JobEnvironment ** workerEnvs, Job ** workerJobs, JobResult ** workerResults,
							int nWorkers)
{
	unsigned char *jobBuffer=0;
	unsigned char *resultBuffer=0;

	managerJob.first(envman); //Initialization, prepare to process the first node.
	
	int n = 0;

	while(! managerJob.isTerminal() ) {
		managerJob.brief();  //brief() is a output function, it does nothing other than printing.

        //pack() and unpack() are used to wrap and unwrap memories.
        //So that information could be passed between master worker and slave worker.
		int buffSize = managerJob.pack(jobBuffer);

		workerJobs[n] -> unpack(jobBuffer, buffSize);

		if(jobBuffer){
			delete [] jobBuffer;
			jobBuffer = 0;
		}

		// workerResults[n] = workerJobs[n] -> process(* workerEnvs[n]);

		bool improved = workerJobs[n] -> process(* workerEnvs[n]);

		// if( workerResults[n] ){
		if( improved ){
	
			workerResults[n] -> extract(* workerJobs[n], * workerEnvs[n]);

			// buffSize = workerResults[n] -> pack(resultBuffer, * workerEnvs[n]);
			buffSize = workerResults[n] -> pack(resultBuffer);

			// delete workerResults[n];
			// workerResults[n] = 0;

			if(resultBuffer) {
			
				managerresult.unpack(resultBuffer, buffSize);

				delete[] resultBuffer;
				resultBuffer = 0;
			}

			managerresult.merge(envman);

		} 

		//Increase the number of possible parents, then check to see if a new set of parents can be found.
        //If number of parents has reached the upper limited, then move on to next node in 'envman'.
        managerJob.next(envman);

		++n;

		n = (n >= nWorkers) ? 0 : n;

	}
	
}

#endif

