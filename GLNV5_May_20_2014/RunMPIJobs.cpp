// RunMPIJobs.cpp
//
// Curtis Luce, Joe Song
// Created 2007
// Modified: September 27, 2008. Joe Song
// Modified: October 3, 2008. Joe Song.  Removed dependency on anything related to GLN.
//    The MPI related functions in this file are now all generic.
// Modified: November 28, 2008. Joe Song.  Changed how the results are generated.  
// Last modified: October 1, 2011. Name changed from "MPIRunJobs.cpp"

#ifdef ENABLE_MPI

#include <mpi.h>

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <new>
using std::bad_alloc;

#include "JobEnvironment.h"
#include "Job.h"
#include "JobResult.h"

void Manager(int size, JobEnvironment & env, Job & job, JobResult & result);

void Worker (JobEnvironment & env, Job & job, JobResult & result);

int MPIRunJobs(JobEnvironment & env, Job & job, JobResult & result)
{
	int rank;
	int size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//	if(size < 3){
	//		cerr << "You must have at least 3 processors!" << endl;
	//		return 1;
	//	}

	if(rank == 0) {

		Manager(size, env, job, result);

	} else {

		Worker(env, job, result);
	}

	//MPI_Finalize();
	return 0;
}

void Manager(int size, JobEnvironment & env, Job & job, JobResult & result)
{
	int unfinishedJobs = 0;

	int init = 0;

	int buffSize = 0;

	unsigned char * jobBuffer = 0;
	unsigned char * resultBuffer = 0;

	MPI_Status status;

	for (int i = 1; i < size; i++) {

		if (init == 0) {

			job.first(env);

			buffSize = job.pack(jobBuffer);
			init = 1;

		} else {

			job.next(env);

			buffSize = job.pack(jobBuffer);

		}

		job.brief();

		MPI_Send(&buffSize, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
		MPI_Send(jobBuffer, buffSize, MPI_UNSIGNED_CHAR, i, 1, MPI_COMM_WORLD);
		delete [] jobBuffer;

		if(! job.isTerminal( ) ) {

			unfinishedJobs++;

		}

	}

	do {
		// Receive a job result
		MPI_Recv(&buffSize, 1, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
		
		if (buffSize > 0) { // Retrieve the job result if it is not empty

			// Allocate result buffer
			try {

				resultBuffer = new unsigned char [buffSize];

			} catch (bad_alloc & ba) {

				cerr << "ERROR: Failed to allocate a buffer for received result in manager!" << endl;
				cerr << "\t" << ba.what() << endl;
			}


			MPI_Recv(resultBuffer, buffSize, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, 1, MPI_COMM_WORLD, &status);

			// result.unpack(resultBuffer, buffSize, env);
			result.unpack(resultBuffer, buffSize);

			// Free result buffer
			delete [] resultBuffer;
			resultBuffer = 0;

			// result.mergeJobResults(env);
			result.merge(env);

			// result.~JobResultRec();	Joe Song September 29, 2008.  
			// The above was fatal in MPI.  Caused the following MPI run time errors depending on number of
			//    processors used:
			//    1. "glibc detected *** double free or corruption"
			//    2. "terminated with signal 11"

		}

		// One job has just finished regardless of its job result being empty or not
		unfinishedJobs--;

		// Generate next job if current job is not terminal.
		job.next(env);

		job.brief();

		buffSize = job.pack(jobBuffer);
	
		MPI_Send(&buffSize, 1, MPI_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD);
		MPI_Send(jobBuffer, buffSize, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, 1, MPI_COMM_WORLD);

		delete [] jobBuffer;
		jobBuffer = 0;

		if(! job.isTerminal( ) ) {
			unfinishedJobs++;
		}

	} while( unfinishedJobs  > 0 );
	
	// env.finalize();

}

void Worker(JobEnvironment & env, Job & job,JobResult & result)
{
	int buffSize = 0;

	unsigned char * jobBuffer = 0;

	// JobResult * result = 0;
	unsigned char * resultBuffer = 0;

	MPI_Status status;

	while(true) {

		MPI_Recv(&buffSize, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);

		// Allocate job buffer
		try {
	
			jobBuffer = new unsigned char [buffSize];
	
		} catch (bad_alloc & ba) {

			cerr << "ERROR: Failed to allocate a buffer for received job in worker!" << endl;
			cerr << "\t" << ba.what() << endl;
		}
	

		MPI_Recv(jobBuffer, buffSize, MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD, &status);

		job.unpack(jobBuffer, buffSize);

		// Free job buffer
		delete [] jobBuffer;
		jobBuffer = 0;

		// job.transfer(env);

		if(job.isTerminal() == false){

			// result = job.process(env);
			bool improved = job.process(env);

			// if(result) {
			if(improved) {

				result.extract(job, env);

				// buffSize = result -> pack(resultBuffer, env);
				buffSize = result.pack(resultBuffer);

				// delete result;
				// result = 0;

				MPI_Send(&buffSize, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
				MPI_Send(resultBuffer, buffSize, MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);

				delete [] resultBuffer;
				resultBuffer = 0;

			} else {
			
				buffSize = 0;
				MPI_Send(&buffSize, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);

			}

		} else {
			break;
		}

	}

}
#endif
