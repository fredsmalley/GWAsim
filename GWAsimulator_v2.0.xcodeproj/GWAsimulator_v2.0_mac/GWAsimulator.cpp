////////////////////////////////////////////////////////////////////
//
//  GWAsimulator.cpp is a C++ program for simulating case-control or
//  population genetic data using SNP chips.
//
//  version 2.0 (11/19/2007)
//  Authors:  Chun Li, Mingyao Li
////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <math.h>
using namespace std;

////////////////////////////////////////////////////////////////////
//  This line is needed if the program is compiled by gcc version 3.4
//  or above.  For earlier versions, comment this line out.
////////////////////////////////////////////////////////////////////
typedef std::vector<bool, std::allocator<bool> > bit_vector;


////////////////////////////////////////////////////////////////////
//
//  Global variables
//
////////////////////////////////////////////////////////////////////
//  Variables for compiler preprocessing.  The numbers are set quite
//  high, and probably won't need to be changed for most situations.
#define MAXPHASE 1000   // maximum number of phased chromosomes in input data
#define MAXSUBJ 20000   // maximum number of subjects to be simulated
#define MAXSNP 100000   // maximum number of SNPs per chromosome
#define MAXI2 256       // maximum number of SNP pairs with 2-way interaction

//////////////////////////////////////////////////
// seed is read in from command line
static long seed;       // seed number for random number generator

/////////////////////////////////////////
// variables for input and simulated data
/////////////////////////////////////////
// Chr[] connects internal coding index (0-22) to chromosome numbers
// The first several are disease chromosomes.
int Chr[23];
int Start[23], End[23]; // start and end positions for regional simulation
int numMarker[23];                 // number of SNPs on each chromosome
bit_vector hap[23][MAXPHASE];      // input phased data
vector<double> freq[23];           // allele frequencies
bit_vector sim_hap[23][MAXSUBJ];   // simulated data

/////////////////////////////////////////////////
// model-related values calculated by the program
/////////////////////////////////////////////////
int GTOTAL;             // number of multi-marker disease genotypes (=3^NUMDL)
double *prob, *summ;    // intermediate values for model calculation
                        // and conditional probability calculation
double *cum_prob_caseF, *cum_prob_caseM, *cum_prob_contF, *cum_prob_contM;
                        // cumulative genotype distributions
int I2loci[MAXI2][2];   // SNP pairs with 2-way interaction effects
double I2effects[MAXI2][4];  // effects of 2-way interactions

////////////////////////////////////////////////////////////////////
//
//  Random number generator of uniform (0, 1)
//  Initialized with the address of a negative number (static long).
//
////////////////////////////////////////////////////////////////////
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy) { 
    if (-(*idum) < 1) *idum=1; 
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }

  k=(*idum)/IQ; 
  *idum=IA*(*idum-k*IQ)-IR*k; 
  if (*idum < 0) *idum += IM;
  j=iy/NDIV; 
  iy=iv[j]; 
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX; 
  else return temp;
}


////////////////////////////////////////////////////////////////////
//
//  Read in the control file.  See the example control file for more
//  details.
//
////////////////////////////////////////////////////////////////////
void controlin(char *CONTFILE,
	       char *DATADIR1, char *DATADIR2,
	       int &NUMCHROM, int &NUMCHROM_X,
	       int &NEEDOUTPUT, char *OUTFORMAT,
	       int &WINDOW_SIZE, int &NUMCASEF, int &NUMCASEM,
	       int &NUMCONTF, int &NUMCONTM,
	       double &PREV, int &NUMDL, int &REGIONAL,
	       int *DLPOS, int *DV, double *GRR, double *GRR2,
	       int &INTER, int &NUMI2)
{
  int i, j, nNonDisChr, flag, needexit=0, fileexist=0, nfield, tmpint;
  double tmpdouble;
  char *pt;
  char tmp[1024], tmp2[1024], testfilename[1024];
  char t1[128], t2[128], t3[128], t4[128], t5[128], t6[128];
  string line;
  ifstream controlfile, testfile;

  controlfile.open(CONTFILE);
  if(controlfile.fail()){
    printf("Cannot open the control file %s\n", CONTFILE);
    exit(1);
  }

  ////
  //// input filename pre- and affix
  getline(controlfile, line);
  sscanf(line.c_str(), "%s%s", DATADIR1, DATADIR2);
  printf("\nInput data files:  %s#%s\n", DATADIR1, DATADIR2);

  ////
  //// number of phased chromosomes/lines for autosomes and X
  getline(controlfile, line);
  sscanf(line.c_str(), "%d%d", &NUMCHROM, &NUMCHROM_X);
  printf("Number of phased chromosomes:   %d autosomes, %d X\n",
	 NUMCHROM, NUMCHROM_X);
  if(NUMCHROM > MAXPHASE | NUMCHROM_X > MAXPHASE) {
    printf("\nThe maximum number of phased chromosomes is set to be %d.\n",
	   MAXPHASE);
    printf("Increase the MAXPHASE variable in the program and recompile.\n\n");
    exit(1);
  }

  ////
  //// NEEDOUTPUT = 1: output simulated data to files, 0: no output
  getline(controlfile, line);
  sscanf(line.c_str(), "%d%s", &NEEDOUTPUT, OUTFORMAT);
  if(! NEEDOUTPUT) {
    puts("Write simulated data to files:  NO");
  } else {
    // change all letters to lower case and check
    pt=OUTFORMAT; while (*pt != '\0') {*pt = tolower(*pt); pt++;}
    printf("Write simulated data to files:  YES (format: %s)\n", OUTFORMAT);
    if(strcmp(OUTFORMAT, "linkage") &
       strcmp(OUTFORMAT, "genotype") &
       strcmp(OUTFORMAT, "phased")) {
      puts("\nOutput format should be linkage or genotype or phased");
      exit(1);
    }

    // detect any files named chr#.dat or chr#.dat.gz
    for(i=0; i<23; i++) {
      sprintf(testfilename, "chr%d.dat", i);
      testfile.open(testfilename);
      if(testfile.is_open()) {
	fileexist = 1;
	testfile.close();
	break;
      }
      sprintf(testfilename, "chr%d.dat.gz", i);
      testfile.open(testfilename);
      if(testfile.is_open()) {
	fileexist = 1;
	testfile.close();
	break;
      }
    }
    if(fileexist) {
      puts("\nWARNING:  Any files named chr#.dat or chr#.dat.gz will be removed.");
      puts("          If you don't want them removed, kill the current job immediately.\n");
    }    
  }

  ////
  //// simulation window size
  getline(controlfile, line);
  sscanf(line.c_str(), "%d", &WINDOW_SIZE);
  printf("Window size for simulation:     %d\n", WINDOW_SIZE);
  if(WINDOW_SIZE <= 1) {
    puts("\nThe simulation window size needs to be >=2\n");
    exit(1);
  }
  if(WINDOW_SIZE >= log(NUMCHROM+.0)/log(2.0)-1) {
    puts("\n-----------------  WARNING  --------------------");
    printf("\nWindow size %d may be too large for %d input chromosomes\n", WINDOW_SIZE, NUMCHROM);
    puts("\n------------------------------------------------\n");
  }

  ////
  //// numbers of subjects to be simulated
  getline(controlfile, line);
  sscanf(line.c_str(), "%d%d%d%ds", &NUMCASEF, &NUMCASEM, &NUMCONTF,&NUMCONTM);

  ////
  //// number of disease loci.  If zero, population sampling.
  getline(controlfile, line);
  sscanf(line.c_str(), "%d%d", &NUMDL, &REGIONAL);
  if(NUMDL<0) {
    puts("\nThe number of disease variants is either positive (for case-control data) or zero (for population data)");
    exit(1);
  }
  if(NUMDL>23) {
    puts("\nThe number of disease variants cannot exceed 23");
    exit(1);
  }

  ////
  //// when NUMDL>0, check REGIONAL setting
  if(NUMDL>0 & REGIONAL != 0 & REGIONAL != 1) {
    puts("Sampling:                       Case-control (?)");
    puts("\nREGIONAL indicator should be either 0 (genome) or 1 (regions).");
    exit(1);
  }

  // output the sampling option (case-only or population sampling)
  printf("Sampling:                       ");
  if(NUMDL>0)
    printf("Case-control (%s)\n",
	   (REGIONAL==0) ? "whole genome" : "regions only");
  else puts("Population (whole genome)");

  // output the number of subjects to be simulated
  printf("Number of subjects:             ");
  if(NUMDL>0) {
    printf("%d cases (%d females, %d males)\n",
	   NUMCASEF+NUMCASEM, NUMCASEF, NUMCASEM);
    printf("%32s%d controls (%d females, %d males)\n",
	   "", NUMCONTF+NUMCONTM, NUMCONTF, NUMCONTM);
  } else {
    // force these variable to zero so that the rest of the program can 
    // be used safely without modification
    NUMCASEF = 0; NUMCASEM = 0;
    printf("%d (%d females, %d males)\n",
	   NUMCONTF+NUMCONTM, NUMCONTF, NUMCONTM);
  }

  // check the number of subjects
  if(NUMCASEF<0 | NUMCASEM<0 | NUMCONTF<0 | NUMCONTM<0) {
    puts("Numbers of subjects must be non-negative numbers");
    exit(1);
  }
  // check MAXSUBJ
  if(NUMCASEF+NUMCASEM+NUMCONTF+NUMCONTM > MAXSUBJ) {
    printf("\nThe number of subjects exceeds %d, the maximum allowed.\n",
	   MAXSUBJ);
    printf("Increase the MAXSUBJ variable in the program and recompile.\n\n");
    exit(1);
  }

  //////////////////////////////////////////////////
  //// If NUMDL>0, read in disease model information
  //////////////////////////////////////////////////
  if(NUMDL>0) {
    // disease model prevalence
    getline(controlfile, line);
    sscanf(line.c_str(), "%s", &tmp);
    PREV=atof(tmp);
    printf("Number of disease variants:     %d\n", NUMDL);
    printf("Disease model prevalence:       %6.4f\n", PREV);
    if(PREV <= 0 | PREV >= 1) {
      puts("\nThe prevalence needs to be between 0 and 1.\n");
      exit(1);
    }

    // chromosome, position, GRR for all disease loci
    printf("Locus  chromosome  position  DV  GRR    GRR2");
    if(REGIONAL==0) printf("\n");
    else printf("    Start     End\n");

    ////
    //// read in information for disease loci, one line for each locus
    for(i=0; i<NUMDL; i++) {
      // reading in a line of parameters for a disease locus
      getline(controlfile, line);
      sscanf(line.c_str(), "%d%d%d%s%s%d%d", &Chr[i], &DLPOS[i], &DV[i],
	     &tmp, &tmp2, &Start[i], &End[i]);
      GRR[i]=atof(tmp);

      // if multiplicative effect, then define GRR2=GRR*GRR
      // if dominance effect, then define GRR2=GRR
      if(!strcmp(tmp2, "M") | !strcmp(tmp2, "m"))
	GRR2[i] = GRR[i]*GRR[i];
      else if(!strcmp(tmp2, "D") | !strcmp(tmp2, "d"))
	GRR2[i] = GRR[i];
      else GRR2[i]=atof(tmp2);

      printf("   %2d          %2d    %6d   %1d %6.3f %6.3f",
	     i+1, Chr[i], DLPOS[i], DV[i], GRR[i], GRR2[i]);
      if(REGIONAL==0) printf("\n");
      else printf("  %6d  %6d\n", Start[i], End[i]);

      // check parameters for each disease locus
      if(Chr[i] < 1 | Chr[i] > 23) {
        puts("\nInvalide chromosome number.");
	needexit=1;
      }
      for(j=0; j<i; j++) {
	if(Chr[i] == Chr[j]) {
	  puts("\nThe program doesn't allow >1 disease variant on the same chromosome.");
	  needexit=1;
	}
      }

      if(DLPOS[i] <= (WINDOW_SIZE-1)/2) {
	puts("\nDisease locus position must be > (WINDOW_SIZE-1)/2.");
	needexit=1;
      }

      if(DV[i] != 0 & DV[i] != 1) {
	puts("\nDisease variant must be either 0 or 1.");
	needexit=1;
      }

      if(GRR[i] < 1 | GRR2[i] < 1) {
	puts("\nGenotypic risk ratio needs to be >= 1.");
	needexit=1;
      }

      if(REGIONAL == 1) {
	if(Start[i] <= 0 | End[i] <= 0) {
	  puts("\nStart and end positions need to be positive integers.");
	  needexit=1;
	} else {
	  if(Start[i] > DLPOS[i] - (WINDOW_SIZE-1)/2) {
	    printf("\nStart position should be before disease locus (<= %d - (WINDOW_SIZE-1)/2)\n", DLPOS[i]);
	    needexit=1;
	  }
	  if(End[i] < DLPOS[i] + (WINDOW_SIZE-1)/2) {
	    printf("\nEnd position should be after disease locus (>= %d + (WINDOW_SIZE-1)/2)\n", DLPOS[i]);
	    needexit=1;
	  }
	}
      }

      // this is the end of a disease locus.  Exit if needed.
      if(needexit) exit(1);
    }

    ////
    //// read in interaction effects information if needed
    NUMI2=0;
    do {
      // reading in a line of parameters
      getline(controlfile, line);
      nfield = sscanf(line.c_str(), "%s%s%s%s%s%s%s",
		      &tmp, &t1,&t2,&t3,&t4,&t5,&t6);
      if(nfield < 1) break;

      // if the first word is "Inter2", process the rest of the line
      pt=tmp; while (*pt != '\0') {*pt = tolower(*pt); pt++;}
      if(!strcmp(tmp, "inter2")) {
	INTER=1;
	I2loci[NUMI2][0]=atoi(t1);
	I2loci[NUMI2][1]=atoi(t2);
	I2effects[NUMI2][0]=atof(t3);
	I2effects[NUMI2][1]=atof(t4);
	I2effects[NUMI2][2]=atof(t5);
	I2effects[NUMI2][3]=atof(t6);

	if(NUMI2==0) {
	  puts("Interaction effects among disease loci");
	  puts("Pair  DL1  DL2  1/1    1/2    2/1    2/2");
	}
	printf("%3d    %2d   %2d %6.3f %6.3f %6.3f %6.3f\n", NUMI2+1,
	       I2loci[NUMI2][0], I2loci[NUMI2][1],
	       I2effects[NUMI2][0], I2effects[NUMI2][1],
	       I2effects[NUMI2][2], I2effects[NUMI2][3]);

	// check for errors
	if(I2loci[NUMI2][0] < 1 | I2loci[NUMI2][0] > NUMDL |
	   I2loci[NUMI2][1] < 1 | I2loci[NUMI2][1] > NUMDL) {
	  printf("Disease loci DL1 and DL2 should be between 1 and %d\n",
		 NUMDL);
	  needexit=1;
	}
	if(I2loci[NUMI2][0] == I2loci[NUMI2][1]) {
	  puts("Disease loci DL1 and DL2 should be different.");
	  needexit=1;
	}
	if(I2effects[NUMI2][0] < 0 | I2effects[NUMI2][1] < 0 |
	   I2effects[NUMI2][2] < 0 | I2effects[NUMI2][3] < 0) {
	  puts("Interactive effects should be non-negative.");
	  needexit=1;
	}

	// this is the end of an interaction line.  Exit if needed.
	if(needexit) exit(1);

	// straighten out the parameters
	// if necessary, swap to have lower locus number first, and
	// also swap the effect of 2/1 and 1/2
	if(I2loci[NUMI2][0] > I2loci[NUMI2][1]) {
	  tmpint = I2loci[NUMI2][0];
	  I2loci[NUMI2][0] = I2loci[NUMI2][1];
	  I2loci[NUMI2][1] = tmpint;
	  tmpdouble = I2effects[NUMI2][1];
	  I2effects[NUMI2][1] = I2effects[NUMI2][2];
	  I2effects[NUMI2][2] = tmpdouble;
	}

	NUMI2++;
	if(NUMI2 > MAXI2) {
	  printf("\nA maximum of %d SNP pairs are allowed to have 2-way interaction effects.\n",
		 MAXI2);
	  puts("Increase the MAXI2 variable in the program and recompile.\n");
	  exit(1);
	}
      }
    } while (1);
  }

  controlfile.close();

  // assign the non-disease chromosomes to the index Chr[]
  nNonDisChr = 0;
  for(i=1; i<=23; i++) {
    flag=0;
    // check whether the chromosome is a disease chromosome
    for(j=0; j<NUMDL; j++) {
      if(Chr[j] == i) flag=1;
    }
    // if not a disease chromosome
    if(flag==0){
      Chr[NUMDL+nNonDisChr]=i;
      nNonDisChr++;
    }
  }
}


////////////////////////////////////////////////////////////////////
//
//  Read in all phased data.
//  Each row in the phase data file is a phased chromosome.  Due to
//  the use of bit_vector to store haplotypes, only 0 and 1 are
//  allowed to denote marker alleles.
//    data_in:  main function for this section
//    read_ped_f:  read in a phased data file, tally #SNPs
//    make_marker_s:  calculate allele frequencies
//
////////////////////////////////////////////////////////////////////
// read in a phased data file, store them, tally #SNPs
//////////////////////////////////////////////////////
void read_ped_f(int nHap,         // number of chromosomes
                char *pedFile,    // phased data filename
                int &numMarker,   // # markers on this chromosome
                bit_vector hap[]  // haplotype matrix row=hap, col=marker
                )
{
  int j, g;
  string line;
  char *allele;

  ifstream ped_stream;
  ped_stream.open(pedFile);
  if(ped_stream.fail()){
    cout<< "Chr file: "<< pedFile << " cannot be opened!\n"; 
    exit(1);
  }

  numMarker=0;
  for(j=0; j<nHap; j++) {
    getline(ped_stream, line);  // read in a line

    // break the line into tokens (i.e. alleles, one at a marker)
    // and store the alleles 
    allele = strtok((char*)line.c_str(), " \t");  // first allele
    while (allele != NULL && strcmp(allele, "\n")){
      if (j==0) numMarker++;  // tally number of markers
      g = atoi(allele);
      hap[j].push_back(g);
      allele = strtok(NULL, " \t");  // next allele
    }
  }
  ped_stream.clear();
  ped_stream.close();
}

///////////////////////////////
// calculate allele frequencies 
///////////////////////////////
void make_marker_s(bit_vector hap[],     // haplo matrix:row=hap,col=marker
		   vector<double> *freq, // frequency vector (MAF)
		   int numMarker,        // number of SNPs
		   int nHap              // number of chromosomes
                   )
{
  int i, j, count;
  double p1, revnHap = 1/(double)nHap;

  for(i=0; i<numMarker; i++) {
    count=0;
    for(j=0; j<nHap; j++) {
      count += (hap[j][i]==0);
    }
    p1 = count*revnHap;
    (*freq).push_back(p1);
  }
}

//////////////////////////////
// Read in all the phased data
//////////////////////////////
void data_in(char DATADIR1[], char DATADIR2[],
	     int NUMCHROM, int NUMCHROM_X,
	     int NUMDL, int REGIONAL, int DLPOS[], int WINDOW_SIZE)
{
  int i, j, nHap;
  char pedFile[1024]=""; // input data filename
  int IsXchr;    // 1: X chromosome, 0: autosomal chromosomes

  puts("\nReading in phased data ...");
  for(i=0; i<23; i++){
    if(NUMDL>0 & REGIONAL==1 & i>=NUMDL) continue;
    if(Chr[i]<23) {nHap=NUMCHROM; IsXchr=0;}
    else {nHap=NUMCHROM_X; IsXchr=1;}

    // reserve space, read in data
    for(j=0; j<nHap; j++) hap[i][j].reserve(MAXSNP);
    sprintf(pedFile, "%s%d%s", DATADIR1, Chr[i], DATADIR2);
    read_ped_f(nHap, pedFile, numMarker[i], hap[i]);

    // error checking
    if(i<NUMDL & DLPOS[i] > numMarker[i]-(WINDOW_SIZE-1)/2) {
      printf("Chromosome %d has %d SNPs.  Its disease locus must be <= %d-(WINDOW_SIZE-1)/2\n", Chr[i], numMarker[i], numMarker[i]);
      exit(1);
    }
    if(i<NUMDL & REGIONAL==1) {
      if(End[i] > numMarker[i]) {
	printf("Chromosome %d end position %d is > its total number of %d SNPs\n", Chr[i], End[i], numMarker[i]);
	exit(1);
      }
    }
    if(numMarker[i] > MAXSNP) {
      printf("\nChromosome %d has %d SNPs, exceeding the maximum %d.\n",
	     Chr[i], numMarker[i], MAXSNP);
      printf("Increase the MAXSNP variable in the program and recompile.\n\n");
      exit(1);
    }

    // calculate allele frequencies
    freq[i].reserve(numMarker[i]);
    make_marker_s(hap[i], &freq[i], numMarker[i], nHap);
  }
}


////////////////////////////////////////////////////////////////////
//
//  Build disease model: find out the coefficients for logit penetrance
//  that satisfy the user's constraits: (a) prevalence, and (b) GRR
//  of genotypes 1 vs 0.
//    modelbuild:  main function for this section
//    pow:  power function a^b
//    index_to_disease_geno:  translate a multi-marker index to genotypes
//    Pgeno:  single-marker genotype frequency (autosomes)
//    PgenoX: single-marker genotype frequency (X, in population)
//    PgenoXmale:  single-marker genotype frequency (X, in males)
//    prevalence:  calculate prevalence
//    cum_disease_geno:  calculate conditional genotype frequencies
//
////////////////////////////////////////////////////////////////////
// power function a^b
/////////////////////
int pow(int a, int b)
{
  int tmp=1;
  for(int i=0; i<b; i++)
    tmp*=a;
  return(tmp);
}

//////////////////////////////////////////////
// translate a multi-marker index to genotypes
//////////////////////////////////////////////
void index_to_disease_geno(int index, int d[], int NUMDL)
{
  int remain = index;
  for(int i=NUMDL-1; i>=0; i--) {
    d[i] = int(remain/pow(3,i));
    remain = remain - d[i] * pow(3,i);
  }
}

///////////////////////////////////
// single-marker genotype frequency
///////////////////////////////////
double Pgeno(int g, double p)
{
  if(g==0) return (1-p)*(1-p);
  else if(g==1) return 2*p*(1-p);
  else return p*p;
}

double PgenoX(int g, double p)
{
  if(g==0) return 0.5*(1-p)*(2-p);
  else if(g==1) return p*(1-p);
  else return 0.5*p*(1+p);
}

double PgenoXmale(int g, double p)
{
  if(g==0) return 1-p;
  else if(g==1) return 0;
  else return p;
}

//////////////////////////////////////////
// calculate prevalence given an intercept
//////////////////////////////////////////
double prevalence(double b0, double prob[], double summ[])
{
  int i;
  double aa;

  aa=0;
  for(i=0; i<GTOTAL; i++)
    aa += prob[i]/(1+exp(-b0-summ[i]));
  return aa;
}

////////////////////////////////////////////////////////////////////
// calculate conditional genotype frequencies for cases and controls
////////////////////////////////////////////////////////////////////
void cum_disease_geno(double beta0, int NUMDL, double p[])
{
  int i, j, g[23];
  double prevtmp;

  // females
  for(i=0; i<GTOTAL; i++) {
    index_to_disease_geno(i, g, NUMDL);
    // prob[i]:  frequency of genotype i (in females)
    prob[i]=1;
    for(j=0; j<NUMDL; j++)
      prob[i] *= Pgeno(g[j], p[j]);
  }
  prevtmp = prevalence(beta0, prob, summ);
  cum_prob_caseF[0] = prob[0]/(1+exp(-beta0-summ[0]))/prevtmp;
  cum_prob_contF[0] = prob[0]/(1+exp(beta0+summ[0]))/(1-prevtmp);
  for(i=1; i<GTOTAL; i++) {
    cum_prob_caseF[i] = cum_prob_caseF[i-1] + 
      prob[i]/(1+exp(-beta0-summ[i]))/prevtmp;
    cum_prob_contF[i] = cum_prob_contF[i-1] +
      prob[i]/(1+exp(beta0+summ[i]))/(1-prevtmp);
  }

  // males
  for(i=0; i<GTOTAL; i++) {
    index_to_disease_geno(i, g, NUMDL);
    // prob[i]:  frequency of genotype i (in males)
    prob[i]=1;
    for(j=0; j<NUMDL; j++) {
      if(Chr[j]<23) prob[i] *= Pgeno(g[j], p[j]);
      else prob[i] *= PgenoXmale(g[j], p[j]);
    }
  }
  prevtmp = prevalence(beta0, prob, summ);
  cum_prob_caseM[0] = prob[0]/(1+exp(-beta0-summ[0]))/prevtmp;
  cum_prob_contM[0] = prob[0]/(1+exp(beta0+summ[0]))/(1-prevtmp);
  for(i=1; i<GTOTAL; i++) {
    cum_prob_caseM[i] = cum_prob_caseM[i-1] + 
      prob[i]/(1+exp(-beta0-summ[i]))/prevtmp;
    cum_prob_contM[i] = cum_prob_contM[i-1] +
      prob[i]/(1+exp(beta0+summ[i]))/(1-prevtmp);
  }

  //  for(i=1; i<GTOTAL; i++) {
  //    printf("%8.5f %8.5f %8.5f %8.5f\n",
  //	   cum_prob_caseF[i], cum_prob_contF[i],
  //	   cum_prob_caseM[i], cum_prob_contM[i]);
  //  }
}

/////////////////////////////////////////////////
// Find out the coefficients for logit penetrance
//   Look at the documentation for details.
/////////////////////////////////////////////////
void modelbuild(double PREV, int NUMDL, int DLPOS[], int DV[],
		double *GRR, double *GRR2, int INTER, int NUMI2)
{
  int i, j, g[23], idx1, idx2;
  double pi, ri1, ri2, expai, p[23], b0hi, b0lo, b0tmp, aa;
  double tmpbeta1, tmpbeta2, eff11, eff12, eff21, eff22;
  double beta0;                  // intercept in the logit penetrance model
  double beta1[23], beta2[23];   // coefficients in the logit penetrance model
  char diseasemodel[65536], tmpstr[1024];  // disease model info string
  FILE *modelout;

  puts("\nConstructing disease models ...");

  GTOTAL=pow(3, NUMDL);
  prob=(double *)malloc(GTOTAL*sizeof(double));
  summ=(double *)malloc(GTOTAL*sizeof(double));
  cum_prob_caseF=(double *)malloc(GTOTAL*sizeof(double));
  cum_prob_caseM=(double *)malloc(GTOTAL*sizeof(double));
  cum_prob_contF=(double *)malloc(GTOTAL*sizeof(double));
  cum_prob_contM=(double *)malloc(GTOTAL*sizeof(double));

  ////
  //// calculate beta for the disease loci
  //// see program manual for derivation of the formula
  for(i=0; i<NUMDL; i++) {
    // p[i] = frequency of allele 1 at disease locus i
    p[i] = freq[i][DLPOS[i]-1];

    // pi = frequency of risk allele at disease locus
    pi = (DV[i]==1) ? p[i] : 1-p[i];
    ri1 = GRR[i];
    ri2 = GRR2[i];

    if(Chr[i]<23) {
      expai = ((1-pi)*(1-pi) + ri1*2*pi*(1-pi) + ri2*pi*pi)/ PREV - 1;
    } else if(Chr[i]==23) {
      expai = (0.5*(1-pi)*(2-pi) + ri1*pi*(1-pi) + 0.5*ri2*pi*(1+pi))/ PREV -1;
    }
    tmpbeta1 = -log(((expai + 1)/ri1 - 1)/expai);
    tmpbeta2 = -log(((expai + 1)/ri2 - 1)/expai);

    // If there are no interaction effects, the beta coefficients can be 
    // calculated in a more meaningful way, depending on whether disease
    // allele is 0 or 1.
    // But if there are interaction effects, it is complicated and I follow
    // the way that allows simple addition of interaction terms.
    beta1[i] = tmpbeta1; beta2[i] = tmpbeta2;
  }

  ////
  //// calculate the components for the logit function
  for(i=0; i<GTOTAL; i++) {
    index_to_disease_geno(i, g, NUMDL);

    //// prob[i]:  frequency of genotype i (in population)
    prob[i]=1;
    for(j=0; j<NUMDL; j++) {
      if(Chr[j]<23) prob[i] *= Pgeno(g[j], p[j]);
      else prob[i] *= PgenoX(g[j], p[j]);
    }

    //// summ[i]:  Sum(beta*g) for genotype i
    summ[i]=0;
    for(j=0; j<NUMDL; j++) {
      if(DV[j]==1) summ[i] += beta1[j]*(g[j]==1) + beta2[j]*(g[j]==2);
      if(DV[j]==0) summ[i] += beta1[j]*(g[j]==1) + beta2[j]*(g[j]==0);
    }
    // additional interaction effect terms
    if(INTER) {
      for(j=0; j<NUMI2; j++) {
	idx1 = I2loci[j][0]-1; idx2 = I2loci[j][1]-1;
	eff11 = log(I2effects[j][0]); eff12 = log(I2effects[j][1]);
	eff21 = log(I2effects[j][2]); eff22 = log(I2effects[j][3]);
	//	printf("%f %f %f %f\n", eff11, eff12, eff21, eff22);
	if(DV[idx1]==1 & DV[idx2]==1)
	  summ[i] += eff11*(g[idx1]==1 & g[idx2]==1) +
	    eff12*(g[idx1]==1 & g[idx2]==2) +
	    eff21*(g[idx1]==2 & g[idx2]==1) +
	    eff22*(g[idx1]==2 & g[idx2]==2);
	if(DV[idx1]==0 & DV[idx2]==1)
	  summ[i] += eff11*(g[idx1]==1 & g[idx2]==1) +
	    eff12*(g[idx1]==1 & g[idx2]==2) +
	    eff21*(g[idx1]==0 & g[idx2]==1) +
	    eff22*(g[idx1]==0 & g[idx2]==2);
	if(DV[idx1]==1 & DV[idx2]==0)
	  summ[i] += eff11*(g[idx1]==1 & g[idx2]==1) +
	    eff12*(g[idx1]==1 & g[idx2]==0) +
	    eff21*(g[idx1]==2 & g[idx2]==1) +
	    eff22*(g[idx1]==2 & g[idx2]==0);
	if(DV[idx1]==0 & DV[idx2]==0)
	  summ[i] += eff11*(g[idx1]==1 & g[idx2]==1) +
	    eff12*(g[idx1]==1 & g[idx2]==0) +
	    eff21*(g[idx1]==0 & g[idx2]==1) +
	    eff22*(g[idx1]==0 & g[idx2]==0);
      }
    }
  }

  ////
  //// search for beta0
  b0hi = 100;  while( prevalence(b0hi, prob, summ) < PREV) b0hi = 2*b0hi;
  b0lo=-10000; while( prevalence(b0lo, prob, summ) > PREV) b0lo = 2*b0lo;
  do {
    b0tmp = (b0hi + b0lo)*.5;
    aa = prevalence(b0tmp, prob, summ);
    if(aa > PREV) b0hi=b0tmp;
    else b0lo=b0tmp;
  } while( fabs(aa/PREV - 1) > .001 );
  beta0 = b0tmp;

  ////
  //// Output disease model to standard output and file diseasemodel.txt
  sprintf(diseasemodel, "Disease model:\n");
  sprintf(tmpstr, "beta0 = %9.4f\n", beta0);
  strcat(diseasemodel, tmpstr);
  sprintf(tmpstr, "Locus  chr   #SNPs  DLpos  DV  DVFreq   GRR     GRR2    beta1    beta2\n");
  strcat(diseasemodel, tmpstr);
  for(i=0; i<NUMDL; i++) {
    sprintf(tmpstr, "   %2d   %2d  %6d %6d   %1d  %6.4f  %6.3f  %6.3f   %6.4f  %6.4f\n",
	    i+1, Chr[i], numMarker[i], DLPOS[i], DV[i],
	    DV[i] ? freq[i][DLPOS[i]-1] : 1-freq[i][DLPOS[i]-1],
	    GRR[i], GRR2[i], beta1[i], beta2[i]);
    strcat(diseasemodel, tmpstr);
  }
  if(INTER) {
    sprintf(tmpstr, "Interaction terms gamma_ij:\nPair  DL1  DL2  gamma11  gamma12  gamma21  gamma22\n");
    strcat(diseasemodel, tmpstr);
    for(j=0; j<NUMI2; j++) {
      idx1 = I2loci[j][0]-1; idx2 = I2loci[j][1]-1;
      eff11 = log(I2effects[j][0]); eff12 = log(I2effects[j][1]);
      eff21 = log(I2effects[j][2]); eff22 = log(I2effects[j][3]);
      sprintf(tmpstr, "%3d    %2d   %2d %6.3f   %6.3f   %6.3f   %6.3f\n",
	       j+1, idx1+1, idx2+1, eff11, eff12, eff21, eff22);
      strcat(diseasemodel, tmpstr);
    }
  }

  puts(diseasemodel);

  if( (modelout = fopen("diseasemodel.txt", "w")) == NULL) {
    printf("Cannot output to diseasemodel.txt.");
  } else {
    fprintf(modelout, diseasemodel);
  }
  fclose(modelout);

  ////
  //// calculate conditional genotype probabilities for cases and controls
  cum_disease_geno(beta0, NUMDL, p);
}


////////////////////////////////////////////////////////////////////
//
//  Simulate disease marker genotypes and flanking marker alleles
//  from -(WINDOW_SIZE-1)/2 to +(WINDOW_SIZE-1)/2
//
////////////////////////////////////////////////////////////////////
void sample_disease_markers(int NUMDL,
			    int *mkpos,  // disease markers positions
			    int WINDOW_SIZE,
			    int NUMCASEF,
			    int NUMCASEM,
			    int NUMCONTF,
			    int NUMCONTM,
			    int NUMCHROM,
			    int NUMCHROM_X
			    )
{
  int step2 = NUMCASEF + NUMCASEM;
  int step3 = step2 + NUMCONTF;
  int n_person = NUMCASEF + NUMCASEM + NUMCONTF + NUMCONTM;
  int n_hap = 2*n_person;
  int n_sample, dd[23], index, allele1, allele2, i, j, k, id;
  double p_rand;

  // simulate genotypes for case females
  for(i=0; i<NUMCASEF; i++){
    p_rand = ran1(&seed);
    if(p_rand < cum_prob_caseF[0]) {index=0;}
    else {
      for(j=1; j<GTOTAL; j++)
	if(p_rand < cum_prob_caseF[j] && p_rand >= cum_prob_caseF[j-1]) {
	  index=j;
	  break;
        }
    }
    index_to_disease_geno(index, dd, NUMDL);

    for(j=0; j<NUMDL; j++){
      if(dd[j]==0) {allele1=0; allele2=0;}
      else if(dd[j]==1) {allele1=0; allele2=1;}
      else {allele1=1; allele2=1;}
      sim_hap[j][2*i][mkpos[j]] = allele1;
      sim_hap[j][2*i+1][mkpos[j]] = allele2;
    }
  }

  // simulate genotypes for case males
  for(i=NUMCASEF; i<step2; i++){
    p_rand = ran1(&seed);
    if(p_rand < cum_prob_caseM[0]) {index=0;}
    else {
      for(j=1; j<GTOTAL; j++)
	if(p_rand < cum_prob_caseM[j] && p_rand >= cum_prob_caseM[j-1]) {
	  index=j;
	  break;
        }
    }
    index_to_disease_geno(index, dd, NUMDL);

    for(j=0; j<NUMDL; j++){
      if(dd[j]==0) {allele1=0; allele2=0;}
      else if(dd[j]==1) {allele1=0; allele2=1;}
      else {allele1=1; allele2=1;}
      sim_hap[j][2*i][mkpos[j]] = allele1;
      sim_hap[j][2*i+1][mkpos[j]] = allele2;
    }
  }

  // simulate genotypes for control females
  for(i=step2; i<step3; i++){
    p_rand=ran1(&seed);
    if(p_rand < cum_prob_contF[0]) {index=0;}
    else {
      for(j=1; j<GTOTAL; j++)
	if(p_rand < cum_prob_contF[j] && p_rand >= cum_prob_contF[j-1]) {
	  index=j;
	  break;
	}
    }
    index_to_disease_geno(index, dd, NUMDL);

    for(j=0; j<NUMDL; j++){
      if(dd[j]==0) {allele1=0; allele2=0;}
      else if(dd[j]==1) {allele1=0; allele2=1;}
      else {allele1=1; allele2=1;}
      sim_hap[j][2*i][mkpos[j]] = allele1;
      sim_hap[j][2*i+1][mkpos[j]] = allele2;
    }
  }

  // simulate genotypes for control males
  for(i=step3; i<n_person; i++){
    p_rand=ran1(&seed);
    if(p_rand < cum_prob_contM[0]) {index=0;}
    else {
      for(j=1; j<GTOTAL; j++)
	if(p_rand < cum_prob_contM[j] && p_rand >= cum_prob_contM[j-1]) {
	  index=j;
	  break;
	}
    }
    index_to_disease_geno(index, dd, NUMDL);

    for(j=0; j<NUMDL; j++){
      if(dd[j]==0) {allele1=0; allele2=0;}
      else if(dd[j]==1) {allele1=0; allele2=1;}
      else {allele1=1; allele2=1;}
      sim_hap[j][2*i][mkpos[j]] = allele1;
      sim_hap[j][2*i+1][mkpos[j]] = allele2;
    }
  }

  // find a haplotype at [n-2,n+2] that match the allele at the disease locus
  for(i=0; i<n_hap; i++){
    for(j=0; j<NUMDL; j++) {
      n_sample = (Chr[j]<23) ? NUMCHROM : NUMCHROM_X;
      // random draw until match
      do {
	id=int(ran1(&seed)*n_sample);
      } while(hap[j][id][mkpos[j]] != sim_hap[j][i][mkpos[j]]);
      for(k=mkpos[j]-(WINDOW_SIZE-1)/2;
	  k<=mkpos[j]+(WINDOW_SIZE-1)/2; k++)
	sim_hap[j][i][k]=hap[j][id][k];
    }
  }
}


/////////////////////////////////////////////////////////////////////
//
//  Simulate non-disease marker genotypes and flanking marker alleles
//  from -(WINDOW_SIZE-1)/2 to +(WINDOW_SIZE-1)/2
//
/////////////////////////////////////////////////////////////////////
void sample_start_marker(bit_vector hap[],     // original haplotypes
                         bit_vector sim_hap[], // simulated haplotypes
                         int markerPos,        // starting position
                         double p,             // starting marker allele freq
			 int WINDOW_SIZE,
                         int simPerson,        // number of people to simulate
                         int n_sample
                         )
{
  int i, j, id;
  int n_hap=2*simPerson;
  double p_rand;

  // simulate the two alleles at the start position
  // The following procedure works for both autosomes and X chromosome
  // freq00 = Pr(0/0) = (1-p)*(1-p)
  // freq01 = Pr(0/0, 0/1, but not 1/0) = (1-p)*(1-p) + p*(1-p)
  // freq10 = Pr(0/0, 0/1, 1/0) = (1-p)*(1-p) + 2*p*(1-p)
  double freq00=(1-p)*(1-p), freq01=1-p, freq10=1-p*p;
  for(i=0; i<simPerson; i++){
    p_rand = ran1(&seed);
    if(p_rand < freq00) {
      sim_hap[2*i][markerPos]=0; sim_hap[2*i+1][markerPos]=0;
    } else if(p_rand < freq01){
      sim_hap[2*i][markerPos]=0; sim_hap[2*i+1][markerPos]=1;
    } else if(p_rand < freq10){
      sim_hap[2*i][markerPos]=1; sim_hap[2*i+1][markerPos]=0;
    } else {
      sim_hap[2*i][markerPos]=1; sim_hap[2*i+1][markerPos]=1;
    }
  }

  // simulate the flanking marker alleles for each haplotype
  // randomly draw a haplotype until the alleles at marker position match
  for(i=0; i<n_hap; i++){
    do{
      id=int(ran1(&seed)*n_sample);
    } while(hap[id][markerPos] != sim_hap[i][markerPos]);
    for(j=markerPos-(WINDOW_SIZE-1)/2; j<=markerPos+(WINDOW_SIZE-1)/2; j++) {
      sim_hap[i][j]=hap[id][j];
    }
  }
}


/////////////////////////////////////////////////////////////////////
//
//  Simulate to the left and right ends of the entire chromosome
//    sample_hap:  main function for this section
//    match_hap:  match checking WINDOW_SIZE-1 alleles between two haplotypes
//
/////////////////////////////////////////////////////////////////////
// match checking WINDOW_SIZE-1 alleles between two haplotypes
//////////////////////////////////////////////////////////////
bool match_hap(bit_vector *sim_hap,  // simulated haplotypes 
	       bit_vector *hap_s,    // original data
	       long mlist,           // start position for match checking
	       int WINDOW_SIZE
	       )
{
  for (int i=0; i<WINDOW_SIZE-1; i++) {
    if ( (*hap_s)[mlist+i] != (*sim_hap)[mlist+i] ) 
      return(0);
  }
  return(1);
}

///////////////////////////////////////////////////////////////
// simulate to the left and right ends of the entire chromosome
///////////////////////////////////////////////////////////////
void sample_hap(bit_vector sim_hap[], // simulated haplotypes 
		bit_vector hap[],     // origianl data
		int n_hap,            // # chromosomes to simulate
		int markerPos,        // algorithm start position
		int startpos,         // start position
		int endpos,           // end position
		int WINDOW_SIZE,
		int n_sample
                )
{
  int i, j, k, id;

  // ->genome end
  for (j=markerPos+(WINDOW_SIZE-1)/2+1; j<=endpos; j++) {
    for (i=0; i<n_hap; i++) {
      // randomly draw a haplotype to match the first WINDOW_SIZE-1 markers
      // match checking from j-4 to j-1 (if WINDOW_SIZE=5)
      do {
	id = int(ran1(&seed)*n_sample);
      } while( !match_hap(&sim_hap[i], &hap[id], j-(WINDOW_SIZE-1),
			  WINDOW_SIZE) );
      sim_hap[i][j] = hap[id][j];
    }
  }

  // genome begining<-
  for (j=markerPos-(WINDOW_SIZE-1)/2-1; j>=startpos; j--) {
    for (i=0; i<n_hap; i++) {
      // randomly draw a haplotype to match the first WINDOW_SIZE-1 markers
      // match checking from j+1 to j+4 (if WINDOW_SIZE=5)
      do {
	id = int(ran1(&seed)*n_sample);
      } while( !match_hap(&sim_hap[i], &hap[id], j+1, WINDOW_SIZE) );
      sim_hap[i][j] = hap[id][j];        
    }
  }
}


////////////////////////////////////////////////
//
//  Simulate data
//
////////////////////////////////////////////////
void simulate(int NUMCASEF, int NUMCASEM, int NUMCONTF, int NUMCONTM,
	      int NUMDL, int REGIONAL, int DLPOS[], int WINDOW_SIZE,
	      int NUMCHROM, int NUMCHROM_X)
{
  int i, j, k, markerPos, countMarker=0, xx;
  int simPerson=NUMCASEF+NUMCASEM+NUMCONTF+NUMCONTM;
  int mkpos[23];

  int numPerson; // 60 for autosomal chromosomes, 45 for X chromosome
  int n_sample;
  int IsXchr;    // 1: X chromosome, 0: autosomal chromosomes

  for(i=0; i<23; i++) {
    for(j=0; j<2*simPerson; j++) sim_hap[i][j].reserve(numMarker[i]);
  }
  for(i=0; i<NUMDL; i++) mkpos[i]=DLPOS[i]-1;

  ////////////////////////////////
  //  Simulate disease chromosomes
  ////////////////////////////////
  if(NUMDL>0) {
    puts("Simulating disease chromosomes ...");

    // simulate genotypes at disease loci and flanking markers
    sample_disease_markers(NUMDL, mkpos, WINDOW_SIZE,
			   NUMCASEF, NUMCASEM, NUMCONTF, NUMCONTM,
			   NUMCHROM, NUMCHROM_X);

    // simulate the rest of the disease chromosomes
    for(i=0; i<NUMDL; i++){
      if(REGIONAL==0) {
	printf("  Simulating chromosome %d (%d SNPs) ...\n",
	       Chr[i], numMarker[i]);
	sample_hap(sim_hap[i], hap[i], 2*simPerson,
		   mkpos[i], 0, numMarker[i]-1, WINDOW_SIZE,
		   Chr[i]<23 ? NUMCHROM : NUMCHROM_X);
      } else {
	printf("  Simulating chromosome %d (%d SNPs), positions %d-%d ...\n",
	       Chr[i], numMarker[i], Start[i], End[i]);
	sample_hap(sim_hap[i], hap[i], 2*simPerson,
		   mkpos[i], Start[i]-1, End[i]-1, WINDOW_SIZE,
		   Chr[i]<23 ? NUMCHROM : NUMCHROM_X);
      }
    }
  }

  //////////////////////////////////
  //  Simulate remaining chromosomes
  //////////////////////////////////
  if(NUMDL==0) 
    puts("\nSimulating the genome ...");
  else if(REGIONAL==0)
    puts("Simulating remaining non-disease chromosomes ...");

  if(NUMDL==0 | REGIONAL==0) {
    for(i=NUMDL; i<23; i++) {
      printf("  Simulating chromosome %d (%d SNPs) ...\n", Chr[i], numMarker[i]);

      // randomly pick a starting marker and then simulate the rest
      markerPos=int(ran1(&seed)*(numMarker[i]-WINDOW_SIZE))+(WINDOW_SIZE-1)/2;
      sample_start_marker(hap[i], sim_hap[i], markerPos, freq[i][markerPos],
			  WINDOW_SIZE, simPerson,
			  Chr[i]<23 ? NUMCHROM : NUMCHROM_X);
      sample_hap(sim_hap[i], hap[i], 2*simPerson,
		 markerPos, 0, numMarker[i]-1, WINDOW_SIZE,
		 Chr[i]<23 ? NUMCHROM : NUMCHROM_X);
    }
  }

  // Chr X males:  second chromosome reset to 0
  // find out which i has Chr[i]==23
  for(i=0; i<23; i++) if(Chr[i]==23) xx=i;
  for(j=NUMCASEF; j<simPerson; j++) {
    if(j<NUMCASEF+NUMCASEM | j>=NUMCASEF+NUMCASEM+NUMCONTF)
      for(k=0; k<numMarker[xx]; k++)
	sim_hap[xx][2*j+1][k]=0;
  }
}


////////////////////////////////////////////////
//
//  Output data to files
//
////////////////////////////////////////////////
void output(int NUMCASEF, int NUMCASEM, int NUMCONTF, int NUMCONTM,
	    char *OUTFORMAT, int NUMDL, int REGIONAL)
{
  int i, j, k, sex, status, startpos, endpos, fileexist=0;
  char outfilename[1024], testfilename[1024], command[1024];
  ofstream outfile;
  ifstream testfile;

  int step2 = NUMCASEF + NUMCASEM;
  int step3 = step2 + NUMCONTF;
  int n_person = NUMCASEF + NUMCASEM + NUMCONTF + NUMCONTM;
  int n_hap = 2*(NUMCASEF + NUMCASEM + NUMCONTF + NUMCONTM);

  // remove old chr#.dat or chr#.dat.gz files
  // This is useful especially for regional simulations to avoid having 
  // files for chromosomes that are not simulated in the last run.
  for(i=0; i<23; i++) {
    sprintf(testfilename, "chr%d.dat", Chr[i]);
    testfile.open(testfilename);
    if(testfile.is_open()) {
      if(fileexist==0) {
	puts("\nRemoving old data files chr#.dat or chr#.dat.gz");
	fileexist=1;
      }
      testfile.close();
      if(remove(testfilename) !=0)
	cout<<"Remove operation failed for " << testfilename << endl;
    }

    sprintf(testfilename, "chr%d.dat.gz", Chr[i]);
    testfile.open(testfilename);
    if(testfile.is_open()) {
      if(fileexist==0) {
	puts("\nRemoving old data files chr#.dat or chr#.dat.gz");
	fileexist=1;
      }
      testfile.close();
      if(remove(testfilename) !=0)
	cout<<"Remove operation failed for " << testfilename << endl;
    }
  }


  // write data to files
  puts("\nWriting to files ...");
  for(i=0; i<23; i++) {
    if(NUMDL==0 | REGIONAL ==0 | i<NUMDL) {
      startpos = 0; endpos = numMarker[i]-1;
      if(NUMDL>0 & REGIONAL==1 & i<NUMDL) {
	startpos=Start[i]-1; endpos=End[i]-1;
      }

      sprintf(outfilename, "chr%d.dat", Chr[i]);
      outfile.open(outfilename);
      if(outfile.fail()){
	cout<< "Output file: "<< outfilename << " cannot be opened!"; 
	exit(1);
      }

      // linkage format
      if(!strcmp(OUTFORMAT, "linkage")) {
	for(j=0; j<n_person; j++) {
	  // sex:  1=male, 2=female
	  sex = (j < NUMCASEF | (j >= step2 & j < step3)) + 1;
	  // affection status:  0=unknown, 1=unaffected, 2=affected
	  // for population sampling, all subjects' status is 0
	  status = (NUMDL>0) ? (j < NUMCASEF + NUMCASEM) + 1 : 0;
	  outfile << j+1 << " 1 0 0 " << sex << " " << status << " ";
	  // alleles:  0 -> "1", 1 -> "2"
	  for(k=startpos; k<=endpos; k++)
	    outfile << sim_hap[i][2*j][k]   +1 << " "
		    << sim_hap[i][2*j+1][k] +1 << " ";
	  outfile << "\n";
	}
      }

      // genotype format
      if(!strcmp(OUTFORMAT, "genotype")) {
	if(Chr[i]<23) {
	  // for autosomes, genotype = 0,1,2
	  for(j=0; j<n_person; j++) {
	    for(k=startpos; k<=endpos; k++)
	      outfile << sim_hap[i][2*j][k] + sim_hap[i][2*j+1][k] << " ";
	    outfile << "\n";
	  }
	} else {
	  // for X: female genotype=0,1,2, male genotype=0,2
	  for(j=0; j<n_person; j++) {
	    if(j < NUMCASEF | (j >= step2 & j < step3)) {
	      for(k=startpos; k<=endpos; k++)
		outfile << sim_hap[i][2*j][k] + sim_hap[i][2*j+1][k] << " ";
	    } else {
	      for(k=startpos; k<=endpos; k++)
		outfile << 2*sim_hap[i][2*j][k] << " ";
	    }
	    outfile << "\n";
	  }
	}
      }

      // phased format
      if(!strcmp(OUTFORMAT, "phased")) {
	for(j=0; j<n_hap; j++) {
	  for(k=startpos; k<=endpos; k++)
	    outfile << sim_hap[i][j][k] << " ";
	  outfile << "\n";
	}
      }

      outfile.close();

      sprintf(command, "gzip -f %s", outfilename, outfilename);
      if (std::system(0)) {
	// A command processor is available.
	std::system(command);
      } else {
	cout << "A command processor is not available.\n";
      }
    }
  }
}


///////////////////////////////////////////////////////////////////
//  Include the data analysis functions from another file
//    Uncomment the following line if you want to test the example
//    analysis function included in the software distribution.
///////////////////////////////////////////////////////////////////
//#include "dataanalysis.cpp"


////////////////////////////////////////////////
//
//  The main function:  The command is like
//     GWAsimulator <control file> <seed>
////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  char* CONTFILE;         // name of the control file

  if(argc != 3) {
    puts("\nSyntax:  GWAsimulator <control file> <seed>\n");
    puts("  Control file contains key parameters for the simulations.");
    puts("\tSee the example control file or the program manual for details.\n");
    puts("  The seed is a positive integer.\n");
    exit(1);
  }

  CONTFILE = argv[1];
  seed = -abs(atol(argv[2]));

  //////////////////////////
  // defined in control file
  //////////////////////////
  char DATADIR1[1024];    // input phased data name path prefix
  char DATADIR2[1024];    // input phased data name path affix
  int NUMCHROM;           // number of phased autosomes in input data
  int NUMCHROM_X;         // number of phased X chromosomes in input data
  int NEEDOUTPUT;         // indicator for writing simulated data to files
  char OUTFORMAT[128];    // output format: HapMap or LINKAGE
  int WINDOW_SIZE = 5;    // window size in simulation (default to 5)
  int NUMCASEF, NUMCASEM; // numbers of cases (females and males)
  int NUMCONTF, NUMCONTM; // numbers of controls (females and males)
  double PREV;            // disease model prevalence
  int NUMDL;              // number of diseaes loci
  int REGIONAL;           // indicator for regional simulation if NUMDL>0
  int DLPOS[23];          // disease loci positions
  int DV[23];             // disease variant/allele:  0 or 1
  double GRR[23];         // genotypic relative risk of genotypes 1 vs 0
  double GRR2[23];        // genotypic relative risk of genotypes 2 vs 0
  int INTER=0;            // indicator for interactive effects
  int NUMI2=0;            // number of 2-way interaction terms

  // read in control file
  controlin(CONTFILE, DATADIR1, DATADIR2, NUMCHROM, NUMCHROM_X,
	    NEEDOUTPUT, OUTFORMAT,
	    WINDOW_SIZE, NUMCASEF, NUMCASEM, NUMCONTF, NUMCONTM,
	    PREV, NUMDL, REGIONAL, DLPOS, DV, GRR, GRR2, INTER, NUMI2);

  // read in phased data
  data_in(DATADIR1, DATADIR2, NUMCHROM, NUMCHROM_X,
	  NUMDL, REGIONAL, DLPOS, WINDOW_SIZE);

  // build disease model given user's input
  if(NUMDL>0) modelbuild(PREV, NUMDL, DLPOS, DV, GRR, GRR2, INTER, NUMI2);

  // start simulation
  simulate(NUMCASEF, NUMCASEM, NUMCONTF, NUMCONTM,
	   NUMDL, REGIONAL, DLPOS, WINDOW_SIZE, NUMCHROM, NUMCHROM_X);

  // output data to files
  if(NEEDOUTPUT) output(NUMCASEF, NUMCASEM, NUMCONTF, NUMCONTM,
			OUTFORMAT, NUMDL, REGIONAL);

  //////////////////////////////////////////////////////////////////////
  // Data analysis:  replace with your programs
  //   Need to include data analysis programs before the main() function
  //   Uncomment the following line if you want to test the example
  //   analysis function included in the software distribution.
  //////////////////////////////////////////////////////////////////////
  //dataanalysis(NUMCASEF, NUMCASEM, NUMCONTF, NUMCONTM);

  // free up memory
  free(cum_prob_caseF); free(cum_prob_caseM); free(cum_prob_contF);
  free(cum_prob_contM); free(prob); free(summ);
  return(0);
}
