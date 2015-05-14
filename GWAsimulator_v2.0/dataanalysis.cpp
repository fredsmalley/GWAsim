/////////////////////////////////////////////////////////////////////
//
//  Data analysis example function for GWAsimulator
//  Calculate chi-squared (df = 1) allelic association test statistic
//
//  version 1.0 (07/16/2007)
//  Authors:  Chun Li, Mingyao Li
/////////////////////////////////////////////////////////////////////

double *chisq;          // a pointer to a sequence of chi-squared statistics

void assoc_chisq(int chr,
		 bit_vector sim_hap[], 
                 double *chisq,
                 long numMarker,
                 int NUMCASEF, int NUMCASEM, int NUMCONTF, int NUMCONTM)
{
  int i, j, n0, n1, m0, m1;
  double numerator, denominator;
  int step2 = NUMCASEF + NUMCASEM;
  int step3 = step2 + NUMCONTF;
  int TOTAL = step3 + NUMCONTM;
  int caseallele = (chr == 23) ? 2*NUMCASEF+NUMCASEM : 2*(NUMCASEF+NUMCASEM);
  int contallele = (chr == 23) ? 2*NUMCONTF+NUMCONTM : 2*(NUMCONTF+NUMCONTM);

  for(j=0; j<numMarker; j++) {
    // count the number of allele 1 at each SNP for cases
    n1=0; 
    for(i=0; i<step2; i++) {
      if(sim_hap[2*i][j]==1) n1++;
      if(chr < 23 | i < NUMCASEF) {
	if(sim_hap[2*i+1][j]==1) n1++;
      }
    }
    n0 = caseallele - n1;

    // count the number of allele 1 at each SNP for controls
    m1=0;
    for(i=step2; i<TOTAL; i++) {
      if(sim_hap[2*i][j]==1) m1++;
      if(chr < 23 | i < step3) {
	if(sim_hap[2*i+1][j]==1) m1++;
      }
    }
    m0 = contallele - m1;

    // calculate the chi-squared test statistic
    numerator = 1.0*(caseallele+contallele)*(n0*m1-n1*m0)*(n0*m1-n1*m0);
    denominator = 1.0*caseallele*contallele*(n0+m0)*(n1+m1);
    chisq[j] = (fabs(denominator)<1e-8) ? 0 : (numerator/denominator);
  }
}


void dataanalysis(int NUMCASEF, int NUMCASEM, int NUMCONTF, int NUMCONTM)
{
  int i, countMarker, NUMMARKERTOTAL;

  NUMMARKERTOTAL=0;
  for(i=0; i<23; i++) NUMMARKERTOTAL += numMarker[i];

  chisq=(double *)malloc(NUMMARKERTOTAL*sizeof(double));

  countMarker=0;
  for(i=0; i<23; i++){
    assoc_chisq(Chr[i], sim_hap[i], chisq+countMarker, numMarker[i],
		NUMCASEF, NUMCASEM, NUMCONTF, NUMCONTM);
    countMarker += numMarker[i];
  }

  //  for(i=900; i<1100; i++) cout << chisq[i] << "\n";
}
