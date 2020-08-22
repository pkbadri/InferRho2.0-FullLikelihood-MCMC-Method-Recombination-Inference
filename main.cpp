/* 
Copyright 2009 Ying Wang (ygwang@uchicago.edu)
Created 2007 at University of California, Davis.
*/
#include "main.h"

void usage()
{
  cout<<"\n  IR - Inference of Recombination Rates and Hotspots"<<endl;
  cout<<"  version: 1.0.0\n"<<endl;
  cout<<"  Usage: IR <inputfile> [options]\n"<<endl;
  cout<<"  Options: -burnin n    Number of burnin iterations. Default is 50000."<<endl;
  cout<<"           -iter n      Number of iterations. Default is 100000."<<endl;
  cout<<"           -thin n      Thinning iterations. Default is 10."<<endl;
  cout<<"           -dTheta x    Mixing parameter for theta. Default is 0.2."<<endl;
  cout<<"           -dRho x      Mixing parameter for background recombination rate. Default is 100."<<endl;
  cout<<"           -dScale x    Mixing parameter for background recombination rate. Default is 0.1."<<endl;
  cout<<"           -nChains n   Number of chains. Default is 4."<<endl;
  cout<<"           -temp x      Temperature of chains. Default is 1.05."<<endl;
  cout<<"           -seed x      Seed. Default seed is set based on system time."<<endl;
  exit(1);
}




void parseOptions(int argc, char* argv[], int& burnin, int& iter, int& thin, double& dTheta, double& dRho, double& dScale, int& nChains, double& temp, long& seed)
{
  burnin = 0;
  iter = 1000000;
  thin = 10;
  dTheta = 0.2;
  dRho = 100;
  dScale = 0.1;
  nChains = 4;
  temp = 1.05;
  seed = time(NULL);

  //Check Inputfile
  ifstream ifs;
  ifs.open(argv[1],ios::in);
  if (!ifs){
  cout<<"The inputfile: "<<argv[1]<<" cannot be opened."<<endl;
  exit(0);
  }
  ifs.close();

  for(int i = 2; i < argc; i += 2)
  {
  string op = argv[i];
  if(op.compare("-burnin") == 0) burnin = atoi(argv[i+1]);
  else if(op.compare("-iter") == 0) iter = atoi(argv[i+1]);
  else if(op.compare("-thin")==0) thin = atoi(argv[i+1]);
  else if(op.compare("-dTheta")==0) dTheta = atof(argv[i+1]);
  else if(op.compare("-dRho")==0) dRho = atof(argv[i+1]);
  else if(op.compare("-dScale")==0) dScale = atof(argv[i+1]); 
  else if(op.compare("-nChains")==0) nChains = atoi(argv[i+1]); 
  else if(op.compare("-temp")==0) temp = atof(argv[i+1]);
  else if(op.compare("-seed")==0) seed = atoi(argv[i+1]);
  else{cout<<"Unrecognized option. Exit."<<endl; exit(0);}
  }


  //Simple Option Checking
  if(burnin < 1000000 || iter < 1000000) {cout<<"warning: small number of burnin and iterations."<<endl;} 
  if(thin < 300) {cout<<"warning: small number of thinning iterations."<<endl;} 
  if(dTheta < 0.05 || dTheta > 1 || dScale < 0.02 || dScale > 1 || dRho < 1 || dRho > 200)
  {cout<<"warning: unexpected mixing parameter values. Need to check the acceptance rate."<<endl;}
  if(nChains < 1) {cout<<"nChains has to be > 0. Exit."<<endl; exit(0);} 
  if(nChains > 5) {cout<<"warning: higher number of chains will slow down the program."<<endl; exit(0);}
  if(temp <= 1)   {cout<<"Incorrect chain temperature. It has to be > 1."<<endl; exit(0);}
  if(temp > 1.5)  {cout<<"warning: higher temperature will make the chain swapping rate low."<<endl;}
  cout<<endl;
} 




int main(int argc, char* argv[])
{   
   if(argc < 2) usage();  
   gsl_rng* gslr;
   gslr = gsl_rng_alloc(gsl_rng_taus2);
 
   clock_t start, end;
   start = clock();

   //Data and Parameters
   int ngraphs;
   char datatype;
   int* nseqs; int* nmarkers;  
   int*** seqmatrix;
   double** phydis;
   double** scaleddis;
   double* length;
   int*** SNPs;
   int burnin;
   int iter;
   int thinning;
   double* stafreq;  
   double deltaTheta;
   double deltaRhoBackg;
   double deltaRhoBackgScale;
   double deltaGamma;	
   int nchains;
   double* beta;
   double deltaTemp;
   double startint;
   double endint;
   long seed;
   
   parseOptions(argc, argv, burnin, iter, thinning, deltaTheta, deltaRhoBackg, deltaRhoBackgScale, nchains, deltaTemp, seed);
   gsl_rng_set(gslr, seed);
   deltaGamma = 100;

   string inputFileName = argv[1];
 
   ifstream ifs;
   ifs.open(argv[1],ios::in);
   if (!ifs){
   cout<<"The inputfile ("<<inputFileName<<") cannot be opened."<<endl;
   exit(0);
   }

   ofstream ofsLog;
   string ofsLogF = inputFileName+".log";
   ofsLog.open(ofsLogF.c_str());
   if (!ofsLog){
   cout<<"Cannot open file "<<inputFileName+".log"<<"."<<endl;
   exit(0);
   }

   ofsLog<<"seed: "<<seed<<endl;
 
   ifs>>ngraphs;  
   nseqs = new int[ngraphs];
   nmarkers = new int[ngraphs];
   length = new double[ngraphs];
   phydis = new double*[ngraphs];
   scaleddis = new double*[ngraphs];
   seqmatrix = new int**[ngraphs];
   SNPs = new int**[ngraphs];
   stafreq = new double[4];
   
   ifs>>datatype;
   ofsLog<<"data type: "<<datatype<<endl; //g(genotype) or h(haplotype)

   for(int tao=0; tao<ngraphs; tao++)
   {
   ifs>>nseqs[tao]; //# individuals for the locus 
   if(datatype=='g' || datatype=='G') nseqs[tao] *= 2;
   ofsLog<<"nseqs: "<<nseqs[tao]<<endl;
   ifs>>nmarkers[tao];   //# markers
   ofsLog<<"nmarkers: "<<nmarkers[tao]<<endl;

   seqmatrix[tao] = new int*[nseqs[tao]];
   
   for(int k=0; k<nseqs[tao]; k++){
  
   seqmatrix[tao][k] = new int[nmarkers[tao]];
   
   ofsLog<<k<<"---";
   for(int j=0; j<nmarkers[tao]; j++)
   { 
   int ch = ifs.get();
   while(ch != (int)'A' && ch != (int)'T' && ch != (int)'G' && ch != (int)'C' && ch != (int)'N')
   ch = ifs.get();
   seqmatrix[tao][k][j]=ch; ofsLog<<(char)seqmatrix[tao][k][j]<<" ";
   }
   ofsLog<<endl;
   } 
   ofsLog<<endl;

   try{SNPs[tao] = new int*[nmarkers[tao]];}
   catch ( bad_alloc ex )  {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}
   for(int j=0; j<nmarkers[tao]; j++)
   {
   try{SNPs[tao][j]=new int[2];}  
   catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}

   int i=0;
   while(i<nseqs[tao] && seqmatrix[tao][i][j] == (int)'N') {i++;} 
   if(i==nseqs[tao]) {SNPs[tao][j][0]='A';} 
   else {SNPs[tao][j][0] = seqmatrix[tao][i][j];}
   int found =0;
   i++;
  
   for(; i<nseqs[tao]; i++) 
   {
   if(seqmatrix[tao][i][j] != (int)'N' && seqmatrix[tao][i][j] != SNPs[tao][j][0]) 
   {SNPs[tao][j][1]=seqmatrix[tao][i][j]; found=1; break;}
   } 
   if(found == 0) {SNPs[tao][j][1]= Util::randNuc(SNPs[tao][j][0], gslr); }
   } 

   for(int j=0; j<nmarkers[tao]; j++)
   {
   ofsLog<<"site"<<j<<": "<<(char)SNPs[tao][j][0]<<"\t"<<(char)SNPs[tao][j][1]<<endl;
   }

   try{phydis[tao] = new double[nmarkers[tao]];}
   catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}
   for(int i=0; i<nmarkers[tao]; i++)
   {
   ifs>>phydis[tao][i]; ofsLog<<phydis[tao][i]<<" "; 
   if(i>0 && phydis[tao][i]<=phydis[tao][i-1])
   {
   phydis[tao][i] = phydis[tao][i-1]+1;
   }
   }
   ofsLog<<endl;

   length[tao] = phydis[tao][nmarkers[tao]-1] - phydis[tao][0];

   try{scaleddis[tao]=new double[nmarkers[tao]];}
   catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}
   for(int i=0; i<nmarkers[tao]; i++)
   {
   scaleddis[tao][i] = (phydis[tao][i]-phydis[tao][0])/length[tao]; 
   ofsLog<<scaleddis[tao][i]<<" "; 
   if(i>0 && scaleddis[tao][i]-scaleddis[tao][i-1] < 1e-6 ) { cout<<"minimum distance is 1e-6 between markers."<<endl; exit(0);}
   }
   ofsLog<<endl;
   ofsLog<<"length: "<<length[tao]<<endl;
   }
   ifs.close();

   startint = phydis[0][0];
   endint = phydis[ngraphs-1][nmarkers[ngraphs-1]-1];
   ofsLog<<"Interval: "<<startint<<"\t"<<endint<<endl;


   for(int i=0; i<4; i++) stafreq[i] = 0.25;//options for changing the stationary distribution 
   ofsLog<<"Stationary distribution of nucleotides:"<<endl;
   ofsLog<<"\t"<<stafreq[0]<<"\t"<<stafreq[1]<<"\t"<<stafreq[2]<<"\t"<<stafreq[3]<<endl;

   ofsLog<<"Burnin, iteration, thinning parameters: "<<endl;
   ofsLog<<"\t"<<burnin<<"\t"<<iter<<"\t"<<thinning<<endl;

   //mixing parameters
   ofsLog<<"Mixing parameters: deltaTheta, deltaRhoBackg, delteRhoBackgScale, deltaGamma: "<<endl;
   ofsLog<<"\t"<<deltaTheta<<"\t"<<deltaRhoBackg<<"\t"<<deltaRhoBackgScale<<"\t"<<deltaGamma<<endl;

   //remember which one is missing data, and initialize missing data
   ofsLog<<"missing data: "<<endl;
   vector<missSnp>* missing = new vector<missSnp>;
   for(int tau=0; tau<ngraphs; tau++)
   {
      for(int k=0; k<nseqs[tau]; k++)
      {
         for(int j=0; j<nmarkers[tau]; j++)
         {
            if(seqmatrix[tau][k][j]=='N')
            {
               missSnp one; one.graph = tau; one.chrid = k; one.locus = j;
               ofsLog<<"site: "<<tau<<", "<<k<<", "<<j<<endl;
               missing->push_back(one);
               seqmatrix[tau][k][j]  =  SNPs[tau][j][0];
            }
         }
      }
   }
   ofsLog<<"Total missing sites: "<<missing->size()<<endl;
//------------------------------finish processing input data---------------------------------




//----------------------initializing before MC^3---------------------------------------------
ofsLog<<"# chains: "<<nchains<<endl;
beta = new double[nchains];
ofsLog<<"delta Temperature: "<<deltaTemp<<endl;

Chain** chains = new Chain*[nchains];
for(int ch = 0; ch < nchains; ch++){
beta[ch] = 1.0/pow(deltaTemp, ch);
ofsLog<<"beta for chain "<<ch<<": "<<beta[ch]<<endl;
chains[ch] = new Chain(beta[ch], ngraphs, nseqs, nmarkers, phydis, scaleddis, seqmatrix, length, startint, endint, stafreq, missing, SNPs, gslr, ch);
chains[ch]->setupMixingPar(deltaTheta, deltaRhoBackg, deltaRhoBackgScale, deltaGamma);
}
ofsLog<<endl;
  
ofsLog<<"Initial rho's: "<<endl;
chains[0]->printRho(ofsLog, 0); 
ofsLog<<"Number of chains: "<<nchains<<endl;  




/*-------------------------------------------------------------------------------------------
------------------------------------MCMCMC---------------------------------------------------
---------------------------------------------------------------------------------------------*/
int nsamples = 0;
//Open output files
//--->posterior distribution of haplotypes
chrHaps** postHap;  
try{postHap = new chrHaps*[ngraphs];}
catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}
for(int i=0; i<ngraphs; i++)
{
try{postHap[i] = new chrHaps[nseqs[i]];}
catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}
}


//--->file for saving the log-prior and log-likelihood
ofstream ofsMon;
string ofsMonName = inputFileName+".monitor";
ofsMon.open(ofsMonName.c_str());
if (!ofsMon){
cout<<"Cannot open file "<<ofsMonName<<"."<<endl;
exit(1);
} 


//header: 
ofsMon<<"iterations\tlog-prior\tlog-likelihood\tsum"<<endl;

  
//--->file for saving recombination hotspots
ofstream ofsHot;
string ofsHotName = inputFileName+".hotsp";
ofsHot.open(ofsHotName.c_str());
if (!ofsHot){
cout<<"Cannot open file "<<inputFileName+".hotsp"<<"."<<endl;
exit(1);
}


//--->file for saving recombination events (# recombinations within each bin) 
ofstream ofsRec;
string ofsRecName = inputFileName+".rec";
ofsRec.open(ofsRecName.c_str());
if (!ofsRec){
cout<<"Cannot open file "<<inputFileName+".rec"<<"."<<endl;
exit(1);
}

  
//--->file for saving rhoBackg
ofstream ofsRhoBackg;
string ofsRhoBackgName = inputFileName+".rhoBackg";
ofsRhoBackg.open(ofsRhoBackgName.c_str());
if (!ofsRhoBackg){
cout<<"Cannot open file "<<ofsRhoBackgName<<"."<<endl;
exit(1);
}


//header: midpoint between SNPs
for(int i=0; i<ngraphs; i++){
for(int j=0; j<nmarkers[i]-1; j++){
ofsRhoBackg<<(phydis[i][j]+phydis[i][j+1])*0.5;
if(i == ngraphs-1  &&  j == nmarkers[i]-2) {}
else {ofsRhoBackg<<"\t";}
}}
ofsRhoBackg<<endl;


//--->file for saving rho
ofstream ofsRho;
string ofsRhoName = inputFileName+".rho";
ofsRho.open(ofsRhoName.c_str());
if (!ofsRho){
cout<<"Cannot open file "<<ofsRhoName<<"."<<endl;
exit(1);
}


//header: SNP locations
for(int i=0; i<ngraphs; i++){
for(int j=0; j<nmarkers[i]-1; j++){
ofsRho<<(phydis[i][j]+phydis[i][j+1])*0.5;
if(i == ngraphs-1 && j == nmarkers[i]-2) {}
else {ofsRho<<"\t";}
}
}
ofsRho<<endl;

//chains[0]->monitor(chains, nchains, 0, ofsMon);
	
//ASSIGN PROBABILITIES FOR EACH OF THE MOVES IN THE METROPOLIS HASTINGS ALGORITHM
double prGen, prTree, prMiss, prHap, prAnc, prRhoBackg, prHotsp, prTheta;
if(missing->empty()) prMiss = 0.00; else prMiss = 0.02;
if(datatype == 'g' || datatype == 'G') prHap = 0.05; else prHap = 0.00;
prAnc = 0.00; prTree = 0.00; prRhoBackg = 0.40; prHotsp = 0.00; prTheta = 0.10;
prGen = 1 - (prTree + prRhoBackg + prHotsp + prTheta + prMiss + prAnc + prHap);

ofsLog<<prGen<<" "<<prTree<<" "<<prMiss<<" "<<prHap<<" "<<prAnc<<" "<<prRhoBackg<<" "<<prHotsp<<" "<<prTheta<<endl<<endl<<endl;
ofsLog.close();


//RUN THE MCMCMC ALGORITHM WITH SWAPPING BETWEEN THE CHAINS
for(int chain = 0; chain < nchains; chain++){chains[chain]->setupProb(prGen, prTree, prMiss, prHap, prAnc, prRhoBackg, prHotsp, prTheta);}


//for omp, added on Jan 23, 2010
double wtime;
int ch;
wtime = omp_get_wtime();
cout<<"Number of processors available = "<<omp_get_num_procs ( )<<endl;
cout<<"Number of threads = "<<omp_get_max_threads ( )<<endl;
omp_set_num_threads(nchains);

//cout<<"iter Pr(gen) Pr(tree) Pr(miss) Pr(hap) Pr(anc) Pr(RhoBgS) Pr(rhoBg) Pr(hotspC) Pr(hotspA) Pr(hotspD) Pr(theta) Pr(swap) theta rhoBGScale TMRCA"<<endl;

int index[nchains], which; for(which = 0;which < nchains; which++){index[which] = which;}
     
for(int num = 1; num <= burnin + iter; num++)//BURNIN PLUS ITERATION LOOP
{  
#pragma omp parallel for private(ch) 
for(int ch = 0; ch < nchains; ch++){
chains[index[ch]]->modifyChain(int(drand48()*SWINT));
}
  
#pragma omp critical  
if(nchains > 1 && num > 1)//SWAP CHAINS WITH CERTAIN PROBABILITY AFTER EACH ITERATION.
{
int chj = gsl_rng_uniform_int(gslr, nchains-1); 
int chk = chj + 1;
chj = index[chj]; chk = index[chk];
chains[chj]->swapChains(chains[chk]);
}


for(which = 0; which <nchains; which++){
for(int iii = 0; iii < nchains; iii++){
if(chains[which]->getbeta() == beta[iii]){index[iii] = which; break;}
}
}

    
//SAMPLE AFTER THE BURNIN PERIOD BASED ON THINNING INTERVALS
if(num > burnin && num%thinning == 0){
nsamples++;    
chains[index[0]]->printRho(ofsRho, num);
chains[index[0]]->printRhoBackg(ofsRhoBackg, num); 
//chains[index[0]]->printRec(ofsRec, num);
chains[index[0]]->printHotspots(ofsHot, num);  
//if(datatype == 'g' || datatype == 'G'){
//chains[index[0]]->sampleHaplotypes(postHap);
//}
}

    
//OUTPUT TO SCREEN AFTER MULTIPLES OF THE THINNING INTERVAL
if(num > 0 && num%thinning == 0)
{
//output acceptance probability and hotspots to screen
//cout<<num<<endl;
//for(int i=0; i<nchains; i++){chains[i]->outputAccProb(); cout<<endl;}
//cout<<chains[0]->getTheta()<<" "<<chains[0]->getRhoBackgScale()<<endl;
//chains[0]->printHotspots();
//chains[0]->monitor(chains, nchains, num, ofsMon);
} 
}

wtime = omp_get_wtime()  -  wtime;
cout<<"wtime: "<<wtime<<endl;


//finish burnin+iteration
ofsRec.close();
ofsRho.close();
ofsRhoBackg.close();
ofsMon.close(); 
ofsHot.close();  

/*----------------------------------finish MCMC---------------------------------------
-------------------------------------------------------------------------------------*/
end = clock();
double elapsed = double(end - start);
double secs = elapsed/CLOCKS_PER_SEC;
cout<<"-----mcmc-----"<<endl;
cout<<"running time: "<<secs<<endl;

//---output the posterior distribution of haplotypes to file----
if(datatype == 'g' || datatype == 'G'){
string hapext(".hap");  
string hapfile = inputFileName+hapext;  
ofstream ofsHap;
     
ofsHap.open(hapfile.c_str());
if (!ofsHap) {cout<<hapfile<<" can't be opened."<<endl; exit(1);}
  
chains[0]->postHaplotypes(postHap, ofsHap);
ofsHap.close();
}//only for genotypic data



    
//---------------------finish writing files--------------------------------------------
//release memory
for(int i = 0; i<ngraphs; i++){    
for(int j = 0; j<nseqs[i];    j++) delete [] seqmatrix[i][j];
for(int j = 0; j<nmarkers[i]; j++) delete [] SNPs[i][j];
    
delete [] seqmatrix[i];
delete [] SNPs[i];
delete [] phydis[i];
delete [] scaleddis[i];
}
    
for(int ch=0; ch<nchains; ch++){delete chains[ch];}
delete chains;
delete [] nseqs;
delete [] nmarkers;   
delete [] length;
delete [] seqmatrix;
delete [] SNPs;
delete [] phydis;
delete [] scaleddis;

gsl_rng_free(gslr);
return 0;
}
