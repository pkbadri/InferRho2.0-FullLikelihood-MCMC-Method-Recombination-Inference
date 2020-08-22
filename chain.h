/*
Copyright 2009 Ying Wang (ygwang@uchicago.edu)
Created 2007 at University of California, Davis.
*/

#include "graph.h"
#define MAXNHS 10000
#define GEN 0
#define TREE 1
#define MISS 2
#define HAP 3
#define ANC 4
#define RHOBGS 5
#define RHOBGR 6
#define HOTSPC 7
#define HOTSPA 8
#define HOTSPD 9
#define THE 10
#define SWAP 11
#define NC 12


typedef struct {
  double X1;  //start position of the hotspot
  double X2;  //end position of the hotspot
  double Z;   //strength of the hotspot
}Hotspot;




typedef struct{
  int graph;
  int chrid;
  int locus;
}missSnp;




typedef struct{
  int* haplotype;
  int count;
}haplotypeCount;




typedef struct{
  vector<haplotypeCount*> haps;
  int chrid;
  int* hap0;
}chrHaps;




class Chain
{
  private:
  double beta; //relate to the temperature of the chain
  int nGraphs;
  int* nSeq;   //pointers point to the input data
  int* nMarkers;
  double** phyDis;
  double** scaledDis;
  int*** seqMatrix;
  double* length;
  double startInterv;
  double endInterv;
  double* statFreq;
  vector<missSnp>* missing;
  int*** SNPs;
  vector<Graph*>* genealogy;//genealogy consists of an array of graphs (ARGs)
  
  

  //The following are the variables and parameters in the recombination rate model
  double** rho;//matrix of rho for ngraphs*(nmarkers-1)
  double** rhoBackg;
  double rhoBackgScale, gamscale;//scale parameter in the Gamma distribution;
  vector<Hotspot*>* hotspots;
  double theta;//theta is constant across the interval
  double **Gamma;
  
  //Fixed parameters
  double lambdaHalfWidth;
  double lambdaSimZ;
  double deltaNode;
  double hotspMuZ;
  double hotspSigmaZ;
  double hotspLambdaX1;
  double hotspLambdaX2;
  double rhoBackgShape, gamshape;  
  double rhoBackgScaleLambda, gamscalelambda;
  double thetaShape;
  double thetaScale;
  
  //Mixing and adjustment parameters in the MCMC  
  double deltaTheta;
  double deltaRhoBackg;
  double deltaRhoBackgScale, deltagamscale;
  double deltaGamma;	
    
  //Prob of proposing each change
  double prGen, prTree, prMiss, prHap, prAnc, prRhoBackg, prHotsp, prTheta, prGamma;

  int* numProp;
  int* numAcc;
  gsl_rng* gslr;

  public:
  double* MEANTRACT;
  double* priorlik;
  double* treelik;
  double f1, f2;
  int chain_number;    
  Chain(double, int, int*, int*, double**, double**, int***, double*, double, double, double*, vector<missSnp>*, int***, gsl_rng*, int);
  ~Chain();

  vector<Graph*>* getGenealogy() {return genealogy;}
  double** getRho() {return rho;}
  double** getRhoBackg() {return rhoBackg;}
  vector<Hotspot*>* getHotspots() {return hotspots;}
  double getTheta() {return theta;}
  double getbeta() {return beta;}
  void setbeta(double bt){beta = bt;}
  double getRhoBackgScale() {return rhoBackgScale;}
  double getRhoBackgShape() {return rhoBackgShape;}

  void setupMixingPar(double, double, double, double);
  void copyChain(Chain&, Chain&);
  void simulateGenealogy();
  void initializeChain();
  void updateRho(double*, double*, vector<Hotspot*>*, double*, int, double*);
  void printRho();
  void readMcmcFile(ifstream&);
  void writeMcmcFile(ofstream&);

  void modifyGenealogy();
  void modifyTree();
  void modifytractlength();
  void modifyMissingData();
  void modifyHaplotype();
  void modifyAncestralStates();
  void modifyRhoBackg();
  void modifyHotspot();
  void modifyTheta();
  
  void modifyChain(int num);

  void changeRhoBackgScale();
  void changegamscale();
  void changeRhoBackgRate();
  void changeRho();	
  void changeGamma();	
  void changeMEANTRACT();
  void changeHotspAreaConst();  
  void changeHotspot(double);
  void changeHotspotGlobal1(double);
  void changeHotspotGlobal2(double);
  void insertHotspot(double);
  void deleteHotspot(double);
  void getProposalProb(double*);
  void changef1();
  void changef2();
  void swapChains(Chain*);
  void exchangeVariablesChains(Chain*, Chain*);

  void copyHotspots(vector<Hotspot*>*, vector<Hotspot*>*);
  void clearHotspots(vector<Hotspot*>*);
  void printHotspots();
  void printHotspots(ofstream&, int);
  void printHotspots(vector<Hotspot*>*);
  double logPriorHotspots(vector<Hotspot*>*);
  void clearGenealogy(vector<Graph*>*);
  void copyGenealogy(vector<Graph*>*, vector<Graph*>*);
  double logPriorGenealogy(vector<Graph*>*, double**);
  void copyRate(double**, double**);
  double logllhGenealogy(vector<Graph*>*, double);
  void recomGenealogy(vector<double>*);
  //void treeLengths(double**);
  double proposeHotspots(vector<Hotspot*>*, vector<Hotspot*>*);

  void setupProb(double, double, double, double, double, double, double, double);
  void outputAccProb();
  double timeMRCA();
  double logPriorRhoBackg(double**, double);
  double logPriorgam(double**, double);

  //void printGenealogy();
  //void printGenealogy(ofstream&);
  //void sampleHaplotypes(chrHaps**);
  bool hapsAreSame(int*, int*, int);
  void printRec(ofstream&, int);
  void printRho(ofstream&, int);
  void printRhoBackg(ofstream&, int);
  void postHaplotypes(chrHaps**, ofstream&);
  void monitor(Chain**, int, int, ofstream&);
};
