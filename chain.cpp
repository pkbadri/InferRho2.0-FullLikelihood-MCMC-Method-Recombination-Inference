#include "chain.h"

//Declare a new chain, allocate memories and initialize the chain
Chain::Chain(double bt, int ng, int* ns, int* nm, double** pd, double** sd, int*** seqm, double* len, double start, double end, double* statf, vector<missSnp>* miss, int*** snps, gsl_rng* gr, int nch):
beta(bt), nGraphs(ng), nSeq(ns), nMarkers(nm), phyDis(pd), scaledDis(sd), seqMatrix(seqm), length(len), startInterv(start), endInterv(end), statFreq(statf), missing(miss), SNPs(snps), gslr(gr), chain_number(nch)
{
  //Allocate Memory
  genealogy =  new   vector<Graph*>;
  rho       =  new double*[nGraphs]; 
  rhoBackg  =  new double*[nGraphs]; 
  Gamma     =  new double*[nGraphs];
  MEANTRACT =  new  double[nGraphs];
  priorlik  =  new  double[nGraphs];
  treelik   =  new  double[nGraphs];
  
  for(int tau = 0; tau < nGraphs; tau++)
  {
  rho[tau]      = new double[nMarkers[tau]-1]; 
  rhoBackg[tau] = new double[nMarkers[tau]-1];
  Gamma[tau]    = new double[nMarkers[tau]-1];
  }

  hotspots = new vector<Hotspot*>;
  lambdaHalfWidth = 1/500.0;
  lambdaSimZ = 1/5000.0;
  deltaNode = 0.4;
  hotspMuZ = 9;
  hotspSigmaZ = 1;
  hotspLambdaX1 = 1/50000.0;
  hotspLambdaX2 =  1/1000.0;
  rhoBackgShape = 100; 
  rhoBackgScaleLambda = 0.01;
  thetaShape = 0.25; 
  thetaScale = 2;
    
  numProp = new int[NC];
  numAcc  = new int[NC]; 	
  initializeChain();    
}




Chain::~Chain()
{
  for(int tau = 0; tau < nGraphs; tau++){delete genealogy->at(tau);}
  genealogy->clear();
  delete genealogy;

  for(int tau=0; tau<nGraphs; tau++)
  {
  delete [] rho[tau];
  delete [] Gamma[tau];
  delete [] rhoBackg[tau];
  }

  delete [] rho;
  delete [] Gamma;
  delete [] rhoBackg;
  delete [] MEANTRACT;

  clearHotspots(hotspots);
  delete hotspots;
  delete [] numProp;
  delete [] numAcc;
} 




void Chain::setupMixingPar(double dt, double drb, double drbs, double dg)
{
  deltaTheta = dt;
  deltaRhoBackg = drb;
  deltaRhoBackgScale = drbs; 
  deltagamscale = drbs;
  deltaGamma = dg;
}




void Chain::simulateGenealogy()
{
  for(int tau = 0; tau<nGraphs; tau++){
  Graph* agraph;
  agraph = new Graph(nSeq[tau], nMarkers[tau], gslr);
  agraph->simuBNTree(scaledDis[tau], seqMatrix[tau], length[tau], MEANTRACT[tau], Gamma[tau], rho[tau], chain_number); //Simulates a binary tree.      
  priorlik[tau] = agraph->getlogprior(Gamma[tau], rho[tau], scaledDis[tau], length[tau]); //Calculates likelihood of a graph.
  treelik[tau]  = agraph->logllh(theta, -1, statFreq); //Calculates likelihood of data given a graph.
  genealogy->push_back(agraph);
  }
}




void Chain::initializeChain()
{
  for(int i = 0; i < NC; i++) {numProp[i] = 0; numAcc[i] = 0;}

  f1 = 1.0; f2 = 1.0;
  for(int tau = 0; tau < nGraphs;       tau++){
  for(int i   = 0; i   < nMarkers[tau]-1; i++){
  rhoBackg[tau][i]  = gsl_ran_flat(gslr, EPSILON, 200);
  }
  MEANTRACT[tau]    = 125.0;
  }

  theta = gsl_ran_flat(gslr, EPSILON, 0.1); 
  rhoBackgScale = gsl_ran_flat(gslr, EPSILON, 5); 

  for(int tau = 0; tau < nGraphs; tau++){updateRho(rho[tau], rhoBackg[tau], hotspots, phyDis[tau], nMarkers[tau], Gamma[tau]);}
  simulateGenealogy(); 
}




void Chain::updateRho(double* rhov, double* rhobgv, vector<Hotspot*>* hsp, double* phyDis_, int nMarkers_, double* gamma)
{
  for(int i = 0; i < nMarkers_ - 1; i++){rhov[i] = rhobgv[i];  gamma[i] = f1*rhobgv[i];}
  
  int nhs = (int)hsp->size();
  for(int i = 0; i < nhs; i++){
  Hotspot* curhs = hsp->at(i);
  
  if(curhs->X1 > phyDis_[nMarkers_ - 1]) break;

  if(curhs->X1 >= phyDis_[0]){        
  int index = -1;
  for(int j = 0; j < nMarkers_ - 1; j++)
  {
  if(phyDis_[j] <= curhs->X1 && phyDis_[j+1] > curhs->X1){index = j; break;}
  }
          
  assert(index >=0 && phyDis_[index + 1] >= curhs->X1);
          
  double hsStart = curhs->X1; double hsEnd = curhs->X2;
  for(int j = index; j < nMarkers_ - 1; j++)
  {
  double l1; double l2;
  if(phyDis_[j]   <   hsStart) l1 = hsStart;     else l1 =   phyDis_[j];
  if(phyDis_[j+1] >   hsEnd  )   l2 = hsEnd;     else l2 = phyDis_[j+1];

  //Overall rate is a weighted sum of rates in hotspots and non-hotspots for both crossing-over and gene-conversion
  rhov[j]   = ((l2 - l1)*(rhobgv[j]  +  curhs->Z)      +  (phyDis_[j+1] - phyDis_[j]  - (l2 - l1))*rhov[j] )/ (phyDis_[j+1]  -  phyDis_[j]);  
  gamma[j]  = ((l2 - l1)*(rhobgv[j]  +  curhs->Z)*f2   +  (phyDis_[j+1] - phyDis_[j]  - (l2 - l1))*gamma[j])/ (phyDis_[j+1]  -  phyDis_[j]);
  if(phyDis_[j+1] > hsEnd) break;
  }
  }

  }


  for(int i = 0; i < nMarkers_ - 1; i++){assert(rhov[i] > 0  &&  gamma[i] > 0  &&  nhs >= 0);}
}




void Chain::printRho()
{
  for(int tau=0; tau<nGraphs; tau++){
  for(int i=0; i<nMarkers[tau]-1; i++){
  cout<<rho[tau][i]<<" ";
  }cout<<endl;
  }cout<<endl;
}




void Chain::modifyGenealogy()
{
   int chg   = 0;
   Graph* G  = genealogy->at(chg);
   double pr = gsl_rng_uniform_pos(gslr);

   //Proposes new waiting times between events
   if(pr  < 0.20){
   double logpriorGp = G->getlogprior(Gamma[chg], rho[chg], scaledDis[chg], -length[chg]); //Changes waiting times and calculates prior likelihood for new graph.
   double llhGp      = G->logllh(theta, -1, statFreq);
   double alpha      = exp(beta*(llhGp - treelik[chg]));
  
   numProp[GEN]++;
   if(gsl_rng_uniform_pos(gslr)  <  alpha){
   G->settotal(G->getTottimeGraph());
   priorlik[chg] = logpriorGp; 
   treelik[chg] = llhGp;
   numAcc[GEN]++;
   }

   else{G->revertwait();} //Revert waiting times if new proposed graph is not accepted.
   }
   

   else{
   Graph* Gp; 
   Gp = new Graph(*G);
   
   double praddco, praddgc, prdelco, prdelgc, prChange;
   
   if(G->getnumCO() == 0) {praddco = 0.32; prdelco = 0;}
   else{
   if(G->getnumCO() == MAXNRECS) {praddco = 0; prdelco = 0.32;} 
   else{praddco = 0.16; prdelco = 0.16;}   
   }


   if(G->getnumGC() == 0) {praddgc = 0.32; prdelgc = 0;}
   else{
   if(G->getnumGC() == MAXNRECS) {praddgc = 0; prdelgc = 0.32;}
   else{praddgc = 0.16; prdelgc = 0.16;}
   }

   
   prChange = 1 - praddco - prdelco - praddgc - prdelgc;
   
   
   //Change the topology of the tree
   double propRatio;
   if(pr < prChange) 
   {
   propRatio = Gp->moveNode(deltaNode, length[chg], MEANTRACT[chg], Gamma[chg], rho[chg], chain_number);
   }
  
   //Add or delete a pair of recombination and coalescence nodes
   else{ 
   if(pr < prChange + praddco){
   propRatio = Gp->addRecCoal(length[chg],  MEANTRACT[chg], Gamma[chg], rho[chg], chain_number);
   if(G->getnumCO()  ==  MAXNRECS-1) propRatio *=2;
   if(G->getnumCO()  ==  0) propRatio *=0.5;
   }
   
   else{
   if(pr < prChange + praddco + praddgc){
   propRatio = Gp->addconvCoal(length[chg],  MEANTRACT[chg], Gamma[chg], rho[chg], chain_number);
   if(G->getnumGC()  ==  MAXNRECS-1) propRatio *=2;
   if(G->getnumGC()  ==  0) propRatio *=0.5;	   
   }  
   
   else{ 
   if(pr < prChange + praddco + praddgc + prdelco){
   propRatio = Gp->deleteRecCoal(length[chg], MEANTRACT[chg], Gamma[chg], rho[chg], chain_number);
   if(G->getnumCO()  ==  MAXNRECS) propRatio *= 0.5;
   if(G->getnumCO()  ==  1) propRatio *= 2;
   }
  
   else{
   propRatio = Gp->deleteconvCoal(length[chg], MEANTRACT[chg], Gamma[chg], rho[chg], chain_number);
   if(G->getnumGC()  ==  MAXNRECS) propRatio *= 0.5;
   if(G->getnumGC()  ==  1) propRatio *= 2;
   }
   } 

   }
   }
   
   numProp[GEN]++;
   
   //Calculate the probability of acceptance
   double alpha, llhGp, logpriorGp;  
   if(propRatio == 0){alpha = 0;}
  
   else{ 
   logpriorGp = Gp->getlogprior(Gamma[chg], rho[chg], scaledDis[chg], length[chg]);
   if(logpriorGp == 0){alpha = 0;}
      
   else{
   llhGp = Gp->logllh(theta, -1, statFreq);
   alpha = propRatio*exp((logpriorGp - priorlik[chg])  +  beta*(llhGp - treelik[chg]));
   }
   }
   
   
   //Accept the new graph based on the value of alpha
   if(gsl_rng_uniform_pos(gslr) < alpha)
   {     
   priorlik[chg] = logpriorGp;
   treelik[chg]  = llhGp;
   genealogy->pop_back(); delete G; genealogy->push_back(Gp);
   numAcc[GEN]++;
   }
     
   else{
   //Delete temporary graph variable
   delete Gp;
   }
   }
}




void Chain::modifyTree()
{   
  int chg = 0;
  Graph* G = genealogy->at(chg);
  if(G->getnumRecs() == 0) return;
  
  Graph* Gp = new Graph(*G);
  
  double alpha, llhGp, logpriorGp;
  double propRatio = Gp->changeBkpointFollow(length[chg], MEANTRACT[chg], Gamma[chg], rho[chg], chain_number);
  
  if(propRatio == 0){alpha = 0;}
  else{
  llhGp      = Gp->logllh(theta, -1, statFreq);
  logpriorGp = Gp->getlogprior(Gamma[chg], rho[chg], scaledDis[chg], length[chg]); 

  double logpriorRatio = logpriorGp  - priorlik[chg];
  double logllhRatio   = beta*(llhGp - treelik[chg]);
  alpha = propRatio*exp(logpriorRatio + logllhRatio);
  }


  numProp[TREE]++;
  if(gsl_rng_uniform_pos(gslr)  <  alpha){
  priorlik[chg] = logpriorGp;
  treelik[chg]  = llhGp;
  genealogy->pop_back(); delete G; genealogy->push_back(Gp);  
  numAcc[TREE]++;
  }

  else{
  //Delete temporary graph variable
  delete Gp;
  }
}




void Chain::modifytractlength()
{
  /*
  int chg = 0;
  Graph* G = genealogy->at(chg);
  if(G->getnumGC() == 0) return;

  Graph* Gp;
  try {Gp = new Graph(*G);}
  catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}

  double alpha, llhGp, logpriorGp;  
  double propRatio = Gp->changetractlength(length[chg], MEANTRACT[chg], Gamma[chg], rho[chg]);

  if(propRatio == 0){alpha = 0;}
  else{
  llhGp      = Gp->logllh(theta, -1, statFreq);
  logpriorGp = Gp->getlogprior(Gamma[chg], rho[chg], scaledDis[chg], length[chg]);

  double logpriorRatio = logpriorGp - priorlik[chg];
  double logllhRatio = beta*(llhGp - treelik[chg]);
  alpha = propRatio*exp(logpriorRatio + logllhRatio);
  }

  numProp[TREE]++;
  if(gsl_rng_uniform_pos(gslr)  <  alpha){
  priorlik[chg] = logpriorGp; treelik[chg] = llhGp;
  genealogy->pop_back(); delete G; genealogy->push_back(Gp);
  numAcc[TREE]++;
  }

  else{delete Gp;}
  */
}




void Chain::modifyMissingData()
{
  assert(!missing->empty());
  int chosen = gsl_rng_uniform_int(gslr, (int)missing->size());
  int tau = 0; //missing->at(chosen).graph;
  
  Graph* Glocal = genealogy->at(tau);
  Node*    node = Glocal->findNode(missing->at(chosen).chrid);
  int mk = missing->at(chosen).locus;

  double llhG = Glocal->logllh(theta, mk, statFreq);

  ALLELE* vec = node->getAncLocs(LEFT);
  int bk = vec[mk];
  if(bk == SNPs[tau][mk][0]){vec[mk] = SNPs[tau][mk][1];}
  else{vec[mk] = SNPs[tau][mk][0];}

  assert(SNPs[tau][mk][0] != SNPs[tau][mk][1]);

  Glocal->updatemavec();
  double llhGp = Glocal->logllh(theta, mk, statFreq);

  double alpha;
  if(llhGp == 0){alpha = 0;}
  else {alpha = exp(beta*(llhGp - llhG));}

  numProp[MISS]++;
  if(gsl_rng_uniform_pos(gslr) < alpha) {treelik[0] = Glocal->logllh(theta, -1, statFreq); numAcc[MISS]++;}
  else {vec[mk] = bk; Glocal->updatemavec();}
}




void Chain::modifyHaplotype()
{
  int chg = 0;
  Graph* G = genealogy->at(chg);

  int chosen = gsl_rng_uniform_int(gslr, (int)(nSeq[chg]*0.5));
  Node* node1 = G->findNode(chosen*2);
  Node* node2 = G->findNode(chosen*2+1);
  assert(node1->getID() == chosen*2 && node2->getID() == chosen*2+1);

  ALLELE* vec1 = node1->getAncLocs(LEFT);
  ALLELE* vec2 = node2->getAncLocs(LEFT);

  int* hap1 = new int[nMarkers[chg]];
  int* hap2 = new int[nMarkers[chg]];
  for(int mk=0; mk<nMarkers[chg]; mk++)
  {
  hap1[mk] = vec1[mk];
  hap2[mk] = vec2[mk];
  }

  double llhG = treelik[chg];
  for(int mk=0; mk<nMarkers[chg]; mk++)
  {
  if(gsl_rng_uniform_pos(gslr) < 0.1)
  {
  int bk = vec1[mk]; vec1[mk] = vec2[mk]; vec2[mk] = bk;
  }
  }

  //Update MA vectors after proposing new haplotypes
  G->updatemavec();
  double llhGp = G->logllh(theta, -1, statFreq);

  double alpha;
  if(llhGp == 0){alpha = 0;}
  else {alpha = exp(beta*(llhGp  -  llhG));}

  numProp[HAP]++;
  if(gsl_rng_uniform_pos(gslr) < alpha) {treelik[chg] = llhGp; numAcc[HAP]++;}
  else{ 
  for(int mk=0; mk<nMarkers[chg]; mk++)
  {
  vec1[mk] = hap1[mk];
  vec2[mk] = hap2[mk];
  } 

  //Update MA vectors after reverting
  G->updatemavec();
  }

  delete [] hap1;
  delete [] hap2;
}




/*
void Chain::modifyAncestralStates()
{
  int chg = gsl_rng_uniform_int(gslr, nGraphs);
  Graph* G = genealogy->at(chg);

  Graph* Gp;
  try {Gp  = new Graph(*G);}
  catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}

  Gp->changeOneAncSeq();

  double llhGp = Gp->logllh(theta, -1, statFreq);
  double alpha = exp(beta*(llhGp  -  treelik[chg]));

  numProp[ANC]++;
  if(gsl_rng_uniform_pos(gslr)<alpha)
  {
  treelik[chg] = llhGp;
  genealogy->pop_back();
  delete G;  
  genealogy->push_back(Gp);  
  numAcc[ANC]++;
  }
  
  else{
  //Delete proposed graph
  delete Gp;
  }

}
*/




void Chain::changeRhoBackgScale()
{
  //Propose a new value of scale parameter
  double rhoBackgScalePrime = fabs(rhoBackgScale + gsl_ran_flat(gslr, -deltaRhoBackgScale, deltaRhoBackgScale));

  //Calculate the priors based on current background rates
  double logpriorPrime = logPriorRhoBackg(rhoBackg, rhoBackgScalePrime);
  double logprior      = logPriorRhoBackg(rhoBackg, rhoBackgScale     );

  //Calculate acceptance probability
  double alpha = exp((logpriorPrime  -  logprior)  +  (rhoBackgScaleLambda*(rhoBackgScale  -  rhoBackgScalePrime)));

  //Update proposal count
  numProp[RHOBGS]++;
  
  //If acceptance
  if(gsl_rng_uniform_pos(gslr) < alpha){numAcc[RHOBGS]++; rhoBackgScale = rhoBackgScalePrime;}
}




/*
void Chain::changegamscale()
{
  
  double gamscaleprime = fabs(gamscale + gsl_ran_flat(gslr, -deltagamscale, deltagamscale));
  double logpriorPrime = logPriorgam(Gamma, gamscaleprime);
  double logprior      = logPriorgam(Gamma, gamscale);

  double alpha = exp((logpriorPrime - logprior) + (gamscalelambda*(gamscale - gamscaleprime))  );

  numProp[RHOBGS]++;
  if(gsl_rng_uniform_pos(gslr) < alpha){numAcc[RHOBGS]++; gamscale = gamscaleprime;}
  
}
*/




void Chain::changeRhoBackgRate()
{
  int chg  = 0;
  Graph* G = genealogy->at(chg);

  //Create an array of proposed rho background rates
  double* rhoBackgPrime; rhoBackgPrime = new double[nMarkers[chg]-1]; 
  for(int i=0; i<nMarkers[chg]-1; i++){  
  if(gsl_rng_uniform_pos(gslr) < 0.2){rhoBackgPrime[i] = fabs(rhoBackg[chg][i]  +  gsl_ran_flat(gslr, -deltaRhoBackg, deltaRhoBackg));}
  else{rhoBackgPrime[i] = rhoBackg[chg][i];}
  }

  //Calculate the prior probabilities for background rho
  double logpriorPrime = 0, logprior = 0;
  for(int i=0; i<nMarkers[chg]-1; i++){
  logpriorPrime +=  Util::logGamma(rhoBackgPrime[i], rhoBackgShape, rhoBackgScale);
  logprior      +=  Util::logGamma(rhoBackg[chg][i], rhoBackgShape, rhoBackgScale);
  }

  //Get proposed values for the parameters rho and gamma based on rhobackgprime
  double* rhoPrime;      rhoPrime = new double[nMarkers[chg]-1];
  double* gammaPrime;  gammaPrime = new double[nMarkers[chg]-1];
  updateRho(rhoPrime, rhoBackgPrime, hotspots, phyDis[chg], nMarkers[chg], gammaPrime);

  //Change rate variables in the graph nodes for the proposals  
  G->setnewprob(length[chg], MEANTRACT[chg], gammaPrime, rhoPrime, chain_number);
  
  //Calculate prior likelihood of graph for the newly proposed rates
  double newgraphprior = G->getlogprior(gammaPrime, rhoPrime, scaledDis[chg], length[chg]);
  
  //Calculate the acceptance probability
  double alpha = exp((logpriorPrime + newgraphprior)   -    (logprior + priorlik[chg]));

  numProp[RHOBGR]++;
  if(gsl_rng_uniform_pos(gslr) < alpha){
  //Update prior likelihood based on new rates
  priorlik[chg] = newgraphprior; 
  
  //Update rhobackg, rho and gamma arrays
  for(int i=0; i<nMarkers[chg]-1; i++){
  rhoBackg[chg][i] =  rhoBackgPrime[i];
  rho[chg][i]      =       rhoPrime[i];
  Gamma[chg][i]    =     gammaPrime[i];
  }

  numAcc[RHOBGR]++;
  }

  else{
  //Revert rate variables to older values in the graph nodes
  G->revertprob();
  }

  //Free memory
  delete [] rhoBackgPrime; delete [] rhoPrime; delete [] gammaPrime;
}




//Propose new values for gene-conversion rates
void Chain::changeGamma(){
int chg = 0;
Graph* G = genealogy->at(chg);

double* gamprime;gamprime = new double[nMarkers[chg]-1];
  
for(int i=0; i<nMarkers[chg]-1; i++){
if(gsl_rng_uniform_pos(gslr) < 0.333333) gamprime[i] = fabs(Gamma[chg][i]  +  gsl_ran_flat(gslr, -deltaGamma, deltaGamma)); 
else gamprime[i] = Gamma[chg][i];
}
   
G->setnewprob1(length[chg], MEANTRACT[chg], gamprime, rho[chg], chain_number);
  
//Calculate the prior densities for conversion rates
double logpriorPrime = 0, logprior = 0;

//Calculate the prior likelihood of the proposed ARG
logpriorPrime += G->getlogprior(gamprime, rho[chg], scaledDis[chg], length[chg]);

//Calculate acceptance probability
double alpha = exp(logpriorPrime  -  priorlik[chg]);
    
numProp[RHOBGR]++;	
if(gsl_rng_uniform_pos(gslr) < alpha){    
for(int i=0; i<nMarkers[chg]-1; i++){Gamma[chg][i] = gamprime[i];}  
priorlik[chg] = logpriorPrime; 
numAcc[RHOBGR]++;	
}
else{G->revertprob1();}

delete [] gamprime;
}




//Propose new values for Rho and Gamma
void Chain::changeRho(){
int chg = 0;
Graph* G = genealogy->at(chg);

double* rhoPrime;rhoPrime = new double[nMarkers[chg]-1];
for(int i=0; i<nMarkers[chg]-1; i++){
if(gsl_rng_uniform_pos(gslr) < 0.333333) rhoPrime[i] = fabs(rho[chg][i]+ gsl_ran_flat(gslr, -deltaRhoBackg, deltaRhoBackg));
else rhoPrime[i] = rho[chg][i];
}
 
G->setnewprob2(length[chg], MEANTRACT[chg], Gamma[chg], rhoPrime, chain_number);

double logpriorPrime = G->getlogprior(Gamma[chg], rhoPrime, scaledDis[chg], length[chg]);	
double alpha = exp(logpriorPrime - priorlik[chg]);

numProp[RHOBGR]++;	
if(gsl_rng_uniform_pos(gslr) < alpha){
for(int i=0; i<nMarkers[chg]-1; i++) {rho[chg][i] = rhoPrime[i];}
priorlik[chg] = logpriorPrime; 
numAcc[RHOBGR]++;
}

else{G->revertprob2();}
	
delete [] rhoPrime;
}




//Propose new value for the ratio of gene-conversion to crossing-over outside hotspots
void Chain::changef1(){
int chg = 0;
Graph* G = genealogy->at(chg);

//Store current f1 to tmp
double tmpf1 = f1; 

//Propose a new value of f1 between 0 and 1000
f1 = fabs(f1 + gsl_ran_flat(gslr, -0.5, 0.5));
f1 = 1000 - fabs(f1 - 1000);

double* gammaPrime;gammaPrime = new double[nMarkers[chg]-1];

//Get new values for gene-conversion rate based on proposed f1
updateRho(rho[chg], rhoBackg[chg], hotspots, phyDis[chg], nMarkers[chg], gammaPrime);

//Change the rate variables in the graph nodes and calculate prior for proposed rates
G->setnewprob1(length[chg], MEANTRACT[chg], gammaPrime, rho[chg], chain_number);
double logpriorPrime = G->getlogprior(gammaPrime, rho[chg], scaledDis[chg], length[chg]);

//Calculate acceptance probability
double alpha = exp(logpriorPrime  -  priorlik[chg]);

//If new f1 is accepted
if(gsl_rng_uniform_pos(gslr) < alpha){
//Update the prior likelihood
priorlik[chg] = logpriorPrime;

//Update gamma values
for(int i=0; i<nMarkers[chg]-1; i++) {Gamma[chg][i] = gammaPrime[i];}
}

else{
//Revert the value of f1 and rate variables in the graph
f1 = tmpf1; G->revertprob1();
}

//Free memory
delete [] gammaPrime;
}




//Propose new value for ratio of gene-conversion to crossing-over within hotspots
void Chain::changef2(){
int chg = 0;
Graph* G = genealogy->at(chg);

//Store current f2 in tmp
double tmpf2 = f2;

//Propose new f2 between 0 and 1000
f2 = fabs(f2 + gsl_ran_flat(gslr, -0.5, 0.5));
f2 = 1000 - fabs(f2 - 1000);

//Get the new values of gamma based on proposed f2
double* gammaPrime;gammaPrime = new double[nMarkers[chg]-1];
updateRho(rho[chg], rhoBackg[chg], hotspots, phyDis[chg], nMarkers[chg], gammaPrime);

//Change rate variables in graph and calculate prior likelihood for new proposal
G->setnewprob1(length[chg], MEANTRACT[chg], gammaPrime, rho[chg], chain_number);
double logpriorPrime = G->getlogprior(gammaPrime, rho[chg], scaledDis[chg], length[chg]);

//Calculated acceptance probability
double alpha = exp(logpriorPrime - priorlik[chg]);

//If accepted
if(gsl_rng_uniform_pos(gslr) < alpha){
priorlik[chg] = logpriorPrime;
for(int i=0; i<nMarkers[chg]-1; i++) {Gamma[chg][i] = gammaPrime[i];}
}

else{
//Revert the value of f2 and rate variables
f2 = tmpf2; G->revertprob1();
}

//Free memory
delete [] gammaPrime;
}




//Propose new value for mean conversion tract length
void Chain::changeMEANTRACT()
{ 
int chg = 0;
Graph* G = genealogy->at(chg);
double newmean = 1.0 + fabs(1.0 - MEANTRACT[chg] - gsl_ran_flat(gslr, -10, 10));
newmean = 1000.0 - fabs(newmean - 1000.0);
assert(newmean >= 1.0 && newmean <= 1000.0);

int flag[1]; flag[0] = 0;
G->setnewprob111(length[chg], newmean, Gamma[chg], chain_number, flag); //Updates rates and lookup table for new mean tract length 

double alpha, logpriorPrime;
if(flag[0] == 1){alpha = -1;}
else{
logpriorPrime = G->getlogprior(Gamma[chg], rho[chg], scaledDis[chg], length[chg]);
alpha = exp(logpriorPrime  -  priorlik[chg]);
numProp[RHOBGR]++;
}

if(gsl_rng_uniform_pos(gslr) < alpha){
MEANTRACT[chg] = newmean;    
priorlik[chg] = logpriorPrime;
numAcc[RHOBGR]++;
}

else
{
G->updatelkup(MEANTRACT[chg], chain_number); G->revertprob111(); //Reverts rates and lookup table to older values
}
}




void Chain::modifyRhoBackg()
{ 
  int nhs = (int)hotspots->size();
  changeMEANTRACT();
  double U = gsl_rng_uniform_pos(gslr);
  if(U < 0.25){changeRhoBackgRate();}
  else{
  if(U < 0.50){changef1();}
  else{
  if(nhs >= 1.00 && U < 0.75){changef2();}
  else{changeRhoBackgScale();}
  }
  }
}




void Chain::changeHotspAreaConst()
{
  vector<Hotspot*>* hotspotsPrime = new vector<Hotspot*>;
  double propRatio      = proposeHotspots(hotspotsPrime, hotspots);
  double logprior       = logPriorHotspots(hotspots);
  double logpriorPrime  = logPriorHotspots(hotspotsPrime);
  double alpha          = exp(logpriorPrime - logprior)*propRatio;

  numProp[HOTSPC]++;
  if(gsl_rng_uniform_pos(gslr)  <  alpha){
  copyHotspots(hotspots, hotspotsPrime);
  numAcc[HOTSPC]++;
  } 

  clearHotspots(hotspotsPrime);
  delete hotspotsPrime;
}




void Chain::changeHotspot(double prChange)
{
  int i, j;
  double logpriorHotsp = logPriorHotspots(hotspots);
  int nhs = (int)hotspots->size();
  assert(nhs > 0);
  
  //Randomly choose a hotspot to modify
  int chosen = gsl_rng_uniform_int(gslr, nhs);

  double lowb, highb;
  if(chosen == 0)      lowb = startInterv; else  lowb = hotspots->at(chosen-1)->X2;
  if(chosen == nhs-1) highb = endInterv;   else highb = hotspots->at(chosen+1)->X1;

  Hotspot* backup = new Hotspot;
  backup->X1  = hotspots->at(chosen)->X1; 
  backup->X2  = hotspots->at(chosen)->X2; 
  backup->Z   = hotspots->at(chosen)->Z;
  double y1   = backup->X1;
  double y2   = backup->X2;

  if(y1 == startInterv || y2 == endInterv){return;}

  double y1prime, y2prime;
  y1prime  = fabs(y1 + gsl_ran_flat(gslr, -100 , 100));
  if(y1prime <= lowb || y1prime >= y2) {y1prime = y1;}
  y2prime  = fabs(y2 + gsl_ran_flat(gslr, -50  ,  50));
  if(y2prime <= y1prime || y2prime >= highb) {y2prime = y2;}

  hotspots->at(chosen)->X1 = y1prime;
  hotspots->at(chosen)->X2 = y2prime;
  hotspots->at(chosen)->Z  = fabs(hotspots->at(chosen)->Z + gsl_ran_flat(gslr, -20, 20));

  double* probs = new double[3];
  getProposalProb(probs);
  double propRatio = probs[0]/prChange;
  delete [] probs;

  //Arrays for proposed rho and gamma 
  double** rhoPrime; double** gammaPrime;
  rhoPrime   = new double*[nGraphs];
  gammaPrime = new double*[nGraphs];

  double logprior = 0; double logpriorPrime = 0;
  for(i=0; i<nGraphs; i++){
  rhoPrime[i]   = new double[nMarkers[i]-1];
  gammaPrime[i] = new double[nMarkers[i]-1];

  //Update rho, gamma and the rate variables in the graph nodes based on new hotspots
  updateRho(rhoPrime[i], rhoBackg[i], hotspots, phyDis[i], nMarkers[i], gammaPrime[i]);

  Graph* G = genealogy->at(i);
  G->setnewprob(length[i], MEANTRACT[i], gammaPrime[i], rhoPrime[i], chain_number);//Set new rate variables
  logpriorPrime += G->getlogprior(gammaPrime[i], rhoPrime[i], scaledDis[i], length[i]);//Calculate new prior
  logprior      += priorlik[i];
  }

  double logpriorHotspPrime = logPriorHotspots(hotspots);
  double alpha = exp((logpriorPrime  -  logprior) + (logpriorHotspPrime  -  logpriorHotsp))*propRatio;

  numProp[HOTSPC]++;
  if(gsl_rng_uniform_pos(gslr)  <  alpha){  
  for(i=0; i<nGraphs; i++){
  priorlik[i] = logpriorPrime;                     
  for(j=0; j<nMarkers[i]-1; j++){rho[i][j] = rhoPrime[i][j]; Gamma[i][j] = gammaPrime[i][j];}                  
  }

  numAcc[HOTSPC]++;
  } 


  else
  {
  //Change back the modified hotspot
  hotspots->at(chosen)->X1 = backup->X1; 
  hotspots->at(chosen)->X2 = backup->X2; 
  hotspots->at(chosen)->Z  = backup->Z;                   

  //Revert the rate variables in the graph nodes
  for(i = 0; i<nGraphs; i++){genealogy->at(i)->revertprob();}
  }


  //Free memory
  delete backup;
  for(i = 0; i<nGraphs; i++){delete [] rhoPrime[i]; delete [] gammaPrime[i];}  
  delete [] rhoPrime; delete [] gammaPrime;
}




void Chain::changeHotspotGlobal1(double prChange)
{
  int i, j;
  double logpriorHotsp = logPriorHotspots(hotspots);
  int nhs = (int)hotspots->size();
  assert(nhs  >  0);

  //Choose a hotspot at random to modify
  int chosen = gsl_rng_uniform_int(gslr, nhs);

  double lowb, highb;
  if(chosen == 0)      lowb = startInterv; else  lowb = hotspots->at(chosen  -  1)->X2;
  if(chosen == nhs-1) highb = endInterv;   else highb = hotspots->at(chosen  +  1)->X1;

  Hotspot* backup = new Hotspot;
  backup->X1 = hotspots->at(chosen)->X1;
  backup->X2 = hotspots->at(chosen)->X2;
  backup->Z  = hotspots->at(chosen)->Z;

  if(backup->X1 == startInterv){return;}

  double logprForw = 0; double logprBackw = 0;
  double y1, y2, B;
  y1 = gsl_rng_uniform_pos(gslr)*(highb - lowb) + lowb;
  if(highb - y1 < y1 - lowb) B = highb - y1; else B = y1 - lowb;
  
  if(highb == endInterv && highb - y1  <  y1 - lowb){ 
  double prExceed = exp(-lambdaHalfWidth*(highb  -  y1));
  
  if(gsl_rng_uniform_pos(gslr) < prExceed) {y2 = highb - y1; logprForw += log(prExceed);}
  else{
  y2 = -1.0/lambdaHalfWidth * log(1   -   gsl_rng_uniform_pos(gslr)*(1  -  exp(-lambdaHalfWidth*B)));
  logprForw += log(1-prExceed) + log(lambdaHalfWidth) - lambdaHalfWidth*y2 - log(1-exp(-lambdaHalfWidth*B));
  }  
  }

  else{
  y2  = -1.0/lambdaHalfWidth * log(1   -   gsl_rng_uniform_pos(gslr)*(1  -  exp(-lambdaHalfWidth*B)));
  logprForw += log(lambdaHalfWidth) - lambdaHalfWidth*y2 -  log(1  -  exp(-lambdaHalfWidth*B));
  }
  
  
  hotspots->at(chosen)->X1 = y1-y2;
  hotspots->at(chosen)->X2 = y1+y2;
  hotspots->at(chosen)->Z = fabs(hotspots->at(chosen)->Z + gsl_ran_flat(gslr, -100, 100));


  y1 = 0.5*(backup->X1+backup->X2);
  y2 = 0.5*(backup->X2-backup->X1);
  if(highb-y1 < y1-lowb) B = highb-y1; else B = y1-lowb;
  if(chosen == nhs-1 && highb-y1 < y1-lowb)
  {
  double prExceed = exp(-lambdaHalfWidth*(highb-y1));
  if(fabs(y2-(highb-y1))  <  1e-8) logprBackw += log(prExceed);
  else logprBackw += log(1-prExceed) + log(lambdaHalfWidth) - lambdaHalfWidth*y2 - log(1-exp(-lambdaHalfWidth*B));
  }  
  else logprBackw += log(lambdaHalfWidth) - lambdaHalfWidth*y2 - log(1-exp(-lambdaHalfWidth*B));

  double propRatio = exp(logprBackw  -  logprForw);

  double* probs = new double[3];
  getProposalProb(probs);
  propRatio *= probs[0]/prChange;
  delete [] probs;

  //Proposed rho and gamma
  double** rhoPrime; double** gammaPrime;
  rhoPrime   = new double*[nGraphs];
  gammaPrime = new double*[nGraphs];
  
 
  double logprior = 0; double logpriorPrime = 0;
  for(i=0; i<nGraphs; i++){
  rhoPrime[i]   = new double[nMarkers[i]-1];
  gammaPrime[i] = new double[nMarkers[i]-1];

  //Update rho, gamma and rate variables in the graph nodes
  updateRho(rhoPrime[i], rhoBackg[i], hotspots, phyDis[i], nMarkers[i], gammaPrime[i]);
  Graph* G = genealogy->at(i);
  G->setnewprob(length[i], MEANTRACT[i], gammaPrime[i], rhoPrime[i], chain_number);  //Change rate variables in graph
  logpriorPrime   += G->getlogprior(gammaPrime[i], rhoPrime[i], scaledDis[i], length[i]);
  logprior        += priorlik[i];
  }
   
  double logpriorHotspPrime = logPriorHotspots(hotspots);
  double alpha = exp((logpriorPrime - logprior)    +    (logpriorHotspPrime - logpriorHotsp))*propRatio;

  numProp[HOTSPC]++;
  if(gsl_rng_uniform_pos(gslr)  <  alpha){
  for(i=0; i<nGraphs; i++){
  //Update prior likelihood based on new rates
  priorlik[i] = logpriorPrime; 
 
  //Update rates
  for(j=0; j<nMarkers[i]-1; j++){rho[i][j] = rhoPrime[i][j]; Gamma[i][j] = gammaPrime[i][j];}
  }
  numAcc[HOTSPC]++;
  }

  else{
  //Revert hotspot to older values
  hotspots->at(chosen)->X1  =  backup->X1;
  hotspots->at(chosen)->X2  =  backup->X2;
  hotspots->at(chosen)->Z   =  backup->Z;
  
  //Revert the rate variables in the graph nodes
  for(i = 0; i<nGraphs; i++){genealogy->at(i)->revertprob();}
  }

   
  //Free memory
  delete backup;
  for(i=0; i<nGraphs; i++){delete [] rhoPrime[i]; delete [] gammaPrime[i];}
  delete [] rhoPrime; delete [] gammaPrime;
}




void Chain::changeHotspotGlobal2(double prChange)
{
  int i, j;
  double logpriorHotsp = logPriorHotspots(hotspots);
  int nhs = (int)hotspots->size();
  assert(nhs  >  0);

  //Choose a hotspot at random to modify
  int chosen = gsl_rng_uniform_int(gslr, nhs);

  double lowb, highb;
  if(chosen == 0)      lowb = startInterv; else  lowb = hotspots->at(chosen  -  1)->X2;
  if(chosen == nhs-1) highb = endInterv;   else highb = hotspots->at(chosen  +  1)->X1;

  Hotspot* backup = new Hotspot;
  backup->X1 = hotspots->at(chosen)->X1;
  backup->X2 = hotspots->at(chosen)->X2;
  backup->Z  = hotspots->at(chosen)->Z;

  if(backup->X2 == endInterv){return;}

  double logprForw = 0; double logprBackw = 0;
  double y1, y2, B;
  y1 = gsl_rng_uniform_pos(gslr)*(highb - lowb) + lowb;
  if(highb - y1  <  y1 - lowb){B = highb - y1;} else{B = y1 - lowb;}
  
  if(lowb == startInterv  &&  B == y1 - lowb){ 
  double prExceed = exp(-lambdaHalfWidth*(y1 - lowb)); //Allow for hotspots to go beyond start
  
  if(gsl_rng_uniform_pos(gslr) < prExceed) {y2 = y1 - lowb; logprForw += log(prExceed);}
  else{
  y2 = -1.0/lambdaHalfWidth * log(1-gsl_rng_uniform_pos(gslr)*(1-exp(-lambdaHalfWidth*B)));
  logprForw += log(1 - prExceed) + log(lambdaHalfWidth) - lambdaHalfWidth*y2 - log(1-exp(-lambdaHalfWidth*B));
  }  
  }

  else{
  y2  = -1.0/lambdaHalfWidth * log(1-gsl_rng_uniform_pos(gslr)*(1-exp(-lambdaHalfWidth*B)));
  logprForw += log(lambdaHalfWidth) - lambdaHalfWidth*y2 -  log(1-exp(-lambdaHalfWidth*B));
  }
    
  hotspots->at(chosen)->X1 = y1  -  y2;
  hotspots->at(chosen)->X2 = y1  +  y2;
  hotspots->at(chosen)->Z = fabs(hotspots->at(chosen)->Z + gsl_ran_flat(gslr, -100, 100));

  y1 = 0.5*(backup->X1+backup->X2);
  y2 = 0.5*(backup->X2-backup->X1);
  if(highb - y1 < y1 - lowb){B = highb - y1;} else{B = y1 - lowb;}
  if(chosen == 0 && B == y1 - lowb)
  {
  double prExceed = exp(-lambdaHalfWidth*(y1 - lowb));
  if(fabs(y2  -  (y1 - lowb)) < 1e-8){logprBackw += log(prExceed);}
  else {logprBackw += log(1-prExceed) + log(lambdaHalfWidth) - lambdaHalfWidth*y2 - log(1-exp(-lambdaHalfWidth*B));}
  }  
  else{logprBackw += log(lambdaHalfWidth) - lambdaHalfWidth*y2 - log(1-exp(-lambdaHalfWidth*B));}

  double propRatio = exp(logprBackw  -  logprForw);

  double* probs = new double[3];
  getProposalProb(probs);
  propRatio *= probs[0]/prChange;
  delete [] probs;

  //Proposed rho and gamma
  double** rhoPrime; double** gammaPrime;
  rhoPrime   = new double*[nGraphs];
  gammaPrime = new double*[nGraphs];
   
  double logprior = 0; double logpriorPrime = 0;
  for(i = 0; i<nGraphs; i++){
  rhoPrime[i]   = new double[nMarkers[i]-1];
  gammaPrime[i] = new double[nMarkers[i]-1];

  //Update rho and gamma 
  updateRho(rhoPrime[i], rhoBackg[i], hotspots, phyDis[i], nMarkers[i], gammaPrime[i]);
  Graph* G = genealogy->at(i);
  G->setnewprob(length[i], MEANTRACT[i], gammaPrime[i], rhoPrime[i], chain_number); //Update rate variables in the graph nodes
  logpriorPrime += G->getlogprior(gammaPrime[i], rhoPrime[i], scaledDis[i], length[i]);
  logprior      += priorlik[i];
  }
   
  double logpriorHotspPrime = logPriorHotspots(hotspots);
  double alpha = exp((logpriorPrime - logprior)    +    (logpriorHotspPrime - logpriorHotsp))*propRatio;

  numProp[HOTSPC]++;
  if(gsl_rng_uniform_pos(gslr)  <  alpha){
  for(i=0; i<nGraphs; i++){
  //Update prior likelihood based on new rates
  priorlik[i] = logpriorPrime; 
 
  //Update recombination rates
  for(j=0; j<nMarkers[i]-1; j++){rho[i][j] = rhoPrime[i][j]; Gamma[i][j] = gammaPrime[i][j];}
  }
  numAcc[HOTSPC]++;
  }

  else{
  //Revert hotspot to older values
  hotspots->at(chosen)->X1  =  backup->X1;
  hotspots->at(chosen)->X2  =  backup->X2;
  hotspots->at(chosen)->Z   =  backup->Z;
  
  //Revert the rate variables in the graph nodes
  for(i = 0; i<nGraphs; i++){genealogy->at(i)->revertprob();}
  }
   
  //Free memory
  delete backup;
  for(i=0; i<nGraphs; i++){delete [] rhoPrime[i]; delete [] gammaPrime[i];}
  delete [] rhoPrime; delete [] gammaPrime;
}




//Inserts a new hotspot
void Chain::insertHotspot(double prInsert)
{
   double logpriorHotsp = logPriorHotspots(hotspots);
   int nhs = (int)hotspots->size();
   int chosen;

   double logprForw = 0, logprBackw = 0;
   if(nhs > 0 && fabs(hotspots->at(nhs-1)->X2 - endInterv) < EPSILON){chosen = gsl_rng_uniform_int(gslr, nhs);}
   else {chosen = gsl_rng_uniform_int(gslr, nhs+1);}
 
   double lowb, highb;
   if(chosen  ==   0) lowb = startInterv;  else lowb  = hotspots->at(chosen-1)->X2;
   if(chosen  == nhs) highb  = endInterv;  else highb =   hotspots->at(chosen)->X1;


   if(chosen == 0 && fabs(highb - startInterv) < 1e-08){return;}

   Hotspot* insertedhs = new Hotspot;

   double y1, y2, B, z;
   y1 = gsl_rng_uniform_pos(gslr)*(highb - lowb) + lowb;
   logprForw += log(1.0/(highb - lowb));

   if(highb - y1 < y1 - lowb) B = highb - y1; else B = y1 - lowb;
   
   y2  = -1.0/lambdaHalfWidth * log(1  -  gsl_rng_uniform_pos(gslr)*(1 - exp(-lambdaHalfWidth*B)));
   logprForw += log(lambdaHalfWidth) - lambdaHalfWidth*y2 - log(1 -  exp(-lambdaHalfWidth*B));

   z = -1.0/lambdaSimZ*log(1   -   gsl_rng_uniform_pos(gslr));
   logprForw += log(lambdaSimZ) - lambdaSimZ*z;

   insertedhs->X1 = y1 - y2;
   insertedhs->X2 = y1 + y2;
   insertedhs->Z = z;

   double propRatio = exp(logprBackw  -  logprForw);
   double jacobian = 2;

   hotspots->insert(hotspots->begin()+chosen, insertedhs);

   double* probs = new double[3];
   getProposalProb(probs);
   propRatio *= probs[2]/prInsert;
   delete [] probs;

   double** rhoPrime;      rhoPrime   = new double*[nGraphs];
   double** gammaPrime;    gammaPrime = new double*[nGraphs];
   
   double logprior = 0;  double logpriorPrime = 0;
   for(int i = 0; i < nGraphs; i++){   
   rhoPrime[i]   = new double[nMarkers[i]-1];
   gammaPrime[i] = new double[nMarkers[i]-1];  

   //Update rho, gamma and rate variables in the graph
   updateRho(rhoPrime[i], rhoBackg[i], hotspots, phyDis[i], nMarkers[i], gammaPrime[i]);
   genealogy->at(i)->setnewprob(length[i], MEANTRACT[i], gammaPrime[i], rhoPrime[i], chain_number);
   logpriorPrime += genealogy->at(i)->getlogprior(gammaPrime[i], rhoPrime[i], scaledDis[i], length[i]);         
   logprior      += priorlik[i];
   }

   double logpriorHotspPrime = logPriorHotspots(hotspots);  
   double alpha = exp((logpriorPrime - logprior)    +    (logpriorHotspPrime  -  logpriorHotsp))*propRatio*jacobian;
   int i, j;

   numProp[HOTSPA]++;
   if(gsl_rng_uniform_pos(gslr) < alpha)   
   {      
   for(i=0; i<nGraphs; i++){         
   //Update prior likelihood variable
   priorlik[i] = logpriorPrime;

   //Update rho and gamma values
   for(j=0; j<nMarkers[i]-1; j++){rho[i][j] = rhoPrime[i][j]; Gamma[i][j] = gammaPrime[i][j];}      
   }   
   numAcc[HOTSPA]++;
   }

   else
   {
   //Delete the proposed hotspot
   Hotspot* hs = *(hotspots->begin()  +  chosen);
   hotspots->erase(hotspots->begin()  +  chosen);
   delete hs;
   
   //Revert the rate variables in the graph nodes
   for(i = 0; i < nGraphs; i++){genealogy->at(i)->revertprob();}
   }

   //Free memory
   for(i = 0; i < nGraphs; i++){delete [] rhoPrime[i]; delete gammaPrime[i];}
   delete [] rhoPrime; delete [] gammaPrime;
}




void Chain::deleteHotspot(double prDelete)
{
  double logpriorHotsp = logPriorHotspots(hotspots);
  int nhs = (int)hotspots->size();
  int chosen;
  Hotspot* removedhs;

  assert(nhs > 0);
  double logprForw = 0, logprBackw = 0;
  if(nhs > 1 && fabs(hotspots->at(nhs-1)->X2 - endInterv) < EPSILON){chosen = gsl_rng_uniform_int(gslr, nhs-1);}
  else{chosen = gsl_rng_uniform_int(gslr, nhs);}


  if(chosen == 0 && fabs(hotspots->at(chosen)->X1 - startInterv) < 1e-08){return;}


  removedhs = hotspots->at(chosen);
  hotspots->erase(hotspots->begin()  +  chosen);
  nhs--;

  double lowb, highb;
  if(chosen ==   0) lowb = startInterv; else lowb = hotspots->at(chosen-1)->X2;
  if(chosen == nhs) highb = endInterv;  else highb = hotspots->at(chosen)->X1;

  double y1, y2, B, z;
  y1 = 0.5*(removedhs->X1+removedhs->X2);
  y2 = 0.5*(removedhs->X2-removedhs->X1);

  logprBackw += log(1.0/(highb-lowb));

  if(highb - y1 < y1 - lowb) B = highb - y1; else B = y1 - lowb;
  logprBackw += log(lambdaHalfWidth) - lambdaHalfWidth*y2 - log(1  -  exp(-lambdaHalfWidth*B));

  z = removedhs->Z;
  logprBackw += log(lambdaSimZ) - lambdaSimZ*z;

  double propRatio = exp(logprBackw   -   logprForw);
  double jacobian  = 0.5;

  double* probs = new double[3];
  getProposalProb(probs);
  propRatio *= probs[1]/prDelete;
  delete [] probs;

  double** rhoPrime; double** gammaPrime;
  rhoPrime   = new double*[nGraphs];
  gammaPrime = new double*[nGraphs];

  double logprior = 0;  double logpriorPrime = 0;
  for(int i=0; i<nGraphs; i++){
  rhoPrime[i]   = new double[nMarkers[i]-1];
  gammaPrime[i] = new double[nMarkers[i]-1];

  //Update rho, gamma and rate variables in the graph nodes
  updateRho(rhoPrime[i], rhoBackg[i], hotspots, phyDis[i], nMarkers[i], gammaPrime[i]);
  genealogy->at(i)->setnewprob(length[i], MEANTRACT[i], gammaPrime[i], rhoPrime[i], chain_number);
    
  logpriorPrime += genealogy->at(i)->getlogprior(gammaPrime[i], rhoPrime[i], scaledDis[i], length[i]);
  logprior      += priorlik[i];
  }

  double logpriorHotspPrime = logPriorHotspots(hotspots);
  double alpha = exp((logpriorPrime  -  logprior)  +  (logpriorHotspPrime   -   logpriorHotsp))*propRatio*jacobian;
  int i, j;

  numProp[HOTSPD]++;
  if(gsl_rng_uniform_pos(gslr) < alpha){
  for(i=0; i<nGraphs; i++){
  //Update prior likelihood variable
  priorlik[i] = logpriorPrime;
  //Update rate variables
  for(j=0; j<nMarkers[i]-1;j++){rho[i][j] = rhoPrime[i][j]; Gamma[i][j] = gammaPrime[i][j];}
  }

  delete removedhs;
  numAcc[HOTSPD]++;
  }

  else{
  //Add back the deleted hotspot
  int nhs = hotspots->size();
  if(nhs == 0 || chosen == nhs) hotspots->push_back(removedhs);
  else hotspots->insert(hotspots->begin()+chosen, removedhs);
  
  //Revert the rate variables in the graph to older values
  for(i=0; i<nGraphs; i++){genealogy->at(i)->revertprob();}
  }
   
  //Free memory
  for(i=0; i<nGraphs; i++){delete [] rhoPrime[i]; delete [] gammaPrime[i];} 
  delete [] rhoPrime; delete [] gammaPrime;
}




void Chain::modifyHotspot()
{
  if(gsl_rng_uniform_pos(gslr) < 0.1){changeHotspAreaConst();}
  else
  {
  double* probs = new double[3];
  getProposalProb(probs);
  double prChange = probs[0], prInsert = probs[1], prDelete = probs[2];

  double U = gsl_rng_uniform_pos(gslr);
  if(U  <  prChange){
  if(gsl_rng_uniform_pos(gslr) < 0.1){
  if(drand48() < 0.5){changeHotspotGlobal1(prChange);} else{changeHotspotGlobal2(prChange);}
  } 
  else{changeHotspot(prChange);}
  }
  else if(U  <  prChange+prInsert){insertHotspot(prInsert);}
  else{deleteHotspot(prDelete);}
  
  delete [] probs;
  }
}




void Chain::getProposalProb(double* probs)
{
   int nhs = (int)hotspots->size();

   if(nhs == 0) {probs[0] = 0; probs[1] = 1; probs[2] = 0;}
   else if(nhs == MAXNHS)
   {
   if(nhs == 1 && fabs(hotspots->at(0)->X2 - endInterv) < EPSILON) 
   {probs[0] = 1;   probs[1] = 0;   probs[2] = 0;}
   else{probs[0] = 0.7; probs[1] = 0;   probs[2] = 0.3;}
   }

   else
   {
   if(nhs == 1 && fabs(hotspots->at(0)->X2 - endInterv) < EPSILON) 
   {probs[0] = 0.2; probs[1] = 0.8; probs[2] = 0;}
   else{probs[0] = 0.5; probs[1] = 0.2; probs[2] = 0.3;}
   }
}




void Chain::modifyTheta()
{
    double dTheta = gsl_ran_flat(gslr, -deltaTheta, deltaTheta);
    double thetaPrime = fabs(theta  +  dTheta);
         
    double logllhPrime = 0; double logllh = 0;
    for(int tau=0; tau<nGraphs; tau++){
    logllhPrime += genealogy->at(tau)->logllh(thetaPrime, -1, statFreq);               
    logllh      += treelik[tau];
    }
   
    double logpriorRatio = Util::logGamma(thetaPrime, thetaShape, thetaScale) - Util::logGamma(theta, thetaShape, thetaScale);
    double alpha = exp(beta*(logllhPrime - logllh)  +  logpriorRatio);

    numProp[THE]++;
    if(gsl_rng_uniform_pos(gslr)  <  alpha){
    theta = thetaPrime;
    treelik[0] = logllhPrime;
    numAcc[THE]++;  
    }
}




void Chain::setupProb(double prgen, double prtree, double prmiss, double prhap, double pranc, double prrhobg, double prhotsp, double prtheta)
{
  prGen = prgen;
  prTree = prtree;
  prMiss = prmiss;
  prHap = prhap;
  prAnc = pranc;
  prRhoBackg = prrhobg;
  prHotsp = prhotsp;
  prTheta = prtheta;
}




void Chain::modifyChain(int num)
{  
  for(int i = 0; i < num; i++){
  double U = gsl_rng_uniform_pos(gslr);
  if(U < prGen) {modifyGenealogy();}
  else if(U < prGen + prTree){modifyTree();}
  else if(U < prGen + prTree + prMiss){modifyMissingData();}
  else if(U < prGen + prTree + prMiss + prHap) {modifyHaplotype();}
  else if(U < prGen + prTree + prMiss + prHap + prAnc + prRhoBackg) {modifyRhoBackg();}
  else if(U < prGen + prTree + prMiss + prHap + prAnc + prRhoBackg + prHotsp) {modifyHotspot();}
  else {modifyTheta();}
  modifyTree();
  modifyHotspot();
  }
}




//Swap chains for MCMCMC according to the temperature
void Chain::swapChains(Chain* chj)
{  
  double logprior    = 0; 
  double logpriorChj = 0; 
  double logpriorRatio = (logprior - logpriorChj) + (logpriorChj - logprior);

  double logllh = 0; double logllhChj = 0;
  for(int i = 0; i < nGraphs; i++){logllh += treelik[i]; logllhChj += chj->treelik[i];}

  double logllhRatio = chj->beta*(logllh  -  logllhChj) + beta*(logllhChj  -  logllh);

  double alpha = exp(logpriorRatio + logllhRatio);

  numProp[SWAP]++;
  if(gsl_rng_uniform_pos(gslr) < alpha){
  exchangeVariablesChains(this, chj);
  numAcc[SWAP]++;
  }
}

  


//Swap beta variables between chains
void Chain::exchangeVariablesChains(Chain* chk, Chain* chj){
  double tmpbeta = chk->getbeta();
  chk->setbeta(chj->getbeta());
  chj->setbeta(tmpbeta);

  for(int i=0; i<NC; i++){
  int tmp;
  tmp = chk->numProp[i];
  chk->numProp[i] = chj->numProp[i];
  chj->numProp[i] = tmp;
  
  tmp = chk->numAcc[i];
  chk->numAcc[i] = chj->numAcc[i];
  chj->numAcc[i] = tmp;
  }

  return;
}




void Chain::copyHotspots(vector<Hotspot*>* newhsvec, vector<Hotspot*>* hsvec)
{
  int nhs = (int)hsvec->size();

  clearHotspots(newhsvec);

  for(int i=0; i<nhs; i++)   
  {
  Hotspot* hs = new Hotspot;      
  hs->X1 = hsvec->at(i)->X1;      
  hs->X2 = hsvec->at(i)->X2;      
  hs->Z = hsvec->at(i)->Z;      
  newhsvec->push_back(hs);   
  }
}




void Chain::clearHotspots(vector<Hotspot*>* hsvec)
{
   vector<Hotspot*>::iterator it;
   for(it=hsvec->begin(); it!=hsvec->end(); it++)
   {
   delete (*it);
   }
   hsvec->clear();
}




void Chain::printHotspots()
{
   cout<<"Hotspots: "<<int(hotspots->size())<<endl;
   vector<Hotspot*>::iterator it;
   for(it=hotspots->begin(); it!=hotspots->end(); it++)
   {cout<<(*it)->X1<<"\t"<<(*it)->X2<<"\t"<<(*it)->Z<<endl;}
   cout<<endl;
}




void Chain::printHotspots(ofstream& ofs, int iter)
{
   ofs<<iter<<"\t"<<int(hotspots->size())<<endl;
   vector<Hotspot*>::iterator it;
   for(it = hotspots->begin(); it != hotspots->end(); it++)
   {
   //assert((*it)->X1 > 0);assert((*it)->X2 > 0);
   //assert((*it)->Z > 0);
   ofs<<(*it)->X1<<"\t"<<(*it)->X2<<"\t"<<(*it)->Z<<endl;
   }
   ofs<<endl;
}




void Chain::printHotspots(vector<Hotspot*>* hsp)
{  
   cout<<hsp->size()<<": ";
   vector<Hotspot*>::iterator it;
   for(it = hsp->begin(); it != hsp->end(); it++)
   cout<<(*it)->X1<<"\t"<<(*it)->X2<<"\t"<<(*it)->Z<<endl;
   cout<<endl;
}     




double Chain::logPriorHotspots(vector<Hotspot*>* hsp)
{
   double logprior  =  0;
   int nhs = (int)hsp->size();

   for(int i = 0; i<nhs; i++){
   double x1Pre;
   if(i == 0){x1Pre = startInterv;} else{x1Pre = hsp->at(i-1)->X2;}
   
   double x1 = hsp->at(i)->X1;
   double x2 = hsp->at(i)->X2;
   if(x1Pre > x1 || x1 > x2){
   cout<<x1Pre<<">"<<x1<<"\tor\t"<<x1<<">"<<x2<<endl;
   exit(1);
   }

   if(i == 0  &&  x1 == startInterv){logprior += log(hotspLambdaX1) - log(hotspLambdaX1 + hotspLambdaX2);}   
   if(i == 0  &&  x1  > startInterv){logprior += log(hotspLambdaX2) - log(hotspLambdaX1 + hotspLambdaX2);}
  
   if(x1 > startInterv){logprior += log(hotspLambdaX1) - hotspLambdaX1*(x1 - x1Pre);}

      
   if(x2 == endInterv){logprior += -hotspLambdaX2*(x2 - x1);}
   else{logprior += log(hotspLambdaX2) - hotspLambdaX2*(x2 - x1);}

   
   logprior += Util::logNormal(hsp->at(i)->Z, hotspMuZ, hotspSigmaZ);
   }

   if(nhs == 0){logprior += log(hotspLambdaX2) - log(hotspLambdaX1 + hotspLambdaX2) - hotspLambdaX1*(endInterv - startInterv);}
   else if(hsp->at(nhs-1)->X2 < endInterv){logprior += -hotspLambdaX1*(endInterv  -  hsp->at(nhs-1)->X2);}

   return logprior;
}




void Chain::clearGenealogy(vector<Graph*>* gen)
{
   vector<Graph*>::iterator it;
   for(it = gen->begin(); it != gen->end(); it++)
   delete (*it);
   gen->clear();
}




void Chain::copyGenealogy(vector<Graph*>* newgen, vector<Graph*>* gen)
{
   vector<Graph*>::iterator it;
   for(it=newgen->begin(); it!=newgen->end(); it++)
   delete (*it);
   newgen->clear();

   for(it=gen->begin(); it!=gen->end(); it++)
   {
   Graph* G = new Graph(*(*it));
   newgen->push_back(G);
   }
}




double Chain::logPriorGenealogy(vector<Graph*>* gen, double** rhov)
{
   double logprior =0;
   vector<Graph*>::iterator it;

   for(int i = 0; i < nGraphs; i++)
   {
   logprior += gen->at(i)->getlogprior(Gamma[i], rhov[i], scaledDis[i], length[i]);
   }

   return logprior;
}




void Chain::copyRate(double** newRho, double** givenRho)
{
   for(int i=0; i<nGraphs; i++)
   {
   for(int j=0; j<nMarkers[i]-1; j++) 
   {
   newRho[i][j]=givenRho[i][j];
   }
   }
}




double Chain::logllhGenealogy(vector<Graph*>* gen, double the)
{
   double logllh =0;
   for(int i=0; i<nGraphs; i++)
   {
   logllh += gen->at(i)->logllh(the, -1, statFreq);
   }
   return logllh;
}




void Chain::recomGenealogy(vector<double>* recom)
{
   for(int tau=0; tau<nGraphs; tau++)
   {
      vector<Node*> nodes = genealogy->at(tau)->getNodesVector();
      int nnodes = (int)nodes.size();

    double startint = phyDis[tau][0];
    double endint = phyDis[tau][nMarkers[tau]-1];
    if(nGraphs==1) assert(startint==startInterv && endint==endInterv); 

      for(int j=0; j<nnodes; j++)
      {
         if(nodes[j]->isRec())
         {
        double loc = (endint-startint)*nodes[j]->getBkpoint()+startint;
            recom->push_back(loc);
         }
      }
   }
   sort(recom->begin(), recom->end(), Util::compareTwoNumbers);
}




/*
void Chain::treeLengths(double** trlen)
{
   for(int i=0; i<nGraphs; i++)
   {
      vector<Node*>** trees;
      try { trees = new vector<Node*>*[nMarkers[i]];}
      catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}
      for(int j=0; j<nMarkers[i]; j++)
      {
         try { trees[j] = new vector<Node*>;}
         catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}
      }

      genealogy->at(i)->getTreesFast(trees);
      for(int j=0; j<nMarkers[i]; j++)
      {
         trlen[i][j] = genealogy->at(i)->totalTreeLengths(trees[j]);
      }

      for(int j=0; j<nMarkers[i]; j++) delete trees[j];
      delete [] trees;
   }
}
*/




double Chain::proposeHotspots(vector<Hotspot*>* hspPrime, vector<Hotspot*>* hsp)
{
   copyHotspots(hspPrime, hsp);

   int nhs = (int)hspPrime->size();
   double propRatio = 1;
   for(int i = 0; i < nhs; i++)
   {
      if(gsl_rng_uniform_pos(gslr) < 0.5) continue;

      double x1 = hspPrime->at(i)->X1;
      double x2 = hspPrime->at(i)->X2;
      double z  =  hspPrime->at(i)->Z;

      for(int tau=0; tau<nGraphs; tau++)
      {
         for(int k=1; k<nMarkers[tau]; k++)
         {
               if(x1  >  phyDis[tau][k-1] && x2  <  phyDis[tau][k])//eligible to modify
               {
               if((i > 0 && hspPrime->at(i-1)->X2 > phyDis[tau][k-1]) || (i<=nhs-2 && hspPrime->at(i+1)->X1 < phyDis[tau][k])) continue;

               double area  = (x2   -   x1)*z;
               double lowb  = phyDis[tau][k-1];
               double highb = phyDis[tau][k];
               double B, y1, y2;

               y1 = gsl_rng_uniform_pos(gslr)*(highb-lowb) + lowb;//uniform(lb, hb)
               if(highb-y1 < y1-lowb) B = highb  -  y1; else B = y1  -  lowb;
               y2  = -1.0/lambdaHalfWidth * log(1    -    gsl_rng_uniform_pos(gslr)*(1-exp(-lambdaHalfWidth*B)));
               double logprForw = log(lambdaHalfWidth) - lambdaHalfWidth*y2 - log(1-exp(-lambdaHalfWidth*B));

               hspPrime->at(i)->X1 = y1   -   y2;
               hspPrime->at(i)->X2 = y1   +   y2;
               hspPrime->at(i)->Z = area/(hspPrime->at(i)->X2    -    hspPrime->at(i)->X1);
                                       
               y1 = 0.5*(hsp->at(i)->X1+hsp->at(i)->X2);
               y2 = 0.5*(hsp->at(i)->X2-hsp->at(i)->X1);
               if(highb-y1 < y1-lowb) B = highb  -  y1; else B = y1  -  lowb;
               double logprBackw = log(lambdaHalfWidth) - lambdaHalfWidth*y2 - log(1-exp(-lambdaHalfWidth*B));

               propRatio *= exp(logprBackw    -    logprForw);
               break;
	       }
         }
      }
   }

   return propRatio;
}




void Chain::outputAccProb()
{
  if(numProp[GEN]  > 0) cout<<(double)numAcc[GEN]/numProp[GEN]<<" "; else cout<<"NA ";
  if(numProp[TREE] > 0) cout<<(double)numAcc[TREE]/numProp[TREE]<<" ";  else cout<<"NA ";
  if(numProp[MISS] > 0) cout<<(double)numAcc[MISS]/numProp[MISS]<<" ";  else cout<<"NA "; 
  if(numProp[HAP]  > 0) cout<<(double)numAcc[HAP]/numProp[HAP]<<" "; else cout<<"NA ";
  if(numProp[ANC]  > 0) cout<<(double)numAcc[ANC]/numProp[ANC]<<" "; else cout<<"NA ";
  if(numProp[RHOBGS]>0) cout<<(double)numAcc[RHOBGS]/numProp[RHOBGS]<<" ";  else cout<<"NA ";
  if(numProp[RHOBGR]>0) cout<<(double)numAcc[RHOBGR]/numProp[RHOBGR]<<" ";  else cout<<"NA ";
  if(numProp[HOTSPC]>0) cout<<(double)numAcc[HOTSPC]/numProp[HOTSPC]<<" ";  else cout<<"NA ";
  if(numProp[HOTSPA]>0) cout<<(double)numAcc[HOTSPA]/numProp[HOTSPA]<<" ";  else cout<<"NA ";
  if(numProp[HOTSPD]>0) cout<<(double)numAcc[HOTSPD]/numProp[HOTSPD]<<" ";  else cout<<"NA ";
  if(numProp[THE]>0) cout<<(double)numAcc[THE]/numProp[THE]<<" ";  else cout<<"NA ";  
  if(numProp[SWAP]>0) cout<<(double)numAcc[SWAP]/numProp[SWAP]<<" ";  else cout<<"NA ";
}




double Chain::timeMRCA()
{
  double tmrca = genealogy->at(0)->getTMRCA();

  for(int tau=1; tau<nGraphs; tau++)
  {
  double age = genealogy->at(tau)->getTMRCA();
  if(age>tmrca) tmrca = age;
  }

  return tmrca;
}




double Chain::logPriorRhoBackg(double** rhobg, double rhobgscale)
{
  double logprior = 0;
  for(int tau=0; tau<nGraphs; tau++){
  for(int i=0; i<nMarkers[tau]-1; i++){
  logprior += Util::logGamma(rhobg[tau][i], rhoBackgShape, rhobgscale);
  }  
  }  
  
  return logprior;
}




double Chain::logPriorgam(double** gama, double gms)
{
  double logprior = 0;
  for(int tau=0; tau<nGraphs; tau++){
  for(int i=0; i<nMarkers[tau]-1; i++){
  logprior += Util::logGamma(gama[tau][i], gamshape, gms);
  }
  }

  return logprior;
}




/*
void Chain::printGenealogy()
{
   for(int tau=0; tau<nGraphs; tau++)
   {
      genealogy->at(tau)->printGraph();
      cout<<endl;
   }
   cout<<endl;
}


void Chain::printGenealogy(ofstream& ofstr)
{
   for(int tau=0; tau<nGraphs; tau++)
   { 
      genealogy->at(tau)->printGraph(ofstr);
      ofstr<<endl;
   } 
   ofstr<<endl; 
}


void Chain::sampleHaplotypes(chrHaps** postHap)
{
  vector<haplotypeCount*>::iterator vi;
  for(int tau=0; tau<nGraphs; tau++)
   {
    int nloci = nMarkers[tau];

    for(int chr=0; chr<nSeq[tau]; chr++)
      {
        Node* atip = genealogy->at(tau)->findNode(chr);
         vector<locAlle*> vec = atip->getAncLocs(LEFT);

         int* tiphap;
         try{ tiphap = new int[nloci]; }
         catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}
         for(int i=0; i<nloci; i++) {tiphap[i]=vec[i]->allele;}

         bool hapexit = false;
         for(vi=postHap[tau][chr].haps.begin(); vi!=postHap[tau][chr].haps.end(); vi++)
         {
          if(hapsAreSame((*vi)->haplotype, tiphap, nloci))
            {
              hapexit = true;
               (*vi)->count++;
               break;
            }
         }//loop over all the entries in postHap 
         if(!hapexit) //add the new haplotype to the vector             
         {
          haplotypeCount* anewhap = new haplotypeCount;
          try{anewhap->haplotype = new int[nloci];}
          catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}
            for(int i=0; i<nloci; i++) anewhap->haplotype[i]=tiphap[i];
            anewhap->count = 1;
            postHap[tau][chr].haps.push_back(anewhap);
          }
         delete tiphap;
      }//loop over all tips
    }//loop over all graphs
}
*/




bool Chain::hapsAreSame(int* hap1, int* hap2, int nloci)
{
    for(int i=0; i<nloci; i++)
    {
    if(hap1[i]  !=  hap2[i]){return false;}
    }
    return true;
}




void Chain::printRec(ofstream& ofsRec, int iter)
{
   vector<double>* recom = new vector<double>;
   recomGenealogy(recom);
   int nrec = (int)recom->size();

   ofsRec<<iter<<"\t";
   for(int i=0; i<nrec; i++)
   {
   ofsRec<<recom->at(i);
   if(i<nrec-1) ofsRec<<"\t";
   }
   ofsRec<<endl;
}




void Chain::printRho(ofstream& ofs, int iter)
{
   ofs<<iter<<"\t"<<f1<<"\t"<<f2<<"\t";
   for(int tau = 0; tau < 1; tau++){
   for(int k   = 0;   k < nMarkers[0]-1; k++){
   ofs<<"("<<rho[0][k]<<" , "<<Gamma[0][k]<<" , "<<MEANTRACT[0]<<" , "<<theta<<")";
   if(tau == nGraphs-1 && k == nMarkers[0]-2) {} 
   else ofs<<"\t";
   }}
   ofs<<endl;
}




void Chain::printRhoBackg(ofstream& ofs, int iter)
{
   ofs<<iter<<"\t";
   for(int tau = 0; tau< nGraphs; tau++)
   {
      for(int k= 0; k  < nMarkers[tau]-1; k++)
      {
      ofs<<rhoBackg[tau][k];
      if(tau==nGraphs-1 && k==nMarkers[tau]-2) {}
      else ofs<<"\t";
      }
   }
   ofs<<endl;
}




void Chain::postHaplotypes(chrHaps** postHap, ofstream& ofsHap)
{
   MapHap**** Hapdist;
   try { Hapdist = new MapHap***[nGraphs];}
   catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}

   for(int l=0; l<nGraphs; l++)
   {
    ofsHap<<"Graph # "<<l<<endl<<endl;
      try {Hapdist[l] = new MapHap**[nSeq[l]];}
      catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}

      for(int tipth=0; tipth<nSeq[l]; tipth++)
      {
        ofsHap<<"Tip "<<tipth<<endl;

          int nbins = (int)postHap[l][tipth].haps.size();
         try { Hapdist[l][tipth] = new MapHap*[nbins]; }
        catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}

         double sum=0;
         for(int j=0; j<nbins; j++) sum += postHap[l][tipth].haps[j]->count;
         for(int j=0; j<nbins; j++)
         {
          try{Hapdist[l][tipth][j] = new MapHap;}
            catch (bad_alloc ex) {cout<<"Exception occurred: "<<ex.what()<<endl; exit(1);}

            Hapdist[l][tipth][j]->key = postHap[l][tipth].haps[j]->haplotype;
            Hapdist[l][tipth][j]->frequency = postHap[l][tipth].haps[j]->count/sum;
         }

         qsort(Hapdist[l][tipth], nbins, sizeof(MapHap*), Util::compQsortFreq);

         double sumfreq =0;
         for(int j=0; j<nbins; j++)
         {
           sumfreq += Hapdist[l][tipth][j]->frequency;
           if(sumfreq<=0.95)  Hapdist[l][tipth][j]->sig =1;
           else Hapdist[l][tipth][j]->sig =0;
         }

          for(int j=0; j<nbins; j++)
          {
          int* ahap = Hapdist[l][tipth][j]->key;
            for(int k=0; k<nMarkers[l]; k++) ofsHap<<(char)ahap[k]<<" "; ofsHap<<":\t";
            ofsHap<<Hapdist[l][tipth][j]->frequency<<"\t"<<Hapdist[l][tipth][j]->sig<<endl;
         }
         ofsHap<<endl;
      }//loop over all tips
  }//loop over all graphs
}




void Chain::monitor(Chain** chains, int nchains, int num, ofstream& ofs)
{
  ofs<<num<<"\t";
  for(int i=0; i<nchains; i++){
  double lll = chains[i]->logllhGenealogy(chains[i]->getGenealogy(), chains[i]->getTheta());
  double lprior = chains[i]->logPriorGenealogy(chains[i]->getGenealogy(), chains[i]->getRho());
  lprior += chains[i]->logPriorRhoBackg(chains[i]->getRhoBackg(), chains[i]->getRhoBackgScale());
  lprior += chains[i]->logPriorHotspots(chains[i]->getHotspots());
  ofs<<lprior<<"\t"<<lll<<"\t"<<lprior+lll<<"\t"<<chains[i]->getHotspots()->size()<<"\t";
  }
  ofs<<endl;
}
