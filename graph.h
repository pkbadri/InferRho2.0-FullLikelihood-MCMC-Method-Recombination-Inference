/*
Copyright 2009 Ying Wang (ygwang@uchicago.edu)
Created 2007 at University of California, Davis.
*/

/*
Two methods for calculating likelihood given the genealogy.
1) Assuming only one mutation per site. It is appropriate for low mutation rate. 
2) Any number of mutations per site is possible. The ancestral states at marker locations need to be specified.
Modified on Nov 12, 2007: implement post-order traversal functions
*/

#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <ctype.h>
#include <time.h>
#include <cstring>
#include "gsl/gsl/gsl_rng.h"
#include "gsl/gsl/gsl_randist.h"
using namespace std;


#define TIP -1 
#define REC 0 
#define COAL 1

#define LEFT 1 
#define RIGHT 2

#define LEFT_LEFT 3
#define LEFT_RIGHT 4

#define EPSILON 0.000001
#define MAXNRECS 100000
#define SQR2PI 2.5066
#define PERMB 0.000001

typedef int ALLELE;

class Node
{
  private:
  int ID;
  Node *leftA;
  Node *rightA;
  Node *leftD;
  Node *rightD;
  double coalTime, recTime;
  double bkpoint, tractlength;	
  int follow;  
  ALLELE* leftAncLocs;  
  ALLELE* rightAncLocs;
  gsl_rng* gslr;
 
  public: 
  double  totancs, prob, lendist, lta, rta, col, cor;
  double  temptime, tacs, tprob, tlendist, tlta, trta, tcol, tcor;
  int ls, le, rs, re, mark;

  Node(int, double, double);
  Node(Node&);
  ~Node();
  void setID(int id) {ID  =  id;}
  void setleftA(Node* lefta){leftA  =  lefta;}
  void setrightA(Node* righta){rightA  =  righta;}
  void setleftD(Node* leftd){leftD  =  leftd;}
  void setrightD(Node* rightd){rightD  =  rightd;}
  void setcoalTime(double ct) {coalTime  =  ct;}
  void setrecTime(double rt)  {recTime  =  rt;}
  void setTime(double time){if(recTime != 0){recTime = time;} else{coalTime = time;}}
  int  getID(){return ID;}
  Node* getleftA() {return leftA;}
  Node* getrightA(){return rightA;}
  Node* getAnc(){if(leftA != NULL){return leftA;} else{return rightA;}}
  Node* getDec(){if(leftD != NULL){return leftD;} else{return rightD;}}
  Node* getleftD(){ return leftD;}
  Node* getrightD(){ return rightD;}
  double getcoalTime(){return coalTime;}
  double getrecTime(){return  recTime;}
  double getTime(){if(recTime != 0){return recTime;} else{return coalTime;}}
  bool isTip();
  bool isCoal(){if(coalTime != 0){return true;} else{return false;}}
  bool isRec(){if(recTime   != 0){return true;} else{return false;}}
  	
  void setBkpoint(double bk) {bkpoint  =  bk;}
  double getBkpoint() {return bkpoint;}
  void   settractlength(double tl) {tractlength = tl;}
  double gettractlength() {return tractlength;}	
  int getFollow(){return follow;}
  void setFollow(int fl) {follow  =  fl;}

  ALLELE* getAncLocs(int dir) {if(dir == LEFT){return leftAncLocs;} else{return rightAncLocs;}}
  void newtract(ALLELE*, double, gsl_rng*, double, int, int);
  void  setconv(ALLELE*, double, gsl_rng*, double, int, int, int);
  void changemean(ALLELE*, double, double, int, int, int);
  void unionAncLocs(ALLELE*, ALLELE*);
  void chunkAncLocs(ALLELE*);
  void convAncLocs(ALLELE*);
  double  updateAncLocs(double, bool, bool, gsl_rng*, double, double*, double*, int);
  void    updateancs();
  double updateAncLocs1(double, gsl_rng*, double, double*, double*);
  
  ALLELE* getLeftAncLocs()  {return leftAncLocs;}
  ALLELE* getRightAncLocs() {return rightAncLocs;}
  void setLeftAncLocs(ALLELE*);
  void setRightAncLocs(ALLELE*);
};




typedef struct{
   Node* node;
   vector<Node*> ancestors;
}nodeAncestors;




typedef struct{
   Node* dec;
   Node* anc;
}decAnc;




typedef struct{
   Node* dec;
   int ancdir;
}decAncdir;




class Graph
{
  private:
  int numTips;    
  int numLoci;

  Node* Root;   
  double TotalTime;
  int numRecs;int numCO;int numGC;
  int numInternNodes; 
  vector<Node*> NodesVector;
  gsl_rng* gslr;

  public:
  Graph(int, int, gsl_rng*);
  Graph(Graph&);
  ~Graph();
  Graph& operator=(Graph&);
  void settotal(double tot){TotalTime = tot;}
  //void changetime(double);
  void copyGraph(Graph&, Graph&);
  void destroy();
  //void resetGraph(vector<Node*>&);  
  //void printGraph(ofstream&);
  //void printGraph();
  int getnumRecs(){return numRecs;} int getnumCO(){return numCO;} int getnumGC(){return numGC;} 
  Node* getRoot(){return Root;}
  double getTMRCA();
  int getNumInternNodes(){ return numInternNodes;}
  vector<Node*> getNodesVector(){return NodesVector;}
  void updatelkup(double, int);
  void revertwait();
  //void traverse(vector<Node*>&, Node*);
  bool containsNode(vector<Node*>&, int);
  void     simuBNTree(double*, int**, double, double, double*, double*, int);
  void createTipNodes(double*, int**, double, double, double*, double*, int);
  void updatemavec();

  void  setnewprob(double, double, double*, double*, int); 
  void  setnewprob1(double, double, double*, double*, int);
  void  setnewprob111(double, double, double*, int, int*);  
  void  setnewprob2(double, double, double*, double*, int); 
  void  revertprob();
  void  revertprob1();
  void  revertprob111();
  void  revertprob2();
 
  //void simuGenealogy(double);
  //void addMutations(double, double, vector<double>&);
  //void forwardTraversal(vector<double>&, double);
  //void mutSequency(int*, vector<double>&, double*, double, double);
  int ranchar();
  int ranchar2(int);
  int ranchar3(int, int);
  double getTottimeGraph();
  
  //void simuGenealogyMarkers(double*, double, double*);
  //void simuHaplotypes(double*, double, double, double*, double*,int**, double);
  void getTreesFast(vector<Node*>**);
  //double totalTreeLengths(vector<Node*>*);
  //int simuSNP(double*);

  //void getConnections(vector<Node*>&, Node*);
  void insertNodeBtTwoNodes(Node*, Node*, Node*);
  bool containsLoop();
  double     addRecCoal(double, double, double*, double*, int);
  double    addconvCoal(double, double, double*, double*, int);
  double  deleteRecCoal(double, double, double*, double*, int);
  double deleteconvCoal(double, double, double*, double*, int);

  void insertNode(Node*, Node*, int);
  double changeBkpointFollow(double, double, double*, double*, int);
  double   changetractlength(double, double, double*, double*);
 
  //double rotate();
  void eligibleRecBranches(vector<decAncdir>&, double);
  void eligibleCoalBranches(vector<decAncdir>&, double);
  void eligiblePairRec(vector<decAncdir>&);
  void eligiblePairConv(vector<decAncdir>&);

  double moveNode(double, double, double, double*, double*, int);
  double changeRec(Node*, double, double, double, double*, double*, int);
  double changeCoal(Node*, double, double, double, double*, double*, int);
  //double changeRoot(double);
 
  Node* findNode(int);
  void changeOneAncSeq();
  double getlogprior(double*, double*, double*, double);
  //double getlogpriorARG(double);
  double logllh(double, double, double*);
  double logllhGraphSingleSite(double, double*, int);
  double logllhGraph(double, double*);
  //double logllhGraph1(double, double*, double);
  //double logllhGraph2(double, double*, double);
  //double logllhTrees(double, vector<Node*>**);
  ALLELE* transHap(Node* dec, Node* cursor); 
};




typedef struct {
  int* key;
  double frequency;
  int sig;
}MapHap;




class Util
{
public:
   static int compQsortFreq(const void* a , const void*);
   static bool compareTwoNodesT(Node*, Node*);
   static bool compareTwoNodesID(Node*, Node*);
   static bool compareFirstConnection(nodeAncestors*, nodeAncestors*);
   static bool compareTwoNumbers(double, double);
   static int randNuc(int, gsl_rng*);
   static double logGamma(double, double, double);
   static double logNormal(double, double, double);
};
