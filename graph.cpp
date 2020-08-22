/*
Copyright  2009   Ying Wang  (ygwang@uchicago.edu).
Created    2007   at University of California, Davis.
*/

#include "graph.h"
double lkup[10][1000000], cum[10][1000000], *scaled;
int snpcount, *pos, nch = 0;

int maximum(int a, int b){if(a > b){return a;} else{return b;}}
int minimum(int a, int b){if(a < b){return a;} else{return b;}}


//Finds the integer closest to the floating number x
int ni(double x){
if(x >= 0){
if(x - int(x) <  0.5) return(int(x));
else return( 1  +  int(x));	
}

else{
if(x - int(x) > -0.5) return(int(x));
else return(-1  +  int(x));	
}
}




//Revert the chain's lookup table to the old tract length                                                                                                                      
void Graph::updatelkup(double oldmean, int chain_number){
int i, k = chain_number;

lkup[k][0] = 1.0; cum[k][0] = 1.0;
for(i = 1; i < 1000000; i++){
lkup[k][i]  = lkup[k][-1+i]*(1.0  - (1.0/oldmean));
cum[k][i] = cum[k][-1+i] + lkup[k][i];
}

}




//Calculates the rate of informative gene-conversion events for a branch
double ancs_tot(double* gamma, ALLELE* vec, double MEANTRACT, int start, int end, int chain_number){
if(end - start < 2){return(0);}

else{
int i, k = chain_number; double tacs;

//Part before pos[0]
if(start > 0){
tacs = gamma[0]*(cum[k][ni(50*MEANTRACT) + pos[start] - 1] - cum[k][pos[start] - 1]   - (cum[k][pos[end-1] + ni(50*MEANTRACT) - 1]   -   cum[k][pos[end-1] - 1]));
}
else{tacs = gamma[0]*(cum[k][ni(50*MEANTRACT) + pos[start] - 1]   -       (cum[k][pos[end-1] + ni(50*MEANTRACT) - 1]   -   cum[k][pos[end-1] - 1]));}

//From pos[0]  to  pos[start]
for(i = 0; i<start; i++){
if(i < start-1){
tacs += gamma[i]*(cum[k][pos[start] - pos[i] - 1]  -  cum[k][pos[start]  - pos[i+1] - 1] -   (cum[k][pos[end-1] - pos[i] - 1]  -  cum[k][pos[end-1] -  pos[i+1] - 1]));
}
else{tacs += gamma[i]*(cum[k][pos[start] - pos[i] - 1]    -    (cum[k][pos[end-1] - pos[i] - 1]   -   cum[k][pos[end-1] - pos[start] - 1]));}
}

//From pos[start]   to   pos[end-1]
int pre = start, j;
for(i = 1 + start; i<end;  i++){
if(vec[i] != 999){
j = i - 1; 
tacs += gamma[j]*cum[k][pos[i] - pos[j] - 1];

while(j > pre){
tacs += gamma[j-1]*(cum[k][pos[i] - pos[j-1] - 1]   -   cum[k][pos[i] - pos[j] - 1]);
--j;}

pre = i;}
}

assert(tacs > 0);return(tacs*PERMB);
}
}




//Calculates the rate of informative gene-conversion events for a branch                                                                                                       
double ancs_tot111(double* gamma, ALLELE* vec, double MEANTRACT, int start, int end, int chain_number){  
if(end - start < 2){return(0);}

else{
int i, k = chain_number; double tacs;

//Part before pos[0]                                                                                                                                                           
if(start > 0){
tacs = gamma[0]*(cum[k][ni(50*MEANTRACT) + pos[start] - 1]    -    cum[k][pos[start] - 1]    -    (cum[k][pos[end-1] + ni(50*MEANTRACT) - 1]   -   cum[k][pos[end-1] - 1]));
}
else{tacs = gamma[0]*(cum[k][ni(50*MEANTRACT) + pos[start] - 1]    -    (cum[k][pos[end-1] + ni(50*MEANTRACT) - 1]   -   cum[k][pos[end-1] - 1]));}

//From pos[0]  to  pos[start]                                                                                                                                                  
for(i = 0; i<start; i++){
if(i < start-1){
tacs += gamma[i]*(cum[k][pos[start] - pos[i] - 1]  -  cum[k][pos[start]  - pos[i+1] - 1]   -   (cum[k][pos[end-1] - pos[i] - 1]  -  cum[k][pos[end-1] -  pos[i+1] - 1]));
}
else{tacs += gamma[i]*(cum[k][pos[start] - pos[i] - 1]  - (cum[k][pos[end-1] - pos[i] - 1] - cum[k][pos[end-1] - pos[start] - 1]));}
}


//From pos[start]   to   pos[end-1]                                                                                                                                            
int pre = start, j;
for(i = 1 + start; i<end;  i++){
if(vec[i] != 999){
j = i - 1;
tacs += gamma[j]*cum[k][pos[i] - pos[j] - 1];

while(j > pre){
tacs += gamma[j-1]*(cum[k][pos[i] - pos[j-1] - 1]   -   cum[k][pos[i] - pos[j] - 1]);
--j;}

pre = i;}
}


assert(tacs > 0);return(tacs*PERMB);
}
}




//Calculates a number proportional to informative gene-conversions from ancestry vector                                                                                        
double ancs_tot1(ALLELE* vec, double MEANTRACT, int start, int end, int chain_number){
if(end - start < 2){return(0);}

else{
int i, k = chain_number;
double tacs = cum[k][ni(50*MEANTRACT) - 1]   -   (cum[k][pos[end-1] - pos[start] + ni(50*MEANTRACT) - 1]     -     cum[k][pos[end-1] - pos[start] - 1]);
int previous = pos[start]; for(i = 1 + start; i<end;  i++){if(vec[i] != 999){tacs   +=   cum[k][pos[i] - previous - 1];  previous = pos[i];}}
assert(tacs > 0);return(tacs);
}

}




//Calculates the rate of informative crossing-over events for a branch
double rho_tot(double* rho, int start, int end){
int i; 
double rholen = 0;for(i = start; i< end - 1; i++){rholen += rho[i]*double(pos[i+1] - pos[i]);}
return(rholen*PERMB);
}




//Creates a node 
Node::Node(int id, double ct, double rt):
ID(id), coalTime(ct), recTime(rt)
{
leftD = NULL; rightD = NULL; leftA = NULL; rightA = NULL; 
leftAncLocs = NULL; rightAncLocs = NULL;
bkpoint = -1; follow = -1; tractlength = -1; 
totancs = -1; prob = -1; lendist = -1; 
lta = -1; rta = -1; 
col = -1; cor = -1; 
ls  = -1; le  = -1; 
rs  = -1; re  = -1;
}




//Copies a node
Node::Node(Node& x)
{
ID          = x.ID;
coalTime    = x.coalTime;
recTime     = x.recTime;
bkpoint     = x.bkpoint;
follow      = x.follow;
tractlength = x.tractlength;	
totancs  = x.totancs;
prob     = x.prob;
lendist  = x.lendist;	
lta  = x.lta;  rta = x.rta;
col  = x.col;  cor = x.cor;
ls   =  x.ls;   le  = x.le;
rs   =  x.rs;   re  = x.re;

if(x.leftAncLocs  !=  NULL){
leftAncLocs  = new ALLELE[snpcount];
for(int i=0; i<snpcount; i++){leftAncLocs[i] = x.leftAncLocs[i];}
}
else{leftAncLocs  = NULL;}

if(x.rightAncLocs != NULL){
rightAncLocs = new ALLELE[snpcount];
for(int i=0; i<snpcount; i++){rightAncLocs[i] = x.rightAncLocs[i];}
}
else{rightAncLocs = NULL;}
   	
leftD = NULL; rightD = NULL; leftA = NULL; rightA = NULL;    
}




Node::~Node()
{
    if(this  !=  NULL)
    {
    leftD = NULL; rightD = NULL; leftA = NULL; rightA = NULL;
    delete [] leftAncLocs; delete [] rightAncLocs;
    leftAncLocs = NULL; rightAncLocs = NULL;
    }
}




bool Node::isTip()
{
   if(coalTime == 0 && recTime == 0) return true;
   else return false;
}




void Node::setLeftAncLocs(ALLELE* x)
{
  if(leftAncLocs != NULL){delete [] leftAncLocs; leftAncLocs = NULL;}

  leftAncLocs = x;assert(leftAncLocs != NULL);

  ls = -1; le = -1;
  for(int i=0; i<snpcount; i++)    {if(leftAncLocs[i] != 999) {ls = i;   break;}}
  for(int i=snpcount-1; i>-1; i--) {if(leftAncLocs[i] != 999) {le = i+1; break;}}
  assert(ls != -1 && le != -1);
}




void Node::setRightAncLocs(ALLELE* x)
{
  if(rightAncLocs != NULL){delete [] rightAncLocs; rightAncLocs = NULL;}

  rightAncLocs = x;assert(rightAncLocs != NULL);

  rs = -1; re = -1;
  for(int i=0; i<snpcount; i++)    {if(rightAncLocs[i] != 999) {rs = i;   break;}}
  for(int i=snpcount-1; i>-1; i--) {if(rightAncLocs[i] != 999) {re = i+1; break;}}
  assert(rs != -1 && re != -1);
}




void Node::chunkAncLocs(ALLELE* locs)
{	
  ALLELE* vec1; ALLELE* vec2;
  vec1 = new ALLELE[snpcount];
  vec2 = new ALLELE[snpcount];

  for(int i = 0; i<snpcount; i++){
  if(scaled[i] > bkpoint){vec1[i] = 999; vec2[i] = locs[i];}
  else {vec1[i] = locs[i]; vec2[i] = 999;}
  }

  delete [] leftAncLocs;delete [] rightAncLocs;leftAncLocs = NULL; rightAncLocs = NULL;

  if(follow == LEFT_LEFT) {leftAncLocs = vec1; rightAncLocs = vec2;}
  else{leftAncLocs = vec2; rightAncLocs = vec1;}

  ls = -1; le = -1; rs = -1; re = -1;
  for(int i=0; i<snpcount; i++)   { if(leftAncLocs[i]  != 999) {ls =   i; break;}}
  for(int i=snpcount-1; i>-1; i--){ if(leftAncLocs[i]  != 999) {le = i+1; break;}}
  for(int i=0; i<snpcount; i++)   { if(rightAncLocs[i] != 999) {rs =   i; break;}}
  for(int i=snpcount-1; i>-1; i--){ if(rightAncLocs[i] != 999) {re = i+1; break;}}
  assert(ls != -1 && le != -1 && rs != -1 && re != -1);
}




void Node::convAncLocs(ALLELE* locs)
{
  ALLELE* vec1; ALLELE* vec2;
  vec1 = new ALLELE[snpcount];
  vec2 = new ALLELE[snpcount];

  for(int i=0; i<snpcount; i++){
  if(scaled[i] - bkpoint > EPSILON  && scaled[i] - bkpoint - tractlength < EPSILON){vec1[i] = 999; vec2[i] = locs[i];}
  else{vec1[i] = locs[i]; vec2[i] = 999;}
  }

  delete [] leftAncLocs;delete [] rightAncLocs;leftAncLocs = NULL; rightAncLocs = NULL;
  
  if(follow == LEFT_LEFT){leftAncLocs = vec1; rightAncLocs = vec2;}
  else {leftAncLocs = vec2; rightAncLocs = vec1;} 

  ls = -1; le = -1; rs = -1; re = -1;
  for(int i=0; i<snpcount; i++)     {if(leftAncLocs[i] != 999) {ls =   i;   break;}}
  for(int i=snpcount-1; i>-1; i--)  {if(leftAncLocs[i] != 999) {le = i+1;   break;}}
  for(int i=0; i<snpcount; i++)     {if(rightAncLocs[i]!= 999) {rs =   i;   break;}}
  for(int i=snpcount-1; i>-1; i--)  {if(rightAncLocs[i]!= 999) {re = i+1;   break;}}
  assert(ls != -1 && le != -1 && rs != -1 && re != -1);
}




void Node::unionAncLocs(ALLELE* locs1, ALLELE* locs2)
{ 
  assert(rightAncLocs == NULL);
  
  ALLELE* leftanclocs = new ALLELE[snpcount];
  for(int i=0; i<snpcount; i++){
  leftanclocs[i]  =  locs1[i];
  if(leftanclocs[i] == 999 && locs2[i] != 999){leftanclocs[i] = locs2[i];}
  }
  
  setLeftAncLocs(leftanclocs);
}




//Set the variables in a gene-conversion node                                                                                                                                 
void Node::newtract(ALLELE* vec, double length, gsl_rng* gslr, double MEANTRACT, int start, int end){
/*

int i, k, x = ni(length*bkpoint);

//Find the lookup row
for(k = 0; k < 9; k++){if(fabs(mn[k] - MEANTRACT) < EPSILON){break;}}   assert(k < 9);

//Find the position of the first ancestral marker after the gene-conversion start
for(i = start; i < end; i++){if(vec[i] != 999 && pos[i] > x){break;}} assert(i < end);

int flag = 0, geom;
while(flag == 0){

if(i > start){
geom = gsl_ran_geometric(gslr, 1.0/MEANTRACT);
lendist = lkup[k][geom - 1]/MEANTRACT;
tractlength = -bkpoint  + double(geom + pos[i] - 1)/length;
}

else{
double z = drand48(), sum = 0;

geom = pos[i] - x;
while(geom < pos[end-1] - x){
sum += lkup[k][geom - 1]/(prob*MEANTRACT);
if(sum > z){break;}
++geom;
}

lendist = lkup[k][geom - 1]/(prob*MEANTRACT);
tractlength = double(geom)/length;
}

assert(scaled[i] - bkpoint - tractlength < EPSILON); 

if(lendist == 0){continue;}else{flag = 1; break;}
}
*/
}




//Set the variables in a gene-conversion node
void Node::setconv(ALLELE* vec, double length, gsl_rng* gslr, double MEANTRACT, int start, int end, int chain_number){
int i, k = chain_number; assert(end - start > 1);
totancs = cum[k][ni(50*MEANTRACT) - 1] -  cum[k][pos[end-1] - pos[start] + ni(50*MEANTRACT) - 1]  +  cum[k][pos[end-1] - pos[start] - 1];
int previous = pos[start]; for(i = 1 + start; i < end; i++){if(vec[i] != 999){totancs += cum[k][pos[i] - previous - 1]; previous = pos[i];}}

int flag = 0, x, j, geom;
while(flag == 0){
x = pos[start] - ni(50*MEANTRACT) +  int(drand48()*(pos[end-1] - pos[start] + ni(50*MEANTRACT)));

for(i = start; i < end; i++){if(vec[i] != 999 && pos[i] > x){break;}} assert(i < end);

if(i > start){prob = lkup[k][pos[i] - x - 1];}
else{prob =  lkup[k][pos[i] - x - 1] - lkup[k][pos[end-1] - x - 1];}

if(prob == 0){continue;}
    
if(drand48() < prob){
setBkpoint(double(x)/length);

if(i > start){
geom = gsl_ran_geometric(gslr, 1.0/MEANTRACT);
lendist = lkup[k][geom - 1]/MEANTRACT; 
settractlength(-bkpoint  + double(geom + pos[i] - 1)/length);
}

else{
double z = drand48();
geom = pos[i] - x;
double sum = 0;
while(geom < pos[end-1] - x){
sum += lkup[k][geom - 1]/(prob*MEANTRACT);
if(sum > z){break;}
++geom; 
}

lendist = lkup[k][geom - 1]/(prob*MEANTRACT); 
settractlength(double(geom)/length);
}

if(lendist == 0){continue;}
flag = 1; break;}
}

}




//Change the variables in the graph based on the new mean tract length
void Node::changemean(ALLELE* vec, double length, double newmean, int start, int end, int chain_number){

//Copy existing variables to temporary ones
tacs = totancs;
tprob = prob;
tlendist = lendist;

//Proposed values for the variables prob, lendist and totancs                                                                                                                                         
totancs = ancs_tot1(vec, newmean, start, end, chain_number);

int x, i;
x = ni(length*bkpoint); for(i = start; i<end; i++){if(vec[i] != 999 && pos[i] > x){break;}} assert(i < end);

if(i > start){
prob    = lkup[chain_number][pos[i] - x - 1];
lendist = lkup[chain_number][ni(length*tractlength) - 1]/(prob*newmean);
}

else{
prob    = lkup[chain_number][pos[i] - x - 1] - lkup[chain_number][pos[end-1] - x - 1];
lendist = lkup[chain_number][ni(length*tractlength) - 1]/(prob*newmean);
}

}




void Graph::setnewprob111(double length, double newmean, double* gamma, int chain_number, int* flag){
int i, k = chain_number;

//Generate a chain's lookup table for the new tract length                                                                                                                  
lkup[k][0] = 1.0; cum[k][0] = 1.0;
for(i = 1; i < 1000000; i++){
lkup[k][i]  = lkup[k][-1+i]*(1.0  - (1.0/newmean));
cum[k][i] = cum[k][-1+i] + lkup[k][i];
}


for(i = 0; i< numTips; i++){
NodesVector[i]->tacs    = NodesVector[i]->totancs;
NodesVector[i]->totancs = ancs_tot111(gamma, NodesVector[i]->getLeftAncLocs(), newmean, NodesVector[i]->ls, NodesVector[i]->le, chain_number);
}


for(i = numTips; i<(int)NodesVector.size(); i++){
if(NodesVector[i]->isRec()){
NodesVector[i]->tlta = NodesVector[i]->lta; NodesVector[i]->trta = NodesVector[i]->rta;
NodesVector[i]->lta  = ancs_tot111(gamma, NodesVector[i]->getLeftAncLocs(),  newmean, NodesVector[i]->ls, NodesVector[i]->le, chain_number);
NodesVector[i]->rta  = ancs_tot111(gamma, NodesVector[i]->getRightAncLocs(), newmean, NodesVector[i]->rs, NodesVector[i]->re, chain_number);

if(NodesVector[i]->gettractlength() > 0){
Node* dec = NodesVector[i]->getDec();
ALLELE* vec; int start, end;
if(!dec->isRec() || (dec->isRec() && dec->getleftA()->getID() == NodesVector[i]->getID())) {vec = dec->getLeftAncLocs(); start = dec->ls; end = dec->le;}
else{vec = dec->getRightAncLocs(); start = dec->rs; end = dec->re;}

assert(vec != NULL  && end - start > 1);
NodesVector[i]->changemean(vec, length, newmean, start, end, chain_number);
if(NodesVector[i]->prob == 0 || NodesVector[i]->lendist == 0){flag[0] = 1;}
}

else{assert(NodesVector[i]->totancs == -1  &&  NodesVector[i]->prob == -1  &&  NodesVector[i]->lendist == -1);}
}

else{
NodesVector[i]->tacs = NodesVector[i]->totancs;
NodesVector[i]->totancs = ancs_tot111(gamma, NodesVector[i]->getLeftAncLocs(), newmean, NodesVector[i]->ls, NodesVector[i]->le, chain_number);
}
}

}




//Revert to older values
void Graph::revertprob111(){
int i;
for(i = 0; i < numTips; i++){NodesVector[i]->totancs = NodesVector[i]->tacs;}

for(i = numTips; i<(int)NodesVector.size(); i++){
if(NodesVector[i]->isRec()){
NodesVector[i]->lta = NodesVector[i]->tlta;
NodesVector[i]->rta = NodesVector[i]->trta;

if(NodesVector[i]->gettractlength() > 0){
NodesVector[i]->totancs = NodesVector[i]->tacs;
NodesVector[i]->prob = NodesVector[i]->tprob;
NodesVector[i]->lendist = NodesVector[i]->tlendist;
}
}

else{NodesVector[i]->totancs = NodesVector[i]->tacs;}
}

}




void Graph::setnewprob(double length, double MEANTRACT, double* gamma, double* rho, int chain_number){
int i;
for(i = 0; i< numTips; i++){
NodesVector[i]->tacs    = NodesVector[i]->totancs;
NodesVector[i]->totancs = ancs_tot(gamma, NodesVector[i]->getLeftAncLocs(), MEANTRACT, NodesVector[i]->ls, NodesVector[i]->le, chain_number);
NodesVector[i]->tcol    = NodesVector[i]->col;
NodesVector[i]->col     = rho_tot(rho, NodesVector[i]->ls, NodesVector[i]->le);
}

for(i = numTips; i<(int)NodesVector.size(); i++){
if(NodesVector[i]->isRec()){
NodesVector[i]->tlta = NodesVector[i]->lta; NodesVector[i]->trta = NodesVector[i]->rta;
NodesVector[i]->lta  = ancs_tot(gamma,  NodesVector[i]->getLeftAncLocs(),  MEANTRACT,  NodesVector[i]->ls, NodesVector[i]->le, chain_number);
NodesVector[i]->rta  = ancs_tot(gamma,  NodesVector[i]->getRightAncLocs(), MEANTRACT,  NodesVector[i]->rs, NodesVector[i]->re, chain_number);
    
NodesVector[i]->tcol = NodesVector[i]->col; NodesVector[i]->tcor = NodesVector[i]->cor;
NodesVector[i]->col  = rho_tot(rho, NodesVector[i]->ls, NodesVector[i]->le);
NodesVector[i]->cor  = rho_tot(rho, NodesVector[i]->rs, NodesVector[i]->re);
}

else{
NodesVector[i]->tacs    = NodesVector[i]->totancs;
NodesVector[i]->totancs = ancs_tot(gamma, NodesVector[i]->getLeftAncLocs(), MEANTRACT, NodesVector[i]->ls, NodesVector[i]->le, chain_number);
NodesVector[i]->tcol    = NodesVector[i]->col;
NodesVector[i]->col     = rho_tot(rho, NodesVector[i]->ls, NodesVector[i]->le);
}
}

}




//Revert to older values                                                                                                                                                       
void Graph::revertprob(){
int i;
for(i = 0;i< numTips;i++){
NodesVector[i]->totancs = NodesVector[i]->tacs;
NodesVector[i]->col = NodesVector[i]->tcol;
}

for(i = numTips; i<(int)NodesVector.size(); i++){
if(NodesVector[i]->isRec()){
NodesVector[i]->lta   =   NodesVector[i]->tlta;
NodesVector[i]->rta   =   NodesVector[i]->trta;
NodesVector[i]->col   =   NodesVector[i]->tcol;
NodesVector[i]->cor   =   NodesVector[i]->tcor;
}

else{
NodesVector[i]->totancs = NodesVector[i]->tacs;
NodesVector[i]->col = NodesVector[i]->tcol;
}
}

}




//Set new values in a graph based on proposed mean tract length
void Graph::setnewprob1(double length, double MEANTRACT, double* gamma, double* rho, int chain_number){
int i;
for(i = 0; i< numTips; i++){
NodesVector[i]->tacs    = NodesVector[i]->totancs; 
NodesVector[i]->totancs = ancs_tot(gamma, NodesVector[i]->getLeftAncLocs(), MEANTRACT, NodesVector[i]->ls, NodesVector[i]->le, chain_number);
}

for(i = numTips; i<(int)NodesVector.size(); i++){
if(NodesVector[i]->isRec()){
NodesVector[i]->tlta = NodesVector[i]->lta; NodesVector[i]->trta = NodesVector[i]->rta;
NodesVector[i]->lta  = ancs_tot(gamma,  NodesVector[i]->getLeftAncLocs(),  MEANTRACT,  NodesVector[i]->ls, NodesVector[i]->le, chain_number);
NodesVector[i]->rta  = ancs_tot(gamma,  NodesVector[i]->getRightAncLocs(), MEANTRACT,  NodesVector[i]->rs, NodesVector[i]->re, chain_number); 
}

else{
NodesVector[i]->tacs    = NodesVector[i]->totancs; 
NodesVector[i]->totancs = ancs_tot(gamma, NodesVector[i]->getLeftAncLocs(), MEANTRACT, NodesVector[i]->ls, NodesVector[i]->le, chain_number);
}
}

}




void Graph::setnewprob2(double length, double MEANTRACT, double* gamma, double* rho, int chain_number){
int i;
for(i = 0; i< numTips; i++){
NodesVector[i]->tcol = NodesVector[i]->col;
NodesVector[i]->col = rho_tot(rho, NodesVector[i]->ls, NodesVector[i]->le);
}

  
for(i = numTips; i<(int)NodesVector.size(); i++){
if(NodesVector[i]->isRec()){
NodesVector[i]->tcol = NodesVector[i]->col; NodesVector[i]->tcor = NodesVector[i]->cor;
NodesVector[i]->col = rho_tot(rho, NodesVector[i]->ls, NodesVector[i]->le);
NodesVector[i]->cor = rho_tot(rho, NodesVector[i]->rs, NodesVector[i]->re);
}

else{
NodesVector[i]->tcol = NodesVector[i]->col;
NodesVector[i]->col  = rho_tot(rho, NodesVector[i]->ls, NodesVector[i]->le);
}
}

}




//Revert to older values
void Graph::revertprob1(){
int i;
for(i = 0;i< numTips; i++){NodesVector[i]->totancs = NodesVector[i]->tacs;}

for(i = numTips; i<(int)NodesVector.size(); i++){
if(NodesVector[i]->isRec()){
NodesVector[i]->lta   =   NodesVector[i]->tlta;
NodesVector[i]->rta   =   NodesVector[i]->trta;
}

else{NodesVector[i]->totancs = NodesVector[i]->tacs;}
}

}




//Revert to older values                                                                                                                                                       
void Graph::revertprob2(){
int i;
for(i = 0;i< numTips;i++){NodesVector[i]->col = NodesVector[i]->tcol;}

for(i = numTips; i<(int)NodesVector.size(); i++){
if(NodesVector[i]->isRec()){
NodesVector[i]->col   =   NodesVector[i]->tcol;
NodesVector[i]->cor   =   NodesVector[i]->tcor;
}

else{NodesVector[i]->col = NodesVector[i]->tcol;}
}

}




void Node::updateancs(){
  if(isRec()){
  Node* dec = getDec(); ALLELE* vec; 
  if(!dec->isRec() || (dec->isRec() && dec->leftA->ID == ID)){vec = dec->leftAncLocs;}
  else{vec = dec->rightAncLocs;}

  //If the node is a gene-conversion node                                                                                                                                      
  if(tractlength > 0){convAncLocs(vec);}

  //If the node is a crossing-over node                                                                                                                                     
  else{chunkAncLocs(vec);}
  }


  else
  {
  Node* dec1 = leftD; Node* dec2 = rightD; ALLELE* vec1; ALLELE* vec2;

  if(dec1->getID() == dec2->getID()){vec1 = dec1->leftAncLocs; vec2 = dec1->rightAncLocs;}
  else{
  if(!dec1->isRec() || (dec1->isRec() && dec1->leftA->ID == ID)){vec1 = dec1->leftAncLocs;}
  else{vec1 = dec1->rightAncLocs;}

  if(!dec2->isRec() || (dec2->isRec() && dec2->leftA->ID == ID)){vec2 = dec2->leftAncLocs;}
  else{vec2 = dec2->rightAncLocs;}
  }

      
  assert(vec1 != NULL && vec2 != NULL);
  unionAncLocs(vec1, vec2);
  }
}



void Graph::updatemavec(){
for(int i=numTips; i<(int)NodesVector.size(); i++){NodesVector[i]->updateancs();}
}



double Node::updateAncLocs(double length, bool changeBkpoint, bool changeFollow, gsl_rng* gslr, double MEANTRACT, double* gamma, double* rho, int chain_number){     
double propRatio = 1;

if(isRec()){        		 
Node* dec = getDec(); ALLELE* vec; int start, end;
if(!dec->isRec() || (dec->isRec() && dec->leftA->ID == ID)){vec = dec->leftAncLocs; start = dec->ls; end = dec->le;}
else{vec = dec->rightAncLocs; start = dec->rs; end = dec->re;}
      
if(end - start < 2){return 0;}

if(mark == -1){return(propRatio);}
    				
//If the node is a gene-conversion node 
if(tractlength > 0){		                
//Check if the MA vector has changed due to move
double flag = 0;	
if(start != minimum(ls, rs) || end != maximum(le, re)){flag = 1;}
else{
for(int i = start; i < end; i++){
int temp = 999;
if(leftAncLocs[i]  != 999){temp  = leftAncLocs[i];}
if(rightAncLocs[i] != 999){assert(temp == 999); temp  = rightAncLocs[i];}
if((vec[i] != 999 && temp == 999) || (vec[i] == 999 && temp != 999)){flag = 1; break;}	
}
}	

      
//Propose  a  new  breakpoint and tractlength		
if(changeBkpoint || flag == 1){
double tal = totancs, pp = prob, ld = lendist;	
setconv(vec, length, gslr, MEANTRACT, start, end, chain_number);
propRatio *= exp(log(pp) + log(totancs) + log(ld) - log(prob) - log(tal) - log(lendist));
}
     			
if(changeFollow){if(gsl_rng_uniform_pos(gslr) < 0.5){follow = LEFT_LEFT;} else{follow = LEFT_RIGHT;}}
      
convAncLocs(vec);
lta = ancs_tot(gamma, leftAncLocs,  MEANTRACT, ls, le, chain_number); col = rho_tot(rho, ls, le);
rta = ancs_tot(gamma, rightAncLocs, MEANTRACT, rs, re, chain_number); cor = rho_tot(rho, rs, re);
}
		 		 		 					
      
//If the node is a crossing-over node 
else{
double lowb, highb;
if(follow == LEFT_LEFT) {highb = scaled[re-1]; lowb = scaled[ls];}
else{highb = scaled[le-1]; lowb = scaled[rs];}
	
assert(highb > lowb);
 	      
//Propose a new breakpoint from the ancestral region in the changed graph 
if(changeBkpoint || scaled[start] != lowb || scaled[end-1] != highb){
bkpoint    = gsl_rng_uniform_pos(gslr)*(scaled[end-1] - scaled[start]) + scaled[start];
propRatio *= (scaled[end-1] - scaled[start])/(highb - lowb);
}

if(changeFollow){if(gsl_rng_uniform_pos(gslr) < 0.5){follow = LEFT_LEFT;} else{follow = LEFT_RIGHT;}}        

assert(bkpoint > scaled[start] && bkpoint < scaled[end-1]);
     
chunkAncLocs(vec);
lta = ancs_tot(gamma, leftAncLocs,  MEANTRACT, ls, le, chain_number); col = rho_tot(rho, ls, le);
rta = ancs_tot(gamma, rightAncLocs, MEANTRACT, rs, re, chain_number); cor = rho_tot(rho, rs, re);
}
}
      
	       
else
{
Node* dec1 = leftD; Node* dec2 = rightD; ALLELE* vec1; ALLELE* vec2;
      
if(dec1->getID() == dec2->getID()){vec1 = dec1->leftAncLocs; vec2 = dec1->rightAncLocs;}		
else{
if(!dec1->isRec() || (dec1->isRec() && dec1->leftA->ID == ID)){vec1 = dec1->leftAncLocs;}
else{vec1 = dec1->rightAncLocs;}

if(!dec2->isRec() || (dec2->isRec() && dec2->leftA->ID == ID)){vec2 = dec2->leftAncLocs;}
else{vec2 = dec2->rightAncLocs;}
}
        
assert(vec1 != NULL && vec2 != NULL);

if(mark == -1){return(propRatio);}

unionAncLocs(vec1, vec2);
totancs = ancs_tot(gamma, leftAncLocs, MEANTRACT, ls, le, chain_number); col = rho_tot(rho, ls, le);
}

return propRatio;
}




double Node::updateAncLocs1(double length, gsl_rng* gslr, double MEANTRACT, double* gamma, double* rho){
/*
double propRatio = 1;  
Node* dec = getDec(); ALLELE* vec; int start, end;
if(!dec->isRec() || (dec->isRec() && dec->leftA->ID == ID)){vec = dec->leftAncLocs; start = dec->ls; end = dec->le;}
else{vec = dec->rightAncLocs; start = dec->rs; end = dec->re;}
  
if(end - start < 2){return 0;}
  
assert(tractlength > 0);
double tal = totancs, pp = prob, ld = lendist;
newtract(vec, length, gslr, MEANTRACT, start, end);
propRatio *= exp(log(pp) + log(totancs) + log(ld) - log(prob) - log(tal) - log(lendist));  
assert(tractlength > 0 && tal == totancs && pp == prob);
  
convAncLocs(vec);
lta = ancs_tot(gamma, leftAncLocs,  MEANTRACT, ls, le); col = rho_tot(rho, ls, le);
rta = ancs_tot(gamma, rightAncLocs, MEANTRACT, rs, re); cor = rho_tot(rho, rs, re);
  
return propRatio;    
*/
}




//------------------------End of functions for Node-------------------------------
//--------------------------------------------------------------------------------
Graph::Graph(int numtips, int numloci, gsl_rng* gr): 
numTips(numtips), numLoci(numloci), gslr(gr)
{ 
  Root = NULL; numRecs = 0; numInternNodes = 0; TotalTime = 0; numCO = 0; numGC = 0;
}




Graph::Graph(Graph& X)
{
   copyGraph(*this, X);
}




void Graph::copyGraph(Graph& newg, Graph& oldg)
{
   if(newg.NodesVector.size()  !=  0) newg.destroy();
   newg.numTips = oldg.numTips;
   newg.numLoci = oldg.numLoci;
   newg.TotalTime = oldg.TotalTime;
   newg.numInternNodes = oldg.numInternNodes;
   newg.numRecs = oldg.numRecs; newg.numCO = oldg.numCO; newg.numGC = oldg.numGC;
   newg.gslr = oldg.gslr;

   sort(oldg.NodesVector.begin(), oldg.NodesVector.end(), Util::compareTwoNodesID);

   vector<Node*>::iterator vi;
   for(vi=oldg.NodesVector.begin(); vi!=oldg.NodesVector.end(); vi++)
   {
   Node* anode;anode = new Node(*(*vi));  
   newg.NodesVector.push_back(anode);
   }

   for(int i=0; i<numTips+numInternNodes; i++)
   {
   if(oldg.NodesVector[i]->getleftA()!=NULL)
   newg.NodesVector[i]->setleftA(newg.NodesVector[oldg.NodesVector[i]->getleftA()->getID()]);
   if(oldg.NodesVector[i]->getrightA()!=NULL)
   newg.NodesVector[i]->setrightA(newg.NodesVector[oldg.NodesVector[i]->getrightA()->getID()]);
   if(oldg.NodesVector[i]->getleftD()!=NULL)
   newg.NodesVector[i]->setleftD(newg.NodesVector[oldg.NodesVector[i]->getleftD()->getID()]);
   if(oldg.NodesVector[i]->getrightD()!=NULL)
   newg.NodesVector[i]->setrightD(newg.NodesVector[oldg.NodesVector[i]->getrightD()->getID()]);
   }
   newg.Root = newg.NodesVector[oldg.Root->getID()];

   sort(newg.NodesVector.begin(), newg.NodesVector.end(), Util::compareTwoNodesT);
   sort(oldg.NodesVector.begin(), oldg.NodesVector.end(), Util::compareTwoNodesT);
}




Graph::~Graph()
{
   destroy(); 
}




Graph& Graph::operator=(Graph& X)
{
   destroy();
   copyGraph(*this, X);
   return *this;
}




void Graph::destroy()
{
   vector<Node*>::iterator it;
   for(it=NodesVector.begin(); it!=NodesVector.end(); it++) delete (*it);
   NodesVector.clear();
}




double Graph::getTMRCA(){return Root->getcoalTime();}




//Create the initial nodes corresponding to the individuals in the sample
void Graph::createTipNodes(double* scaleddis, int** seqmatrix, double length, double MEANTRACT, double* gamma, double* rho, int chain_number)
{    
    int i;  
    if(nch == 0){
    nch = 1;
    snpcount = numLoci; scaled = new double[snpcount]; pos = new int[snpcount];
    for(i = 0; i<snpcount; i++){scaled[i] = scaleddis[i]; pos[i] = ni(length*scaled[i]);}
    }


    //Initialize the lookup arrays
    lkup[chain_number][0] = 1; cum[chain_number][0] = 1;
    for(i = 1;i<1000000;i++){
    lkup[chain_number][i] = lkup[chain_number][i-1]*(1  - (1.0/MEANTRACT));  
    cum[chain_number][i]  =  cum[chain_number][i-1] + lkup[chain_number][i]; 
    }
    

    for(i=0; i<numTips; i++){
    Node* anode; anode = new Node(i, 0, 0);
    ALLELE* leftanc = new ALLELE[numLoci];
    for(int j=0; j<numLoci; j++){leftanc[j] = seqmatrix[i][j];}
    
    anode->setLeftAncLocs(leftanc);
    anode->totancs =  ancs_tot(gamma, leftanc, MEANTRACT, anode->ls, anode->le, chain_number); 
    anode->col     =  rho_tot(rho, anode->ls, anode->le);
    NodesVector.push_back(anode); 
    }    
}




//Simulates a binary tree
void Graph::simuBNTree(double* scaleddis, int** seqmatrix, double length, double MEANTRACT, double* gamma, double* rho, int chain_number)
{   
   createTipNodes(scaleddis, seqmatrix, length, MEANTRACT, gamma, rho, chain_number);

   TotalTime = 0;int nlins = numTips;
   
   double wtime = 0;  double tc = -1; vector<int> availNodes;
   for(int i = 0; i < numTips; i++) availNodes.push_back(i);
   int id = numTips;

   while( nlins > 1){
    
   //Waiting time to next coalescence event = exponential random variable
   double wtimePreG = wtime;
   tc = -2*log(1 - gsl_rng_uniform_pos(gslr))/(nlins*(nlins-1)); tc += wtime;
   wtime = tc;

   //Total time across all branches   
   TotalTime += nlins*(wtime - wtimePreG);

   //Create a new coalescence node and pick 2 existing nodes as its descendants
   int index; Node* newnode; Node* found; int size;
   newnode = new Node(id, wtime, 0);
   
   //Randomly choose 2 existing nodes as descendants
   for(int i = 0; i < 2; i++)
   {
   size = int(availNodes.size());
   assert(size > 0);
   index = gsl_rng_uniform_int(gslr, size);
   found = NodesVector[availNodes[index]];
   found->setleftA(newnode);
   availNodes.erase(availNodes.begin()+index);
   if(i == 0) newnode->setleftD(found);
   else newnode->setrightD(found);
   }

   //Update ancestry vector
   ALLELE* cur = new ALLELE[numLoci];
   ALLELE* ancs = newnode->getleftD()->getAncLocs(LEFT);
   for(int i=0; i<numLoci; i++) {cur[i] = ancs[i];}
   newnode->setLeftAncLocs(cur);
   newnode->totancs  =  ancs_tot(gamma, cur, MEANTRACT, newnode->ls, newnode->le, chain_number);//Rate of informative gene-conversions
   newnode->col      =  rho_tot(rho, newnode->ls, newnode->le);//Rate of informative cross-overs
      
   //Add the new node to lists
   availNodes.push_back(id);
   NodesVector.push_back(newnode);

   //Update counts
   nlins--; id++;
   }
  
   //Assign root and count the number of internal nodes 
   Root = NodesVector[id-1];
   numInternNodes = id - numTips;
}




//Randomly change the topology of the graph
double Graph::moveNode(double delta, double length, double MEANTRACT, double* gamma, double* rho, int chain_number)
{
   double propRatio = 1;
   int nnodes;
   	
   //Pick the number of nodes to move randomly	
   double u = gsl_rng_uniform_pos(gslr);
   if(u < 0.9) nnodes = 1;
   else if(u < 0.98) nnodes = 1;
   else nnodes = 1;
	
   //Choose a random node other than the tip nodes.
   for(int i=0; i<nnodes; i++)
   {
   int chosen = numTips + gsl_rng_uniform_int(gslr, numInternNodes);
   Node* moved = NodesVector[chosen];

   if(moved->isRec()) {propRatio *= changeRec(moved, delta, length, MEANTRACT, gamma, rho, chain_number);}
   else{propRatio *= changeCoal(moved, delta, length, MEANTRACT, gamma, rho, chain_number);}
   
   if(propRatio == 0){return 0;}
   }
   
   return propRatio;
}




//Move a recombination node to a new location in the graph
double Graph::changeRec(Node* moved, double delta, double length, double MEANTRACT, double* gamma, double* rho, int chain_number)
{
    double logprBackw = 0, logprForw = 0;
    Node* attached; Node* dec; Node* anc;
        
    //Randomly choose one of the 2 ancestors to remain attached to the recombination node
    if(gsl_rng_uniform_pos(gslr) < 0.5) {attached = moved->getleftA(); anc = moved->getrightA(); moved->setrightA(NULL);}    
    else {attached = moved->getrightA(); anc = moved->getleftA(); moved->setleftA(NULL);}
    
    dec = moved->getDec(); moved->setleftD(NULL); moved->setrightD(NULL);
    
    if(dec->getleftA() != NULL && dec->getleftA()->getID() == moved->getID()) dec->setleftA(anc);    
    else dec->setrightA(anc);
    if(anc->getleftD() != NULL && anc->getleftD()->getID() == moved->getID()) anc->setleftD(dec);
    else anc->setrightD(dec);

    //Get the current time of the node to be moved
    int nbr = 0;double oldtime = moved->getTime();
    	
    //Count the number of eligible branches at old time.
    for(int i=0; i<(int)NodesVector.size(); i++)
    {
    if(NodesVector[i]->getTime() < oldtime)
    {
    if(NodesVector[i]->getleftA()  != NULL  &&  NodesVector[i]->getleftA()->getTime()  > oldtime   && (NodesVector[i]->le - NodesVector[i]->ls)  > 1) nbr++;
    if(NodesVector[i]->getrightA() != NULL  &&  NodesVector[i]->getrightA()->getTime() > oldtime   && (NodesVector[i]->re - NodesVector[i]->rs)  > 1) nbr++;
    }    
    }
    
    assert(nbr > 0);

    logprBackw += log(1.0/nbr);
		
    //Propose a new time for the node to be moved
    double newtime;
    if(gsl_rng_uniform_pos(gslr) < 0.1) newtime = gsl_rng_uniform_pos(gslr)*attached->getTime();
    else 
    {  
    newtime = fabs(oldtime + gsl_ran_flat(gslr, -delta, delta));
    if(newtime > attached->getTime()) newtime = oldtime;
    }
    assert(newtime > 0 && newtime < attached->getTime());
   
    
    moved->setrecTime(newtime);
    sort(NodesVector.begin(), NodesVector.end(), Util::compareTwoNodesT);
   

    //Count the number of branches at the new time    
    vector<decAnc> branches;
    for(int i=0; i<(int)NodesVector.size(); i++){
    if(NodesVector[i]->getTime() < newtime){
    if(NodesVector[i]->getleftA()  != NULL   &&  NodesVector[i]->getleftA()->getTime()  > newtime   &&  (NodesVector[i]->le - NodesVector[i]->ls)  > 1)      
    {decAnc abr; abr.dec = NodesVector[i]; abr.anc = NodesVector[i]->getleftA();  branches.push_back(abr);}
    if(NodesVector[i]->getrightA() != NULL   &&  NodesVector[i]->getrightA()->getTime() > newtime   &&  (NodesVector[i]->re - NodesVector[i]->rs)  > 1)
    {decAnc abr; abr.dec = NodesVector[i]; abr.anc = NodesVector[i]->getrightA(); branches.push_back(abr);}
    }
    }
    if((int)branches.size() == 0){return 0;}
    
    //Randomly choose a new branch to insert a new node
    int chosen = gsl_rng_uniform_int(gslr, (int)branches.size());
    Node* chbr  = branches[chosen].anc;
    Node* chbrD = branches[chosen].dec;
    logprForw += log(1.0/(int)branches.size());

    insertNodeBtTwoNodes(chbrD, chbr, moved);

    double min;if(oldtime < newtime){min = oldtime - EPSILON;}else{min = newtime - EPSILON;}
    for(int i = numTips; i<(int)NodesVector.size(); i++){if(NodesVector[i]->getTime() < min){NodesVector[i]->mark = -1;}else{NodesVector[i]->mark = 1;}}

    double propRatio = exp(logprBackw - logprForw);
    for(int i=numTips; i<(int)NodesVector.size(); i++)
    {
    propRatio *= NodesVector[i]->updateAncLocs(length, false, false, gslr, MEANTRACT, gamma, rho, chain_number);
    if(propRatio == 0) return 0;
    }

    return propRatio;
}




//Move a coalescence node
double Graph::changeCoal(Node* moved, double delta, double length, double MEANTRACT, double* gamma, double* rho, int chain_number)
{
    double logprBackw = 0, logprForw = 0;
   
    //Create tmproot because coalescence time of the new node can exceed the root time
    Node* tmproot;
    tmproot = new Node(-1, Root->getTime()+1, 0);
    tmproot->setleftD(Root);
    Root->setleftA(tmproot);

    //Randomly choose one of the 2 descendants to remain attached to the coalescence node
    Node* attached; Node* dec; Node* anc;
    if(gsl_rng_uniform_pos(gslr) < 0.5) {attached = moved->getleftD(); dec = moved->getrightD(); moved->setrightD(NULL);}
    else {attached = moved->getrightD(); dec = moved->getleftD(); moved->setleftD(NULL);}    
        
    //Link the non-attached node to the ancestor of the node to be moved
    anc = moved->getAnc();moved->setleftA(NULL); moved->setrightA(NULL);
    if(dec->getleftA() != NULL && dec->getleftA()->getID() == moved->getID()) dec->setleftA(anc);
    else dec->setrightA(anc);    
    if(anc->getleftD() != NULL && anc->getleftD()->getID() == moved->getID()) anc->setleftD(dec);
    else anc->setrightD(dec);
	
    int nbr = 0;double oldtime = moved->getTime();
    for(int i=numTips; i<(int)NodesVector.size(); i++)
    {
        if(NodesVector[i]->getTime() > oldtime)
        {
        if(NodesVector[i]->getleftD()!=NULL  &&  NodesVector[i]->getleftD()->getTime() < oldtime) nbr++;
        if(NodesVector[i]->getrightD()!=NULL && NodesVector[i]->getrightD()->getTime() < oldtime) nbr++;
        }
    }
    if(tmproot->getDec()->getTime() < oldtime) nbr++;
    logprBackw += log(1.0/nbr);


    double lamda = 1.0/delta;
    double dtPrime = -log(1-gsl_rng_uniform_pos(gslr))/lamda;
    double newtime = attached->getTime() + dtPrime;assert(newtime > attached->getTime());
    double dt = oldtime - attached->getTime();
    moved->setcoalTime(newtime);
    sort(NodesVector.begin(), NodesVector.end(), Util::compareTwoNodesT);
    tmproot->setcoalTime(tmproot->getcoalTime() + newtime); 

    logprBackw += log(lamda) - lamda*dt;
    logprForw += log(lamda) - lamda*dtPrime;

    vector<decAnc> branches;
    for(int i = numTips; i<(int)NodesVector.size(); i++)
    {
        if(NodesVector[i]->getTime()  >  newtime)
        {
        if(NodesVector[i]->getleftD()!=NULL && NodesVector[i]->getleftD()->getTime() < newtime )
        {decAnc abr; abr.dec=NodesVector[i]->getleftD(); abr.anc=NodesVector[i]; branches.push_back(abr);}
        if(NodesVector[i]->getrightD()!=NULL && NodesVector[i]->getrightD()->getTime() < newtime)
        {decAnc abr; abr.dec=NodesVector[i]->getrightD(); abr.anc=NodesVector[i]; branches.push_back(abr);}        
        }
    }

    if(tmproot->getDec()->getTime()  <  newtime)    
    {decAnc abr; abr.dec=tmproot->getDec(); abr.anc=tmproot; branches.push_back(abr);}

    int chosen = gsl_rng_uniform_int(gslr, (int)branches.size());  
    Node* chbr = branches[chosen].anc;
    Node* chbrD = branches[chosen].dec;
    logprForw += log(1.0/(int)branches.size());

    insertNodeBtTwoNodes(chbrD, chbr, moved);
    Root = NodesVector[numInternNodes + numTips - 1]; Root->setleftA(NULL); Root->setrightA(NULL);
    if(tmproot != NULL) {delete tmproot;}

    double propRatio = exp(logprBackw - logprForw);

    double min;if(oldtime < newtime){min = oldtime - EPSILON;}else{min = newtime - EPSILON;}
    for(int i = numTips; i<(int)NodesVector.size(); i++){if(NodesVector[i]->getTime() < min){NodesVector[i]->mark = -1;}else{NodesVector[i]->mark = 1;}}

    //Update the ancestry vectors for the whole graph due to the change in topology
    for(int i=numTips; i<(int)NodesVector.size(); i++)
    {
    propRatio *= NodesVector[i]->updateAncLocs(length, false, false, gslr, MEANTRACT, gamma, rho, chain_number);
    if(propRatio == 0){return 0;}
    }

    return propRatio;
}




void Graph::insertNodeBtTwoNodes( Node* node1, Node* node2, Node* insertNode)
{
   assert(node1 != NULL && node2 != NULL);
  
   Node* parent; Node* child; 
   if(node1->getTime() < node2->getTime()) {child = node1; parent = node2;}
   else {child = node2; parent = node1;}

   if(child->getleftA()!=NULL && child->getleftA()->getID() == parent->getID()){ child->setleftA(insertNode);}
   else {child->setrightA(insertNode);}
   if(insertNode->isCoal())
   {if(insertNode->getleftD() == NULL) insertNode->setleftD(child); else insertNode->setrightD(child);}
   else{ insertNode->setleftD(child);}

   if(parent->getleftD()!=NULL && parent->getleftD()->getID() == child->getID()) { parent->setleftD(insertNode);}
   else { parent->setrightD(insertNode);}
   if(insertNode->isCoal()) { insertNode->setleftA(parent);}
   else{ if(insertNode->getleftA()==NULL) insertNode->setleftA(parent); else insertNode->setrightA(parent);}
}




void Graph::revertwait(){
vector<Node*>::iterator vi;
for(vi=NodesVector.begin()+numTips; vi!=NodesVector.end(); vi++){(*vi)->setTime((*vi)->temptime);}
}




//Get the coalescent prior probability of a graph given the value of rho and gamma
double Graph::getlogprior(double* gamma, double* rho, double* scaleddis, double length)
{ 
  double logpr = 0;
  vector<Node*>::iterator vi;
  vector<decAnc>::iterator ii;

  vector<decAnc> branches;
  for(vi=NodesVector.begin(); vi!=NodesVector.begin()+numTips; vi++)
  {decAnc abr; abr.dec=(*vi); abr.anc=(*vi)->getAnc(); branches.push_back(abr);}

  int lin = numTips;
  for(vi=NodesVector.begin()+numTips; vi!=NodesVector.end(); vi++)
  {	  	  	  		
      double rholen = 0, gammalen = 0, qty1, qty2;
      for(ii=branches.begin(); ii!=branches.end(); ii++){
       
      //Coalescence between the products of recombination 
      if((*ii).dec->isRec() && (*ii).dec->getleftA()->getID()  ==  (*ii).dec->getrightA()->getID())//isloop
      {
      for(int k=0; k<2; k++){
      if(k == 0){qty1 = (*ii).dec->lta; qty2 = (*ii).dec->col;}
      else      {qty1 = (*ii).dec->rta; qty2 = (*ii).dec->cor;}
	
      assert(qty1 >= 0); assert(qty2 >= 0);      		
      gammalen += qty1;  rholen += qty2;	
      }
			   
      ii++; assert((*ii).dec->getID() == (*(ii-1)).dec->getID());    
      }
          		  
      else{
      if(!(*ii).dec->isRec() || (*ii).dec->getleftA()->getID() == (*ii).anc->getID()){
      qty2 = (*ii).dec->col;
      if(!(*ii).dec->isRec()){qty1 = (*ii).dec->totancs;}else{qty1 = (*ii).dec->lta;}	        
      }
      else {qty1 = (*ii).dec->rta; qty2 = (*ii).dec->cor;}
          
      assert(qty1 >= 0); assert(qty2 >= 0);      		
      gammalen += qty1;  rholen += qty2;	
      }
      }

      assert(lin == (int)branches.size());
            
      if(length < 0){
      double newtime = (*(vi - 1))->getTime()  +  (-log(drand48())/(0.5*lin*(lin - 1)  + 0.5*gammalen + 0.5*rholen));
      (*vi)->temptime = (*vi)->getTime();
      (*vi)->setTime(newtime);
      }

	  
      //If a coalescence node
      if((*vi)->isCoal())
      {	 
      double intervTime = ((*vi)->getcoalTime()  -  (*(vi-1))->getTime());
      logpr += -(0.5*gammalen + 0.5*rholen + (double)lin*(lin-1)/(2))*intervTime;

      int nerase = 0;
      for(int i=0; i<2; i++)
      {
      for(ii=branches.begin(); ii!=branches.end(); ii++)
      {if((*ii).anc->getID() == (*vi)->getID()) {branches.erase(ii); nerase++; break;}}
      }
         
      decAnc abr; abr.dec = (*vi); abr.anc = (*vi)->getAnc(); branches.push_back(abr);
      if(nerase != 2){/*printGraph();*/ cout<<(*vi)->getID()<<"\t"<<(*vi)->getTime()<<endl; exit(1);}
      lin--;
      }
	  
	  	  
      //If a recombination node
      else
      {
      assert(gammalen > 0  &&  rholen > 0);
      double intervTime = (*vi)->getrecTime()  -  (*(vi-1))->getTime();
   
      //If a gene-conversion node add the probability of initiation point and tract length
      if((*vi)->gettractlength() > 0){
      logpr += -(0.5*gammalen + 0.5*rholen +  (double)lin*(lin - 1)/(2))*intervTime + log(0.5*gammalen);
      
      //Find the rate at gene-conversion start point
      int i; double bkp = (*vi)->getBkpoint(); 
      for(i=0; i<snpcount; i++){if(scaled[i] - bkp > EPSILON){break;}} i = maximum(i - 1, 0);
      logpr += log((*vi)->prob) + log((*vi)->lendist) + log(gamma[i]*PERMB) - log(gammalen); 
      }
		  		  		  
      else{//If a crossing-over node get the recombination rate for bkp 
      logpr += -(0.5*gammalen + 0.5*rholen + (double)lin*(lin - 1)/(2))*intervTime + log(0.5*rholen);
      double bkp = (*vi)->getBkpoint(); double rholeni = -1; double scaledLen = -1;
      for(int i=0; i<numLoci-1; i++)
      {
      if(bkp > scaleddis[i] && bkp < scaleddis[i+1]) 
      {rholeni = rho[i]*(scaleddis[i+1] - scaleddis[i])*fabs(length)*PERMB; scaledLen = scaleddis[i+1] - scaleddis[i]; break;}
      }
      assert(rholeni > 0 && scaledLen > 0);
      if(rholeni/rholen - 1 > 1e-7){cout<<"rholeni = "<<rholeni<<"\trholen = "<<rholen<<endl;}
      logpr += log(rholeni) - log(rholen) - log(scaledLen);
      }
		  		  		  
      //Update branches 
      for(ii=branches.begin(); ii!=branches.end(); ii++){if((*ii).anc->getID()  ==  (*vi)->getID()) {branches.erase(ii); break;}}
      decAnc abr1; abr1.dec = (*vi); abr1.anc =  (*vi)->getleftA();  branches.push_back(abr1);
      decAnc abr2; abr2.dec = (*vi); abr2.anc = (*vi)->getrightA();  branches.push_back(abr2);
      lin++;
      }

      if(vi < NodesVector.end()-1 && lin == 1) {return 0;}
  }
	
  assert(lin == 1);
  assert((int)branches.size() == 1 && branches[0].dec->getID() == Root->getID());
  assert(Root->le - Root->ls == snpcount);

  return logpr;
}




bool Graph::containsLoop()
{
   for(int i=numTips; i<numTips+numInternNodes; i++){
   if(NodesVector[i]->isCoal() && NodesVector[i]->getleftD()->getID() == NodesVector[i]->getrightD()->getID()) return true;
   }

   return false;
}




void Graph::eligibleRecBranches(vector<decAncdir>& branches, double rectime)
{
  vector<Node*>::iterator vi;
  branches.clear();

  for(vi=NodesVector.begin(); vi!=NodesVector.end()-1; vi++)
  {
      if((*vi)->getTime() < rectime)
      {
           Node* lefta = (*vi)->getleftA();
           if(lefta!=NULL  && lefta->getTime() > rectime   && ((*vi)->le - (*vi)->ls) > 1)
           {
           decAncdir abr; abr.dec = (*vi); abr.ancdir = LEFT;
           branches.push_back(abr);
           }

           Node* righta = (*vi)->getrightA();
           if(righta!=NULL && righta->getTime() > rectime  && ((*vi)->re - (*vi)->rs) > 1)
           {
           decAncdir abr; abr.dec = (*vi); abr.ancdir = RIGHT;
           branches.push_back(abr);
           }
      }
  }
}




void Graph::eligibleCoalBranches(vector<decAncdir>& branches, double coaltime)
{
  vector<Node*>::iterator vi;
  branches.clear();

  for(vi=NodesVector.begin(); vi!=NodesVector.end(); vi++)
  {
      if((*vi)->getTime()<coaltime)
      {
           Node* lefta = (*vi)->getleftA();
           if(lefta!=NULL && lefta->getTime() > coaltime )
           {
                decAncdir abr; abr.dec=(*vi); abr.ancdir=LEFT;
                branches.push_back(abr);
           }
           Node* righta = (*vi)->getrightA();
           if(righta!=NULL && righta->getTime() > coaltime)
           {
                decAncdir abr; abr.dec=(*vi); abr.ancdir=RIGHT;
                branches.push_back(abr);
           }
      }
  }
}




void Graph::eligiblePairRec(vector<decAncdir>& branches)
{
        vector<Node*>::iterator vi;
        branches.clear();assert(branches.size()==0);

        for(vi=NodesVector.begin()+numTips; vi!=NodesVector.end(); vi++){
        if((*vi)->isRec() && (*vi)->gettractlength() == -1){
        Node* lefta = (*vi)->getleftA();
        if(lefta->isCoal() && lefta->getID()  != Root->getID())
        {
        decAncdir abr; abr.dec=(*vi); abr.ancdir = LEFT;
        branches.push_back(abr);
        }

        Node* righta = (*vi)->getrightA();
        if(righta->isCoal() && righta->getID() != Root->getID())
        {
        decAncdir abr; abr.dec=(*vi); abr.ancdir = RIGHT;
        branches.push_back(abr);
        }
        }        
        }
}




void Graph::eligiblePairConv(vector<decAncdir>& branches)
{
	vector<Node*>::iterator vi;
	branches.clear();assert(branches.size() == 0);
	
	for(vi=NodesVector.begin()+numTips; vi!=NodesVector.end(); vi++)
	{
        if((*vi)->isRec() && (*vi)->gettractlength() > 0)
        {
	Node* lefta = (*vi)->getleftA();
	if(lefta->isCoal() && lefta->getID()   != Root->getID())
        {
        decAncdir abr; abr.dec=(*vi); abr.ancdir = LEFT;
        branches.push_back(abr);
	}
	Node* righta = (*vi)->getrightA();
	if(righta->isCoal() && righta->getID() != Root->getID() )
	{
        decAncdir abr; abr.dec=(*vi); abr.ancdir=RIGHT;
        branches.push_back(abr);
	}
	}
	}
}




void Graph::insertNode( Node* insertNode, Node* decNode, int ancDir)
{
   Node* parent; Node* child;
   child = decNode;
   if(ancDir == LEFT){parent = child->getleftA();}
   else{parent = child->getrightA();}

   if(ancDir == LEFT){child->setleftA(insertNode);}
   else {child->setrightA(insertNode);}
   if(insertNode->isCoal()){
   if(insertNode->getleftD()==NULL ) insertNode->setleftD(child); else insertNode->setrightD(child);
   }
   else{ insertNode->setleftD(child);}

   if(parent->getleftD()!=NULL && parent->getleftD()->getID() == child->getID()) { parent->setleftD(insertNode);}
   else { parent->setrightD(insertNode);}
   if(insertNode->isCoal()) {insertNode->setleftA(parent);}
   else{if(ancDir == LEFT) insertNode->setleftA(parent); else insertNode->setrightA(parent);}
}




//Propose new breakpoint and/ tractlength
double Graph::changeBkpointFollow(double length, double MEANTRACT, double* gamma, double* rho, int chain_number)
{
   double propRatio = 1;
   int chosen = gsl_rng_uniform_int(gslr, numRecs); //Pick a recombination node at random
   int count = 0; int signal = 0;
   int v = int(NodesVector.size());   

   for(int i=numTips; i<v; i++){NodesVector[i]->mark = -1;}

   for(int i=numTips; i<v; i++){
   if(NodesVector[i]->isRec() && signal == 0){
   if(count == chosen){ 
   
   NodesVector[i]->mark = 1;
   for(int j=numTips; j<v; j++){
   if(NodesVector[j]->getleftD()   !=  NULL   &&    NodesVector[j]->getleftD()->mark  == 1){NodesVector[j]->mark = 1;}
   if(NodesVector[j]->getrightD()  !=  NULL   &&    NodesVector[j]->getrightD()->mark == 1){NodesVector[j]->mark = 1;}
   }

   propRatio *= NodesVector[i]->updateAncLocs(length, true, true, gslr, MEANTRACT, gamma, rho, chain_number);
   signal = 1;
   }   
   count++;
   }
    
   if(signal == 1){propRatio *= NodesVector[i]->updateAncLocs(length, false, false, gslr, MEANTRACT, gamma, rho, chain_number);}
   if(propRatio == 0){return 0;}
   }

   return propRatio;
}




//Propose new tractlength alone                                                                                                                                                
double Graph::changetractlength(double length, double MEANTRACT, double* gamma, double* rho)
{ 
  /*

  int i, j, times;
  int v = int(NodesVector.size()); 

  double propRatio = 1;
  for(times = 0; times < 2; times++){

  int chosen = gsl_rng_uniform_int(gslr, numGC);
  int count = 0; 
  for(i=numTips; i<v; i++){  
  if(NodesVector[i]->gettractlength() > 0){
  if(count == chosen){propRatio *= NodesVector[i]->updateAncLocs1(length, gslr, MEANTRACT, gamma, rho);break;}
  count++;
  }
  }

  if(propRatio == 0){return 0;}  

  for(j = i+1; j < v; j++){
  propRatio *= NodesVector[j]->updateAncLocs(length, false, false, gslr, MEANTRACT, gamma, rho);
  if(propRatio == 0){return 0;}
  }

  }

  return propRatio;
  */
}




void Graph::changeOneAncSeq()
{
	vector<Node*>** trees;
	
	trees = new vector<Node*>*[numLoci];
	
	for(int i=0; i<numLoci; i++){ trees[i] = new vector<Node*>; }
	
	getTreesFast(trees);
	
	int nnodes = trees[0]->size();
	assert(nnodes = 3*(numTips-1));
	
	for(int mk=0; mk<numLoci; mk++)
	{
	    if(gsl_rng_uniform_pos(gslr)<0.05) //change an allele of the marker tree
	    {
	    vector<Node*>* markerTree = trees[mk];
			
	    //choose an internal node
	    int ch = gsl_rng_uniform_int(gslr, nnodes); 
	    Node* chosen = *(markerTree->begin()+ch);
			
	                if(!chosen->isTip())
	                {
			ALLELE* hap = chosen->getAncLocs(LEFT);   
			hap[mk] = ranchar2(hap[mk]);
			}
	     }   
	}
	

	for(int i=0; i<numLoci; i++) delete trees[i];
	delete [] trees;
}




int Graph::ranchar()
{
   double p = gsl_rng_uniform_pos(gslr);
   if(p < 0.25) return'A';
   else if(p >= 0.25  && p < 0.50) return 'T';
   else if(p >= 0.50  && p < 0.75) return 'G';
   else return 'C';
}




int Graph::ranchar2(int allele)
{
   double p = gsl_rng_uniform_pos(gslr);

   if(allele == 'A') 
   {
      if(p<1.0/3.0) return 'C'; else if(p<2.0/3.0) return 'G'; else return 'T';
   }
   else if(allele == 'C') 
   {
      if(p<1.0/3.0) return 'A'; else if(p<2.0/3.0) return 'G'; else return 'T';
   }
   else if(allele == 'G') 
   {
      if(p<1.0/3.0) return 'A'; else if(p<2.0/3.0) return 'C'; else return 'T';
   }
   else if(allele == 'T') 
   {
      if(p<1.0/3.0) return 'A'; else if(p<2.0/3.0) return 'C'; else return 'G';
   }
   else {cout<<"allele doesn't exist!"<<endl; exit(1);}
}




int Graph::ranchar3(int allele1, int allele2)
{
   assert(allele1=='A' || allele1=='C' ||allele1=='G'||allele1=='T');
   assert(allele2=='A' || allele2=='C' ||allele2=='G'||allele2=='T');

   double p = gsl_rng_uniform_pos(gslr);

   if((allele1=='A' && allele2=='C') || (allele2=='A' && allele1=='C')) {if(p<0.5) return 'G'; else return 'T';}
   else if((allele1=='A' && allele2=='G') || (allele2=='A' && allele1=='G')) {if(p<0.5) return 'C'; else return 'T';}
   else if((allele1=='A' && allele2=='T') || (allele2=='A' && allele1=='T')) {if(p<0.5) return 'C'; else return 'G';}
   else if((allele1=='C' && allele2=='G') || (allele2=='C' && allele1=='G')) {if(p<0.5) return 'A'; else return 'T';}
   else if((allele1=='C' && allele2=='T') || (allele2=='C' && allele1=='T')) {if(p<0.5) return 'A'; else return 'G';}
   else if((allele1=='G' && allele2=='T') || (allele2=='G' && allele1=='T')) {if(p<0.5) return 'A'; else return 'C';}
   else {cout<<"allele doesn't exist!"<<endl; exit(1);}
}




//Add a pair of crossing-over and coalescence node
double Graph::addRecCoal(double length, double MEANTRACT, double* gamma, double* rho, int chain_number)
{
  double logprForw  = 0;
  double logprBackw = 0;
  vector<Node*>::iterator vi;

  double rectime, coaltime;
  double t1 = gsl_rng_uniform_pos(gslr)*Root->getTime();
  double t2 = gsl_rng_uniform_pos(gslr)*Root->getTime();
  
  while (fabs(t1 - t2) < EPSILON){
  t1 = gsl_rng_uniform_pos(gslr)*Root->getTime();
  t2 = gsl_rng_uniform_pos(gslr)*Root->getTime();
  }

  if(t1 < t2) {rectime = t1; coaltime = t2;}
  else {rectime = t2; coaltime = t1;}

  //Proposal probability for time
  logprForw += log(pow(1.0/Root->getTime(), 2)*2);

  //Gets a vector of eligible branches at rectime
  vector<decAncdir> branches;
  eligibleRecBranches(branches, rectime);
  if(branches.empty()) return 0;  

  int nbr = (int)branches.size();
  int ch1 = gsl_rng_uniform_int(gslr, nbr); 
  Node* brec = branches[ch1].dec;  int brecAdir = branches[ch1].ancdir;
  Node* brecA; if(brecAdir == LEFT) brecA = brec->getleftA(); else brecA = brec->getrightA();
  logprForw += log(1.0/nbr);

  eligibleCoalBranches(branches, coaltime);
  int nbr2 = (int)branches.size(); assert(nbr2 >= 1);
  int ch2 = gsl_rng_uniform_int(gslr, nbr2);
  Node* bcoal = branches[ch2].dec; int bcoalAdir = branches[ch2].ancdir;
  Node* bcoalA; if(bcoalAdir == LEFT) bcoalA=bcoal->getleftA(); else bcoalA = bcoal->getrightA();
  logprForw += log(1.0/nbr2);

  assert(rectime  < coaltime);
  assert(rectime  > brec->getTime()  && rectime  < brecA->getTime());
  assert(coaltime > bcoal->getTime() && coaltime < bcoalA->getTime());

  Node* recnode; Node* coalnode;
  recnode  = new Node(numInternNodes+numTips, 0, rectime);
  coalnode = new Node(numInternNodes+numTips+1, coaltime, 0);
  
  insertNode(recnode, brec, brecAdir);
  if(bcoal->getID() == brec->getID() && bcoalA->getID() == brecA->getID())//loop
  {
  if(recnode->getleftA() != NULL){insertNode(coalnode, recnode, LEFT);}
  else{insertNode(coalnode, recnode, RIGHT);}
  }
  else{insertNode(coalnode, bcoal, bcoalAdir);}

  if(recnode->getleftA()  == NULL){recnode->setleftA(coalnode);}
  else{recnode->setrightA(coalnode);}
  if(coalnode->getleftD() == NULL){coalnode->setleftD(recnode);}
  else{coalnode->setrightD(recnode);}

  NodesVector.push_back(recnode);
  NodesVector.push_back(coalnode);
  sort(NodesVector.begin(), NodesVector.end(), Util::compareTwoNodesT);
  numInternNodes += 2; numRecs++; numCO++;

  ALLELE* vec; int start, end;
  if(!brec->isRec() || brec->getleftA()->getID() == recnode->getID()){vec = brec->getAncLocs(LEFT); start = brec->ls; end = brec->le;}
  else {vec = brec->getAncLocs(RIGHT); start = brec->rs; end = brec->re;}
  assert(end - start  >  1);
  recnode->setBkpoint(gsl_rng_uniform_pos(gslr)*(scaled[end-1] - scaled[start])  +  scaled[start]);
  assert(recnode->getBkpoint() > scaled[start] && recnode->getBkpoint() < scaled[end-1]);
  
  //Sets the direction of follow	
  if(gsl_rng_uniform_pos(gslr) < 0.5){recnode->setFollow(LEFT_LEFT);} else{recnode->setFollow(LEFT_RIGHT);}

  //Proposal term for follow and breakpoint
  logprForw += log(0.5) - log(scaled[end-1] - scaled[start]);
  recnode->chunkAncLocs(vec);
  recnode->lta = ancs_tot(gamma, recnode->getLeftAncLocs(),  MEANTRACT, recnode->ls, recnode->le, chain_number); recnode->col = rho_tot(rho, recnode->ls, recnode->le);
  recnode->rta = ancs_tot(gamma, recnode->getRightAncLocs(), MEANTRACT, recnode->rs, recnode->re, chain_number); recnode->cor = rho_tot(rho, recnode->rs, recnode->re);
  
  //Mark nodes ancestral to modified recombination node  
  for(int i = numTips; i<(int)NodesVector.size(); i++){NodesVector[i]->mark = -1;}
  recnode->mark = 1;
  for(int i=numTips; i<(int)NodesVector.size(); i++){
  if(NodesVector[i]->getleftD()!=NULL  &&  NodesVector[i]->getleftD()->mark  == 1){NodesVector[i]->mark = 1;}
  if(NodesVector[i]->getrightD()!=NULL &&  NodesVector[i]->getrightD()->mark == 1){NodesVector[i]->mark = 1;}
  }
  
  double propRatio = 1;
  for(int i=numTips; i<(int)NodesVector.size(); i++)
  {
  if(NodesVector[i]->getTime() > rectime)
  {
  propRatio *= NodesVector[i]->updateAncLocs(length, false, false, gslr, MEANTRACT, gamma, rho, chain_number);
  if(propRatio == 0) {return 0;}
  }
  }

  eligiblePairRec(branches);
  assert(!branches.empty());
  logprBackw += log(1.0/(double)branches.size());
  
  return exp(logprBackw - logprForw)*propRatio;
}




//Add a pair of gene-conversion and coalescence nodes 
double Graph::addconvCoal(double length, double MEANTRACT, double* gamma, double* rho, int chain_number)
{
  double logprForw  = 0;
  double logprBackw = 0;
  vector<Node*>::iterator vi;
	
  double rectime, coaltime;
  double t1 = gsl_rng_uniform_pos(gslr)*Root->getTime();
  double t2 = gsl_rng_uniform_pos(gslr)*Root->getTime();
  
  while (fabs(t1 - t2) < EPSILON){
  t1 = gsl_rng_uniform_pos(gslr)*Root->getTime();
  t2 = gsl_rng_uniform_pos(gslr)*Root->getTime();
  }

  if(t1 < t2) {rectime = t1; coaltime = t2;}
  else {rectime = t2; coaltime = t1;}
  logprForw += log(pow(1.0/Root->getTime(), 2)*2);
	
  //Get the vector of eligible branches at rectime	
  vector<decAncdir> branches;
  eligibleRecBranches(branches, rectime);
  if(branches.empty()) return 0;  
	
  int nbr = (int)branches.size();
  int ch1 = gsl_rng_uniform_int(gslr, nbr); 
  Node* brec = branches[ch1].dec;  int brecAdir = branches[ch1].ancdir;
  Node* brecA; if(brecAdir == LEFT) brecA = brec->getleftA(); else brecA = brec->getrightA();
  logprForw += log(1.0/nbr);
		
  eligibleCoalBranches(branches, coaltime);
  int nbr2 = (int)branches.size(); assert(nbr2 >= 1);
  int ch2 = gsl_rng_uniform_int(gslr, nbr2);
  Node* bcoal = branches[ch2].dec; int bcoalAdir = branches[ch2].ancdir;
  Node* bcoalA; if(bcoalAdir == LEFT) bcoalA = bcoal->getleftA(); else bcoalA = bcoal->getrightA();
  logprForw += log(1.0/nbr2);
	
  assert(rectime  < coaltime);
  assert(rectime  > brec->getTime()  && rectime<brecA->getTime());
  assert(coaltime > bcoal->getTime() && coaltime<bcoalA->getTime());
	
	
  Node* recnode; Node* coalnode;
  recnode = new Node(numInternNodes+numTips,    0, rectime);
  coalnode= new Node(numInternNodes+numTips+1, coaltime, 0);
  	
  insertNode(recnode, brec, brecAdir);
  if(bcoal->getID()==brec->getID() && bcoalA->getID()==brecA->getID())//loop
  {
  if(recnode->getleftA()!=NULL){insertNode(coalnode, recnode, LEFT);}
  else{insertNode(coalnode, recnode, RIGHT);}
  }
  else{insertNode(coalnode, bcoal, bcoalAdir);}
	
		
  if(recnode->getleftA()  == NULL) recnode->setleftA(coalnode);
  else recnode->setrightA(coalnode);
  if(coalnode->getleftD() == NULL) coalnode->setleftD(recnode);
  else coalnode->setrightD(recnode);
	
	
  NodesVector.push_back(recnode);
  NodesVector.push_back(coalnode);
  sort(NodesVector.begin(), NodesVector.end(), Util::compareTwoNodesT);
  numInternNodes += 2; numRecs++; numGC++;
		
  ALLELE* vec; int start, end;
  if(!brec->isRec() || brec->getleftA()->getID() == recnode->getID()){vec = brec->getAncLocs(LEFT); start = brec->ls; end = brec->le;}
  else {vec = brec->getAncLocs(RIGHT); start = brec->rs; end = brec->re;}
  assert(end - start  >  1);
		
  //Choose the start of conversion tract and tract length such that tract always includes atleast one marker position
  recnode->setconv(vec, length, gslr, MEANTRACT, start, end, chain_number);
	
  //Randomly distribute the 2 products of conversion
  if(gsl_rng_uniform_pos(gslr) < 0.5){recnode->setFollow(LEFT_LEFT);} else{recnode->setFollow(LEFT_RIGHT);}	
  
  //Probability of proposing start point, tract length and follow in the added gene-conversion node
  logprForw += log(recnode->prob) + log(recnode->lendist) - log(recnode->totancs) + log(0.5);
	
  recnode->convAncLocs(vec);
  recnode->lta = ancs_tot(gamma, recnode->getLeftAncLocs(),  MEANTRACT, recnode->ls, recnode->le, chain_number); recnode->col = rho_tot(rho, recnode->ls, recnode->le);
  recnode->rta = ancs_tot(gamma, recnode->getRightAncLocs(), MEANTRACT, recnode->rs, recnode->re, chain_number); recnode->cor = rho_tot(rho, recnode->rs, recnode->re);
  
  //Mark all nodes ancestral to the added gene-conversion node 
  for(int i=numTips; i<(int)NodesVector.size(); i++){NodesVector[i]->mark = -1;}
  recnode->mark = 1;
  for(int i=numTips; i<(int)NodesVector.size(); i++){
  if( NodesVector[i]->getleftD()!=NULL &&   NodesVector[i]->getleftD()->mark == 1){NodesVector[i]->mark = 1;}
  if(NodesVector[i]->getrightD()!=NULL &&  NodesVector[i]->getrightD()->mark == 1){NodesVector[i]->mark = 1;}
  }
	
  double propRatio = 1;
  for(int i=numTips; i<(int)NodesVector.size(); i++){
  if(NodesVector[i]->getTime() > rectime){
  propRatio *= NodesVector[i]->updateAncLocs(length, false, false, gslr, MEANTRACT, gamma, rho, chain_number);
  if(propRatio == 0) {return 0;}
  }
  }
	
  eligiblePairConv(branches);
  assert(!branches.empty());
  logprBackw += log(1.0/(double)branches.size());
        	
  return exp(logprBackw - logprForw)*propRatio;
}
 



//Delete a pair of crossing-over and coalescence nodes
double Graph::deleteRecCoal(double length, double MEANTRACT, double* gamma, double* rho, int chain_number)
{
  double logprForw = 0, logprBackw = 0;
  vector<Node*>::iterator vi;

  vector<decAncdir> branches;
  eligiblePairRec(branches);
  if(branches.empty()) return 0;

  int chosen = gsl_rng_uniform_int(gslr, (int)branches.size());
  logprForw += log(1.0/(double)branches.size());
  Node* rec = branches[chosen].dec;
  
  //Mark all the ancestors of the recombination node that is deleted  
  for(int i=numTips; i<(int)NodesVector.size(); i++){NodesVector[i]->mark = -1;}
  rec->mark = 1;
  for(int i=numTips; i<(int)NodesVector.size(); i++){
  if(NodesVector[i]->getleftD()!=NULL &&  NodesVector[i]->getleftD()->mark == 1){NodesVector[i]->mark = 1;}
  if(NodesVector[i]->getrightD()!=NULL && NodesVector[i]->getrightD()->mark == 1){NodesVector[i]->mark = 1;}
  }
  
  Node* coal;
  int coalDir = branches[chosen].ancdir;
  if(coalDir == LEFT) coal = rec->getleftA(); else coal = rec->getrightA();
  assert(rec->isRec() && coal->isCoal() && rec->gettractlength()== -1);

  double rectime = rec->getTime();
  double coaltime = coal->getTime();
	
  if(rec->getleftA()->getID() ==  coal->getID()) rec->setleftA(NULL);
  else rec->setrightA(NULL);
  if(coal->getleftD()->getID() == rec->getID()) coal->setleftD(NULL);
  else coal->setrightD(NULL);
	
  Node* brec; Node* brecA; Node* bcoal; Node* bcoalA;
  if(rec->getAnc()->getID() != coal->getID()) //not a loop
  {brec = rec->getDec(); brecA = rec->getAnc(); bcoal = coal->getDec(); bcoalA = coal->getAnc();}
  else{brec = rec->getDec(); brecA = coal->getAnc(); bcoal = brec; bcoalA = brecA;}

  ALLELE* vec;int start, end;
  if(!brec->isRec() || brec->getleftA()->getID()  ==  rec->getID()){vec = brec->getAncLocs(LEFT); start = brec->ls; end = brec->le;}
  else{vec = brec->getAncLocs(RIGHT); start = brec->rs; end = brec->re;}
  logprBackw += log(0.5) - log(scaled[end-1]  -  scaled[start]);
  

  //Update ancestry details
  Node* anc; Node* dec; Node* cursor;
  for(int i=0; i<2; i++)
  {
     if(i == 0){cursor =  rec;}
     else{cursor = coal;}

     if(cursor->getAnc()!=NULL)
     {
     anc = cursor->getAnc(); dec=cursor->getDec();
     if(anc->getleftD()!=NULL && anc->getleftD()->getID()==cursor->getID()) anc->setleftD(dec);
     else anc->setrightD(dec);
     if(dec->getleftA()!=NULL && dec->getleftA()->getID()==cursor->getID()) dec->setleftA(anc);
     else dec->setrightA(anc);
     }
     numInternNodes --;
  }


  int nerase=0;
  for(vi=NodesVector.begin(); vi!=NodesVector.end(); vi++)
  {if((*vi)->getID()==rec->getID())  {NodesVector.erase(vi); nerase++; break;}}
  for(vi=NodesVector.begin(); vi!=NodesVector.end(); vi++)
  {if((*vi)->getID()==coal->getID())  {NodesVector.erase(vi); nerase++; break;}}
  
  assert(nerase == 2);

  numRecs--; numCO--;

  delete rec; delete coal;

  sort(NodesVector.begin(), NodesVector.end(), Util::compareTwoNodesT);
  for(int i=numTips; i<numInternNodes+numTips; i++){NodesVector[i]->setID(i);}

   
  //Update the ancestry vectors of all nodes after the added node
  double propRatio = 1;
  for(int i=numTips; i<(int)NodesVector.size(); i++)
  {
  if(NodesVector[i]->getTime() > rectime)
  {
  propRatio *= NodesVector[i]->updateAncLocs(length, false, false, gslr, MEANTRACT, gamma, rho, chain_number);
  if(propRatio == 0) return 0;
  }
  }

  eligibleRecBranches(branches, rectime);
  int nbr = (int)branches.size();
  assert(nbr > 0);
  logprBackw += log(1.0/nbr);

  eligibleCoalBranches(branches, coaltime);
  nbr = (int)branches.size();
  assert(nbr > 0);
  logprBackw += log(1.0/nbr);

  logprBackw += log(pow(1.0/Root->getTime(), 2)*2);
  
  return exp(logprBackw - logprForw)*propRatio;
}




//Delete a pair of gene-conversion and coalescence nodes
double Graph::deleteconvCoal(double length, double MEANTRACT, double* gamma, double* rho, int chain_number)
{	
    double logprForw = 0, logprBackw = 0;
    vector<Node*>::iterator vi;
	
    vector<decAncdir> branches;
    eligiblePairConv(branches);
    if(branches.empty()) return 0;
	
    int chosen = gsl_rng_uniform_int(gslr, (int)branches.size());
    logprForw += log(1.0/(double)branches.size());
	
    Node* rec = branches[chosen].dec;

    //Mark all the ancestors of the conversion node to be deleted  
    for(int i=numTips; i<(int)NodesVector.size(); i++){NodesVector[i]->mark = -1;}
    rec->mark = 1;
    for(int i=numTips; i<(int)NodesVector.size(); i++){
    if(NodesVector[i]->getleftD()!=NULL  &&  NodesVector[i]->getleftD()->mark == 1){NodesVector[i]->mark = 1;}
    if(NodesVector[i]->getrightD()!=NULL && NodesVector[i]->getrightD()->mark == 1){NodesVector[i]->mark = 1;}
    }
    
    Node* coal;
    int coalDir = branches[chosen].ancdir;
    if(coalDir == LEFT) coal = rec->getleftA(); else coal = rec->getrightA();
    assert(rec->isRec() && coal->isCoal() && rec->gettractlength() > 0);
	
    double rectime = rec->getTime();
    double coaltime = coal->getTime();
	
    if(rec->getleftA()->getID()  == coal->getID()) rec->setleftA(NULL);
    else rec->setrightA(NULL);
    if(coal->getleftD()->getID() == rec->getID()) coal->setleftD(NULL);
    else coal->setrightD(NULL);
	
		
    Node* brec; Node* brecA; Node* bcoal; Node* bcoalA;
    if(rec->getAnc()->getID()!=coal->getID())//not a loop
    {brec = rec->getDec(); brecA = rec->getAnc(); bcoal = coal->getDec(); bcoalA = coal->getAnc();}
    else{brec = rec->getDec(); brecA = coal->getAnc(); bcoal = brec; bcoalA = brecA;}
	
    //ALLELE* vec;
    //if(!brec->isRec() || brec->getleftA()->getID() == rec->getID()) vec = brec->getAncLocs(LEFT);
    //else vec = brec->getAncLocs(RIGHT);
    logprBackw += log(rec->prob) + log(rec->lendist) - log(rec->totancs) + log(0.5);
	
    //Change ancestry details	
    Node* anc; Node* dec; Node* cursor;
    for(int i=0; i<2; i++)
    {
    if(i == 0){cursor = rec;}
    else{cursor = coal;}
		
    if(cursor->getAnc()  !=  NULL){
    anc = cursor->getAnc(); dec = cursor->getDec();			
    if(anc->getleftD() != NULL && anc->getleftD()->getID() == cursor->getID()){anc->setleftD(dec);}
    else{anc->setrightD(dec);}
    if(dec->getleftA() != NULL && dec->getleftA()->getID() == cursor->getID()){dec->setleftA(anc);}
    else{dec->setrightA(anc);}
    }
		
    numInternNodes --;
    }
	

    int nerase=0;
    for(vi=NodesVector.begin(); vi!=NodesVector.end(); vi++)
    {if((*vi)->getID() == rec->getID())   {NodesVector.erase(vi); nerase++; break;}}
    for(vi=NodesVector.begin(); vi!=NodesVector.end(); vi++)
    {if((*vi)->getID() == coal->getID())  {NodesVector.erase(vi); nerase++; break;}}
    assert(nerase  ==  2);
	
    numRecs--; numGC--;
	
    delete rec; delete coal;
	
    sort(NodesVector.begin(), NodesVector.end(), Util::compareTwoNodesT);
    for(int i=numTips; i<numInternNodes+numTips; i++){NodesVector[i]->setID(i);}
	
    double propRatio = 1;
    for(int i=numTips; i<(int)NodesVector.size(); i++){
    if(NodesVector[i]->getTime() > rectime){
    propRatio *= NodesVector[i]->updateAncLocs(length, false, false, gslr, MEANTRACT, gamma, rho, chain_number);
    if(propRatio == 0) return 0;
    }
    }
	
    eligibleRecBranches(branches, rectime);
    int nbr = (int)branches.size(); assert(nbr > 0);
    logprBackw += log(1.0/nbr);
		
    eligibleCoalBranches(branches, coaltime);
    nbr = (int)branches.size(); assert(nbr > 0);
    logprBackw += log(1.0/nbr);
	
    logprBackw += log(pow(1.0/Root->getTime(), 2)*2);
        	
    return exp(logprBackw - logprForw)*propRatio;
}




double Graph::getTottimeGraph()
{
    double tottime = 0;
    int nlins = numTips;
    for(int i=numTips; i<(int)NodesVector.size(); i++)
    {
    tottime += nlins*(NodesVector[i]->getTime() - NodesVector[i-1]->getTime());
    if(NodesVector[i]->isRec()) nlins++;
    else nlins--;
    }

    assert(nlins == 1);
    return tottime;
}




Node* Graph::findNode(int id)
{
    for(int i=0; i<(int)NodesVector.size(); i++)
    {
    if(NodesVector[i]->getID() == id){return NodesVector[i];}
    }
    return NULL;
}




double Graph::logllh(double theta, double mkloc, double* staFreq)
{
    double logllh = 0;
    if(mkloc < 0){logllh = logllhGraph(theta, staFreq);}
    else{logllh = logllhGraphSingleSite(theta, staFreq, int(mkloc));}

    return logllh;
}




double Graph::logllhGraphSingleSite(double theta, double* staFreq, int mk)
{
  double logllh = 0;
  vector<Node*>::iterator vi;
  double brlen = 0;
  
  int nlineages = numTips;

  for(vi=NodesVector.begin()+numTips; vi!=NodesVector.end(); vi++)
    {
      Node* cursor = (*vi);
      if((*vi)->isRec() || ((*vi)->getleftD()->getID()==(*vi)->getrightD()->getID())) continue;

      ALLELE* hap = cursor->getAncLocs(LEFT);

      if(hap[mk]==999) continue;

      ALLELE* hapdLeft = transHap(cursor->getleftD(), cursor);
      ALLELE* hapdRight = transHap(cursor->getrightD(), cursor);

      if(hapdLeft[mk]==999 || hapdRight[mk]==999) continue;

      Node* decTrees[2];
      decTrees[0] = cursor->getleftD();
      decTrees[1] = cursor->getrightD();

      int allele = hap[mk];
      assert(allele!=999);

      ALLELE* hap1;
      ALLELE* hap2;

      for(int i=0; i<2; i++)
	{
	  while(decTrees[i]!=NULL)
	    {
	      if(decTrees[i]->isRec()) decTrees[i]=decTrees[i]->getDec();
	      else if(decTrees[i]->isCoal()) 
		{
		  Node* d1 = decTrees[i]->getleftD();
		  Node* d2 = decTrees[i]->getrightD();
		  if(d1->getID()==d2->getID()) {hap1=d1->getAncLocs(LEFT); hap2=d1->getAncLocs(RIGHT);}
		  else { hap1=transHap(d1, decTrees[i]); hap2=transHap(d2, decTrees[i]);}

		  if(hap1[mk]!=999 && hap2[mk]!=999) break;
		  else if(hap1[mk]!=999 && hap2[mk]==999) decTrees[i] = d1;
		  else if(hap1[mk]==999 && hap2[mk]!=999) decTrees[i] = d2;
		  else{ cerr<<"error!"<<endl; exit(-1);}
		}
	      else {assert(decTrees[i]->isTip()); break;}
	    }//while loop
	}//finish finding two descendants

      int alle[2];
      alle[0]=-1; alle[1]=-1;
      for(int i=0; i<2; i++)
	{
	  ALLELE* hapd = decTrees[i]->getAncLocs(LEFT);
	  alle[i] = hapd[mk];
	  assert(alle[i]!=999);

	  double lamda = theta*0.5*(cursor->getTime()-decTrees[i]->getTime());
	  double prnomut = exp(-lamda);

	  double staf;
	  if(alle[i]==(int)'A') staf=staFreq[0];
	  else if(alle[i]==(int)'C') staf=staFreq[1];
	  else if(alle[i]==(int)'G') staf=staFreq[2];
	  else if(alle[i]==(int)'T') staf=staFreq[3];
	  else { cerr<<"Nucleotide ("<<(char)alle[i]<<") does not exist."<<endl; exit(0); }

	  if(allele!=alle[i]) logllh += log((1-prnomut)*staf);
	  else logllh += log(prnomut + (1-prnomut)*staf);

	  brlen += cursor->getTime()-decTrees[i]->getTime();

	  if(i==0)
	    {
	      nlineages--;
	      if(nlineages==1)//root of the current marker tree
		{
		  double staf;
		  if(allele==(int)'A') staf=staFreq[0];
		  else if(allele==(int)'C') staf=staFreq[1];
		  else if(allele==(int)'G') staf=staFreq[2];
		  else if(allele==(int)'T') staf=staFreq[3];
		  else { cerr<<"Nucleotide ("<<(char)alle[i]<<") does not exist."<<endl; exit(0); }
		  logllh += log(staf);
		}
	    }
	}
    }//loop over all internal nodes    

  assert(nlineages==1);
  logllh -= log(1-exp(-theta*0.5*brlen));

  return logllh;
}




/*-----------------------------------------------------------------------------------------------
  Return the haplotype that is transmitted to cursor from the given dec
-----------------------------------------------------------------------------------------------*/
ALLELE* Graph::transHap(Node* dec, Node* cursor)
{
  ALLELE* hap;
  if(!dec->isRec() || dec->getleftA()->getID()  == cursor->getID()){hap = dec->getAncLocs(LEFT);}
  else{hap = dec->getAncLocs(RIGHT);}

  return hap;
}




void Graph::getTreesFast(vector<Node*>** trees)
{
	vector<Node*>::iterator vi;
	vector<ALLELE>::iterator it;
	
	for(int i=0; i<numLoci; i++){trees[i]->clear();}
	
	for(vi=NodesVector.begin()+numTips; vi!=NodesVector.end(); vi++)
	{
		Node* cursor = (*vi);
		if((*vi)->isRec() || ((*vi)->getleftD()->getID() == (*vi)->getrightD()->getID())) continue;
		
		ALLELE* hap = cursor->getAncLocs(LEFT);
		
		vector<ALLELE> coalescedLocs;
		ALLELE* hapdLeft = transHap(cursor->getleftD(), cursor);
		ALLELE* hapdRight = transHap(cursor->getrightD(), cursor);
		
		for(int i=0; i<numLoci; i++)   
		{
		if(hap[i] != 999 && hapdLeft[i] !=999 && hapdRight[i] != 999){coalescedLocs.push_back(i);}
		}
		

		for(it=coalescedLocs.begin(); it!=coalescedLocs.end(); it++)
		{
			Node* decTrees[2];
			decTrees[0] = cursor->getleftD();
			decTrees[1] = cursor->getrightD();
			
			int mk = *it;
			
			for(int i=0; i<2; i++)
			{
				while(decTrees[i]!=NULL)
				{
					if(decTrees[i]->isRec()) decTrees[i]=decTrees[i]->getDec();
					else if(decTrees[i]->isCoal()) 
					{
						Node* d1 = decTrees[i]->getleftD();
						Node* d2 = decTrees[i]->getrightD();
						ALLELE* hap1;
						ALLELE* hap2;
						if(d1->getID()  ==  d2->getID()) {hap1=d1->getAncLocs(LEFT); hap2=d1->getAncLocs(RIGHT);}
						else {hap1=transHap(d1, decTrees[i]); hap2=transHap(d2, decTrees[i]);}
						
						if(hap1[mk]!=999 && hap2[mk]!=999) {break;}
						else if (hap1[mk]!=999 && hap2[mk]==999) {decTrees[i] = d1;}   
						else if(hap1[mk]==999 && hap2[mk]!=999) {decTrees[i] = d2;}
						else{cerr<<"error!"<<endl; exit(-1);}
					}
					else {assert(decTrees[i]->isTip()); break;}   
				}
			}
			trees[mk]->push_back(decTrees[0]);
			trees[mk]->push_back(decTrees[1]);
			trees[mk]->push_back(cursor);
		} 
	}    
}




double Graph::logllhGraph(double theta, double* staFreq)
{
	double logllh = 0;
	vector<Node*>::iterator vi;
	vector<ALLELE>::iterator it;
	
	int nlineages[numLoci];
	for(int i=0; i<numLoci; i++){nlineages[i] = numTips;}
	
	double brlen[numLoci]; 
	for(int i=0; i<numLoci; i++){brlen[i] = 0;}  
	
	for(vi=NodesVector.begin()+numTips; vi!=NodesVector.end(); vi++)
	{
		Node* cursor = (*vi);
		if((*vi)->isRec() || ((*vi)->getleftD()->getID() == (*vi)->getrightD()->getID())){continue;}
		
		ALLELE* hap = cursor->getAncLocs(LEFT);
		
		vector<ALLELE> coalescedLocs;
		ALLELE* hapdLeft =  transHap(cursor->getleftD(),  cursor);
		ALLELE* hapdRight = transHap(cursor->getrightD(), cursor);
		
		assert(hapdLeft != NULL && hapdRight != NULL);
		
                //Find the vector of coalesced loci
		for(int i=0; i<numLoci; i++)
		{
		if(hap[i] != 999 && hapdLeft[i] != 999 && hapdRight[i] != 999){ coalescedLocs.push_back(i); }
		}
		
		for(it=coalescedLocs.begin(); it!=coalescedLocs.end(); it++)
		{
			Node* decTrees[2];
			decTrees[0] = cursor->getleftD();
			decTrees[1] = cursor->getrightD();
			
			int mk = *it; int allele = hap[mk];
			assert(allele == (int)'A' || allele == (int)'T' || allele == (int)'G' || allele == (int)'C');
			
			for(int i=0; i<2; i++)
			{
				while(decTrees[i]!=NULL)
				{
				        if(decTrees[i]->isRec()){decTrees[i] = decTrees[i]->getDec();}
					else if(decTrees[i]->isCoal()) 
					{
					Node* d1 = decTrees[i]->getleftD(); 
					Node* d2 = decTrees[i]->getrightD();
					ALLELE* hap1;
					ALLELE* hap2;
					if(d1->getID() == d2->getID()) {hap1 = d1->getAncLocs(LEFT); hap2 = d1->getAncLocs(RIGHT);}
					else {hap1 = transHap(d1, decTrees[i]); hap2 = transHap(d2, decTrees[i]);}
						
					if(hap1 != NULL && hap1[mk] !=999   && ((hap2 == NULL) || (hap2 != NULL  && hap2[mk]== 999)))  decTrees[i] = d1;
					else if(hap2!=NULL && hap2[mk]!=999 && ((hap1 == NULL) || (hap1 != NULL  && hap1[mk]== 999)))  decTrees[i] = d2;
					else if(hap1!=NULL && hap2 != NULL && hap1[mk] != 999 && hap2[mk] != 999) break;
					else{cerr<<"error!"<<endl; exit(-1);}
					}
					else {assert(decTrees[i]->isTip()); break;}
				 }//while loop
			}//finish finding two descendants
			
			int alle[2];
			for(int i=0; i<2; i++){
			ALLELE* hapd = decTrees[i]->getAncLocs(LEFT);
			alle[i] = hapd[mk];
			assert(alle[i] == (int)'A' || alle[i] == (int)'T' || alle[i] == (int)'G' || alle[i] == (int)'C');
				
			double lamda = theta*0.5*(cursor->getTime()  -  decTrees[i]->getTime());
			double prnomut = exp(-lamda);
				
			double staf;
				     if(alle[i] == (int) 'A') staf = staFreq[0];
				else if(alle[i] == (int) 'C') staf = staFreq[1];
				else if(alle[i] == (int) 'G') staf = staFreq[2];
				else if(alle[i] == (int) 'T') staf = staFreq[3];
				else {cerr<<"Nucleotide ("<<(char)alle[i]<<") does not exist."<<endl; exit(0); }
				
				if(allele != alle[i]) logllh += log((1-prnomut)*staf);
				else logllh += log(prnomut + (1-prnomut)*staf);
				
				brlen[mk] += cursor->getTime()  -  decTrees[i]->getTime(); 
				
				if(i  ==  0)
				{
					nlineages[mk]--;
					if(nlineages[mk]  ==  1)
					{
						double staf;
						     if(allele == (int)'A') staf = staFreq[0];
						else if(allele == (int)'C') staf = staFreq[1];
						else if(allele == (int)'G') staf = staFreq[2];
						else if(allele == (int)'T') staf = staFreq[3];
						else {cerr<<"Nucleotide ("<<(char)alle[i]<<") does not exist."<<endl; exit(0); }
						logllh += log(staf);
					}
				}
			}//two descendants 
		}//loop over all coalescence loci

 
    }//loop over all internal nodes    
	
	
    for(int i=0; i<numLoci; i++){logllh -= log(1-exp(-theta*0.5*brlen[i]));}
    
    return logllh;
}




int Util::compQsortFreq(const void* a , const void* b)
{
  double tmp = (*(MapHap**)a)->frequency - (*(MapHap**)b)->frequency;
  if(tmp<0) return 1;
  else if(tmp>0) return -1;
  else return 0;
}




bool Util::compareTwoNodesT(Node* a, Node* b)
{
   return a->getTime() < b->getTime();
}




bool Util::compareTwoNodesID(Node* a, Node* b)
{
   return a->getID() < b->getID();
}




bool Util::compareFirstConnection(nodeAncestors* a, nodeAncestors* b)
{
  return a->ancestors[0]->getTime() < b->ancestors[0]->getTime();
}




bool Util::compareTwoNumbers(double a, double b)
{
  return a<b;
}




int Util::randNuc(int oldType, gsl_rng* gslr)
{
   double u = gsl_rng_uniform_pos(gslr);
   double x = 1.0/3.0; double y = 2.0/3.0;

   switch (oldType)
   {
   case 'A': if(u < x) return 'T'; else if(u > x  &&  u < y) return 'G'; else return 'C';
   case 'T': if(u < x) return 'G'; else if(u > x  &&  u < y) return 'C'; else return 'A';
   case 'G': if(u < x) return 'C'; else if(u > x  &&  u < y) return 'A'; else return 'T';
   default : if(u < x) return 'A'; else if(u > x  &&  u < y) return 'T'; else return 'G';
   }
}




double Util::logGamma(double x, double shape, double scale)
{
   double logg = 0;
   for(int i=(int)(shape-1); i>=2; i--) logg += log(i);
   double logpr = (shape-1)*log(x) - shape*log(scale) - x/scale - logg;
   return logpr;
}




double Util::logNormal(double x, double mu, double sigma)
{
   double res = -log(SQR2PI*sigma*x) - pow(log(x)-mu,2)/(2*sigma*sigma);
   return res;
}
