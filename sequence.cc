/*---------------------------------------------------------------
PRISM: A Prime-Encoding approach for frequent Sequence Mining. 
Authors: Karam Gouda, Mosab Hassan and Mohammad J. Zaki
-----------------------------------------------------------------*/
#include <errno.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
//#include <sys/mode.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <sys/types.h>
#include <strings.h>
#include <sys/mman.h>
#include <malloc.h>
#include <strings.h>
#include "Eqclass.h"
#include "Itemset.h"
#include "Lists.h"
#include "extl2.h"
#include "partition.h"

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

#define NONMAXFLG -2
#define MBYTE (1024*1024)

//#define GETBIT(a,b) ((a >> b) & 01)
//#define SETBIT(a,v,b)  ((v != 0) ? (a | (01 << b)): (a & ~(01 << b)))

#define CUSTID(a, b) ((a)[2*(b)])
#define TID(a, b) ((a)[2*(b)+1])
#define LJOIN 0
#define EJOIN 1
#define MJOIN 2

long MEMUSED = 0;
long AVAILMEM = 32*MBYTE;
char fname[300];
char constf[300];               // constraint file
char dataf[300];
char idxf[300];
char conf[300];
char it2f[300];
char seqf[300];

double L2ISECTTIME=0, EXTL1TIME=0, EXTL2TIME=0; 

int *NumLargeItemset;
int maxitemsup;
Itemset *item1, *item2; // for use in reading external dbase

Array *interval, *interval2;
GArray *interv1, *interv2, *interv3;
GArray *inter, *inter2, *inter3;

int maxeqsize = 1;
EqGrNode** eqgraph;
double MAXIMAL_THRESHOLD = 1.0;
int ext_l2_pass = 0;
int use_clique = 0;
char use_constraint =0;
int use_hash = 0;
int num_intersect=0;
int recursive = 0;
FILE *out;
int maxiter = 2;
int min_gap = 0;
int max_gap = (int)pow(2,sizeof(int)*8); //infinity
int window = 0;
char memtrace = 0;
char use_newformat = 1;
int use_ascending = -2;
char use_isetonly = 0;
int L2pruning = 0;
char print_seq = 0;

ofstream mout;
int DBASE_NUM_TRANS;
int DBASE_MAXITEM;
float DBASE_AVG_TRANS_SZ;
float DBASE_AVG_CUST_SZ;
int DBASE_TOT_TRANS;

double MINSUP_PER;
int MINSUPPORT=-1;

int FreqArraySz = 100;
FreqIt **FreqArray;
int FreqArrayPos = 0;

//For using Prime Encoding------------------------------
int N, u2;
int it1_pos_arr[9]={0};
int it2_pos_arr[9]={0};
unsigned char Sup_Array[256]={0,1,1,1,1,1,1,1,1};
unsigned char Prime[8]={19,2,3,5,7,11,13,17};
//Ref_Array to give the order of its indices (prime numbers) inside P_Star
unsigned char Ref_Array[2000] ={0,0,1,2,0,3,0,4,0,0,0,5,0,6,0,0,0,7,0,8};
//unsigned char *GCD_Matrix[256]={0};
unsigned char GCD_Matrix[256][256]={0};
//unsigned char P_Factors[9][256]={{0},{0,1,2,3,4,5,6,7,8}, {0},{0},{0},{0},{0},{0},{0}}; 
unsigned char P_Rem[9][256]={{0},{0}, {0},{0},{0},{0},{0},{0},{0}}; 
//GArray *Common_Factors[256][256]={0}; //used in intersecting compressed sids
unsigned char *Common_Factors[256][256]={0};
//unsigned char **Common_Factors[256]={0}; 
//unsigned char FFGreater[9]={0};  //only for sequence mining
//AspArray contains the products of all prime factors greater than the index
unsigned char AspArray[256]={0}; 
//unsigned char Temp[257][257]={0};
//________________________________________________________


void process_cluster1(Eqclass *cluster, Lists<Eqclass *> *LargeL, int iter);

void add_freq(Itemset *it, int templ)
{
   FreqIt *freq = new FreqIt(it->itemset()->array(), it->size(), templ);
   if (FreqArrayPos+1 >= FreqArraySz){
      FreqArraySz = (int)(1.5*FreqArraySz);
      FreqArray = (FreqIt **)realloc(FreqArray, FreqArraySz*sizeof(FreqIt*));
      if (FreqArray == NULL){
         perror("no mmeory fro FREqArray ");
         exit(-1);
      }
   }
   FreqArray[FreqArrayPos++] = freq;
}

void print_freqary()
{
   int j=0;
   cout << "FREQARRAY " << FreqArrayPos << ":" << endl;
   for (j=0; j < FreqArrayPos; j++){
      cout << *FreqArray[j];
   }
   cout << "!!!!!!!!!!!!!!!!!!!!" << endl;
}

void parse_args(int argc, char **argv)
{
   extern char * optarg;
   int c;

   if (argc < 2)
      cout << "usage: seq -i<infile> -s<support> -e 1\n";
   else{
      while ((c=getopt(argc,argv,"a:bce:fhi:l:m:ors:t:u:v:w:x:"))!=-1){
         switch(c){
         case 'a':
            //if val = -1 then do ascending generation
            //else only generate the eqclass given by the value
            use_ascending = atoi(optarg);
            break;
         case 'b':
            use_isetonly = 1;
            break;
         case 'c': //use constraints file
            use_constraint = 1;
            sprintf(constf,"%s.constraint", fname);            
            break;
         case 'e': //calculate L2 from inverted dbase
            num_partitions = atoi(optarg);
            ext_l2_pass = 1;
            break;
         case 'f':
            use_newformat = 0;
            break;
         case 'h': //use hashing to prune candidates
            use_hash = 1;
            break;
         case 'i': //input file
            sprintf(fname, "%s", optarg);
            sprintf(dataf,"%s.tpose", optarg);
            sprintf(idxf,"%s.idx", optarg);
            sprintf(conf,"%s.conf", optarg);
            sprintf(it2f,"%s.2it", optarg);
            sprintf(seqf,"%s.2seq", optarg);
            break;
         case 'l': //min-gap between items (not implemented)
            min_gap = atoi(optarg);
            break;
         case 'm': //amount of mem available
            AVAILMEM = (long) atof(optarg)*MBYTE;
            break;
         case 'o': //print sequences
            print_seq = 1;
            break;
         case 'r': //use recursive algorithm (doesn't work with pruning)
            recursive = 1;
            break;
         case 's': //min support
            MINSUP_PER = atof(optarg);
            break;
         case 't': //not used
            MAXIMAL_THRESHOLD = atof(optarg);
            break;
         case 'u': //max-gap between items (not implemented)
            max_gap = atoi(optarg);
            break;
         case 'v':
            u2 = atoi(optarg);
            break;
         case 'w': //sliding window (not implemented)
            window = atoi(optarg);
            break;
         case 'x':
            memtrace = 1;
            mout.open(optarg, ios::app);
            break;
         }
      }
   }
   c= open(conf, O_RDONLY);
   if (c < 0){
      perror("ERROR: invalid conf file\n");
      exit(errno);
   }
   read(c,(char *)&DBASE_NUM_TRANS,ITSZ);
   if (MINSUPPORT == -1)
      MINSUPPORT = (int) (MINSUP_PER*DBASE_NUM_TRANS+0.5);
   //ensure that support is at least 2
   if (MINSUPPORT < 2) MINSUPPORT = 2;
   //cout << "MINSUPPORT " << MINSUPPORT << " " << DBASE_NUM_TRANS << endl;
   read(c,(char *)&DBASE_MAXITEM,ITSZ);
   read(c,(char *)&DBASE_AVG_CUST_SZ,sizeof(float));
   read(c,(char *)&DBASE_AVG_TRANS_SZ,sizeof(float));
   read(c,(char *)&DBASE_TOT_TRANS,ITSZ);
   cout << "CONF " << DBASE_NUM_TRANS << " " << DBASE_MAXITEM << " "
        << DBASE_AVG_CUST_SZ << " " << DBASE_AVG_TRANS_SZ << " "
        << DBASE_TOT_TRANS << endl;
   close(c);
}

int choose(int n, int k)
{
   int i;
   int val = 1;

   if (k >= 0 && k <= n){
      for (i=n; i > n-k; i--)
         val *= i;
      for (i=2; i <= k; i++)
         val /= i;
   }

   return val;
}

void get_2OldF_intersect(Itemset *ljoin, Itemset *ejoin,
                         int *it1, int *it2, int sup1, int sup2)
{
   int i,j,k,l,p,f;
   int nval1, nval2;
   int l_increment_flag, e_increment_flag;
   unsigned char l_temp1, e_temp1;
   int temp, l_temp=1, e_temp=1;
   int rm, newindex,oldindex,
       l_rm, l_newindex,l_oldindex=0, e_rm, e_newindex,e_oldindex=0;
   int olpos, oepos;  
   int a, rval;
   int itid, jtid;

   num_intersect++;

   //start the first position with dummy value
   if (ljoin) ljoin->ival2()->optadd(0);
   if (ejoin) ejoin->ival2()->optadd(0);

   f=0;//a counter on indexes of sids of the prefix
   for (i=0,j=0; i < sup1 && j < sup2;){
      itid = it1[i];
      jtid = it2[j];
      if (itid > jtid){ 
         nval2 = it2[j+1]; //nval2 seq length
         j += (nval2+2);   //to move to the next seq id in it2
      }
      else if (itid < jtid){
         nval1 = it1[i+1];  //nval1 seq length
         i += (nval1+2);   //to move to the next seq id in it1
         f++;
      }
      else{
         i++;
         j++;
         nval1 = it1[i++]; // seq size in it1
         nval2 = it2[j++]; // seq size in it2
         rval = it1[i]; //get reference value to start encoding from it

         if (ljoin){
            l_increment_flag = 0;
            if ( rval < it2[j+nval2-1]){
               l_increment_flag=1; 
               olpos = ljoin->ival1size(); //get current size of tid array
               temp = 1;
               oldindex = 0;
	       ljoin->increment_support();
               for (l=j; l < j+nval2; l++){
                  if ( rval < it2[l]){
                     a=it2[l] - rval  + 1;
                     newindex = (a-1)/8;
                     rm = a%8;
                     if(newindex==oldindex)  temp=temp*Prime[rm];
                     else{
                        ljoin->ival1()->optadd(Ref_Array[temp%1966]);
	                oldindex++;
                        for(p=oldindex;p<newindex;p++){              
                           ljoin->ival1()->optadd(0);
	                   oldindex++;
	                }
	                temp=Prime[rm];
                     }
                  }              
               }
               ljoin->ival1()->optadd(Ref_Array[temp%1966]);
            }

            if (l_increment_flag){ //i.e., correspondent sids are matching.
               //Add to sid array the starting position (or the size) 
               //of the current sequence block at the tid array
               l_newindex=(f+1-1)/8;
	       l_rm=(f+1)%8;           
               if(l_newindex==l_oldindex){  
                          l_temp=l_temp*Prime[l_rm];
                          //ljoin->ival2()->optadd(olpos);
                          //if ((ljoin->ival1size()-olpos)==1) 
                            //                ljoin->increment_ones();
                          ljoin->ival2()->optadd(ljoin->ival1size()-olpos);
               }
	       else{
                 if (l_temp == 1) l_temp1 = 0; 
                 else l_temp1 = Ref_Array[l_temp%1966];
                 //Put the encoded sids in the right place of sids array
                 ljoin->ival2()->setitem
                     (ljoin->ival2size()-Sup_Array[l_temp1]-1, l_temp1);
	         l_oldindex++;
	         for(p=l_oldindex;p<l_newindex;p++){
                    ljoin->ival2()->optadd(0);
	            l_oldindex++;
	         }
                 //Add dummy value for the next encoded sids
                 ljoin->ival2()->optadd(0);
                 //ljoin->ival2()->optadd(olpos);
                 //if ((ljoin->ival1size()-olpos)==1) 
                   //                         ljoin->increment_ones();
                 ljoin->ival2()->optadd(ljoin->ival1size()-olpos);
	         l_temp=Prime[l_rm];
	       }
            }
         }
         if (ejoin){
            e_increment_flag = 0;
            temp = 1;
            oldindex = 0;
            oepos = ejoin->ival1size();
            for (k=i, l=j; k < i+nval1 && l < j+nval2;){
               if (it1[k] < it2[l]) k++;
               else if (it1[k] > it2[l]) l++;
               else{
                   a=it2[l] - rval + 1;
                   newindex = (a-1)/8;
                   rm = a%8;
                   if(newindex==oldindex)  temp=temp*Prime[rm];
                   else{
                       ejoin->ival1()->optadd(Ref_Array[temp%1966]);
	               oldindex++;
                       for(int p=oldindex;p<newindex;p++){               
                          ejoin->ival1()->optadd(0);
	                  oldindex++;
	               }
	               temp=Prime[rm];
                   }
                   e_increment_flag=1;
                   k++;
                   l++;
               }
            }
            if (e_increment_flag){
               ejoin->ival1()->optadd(Ref_Array[temp%1966]);
	       ejoin->increment_support();
               e_newindex=(f+1-1)/8;
	       e_rm=(f+1)%8;
               if(e_newindex==e_oldindex){  
                          e_temp=e_temp*Prime[e_rm];
                          //ejoin->ival2()->optadd(oepos);
                          ejoin->ival2()->optadd(ejoin->ival1size()-oepos);
               }
	       else{
                 if (e_temp == 1) e_temp1 = 0; 
                 else e_temp1 = Ref_Array[e_temp%1966];
                 //Put the commpressed sids in the right place of the 2nd array
                 ejoin->ival2()->setitem
                     (ejoin->ival2size()-Sup_Array[e_temp1]-1, e_temp1);
	         e_oldindex++;
	         for(p=e_oldindex;p<e_newindex;p++){
                    ejoin ->ival2()->optadd(0);
	            e_oldindex++;
	         }
                 //Add dummy value for the next compressed sids
                 ejoin ->ival2()->optadd(0);
                 //ejoin->ival2()->optadd(oepos);
                 ejoin->ival2()->optadd(ejoin->ival1size()-oepos);
	         e_temp=Prime[e_rm];
	       }
            }
         }
         i += nval1;
         j += nval2;
         f++;
      }
   }
   if (ljoin){
      l_temp1= Ref_Array[l_temp%1966];
      ljoin->ival2()->setitem
                     (ljoin->ival2size()-Sup_Array[l_temp1]-1, l_temp1);
   }
   if (ejoin){
      e_temp1= Ref_Array[e_temp%1966];
      ejoin->ival2()->setitem
                     (ejoin->ival2size()-Sup_Array[e_temp1]-1, e_temp1);
   }

}


void get_2NewF_intersect(Itemset *ljoin, Itemset *ejoin,
                         int *it1, int *it2, int sup1, int sup2)
{
   int i,j,k,l,p,f;
   int nval1, nval2;
   int l_increment_flag, e_increment_flag;
   unsigned char l_temp1, e_temp1;
   int temp, l_temp=1, e_temp=1;
   int rm, newindex,oldindex,
       l_rm, l_newindex,l_oldindex=0, e_rm, e_newindex,e_oldindex=0;
   int olpos, oepos;  
   int a, rval;
   int itid, jtid;

   num_intersect++;

   //Start the first position with dummy value
   if (ljoin) ljoin->ival2()->optadd(0);
   if (ejoin) ejoin->ival2()->optadd(0);

   f=0;//A counter on indexes of the prefix's sids
   for (i=0,j=0; i < sup1 && j < sup2;){
      itid = it1[i];
      jtid = it2[j];
      if (itid > jtid){ 
         nval2 = j;
         while(it2[j] == it2[nval2] && nval2 < sup2) nval2 += 2;
         j = nval2;
      }
      else if (itid < jtid){
         nval1 = i;
         while(it1[i] == it1[nval1] && nval1 < sup1) nval1 += 2;
         i = nval1; 
         f++;
      }
      else{
         nval1 = i;
         nval2 = j;
         rval = it1[i+1]; //Get a reference value to start encoding from it
         while(it1[i] == it1[nval1] && nval1 < sup1) nval1 += 2;
         while(it2[j] == it2[nval2] && nval2 < sup2) nval2 += 2;
         if (ljoin){
            l_increment_flag = 0;
            if ( rval < it2[nval2-1]){
               l_increment_flag=1; 
               olpos = ljoin->ival1size(); //Get current size of tid array
               temp = 1;
               oldindex = 0;
	       ljoin->increment_support();
               for (l=j; l < nval2; l+=2){
                  if ( rval < it2[l+1]){
                     a=it2[l+1] - rval  + 1;
                     newindex = (a-1)/8;
                     rm = a%8;
                     if(newindex==oldindex)  temp=temp*Prime[rm];
                     else{
                        ljoin->ival1()->optadd(Ref_Array[temp%1966]);
	                oldindex++;
                        for(p=oldindex;p<newindex;p++){              
                           ljoin->ival1()->optadd(0);
	                   oldindex++;
	                }
	                temp=Prime[rm];
                     }
                  }              
               }
               ljoin->ival1()->optadd(Ref_Array[temp%1966]);
            }

            if (l_increment_flag){ //i.e., corresponding sids are matching.
               //Add to sid array the starting position (or the size) 
               //of the current sequence block at the tid array
               l_newindex=(f+1-1)/8;
	       l_rm=(f+1)%8;           
               if(l_newindex==l_oldindex){  
                          l_temp=l_temp*Prime[l_rm];
                          //ljoin->ival2()->optadd(olpos);
                          //if ((ljoin->ival1size()-olpos)==1) 
                            //                ljoin->increment_ones();
                          ljoin->ival2()->optadd(ljoin->ival1size()-olpos);
               }
	       else{
                 if (l_temp == 1) l_temp1 = 0; 
                 else l_temp1 = Ref_Array[l_temp%1966];
                 //Put the encoded sids in the right place of sids array
                 ljoin->ival2()->setitem
                     (ljoin->ival2size()-Sup_Array[l_temp1]-1, l_temp1);
	         l_oldindex++;
	         for(p=l_oldindex;p<l_newindex;p++){
                    ljoin->ival2()->optadd(0);
	            l_oldindex++;
	         }
                 //Add dummy value for the next encoded sids
                 ljoin->ival2()->optadd(0);
                 //ljoin->ival2()->optadd(olpos);
                 //if ((ljoin->ival1size()-olpos)==1) 
                   //                         ljoin->increment_ones();
                 ljoin->ival2()->optadd(ljoin->ival1size()-olpos);
	         l_temp=Prime[l_rm];
	       }
            }
         }
         if (ejoin){
            e_increment_flag = 0;
            temp = 1;
            oldindex = 0;
            oepos = ejoin->ival1size();
            for (k=i, l=j; k < nval1 && l < nval2;){
               if (it1[k+1] < it2[l+1]) k+=2;
               else if (it1[k+1] > it2[l+1]) l+=2;
               else{
                   a=it2[l+1] - rval + 1;
                   newindex = (a-1)/8;
                   rm = a%8;
                   if(newindex==oldindex)  temp=temp*Prime[rm];
                   else{
                       ejoin->ival1()->optadd(Ref_Array[temp%1966]);
	               oldindex++;
                       for(int p=oldindex;p<newindex;p++){               
                          ejoin->ival1()->optadd(0);
	                  oldindex++;
	               }
	               temp=Prime[rm];
                   }
                   e_increment_flag=1;
                   k+=2;
                   l+=2;
               }
            }
            if (e_increment_flag){
               ejoin->ival1()->optadd(Ref_Array[temp%1966]);
	       ejoin->increment_support();
               e_newindex=(f+1-1)/8;
	       e_rm=(f+1)%8;
               if(e_newindex==e_oldindex){  
                          e_temp=e_temp*Prime[e_rm];
                          //ejoin->ival2()->optadd(oepos);
                          ejoin->ival2()->optadd(ejoin->ival1size()-oepos);
               }
	       else{
                 if (e_temp == 1) e_temp1 = 0; 
                 else e_temp1 = Ref_Array[e_temp%1966];
                 //Put the encoded sids in the right place of the 2nd array
                 ejoin->ival2()->setitem
                     (ejoin->ival2size()-Sup_Array[e_temp1]-1, e_temp1);
	         e_oldindex++;
	         for(p=e_oldindex;p<e_newindex;p++){
                    ejoin ->ival2()->optadd(0);
	            e_oldindex++;
	         }
                 //Add dummy value for the next encoded sids
                 ejoin ->ival2()->optadd(0);
                 //ejoin->ival2()->optadd(oepos);
                 ejoin->ival2()->optadd(ejoin->ival1size()-oepos);
	         e_temp=Prime[e_rm];
	       }
            }
         }
         i = nval1;
         j = nval2;
         f++;
      }
   }
   if (ljoin){
      l_temp1= Ref_Array[l_temp%1966];
      ljoin->ival2()->setitem
                     (ljoin->ival2size()-Sup_Array[l_temp1]-1, l_temp1);
   }
   if (ejoin){
      e_temp1= Ref_Array[e_temp%1966];
      ejoin->ival2()->setitem
                     (ejoin->ival2size()-Sup_Array[e_temp1]-1, e_temp1);
   }

}

int interval_comp(int a, int b, int c, int d)
{
   if (a < c) return -1;
   else if (a > c) return 1;
   else{
      if (b < d) return -1;
      else if (b > d) return 1;
      else return 0;
   }
}

void make_itemset(Itemset *it, GArray *ary, GArray *ary2)
{
   int i;
   for (i=0; i < ary->size(); i++){
      it->ival1()->optadd((*ary)[i]);
   }
   for (i=0; i < ary2->size(); i++){
      it->ival2()->optadd((*ary2)[i]);
   }
}

/*
void get_tmpEqNewF_intersect(Itemset *&ljoin, int &lcnt, Itemset *it1, 
                               int iter, int x, int list)
{
   int i,j,l, m, f, ii, lz;
   int a, c,ysup,sp;
   int nval1;
   int flag1;
   int olpos;
   int lones=0;
   int Mlpos;
   Itemset *templ=ljoin;
   Itemset *temp=ljoin;
   int it1post=0;
   unsigned char N_gcd_Val;
   num_intersect++;
   lcnt = 0;
   int dc1 = it1->support()-MINSUPPORT;
   int df1=it1->ones();   
   
   //if (df1 > dc1) templ=NULL;
if (df1 <= dc1){
   ljoin = new Itemset(iter, it1->ival1size()-df1, it1->ival2size()-df1);
   for (i=0; i < it1->ival2size(); ){
      if (templ==NULL) break;
      a = it1->ival2(i);
      if(a==0){
          if (list==2) ljoin->ival2()->optadd(0);
          i+=1; 
          continue;
      }
      c=Sup_Array[a];
      Mlpos=ljoin->ival2size();
      ljoin->ival2()->optadd(a);
      lz=0;
      for(int h=1; h<=c; h++){
            if (templ==NULL) break;
            //ii=i+h;
            nval1=it1->ival2(i+h);
            if (nval1 == 1 ){ 
                   flag1=it1->ival1(it1post);
                   if (Sup_Array[flag1]!=1){
                        N_gcd_Val=P_Rem[1][flag1]; 
                        //lary2->optadd(lary->size());
                        ljoin->ival2()->optadd(1);  
                        ljoin->ival1()->optadd(N_gcd_Val); 
                        if (Sup_Array[N_gcd_Val]==1) lones++;
                        lcnt++; lz++;
                   }
                   else ljoin->ival2()->setitem(Mlpos, P_Rem[lz+1][ljoin->ival2()->item(Mlpos)]); 
            }
            else{ 
	       ysup=0;
               f=it1post+nval1;
               for(l=it1post;l<(it1post+nval1);l++){
                  if(it1->ival1(l)!=0){
                     N_gcd_Val=P_Rem[1][it1->ival1(l)];
                     m = l+1;
                     if (m < f){
                        //lary2->optadd(lary->size());
                        ljoin->ival2()->optadd(f-l);                     
                        ljoin->ival1()->optadd(N_gcd_Val); 
                        ysup++;
                        for(; m < f; m++) ljoin->ival1()->optadd(it1->ival1(m));           
                     }
                     else if (N_gcd_Val != 0){ 
                           //lary2->optadd(lary->size());
                           ljoin->ival2()->optadd(1); 
                           ljoin->ival1()->optadd(N_gcd_Val); 
                           if (Sup_Array[N_gcd_Val]==1) lones++;
                           ysup++; 
                     }
                     break;
                  }
               }
	       if (ysup>0){lz++; lcnt++;}
	       else{
                  ljoin->ival2()->setitem(Mlpos, P_Rem[lz+1][ljoin->ival2()->item(Mlpos)]); 
                  if ((++df1) > dc1) templ=NULL;
               }
            }
            it1post+=nval1;
      }   
      i+=c+1; 
   }
   if (lcnt >= MINSUPPORT){
      //ljoin=templ;
      ljoin->set_ones(lones);
      ljoin->set_support(lcnt);
   }
   else delete ljoin;
}
}
*/
void get_tmpEqNewF_intersect(Itemset *&ljoin, int &lcnt, Itemset *it1, 
                               int iter, int x, int list)
{
   int i,j,l, m, h, f, lz;
   int a, c,ysup,sp;
   int nval1;
   int flag1;
   GArray *lary;
   GArray *lary2;
   int olpos;
   int lones=0;
   int Mlpos;
   Itemset *templ=ljoin;
   int it1post=0;
   unsigned char N_gcd_Val;
   num_intersect++;
   lcnt = 0;
   int dc1 = it1->support()-MINSUPPORT;
   int df1=it1->ones();   
   
   //if (df1 > dc1) templ=NULL;
if (df1 <= dc1){

   lary = interv1;
   lary2 = inter;
   lary->reset();   
   lary2->reset();

   for (i=0; i < it1->ival2size(); ){
      if (templ==NULL) break;
      a = it1->ival2(i);
      if(a==0){
          if (list==2) lary2->optadd(0);
          i+=1; 
          continue;
      }
      c=Sup_Array[a];
      Mlpos=lary2->size();
      lary2->optadd(a);
      lz=0;
      for(h=1; h<=c; h++){
            nval1=it1->ival2(i+h);
            if (nval1 == 1 ){ 
                   flag1=it1->ival1(it1post);
                   if (Sup_Array[flag1]!=1){
                        N_gcd_Val=P_Rem[1][flag1]; 
                        //lary2->optadd(lary->size());
                        lary2->optadd(1);  
                        lary->optadd(N_gcd_Val); 
                        if (Sup_Array[N_gcd_Val]==1) lones++;
                        lcnt++; lz++;
                   }
                   else lary2->setitem(Mlpos, P_Rem[lz+1][lary2->item(Mlpos)]); 
            }
            else{ 
	       ysup=0;
               f=it1post+nval1;
               for(l=it1post;l<f;l++){
                  if(it1->ival1(l)!=0){
                     N_gcd_Val=P_Rem[1][it1->ival1(l)];
                     m = l+1;
                     if (m < f){
                        //lary2->optadd(lary->size());
                        lary2->optadd(f-l);                     
                        lary->optadd(N_gcd_Val); 
                        ysup++;
                        for(; m < f; m++) lary->optadd(it1->ival1(m));           
                     }
                     else if (N_gcd_Val != 0){ 
                           //lary2->optadd(lary->size());
                           lary2->optadd(1); 
                           lary->optadd(N_gcd_Val); 
                           if (Sup_Array[N_gcd_Val]==1) lones++;
                           ysup++; 
                     }
                     break;
                  }
               }
	       if (ysup>0){lz++; lcnt++;}
	       else{
                  lary2->setitem(Mlpos, P_Rem[lz+1][lary2->item(Mlpos)]); 
                  if ((++df1) > dc1) templ=NULL;
               }
            }
            it1post+=nval1;
      }   
      i+=c+1; 
   }

   if (lcnt >= MINSUPPORT){
      ljoin = new Itemset(iter, lary->size(), lary2->size());
      ljoin->set_ones(lones);
      ljoin->set_support(lcnt);
      make_itemset(ljoin, lary, lary2);
   }
}
}


void get_tmpNewF_intersect(Itemset *&ljoin, Itemset *&ejoin, Itemset *&mjoin,
                       int &lcnt, int &ecnt, int &mcnt,
                       Itemset *it1, Itemset *it2, int iter, int x, int list)
{

   int i,j,h,k,l, m, f, d, ii, jj, lz, mz, ez;
   //char ff;
   int a, b, c, ysup, sp;
   int nval1, nval2;
   int flag1, lflge, flag, P0;
   GArray *lary, *eary, *mary;
   GArray *lary2, *eary2, *mary2;
   int olpos, oepos, ompos;
   int lones=0, mones=0, eones=0;
   int nlmin=0, nemin=0, nmmin=0;
   int Mlpos, Mmpos, Mepos;
   Itemset *templ=ljoin, *tempe=ejoin, *tempm = mjoin;
   int it1positioner=0, it2positioner=0; 
   int intersct_fact1, intersct_fact2;
   //int count, r1, r2;
   //GArray *Common;
   unsigned char *Common;
   int it1post, it2post;
   int ll, sum;
   unsigned char N_Val, N_gcd_Val;

   //int s1=0, s2=0;
   num_intersect++;
   lcnt = ecnt = mcnt = 0;

   int dc1 = it1->support()-MINSUPPORT;
   int dc2 = it2->support()-MINSUPPORT;
   int df1=0;
   int df2=0;
   

   //Preproccessing loop to get initial values to df1 and df2
   for (i=0,j=0; i < it1->ival2size() && j < it2->ival2size(); ){
      if (df1 > dc1 || df2 > dc2){ 
         templ=NULL; tempm=NULL; tempe=NULL;
         break;
      }
      a = it1->ival2(i);
      b = it2->ival2(j);

/////////////////////////////
      if (a==0){ 
         i++;
         j+=Sup_Array[b]+1;
         df2+=Sup_Array[b];
         continue;
      }
      if (b==0){
         i+=Sup_Array[a]+1; 
         j++;
         df1+=Sup_Array[a];
         continue;
      }      
//////////////////////////
      N_Val=GCD_Matrix[a][b];
      c=Sup_Array[N_Val];


//      //In case you want to build ljoin, ejoin, and mjoin before computations
//      Common=Common_Factors[a][b];
//      d=0;
//      for(h=1; h<=c; h++){
//         s1+=it1->ival2(i+Common->item(d++));
//         s2+=it2->ival2(j+Common->item(d++));
//      }


      df1+=Sup_Array[a]-c;
      df2+=Sup_Array[b]-c;
      i+=Sup_Array[a]+1; 
      j+=Sup_Array[b]+1;
   }

if (templ!=NULL || tempm!=NULL || tempe!=NULL){

   lary = interv1;
   eary = interv2;
   mary = interv3;
   lary->reset();
   eary->reset();
   mary->reset();

   lary2 = inter;
   eary2 = inter2;
   mary2 = inter3;
   lary2->reset();
   eary2->reset();
   mary2->reset();

   for (i=0,j=0; i < it1->ival2size() && j < it2->ival2size(); ){

      if (templ==NULL && tempm==NULL && tempe==NULL) break;
      a = it1->ival2(i);
      b = it2->ival2(j);
      
      if (b!=0){
         sum=0;
         for (ll=2; ll<=Sup_Array[b]+1; ll++){
             sum+= it2->ival2(j+ll-1);
             it2_pos_arr[ll]=sum;
         } 
      }    
      //In case of list==1 we can get rid of zeros in it1 that appeared in the 
      //previous iteration, otherwise we can not
      if(a==0){
          if (list==2){
            lary2->optadd(0);
            eary2->optadd(0);
            mary2->optadd(0);
          }
          i+=1; 
          j+=Sup_Array[b]+1;
          if (b!=0) it2positioner+=it2_pos_arr[Sup_Array[b]+1];
          continue;
      }

      sum=0;
      for (ll=2; ll<=Sup_Array[a]+1; ll++){
          sum+= it1->ival2(i+ll-1);
          it1_pos_arr[ll]=sum;
      }

      if(b==0){ 
          mary2->optadd(0);          
          eary2->optadd(0);
          lary2->optadd(0);
          j+=1; 
          i+=Sup_Array[a]+1;
          it1positioner+=it1_pos_arr[Sup_Array[a]+1];
          continue;
      }
      N_Val=GCD_Matrix[a][b];
      c=Sup_Array[N_Val];      
      if(c==0){
          lary2->optadd(0);
          mary2->optadd(0);
          eary2->optadd(0); 
          i+=Sup_Array[a]+1; 
          j+=Sup_Array[b]+1;
          it1positioner+=it1_pos_arr[Sup_Array[a]+1];
          it2positioner+=it2_pos_arr[Sup_Array[b]+1];
          continue;
      }
      
      //Get the exact position for values that would be modified later
      Mlpos=lary2->size();
      Mmpos=mary2->size();
      Mepos=eary2->size();

      //Add the value, gcd(a,b), that would be modified later
      lary2->optadd(N_Val);
      mary2->optadd(N_Val);
      eary2->optadd(N_Val);

      //if (a<=b){ ff=1; Common=Common_Factors[a][b];}
      //else {ff=0; Common=Common_Factors[b][a];}
      Common=Common_Factors[a][b];
      
      d=0;//A counter on the third dimension of Common matrix for element (a,b) 
      lz=0; mz=0;ez=0;//Variables that help in the modification of N_Val
      for(h=1; h<=c; h++){

        //if (ff){
           intersct_fact1=Common[d++];
           intersct_fact2=Common[d++];
           //intersct_fact1=Common->item(d++);
           //intersct_fact2=Common->item(d++);
        //}
        // else{
          //  intersct_fact2=Common[d++];
          //  intersct_fact1=Common[d++];
          // intersct_fact2=Common->item(d++);
          // intersct_fact1=Common->item(d++);
         //}

         ii=i+intersct_fact1;
         jj=j+intersct_fact2;

/*
         //This part is used in case of using start position of each block
         //instead of using block size

         r1=it1->ival2(ii);//starting index value of block1 in the tid array
         r2=it2->ival2(jj);//starting index value of block2 in the tid array
         //to determine the current block1 and block2 sizes         
         if (h==c){
                if (intersct_fact1==Sup_Array[a]){
                     count=ii;
                     while(count < it1->ival2size()-1){
                          if (it1->ival2(count+1) != 0) break;
                          count++;
                     }
                     if (count==it1->ival2size()-1) 
                          nval1=it1->ival1size()-r1;
                     else nval1=it1->ival2(count+2) - r1;//block1 size 
                }
                else nval1=it1->ival2(ii+1) - r1;//block1 size 

                if (intersct_fact2==Sup_Array[b]){
                     count=jj;
                     while(count < it2->ival2size()-1){
                          if (it2->ival2(count+1) != 0) break;
                          count++;
                     }
                     if (count==it2->ival2size()-1) 
                          nval2=it2->ival1size()-r2;
                     else nval2=it2->ival2(count+2) - r2;//block2 size 
                }
                else nval2=it2->ival2(jj+1) - r2;//block2 size 
         }
         else{
            nval1=it1->ival2(ii+1) - r1;//block1 size 
            nval2=it2->ival2(jj+1) - r2;//block2 size 
         }

*/
         it1post=it1positioner+it1_pos_arr[intersct_fact1];
         it2post=it2positioner+it2_pos_arr[intersct_fact2];
         nval1=it1->ival2(ii);
         nval2=it2->ival2(jj);      
         P0=min(nval1,nval2);
         flag1=0;   
         if (templ){
            if (nval2 == 1 ){ 
                   flag1=1;
                   //sp=AspArray[it1->ival1(r1)];
	           //N_gcd_Val=GCD_Matrix[sp][it2->ival1(r2)];
                   sp=AspArray[it1->ival1(it1post)];
	           N_gcd_Val=GCD_Matrix[sp][it2->ival1(it2post)];
                   //N_gcd_Val=Temp[it1->ival1(it1post)][it2->ival1(it2post)]; 
                   if (N_gcd_Val != 0){ 
                        //lary2->optadd(lary->size());
                        lary2->optadd(1);  
                        lary->optadd(N_gcd_Val);  
                        if (Sup_Array[N_gcd_Val]==1) lones++;
                        lcnt++; lz++;
                   }
                   else{
                     //Diminish the value at Mlpos to the exact value
                     lary2->setitem(Mlpos, P_Rem[lz+1][lary2->item(Mlpos)]);  
                     nlmin++;
                     if ((nlmin + df1) > dc1 ||(nlmin + df2) > dc2 ) templ=NULL;
                   } 
            }
            else{ 
	    ysup=0;
            //f = r2+nval2;
            f = it2post+nval2;
            for(l=it1post, k=it2post;l<(it1post+P0);l++, k++){
            //for(l=r1, k=r2;l<(r1+P0);l++, k++){
	       if(it1->ival1(l)!=0){
                  sp=AspArray[it1->ival1(l)];
	          N_gcd_Val=GCD_Matrix[sp][it2->ival1(k)]; 
                  //N_gcd_Val=Temp[it1->ival1(l)][it2->ival1(k)];        
                  m = k+1;
                  if (m < f){
                     //lary2->optadd(lary->size());
                     lary2->optadd(f-k);                    
                     lary->optadd(N_gcd_Val); 
                     ysup++;
                     for(; m < f; m++) lary->optadd(it2->ival1(m)); 
                  }
                  else if (N_gcd_Val != 0){ 
                           //lary2->optadd(lary->size());
                           lary2->optadd(1); 
                           lary->optadd(N_gcd_Val); 
                           if (Sup_Array[N_gcd_Val]==1) lones++;
                           ysup++; 
                  }
                  break;
               }
            }
	    if (ysup>0){lz++; lcnt++;}
	    else{
               //here you should record the modified N_Val 
               lary2->setitem(Mlpos, P_Rem[lz+1][lary2->item(Mlpos)]); 
               nlmin++;
               if ((nlmin + df1) > dc1 ||(nlmin + df2) > dc2 ) templ=NULL;
            }
           }
         }  
         if (tempm){
            if (nval1 == 1){
               flag1=1;
               //sp=AspArray[it2->ival1(r2)];
	       //N_gcd_Val=GCD_Matrix[sp][it1->ival1(r1)];
               sp=AspArray[it2->ival1(it2post)];
	       N_gcd_Val=GCD_Matrix[sp][it1->ival1(it1post)];
               //N_gcd_Val=Temp[it2->ival1(it2post)][it1->ival1(it1post)]; 
               if (N_gcd_Val != 0){
                           //mary2->optadd(mary->size());
                           mary2->optadd(1);  
                           mary->optadd(N_gcd_Val); 
                           if (Sup_Array[N_gcd_Val]==1) mones++; 
                           mcnt++;mz++;
               }
               else{ 
                  mary2->setitem(Mmpos, P_Rem[mz+1][mary2->item(Mmpos)]);
                  nmmin++;
                  if ((nmmin + df1) > dc1 ||(nmmin + df2) > dc2 ) tempm=NULL;
               } 
            }
            else{
	    ysup=0;
            //f = r1+nval1;
            f = it1post+nval1;
            for(k=it1post, l=it2post;l<(it2post+P0);l++, k++){
            //for(k=r1, l=r2;l<(r2+P0);l++, k++){
               if(it2->ival1(l)!=0){
                  sp=AspArray[it2->ival1(l)];
	          N_gcd_Val=GCD_Matrix[sp][it1->ival1(k)];
                  //N_gcd_Val=Temp[it2->ival1(l)][it1->ival1(k)];                  
                  m = k+1;
                  if (m < f){
                     //mary2->optadd(mary->size());
                     mary2->optadd(f-k);                    
                     mary->optadd(N_gcd_Val); 
                     ysup++;
                     for(; m < f; m++)  mary->optadd(it1->ival1(m));  
                  }
                  else if (N_gcd_Val != 0){ 
                           //mary2->optadd(mary->size());
                           mary2->optadd(1); 
                           mary->optadd(N_gcd_Val); 
                           if (Sup_Array[N_gcd_Val]==1) mones++; 
                           ysup++; 
                  }
                  break;
               }
            }
	    if (ysup>0){mz++; mcnt++;}
	    else{
                mary2->setitem(Mmpos, P_Rem[mz+1][mary2->item(Mmpos)]);
                nmmin++;
                if ((nmmin + df1) > dc1 ||(nmmin + df2) > dc2 ) tempm=NULL;
            }
            }
         }
         if (tempe){      
            if (flag1){
	           N_gcd_Val=GCD_Matrix[it1->ival1(it1post)][it2->ival1(it2post)];
                   //N_gcd_Val=GCD_Matrix[it1->ival1(r1)][it2->ival1(r2)];
                   if (N_gcd_Val != 0){ 
                           //eary2->optadd(eary->size());
                           eary2->optadd(1);  
                           eary->optadd(N_gcd_Val);
                       //The next if statement is closed 
                       //since each pattern here will not join itself later 
                           //if (Sup_Array[N_gcd_Val]==1) eones++;
                           ecnt++;ez++;
                   }
                   else{
                     eary2->setitem(Mepos, P_Rem[ez+1][eary2->item(Mepos)]);
                     nemin++;
                     if ((nemin + df1) > dc1 ||(nemin + df2) > dc2 ) tempe=NULL;
                   } 
            }
            else{ 
	    ysup=0;
            oepos = eary->size();
            flag = 1;
            for(l=it1post, k=it2post; l < (it1post+P0); l++,k++){
	    //for(l=r1, k=r2; l < (r1+P0); l++,k++){
                //if statement to get rid of zeros that appeared in the 
                //previous iteration
                if (flag) 
                      if (x==2) { 
                           if (it2->ival1(k)==0) continue; 
                      }
                      else if (it1->ival1(l)==0) continue;
                flag = 0;
                N_gcd_Val=GCD_Matrix[it1->ival1(l)][it2->ival1(k)];
	        eary->optadd(N_gcd_Val);
                ysup=ysup+Sup_Array[N_gcd_Val];
	    }
	    if (ysup > 0){
               //eary2->optadd(oepos);
               eary2->optadd(eary->size()-oepos);
               //if statement is closed since each pattern here will not join itself later 
               //if (ysup==1) eones++;
               ecnt++; 
               ez++;
            }
	    else{
               eary->set_size(oepos); 
               eary2->setitem(Mepos, P_Rem[ez+1][eary2->item(Mepos)]);
               nemin++;
               if ((nemin + df1) > dc1 ||(nemin + df2) > dc2 ) tempe=NULL;
            }
            }
         } 
      }
      i+=Sup_Array[a]+1; 
      j+=Sup_Array[b]+1;
      it1positioner+=it1_pos_arr[Sup_Array[a]+1];
      it2positioner+=it2_pos_arr[Sup_Array[b]+1];
   }
   if (ljoin && lcnt >= MINSUPPORT){
      ljoin = new Itemset(iter, lary->size(), lary2->size());
      ljoin->set_ones(lones);
      ljoin->set_support(lcnt);
      make_itemset(ljoin, lary, lary2);
   }
   if (ejoin && ecnt >= MINSUPPORT){
      ejoin = new Itemset(iter, eary->size(), eary2->size());
      //ejoin->set_ones(eones); //it is closed since each ejoin will not join itself later 
      ejoin->set_support(ecnt);
      make_itemset(ejoin, eary, eary2);
   }
   if (mjoin && mcnt >= MINSUPPORT){
      mjoin = new Itemset(iter, mary->size(), mary2->size());
      mjoin->set_ones(mones);
      mjoin->set_support(mcnt);
      make_itemset(mjoin, mary, mary2);
   }
}
}


void fill_seq_template(Eqclass *EQ, Eqclass *parent, int LR)
{
   if (LR == 1){
      EQ->set_templ(parent->templ()*2+1);
      EQ->set_templ2(parent->templ()*2);
   }
   else if (LR == 2){
      EQ->set_templ(parent->templ2()*2+1);
      EQ->set_templ2(parent->templ2()*2);
   }
}

int get_valid_el(int it, char *ibvec, char *sbvec)
{
   int i, j;
   int i1, i2;
   int rval = 0;

   // cout << "[" << it << "] = ";
   for (i=0; i < eqgraph[it]->seqnum_elements(); i++){
      sbvec[i] = 0;
      //cout << " " << eqgraph[it]->seqget_element(i);
   }
   //cout << " --- ";
   for (i=0; i < eqgraph[it]->num_elements(); i++){
      ibvec[i] = 0;
      //cout << " " << eqgraph[it]->get_element(i);
   }
   //cout << endl;

   for (i=0; i < eqgraph[it]->seqnum_elements(); i++){
      i1 = eqgraph[it]->seqget_element(i);
      for (j=i; j < eqgraph[it]->seqnum_elements(); j++){
         i2 = eqgraph[it]->seqget_element(j);
         if (eqgraph[i1] && eqgraph[i1]->seqfind(i2)){
            sbvec[i] = 1;
            sbvec[j] = 1;
            rval = 1;
         }
         if (j > i){
            if ((eqgraph[i1] && eqgraph[i1]->find(i2)) ||
                (eqgraph[i2] && eqgraph[i2]->seqfind(i1))){
               sbvec[i] = 1;
               sbvec[j] = 1;
               rval = 1;
            }
         }
      }
   }

   for (i=0; i < eqgraph[it]->num_elements(); i++){
      i1 = eqgraph[it]->get_element(i);
      if (eqgraph[i1]){
         for (j=i+1; j < eqgraph[it]->num_elements(); j++){
            i2 = eqgraph[it]->get_element(j);
            if (eqgraph[i1]->find(i2)){
               ibvec[i] = 1;
               ibvec[j] = 1;
               rval = 1;
            }
         }
         for (j=0; j < eqgraph[it]->seqnum_elements(); j++){
            i2 = eqgraph[it]->seqget_element(j);
            if (eqgraph[i1]->seqfind(i2)){
               ibvec[i] = 1;
               sbvec[j] = 1;
               rval =1;
            }
         }
      }
   }

   for (i=0; i < eqgraph[it]->seqnum_elements(); i++)
      if (!sbvec[i]) L2pruning++;
   for (i=0; i < eqgraph[it]->num_elements(); i++)
      if (!ibvec[i]) L2pruning++;

   return rval;
}

//construct the next set of eqclasses from external disk
Eqclass* get_ext_eqclass(int it)
{
   double t1, t2;
   int temp, nval1, ejoin_size = 0, ljoin_size = 0;

   seconds(t1);
   //cout << "MEMEXT " << it << " " << MEMUSED << endl;
   int i, k, it2, supsz, supsz2;
   Itemset *ljoin = NULL;
   Itemset *ejoin = NULL;

   char *ibvec, *sbvec;
   ibvec = sbvec = NULL;
   if (eqgraph[it]->num_elements() > 0)
      ibvec = new char[eqgraph[it]->num_elements()];
   if (eqgraph[it]->seqnum_elements() > 0)
      sbvec = new char[eqgraph[it]->seqnum_elements()];

   if (!get_valid_el(it, ibvec, sbvec)) return NULL;

   Eqclass *L2 = new Eqclass(1, EQCTYP1);
   if (L2 == NULL)
   {
      perror("memory exceeded : ext_class ");
      exit (errno);
   }
   //init seq pattern templates
   L2->set_templ(1);
   L2->set_templ2(0);

   interval->reset();
   interval2->reset();

   supsz = partition_get_idxsup(it);
   partition_read_item(interval->array(), it);

   /////////////////////////////////////////////////////////////////
   ///to get number of cust. in item it(N)
   N=0;
   for(int m=0;m<supsz-1;m+=2)
    {
     if(interval->array()[m]!=interval->array()[m+2])
       N++;
    }
   /////////////////////////////////////////////////////////////////

   int tmpit;
   for (i=0, k=0; i < eqgraph[it]->num_elements() ||
           k < eqgraph[it]->seqnum_elements();){
      ljoin = NULL;
      ejoin = NULL;
      it2 = DBASE_MAXITEM+1;
      tmpit = DBASE_MAXITEM+1;
      if (i < eqgraph[it]->num_elements() && ibvec[i])
         it2 = eqgraph[it]->get_element(i);
      if (k < eqgraph[it]->seqnum_elements() && sbvec[k])
         tmpit = eqgraph[it]->seqget_element(k);
      if (it2 == tmpit){
         ejoin=(Itemset*)1;
         ljoin=(Itemset*)1;
         k++;
         i++;
         if (it2 == DBASE_MAXITEM+1) continue;
      }
      else if (it2 < tmpit){
         ejoin=(Itemset*)1;
         i++;
      }
      else{
         ljoin=(Itemset*)1;
         k++;
         it2 = tmpit;
      }
      //cout << "JOIN " << it << " " << it2 << " " << ejoin << " " << ljoin << endl << flush;
      supsz2 = partition_get_idxsup(it2);
      partition_read_item(interval2->array(), it2);
      if (ejoin){
         ejoin = new Itemset(2,N*u2, N+N/8+1);
         if (ejoin == NULL){
            perror("memory exceeded");
            exit(errno);
         }
      }
      else ejoin = NULL;

      if (ljoin){ 
         //ljoin = new Itemset(2,N*u2, 2*N);
         ljoin = new Itemset(2,N*u2, N+N/8+1);
         if (ljoin == NULL){
            perror("memory exceeded");
            exit(errno);
         }
      }
      else ljoin = NULL;
      //cout << "ljoin " << ljoin << " " << ejoin << " " <<
      //supsz << " " << supsz2 << " " << it << " " << it2 << endl;

      get_2NewF_intersect(ljoin, ejoin, interval->array(), interval2->array(),
                      supsz, supsz2);

      if (ljoin){
         if (ljoin->support() >= MINSUPPORT && !use_isetonly){             
            //ljoin->reallocival1();
            //ljoin->reallocival2();           
            ljoin->add_item(it);
            ljoin->add_item(it2);
            L2->prepend(ljoin);
            //cout << "LARGE ";
            //ljoin->print_seq(L2->templ());
            //NumLargeItemset[1]++;
         }
         else{
            //cout << "DELETED ";
            //ljoin->print_seq(L2->templ());
            delete ljoin;
         }
      }
      if (ejoin){
         if (ejoin->support() >= MINSUPPORT){
            //ejoin->reallocival1();
            //ejoin->reallocival2();
            ejoin->add_item(it);
            ejoin->add_item(it2);
            L2->prepend2(ejoin);
            //cout << "LARGE ";
            //ejoin->print_seq(L2->templ2());
            //NumLargeItemset[1]++;
         }
         else{
            //cout << "DELETED ";
            //ejoin->print_seq(L2->templ2());
            delete ejoin;
         }
      }
   }

   //cout << "MEMEXTEND " << it << " " << MEMUSED << endl;
   seconds(t2);
   L2ISECTTIME += t2-t1;
   return L2;
}

void delete_eq_list(Lists<Eqclass *> *eqlist)
{
   ListNodes<Eqclass *> *eqhd = eqlist->head();

   for (;eqhd; eqhd=eqhd->next()){
      delete eqhd->item()->list();
      eqhd->item()->set_list(NULL);
      delete eqhd->item();
   }
   delete eqlist;
}

void fill_join(Itemset *join, Itemset *hdr1, Itemset *hdr2)
{
   int i;
   for (i=0; i < hdr1->size(); i++){
      join->add_item((*hdr1)[i]);
   }
   join->add_item((*hdr2)[hdr2->size()-1]);
}

Itemset *prune_decision(Itemset *it1, Itemset *it2, int ptempl, int jflg, int LR)
{
   int l1 = (*it1)[it1->size()-1];
   int l2 = (*it2)[it2->size()-1];
   if (use_hash && (it1->size() > 2)){
      int i,j,k;
      unsigned int bit, ttpl;
      int nsz;
      FreqIt fit(it1->size(), 0);

      //skip the first two subsets (or omit the last two elements)
      nsz = it1->size()-2;

      //cout << "PTEMPL " << ptempl << " " << jflg << endl;
         //cout << *it1;
      //cout << *it2;
      for (i=nsz; i >= 0; i--){
         k=0;
         //form new subset template
         if (i == 0) ttpl = SETBIT(ptempl,0,nsz+1);
         else{
            ttpl = 0;
            for (j=0; j < i-1; j++){
               bit = GETBIT(ptempl,nsz-j+1);
               ttpl = SETBIT(ttpl, bit, nsz-j);
            }
            bit = GETBIT(ptempl, nsz-j+1);
            bit = bit || GETBIT(ptempl, nsz-j);
            ttpl = SETBIT(ttpl, bit, nsz-j);
            j+=2;
            for (; j < nsz+2; j++){
               bit = GETBIT(ptempl,nsz-j+1);
               ttpl = SETBIT(ttpl, bit, nsz-j+1);
            }
         }
         //form new subset by omitting the i-th item
         for (j=0; j < nsz+1; j++){
            if (j != i){
               fit.seq[k++] = (*it1)[j];
            }
         }
         fit.seq[k++] = l1;
         fit.seq[k++] = l2;
         fit.templ = ttpl;

         //if ((*it1)[0] == 101)
         //   cout << "SEARCH " << fit;
         if (fit.seq[0] == (*it1)[0] && !recursive){
            //elements should be in current class
            if (FreqArrayPos > 0){
               if (!EqGrNode::bsearch(0,FreqArrayPos-1,FreqArray,
                                      fit, recursive)){
                  //if ((*it1)[0] == 101) cout << "NOTF" <<endl;
                  //print_freqary();
                  return NULL;
               }
               //if ((*it1)[0] == 101) cout << "FOUND" <<endl;
            }
            else return NULL;
         }
         else if (fit.seq[0] > (*it1)[0]){
            // class must already have been processed, otherwise we can't prune
            if (!eqgraph[fit.seq[0]]->find_freqarray(fit, recursive)){
               //if ((*it1)[0] == 101) cout << "NOTF" <<endl;
               return NULL;
            }
            //if ((*it1)[0] == 101) cout << "FOUND" <<endl;
         }
      }
   }
   else if (it1->size() == 2){
      if (eqgraph[l1]){
         if (jflg == LJOIN || jflg == MJOIN){
            if (!eqgraph[l1]->seqfind(l2))
               return NULL;
         }
         else{
            if (!eqgraph[l1]->find(l2))
               return NULL;
         }
      }
      else return NULL;
      //cout << "FOUND " << endl;
   }
   return (Itemset *)1;
}



void insert_freqarray(Lists<Eqclass *> *LargeL)
{
   //insert frequent itemsets into hash table
   ListNodes<Eqclass *> *chd;
   ListNodes<Itemset *> *hdr1, *hdr2;
   Eqclass *cluster;

   chd = LargeL->head();
   for (; chd; chd = chd->next()){
      cluster = chd->item();
      hdr1 = cluster->list()->head();
      for (; hdr1; hdr1=hdr1->next()){
         add_freq(hdr1->item(), cluster->templ());
         //hdr1->item()->print_seq(cluster->templ());
      }
      hdr2 = cluster->list2()->head();
      for (; hdr2; hdr2=hdr2->next()){
         add_freq(hdr2->item(), cluster->templ2());
         //hdr2->item()->print_seq(cluster->templ2());
      }
   }
}

void process_cluster_list1(ListNodes<Itemset *> *hdr1,
                           Lists<Itemset *> *cluster1,
                           Lists<Itemset *> *cluster2, Lists<Eqclass *> *LargeL,
                           int iter, int eqtype, Eqclass *parent)
{
   ListNodes<Itemset *> *hdr2;
   Eqclass *EQ = new Eqclass(iter-1,eqtype);
   if (EQ == NULL){
      perror("memory exceeded");
      exit(errno);
   }
   fill_seq_template(EQ, parent, 2);
   //int first;
   Itemset *ljoin, *ejoin, *mjoin;
   int lsup, esup, msup;
   //cout << "BEG CLUSERT 1 : " << MEMUSED << endl;

   //first = 1;
   hdr2 = cluster2->head();
   for (; hdr2; hdr2=hdr2->next()){
      //ljoin = (Itemset *)1;
      ljoin = prune_decision(hdr1->item(), hdr2->item(), EQ->templ(), LJOIN,2);
      ejoin = NULL;
      mjoin = NULL;
      lsup = esup = msup = 0;
      //cout << "process 1 0 0" << endl;
      if (ljoin || ejoin || mjoin) 
         get_tmpNewF_intersect(ljoin, ejoin, mjoin, lsup, esup, msup,
                           hdr1->item(), hdr2->item(), iter, 2, 1);

      if (lsup >= MINSUPPORT){
         NumLargeItemset[iter-1]++;
         fill_join(ljoin, hdr1->item(), hdr2->item());
         //cout << "LARGE ";
         if (print_seq) ljoin->print_seq(EQ->templ());
         EQ->append(ljoin);
         //EQ->list()->sortedAscend(ljoin, Itemset::scmp);
         //EQ->list()->sortedAscend(ljoin, Itemset::supportcmp);
      }
   }

   hdr2 = cluster1->head();
   for (; hdr2 != hdr1; hdr2=hdr2->next()){
      //ejoin = (Itemset *)1;
      ejoin = prune_decision(hdr1->item(), hdr2->item(), EQ->templ2(), EJOIN,2);
      ljoin = NULL;
      mjoin = NULL;
      lsup = esup = msup = 0;
      //cout << "process 0 1 0" << endl;
      if (ljoin || ejoin || mjoin)
         get_tmpNewF_intersect(ljoin, ejoin, mjoin, lsup, esup, msup,
                           hdr1->item(), hdr2->item(), iter, 1, 1);

      //cout << "AFT JOIN " << MEMUSED << endl;
      if (esup >= MINSUPPORT){
         NumLargeItemset[iter-1]++;
         fill_join(ejoin, hdr1->item(), hdr2->item());
         //cout << "LARGE ";
         if (print_seq) ejoin->print_seq(EQ->templ2());
         EQ->append2(ejoin);
         //EQ->list2()->sortedAscend(ejoin, Itemset::scmp);
         //EQ->list2()->sortedAscend(ejoin, Itemset::supportcmp);
      }
   }

   if (EQ){
      if ((EQ->list()->size() > 0) || (EQ->list2()->size() > 0)){
         if (recursive){
            //if (use_hash) insert_freqarray(EQ);
            process_cluster1(EQ, NULL, iter+1);
            delete EQ;
         }
         else LargeL->append(EQ);
      }
      else{
         //   if (use_hash && EQ->list2()->size() == 1)
         //      add_freq(EQ->list2()->head()->item(), EQ->templ2());
         delete EQ;
         EQ = NULL;
      }
   }
   //cout << "END CLUSTER1 : " << MEMUSED << endl;
}

void process_cluster_list2(ListNodes<Itemset *> *hdr1, int i, Eqclass ** EQ,
                           Lists<Itemset *> *cluster, Lists<Eqclass *> *LargeL,
                           int iter, int eqtype, Eqclass *parent)
{
   int j;

   ListNodes<Itemset *> *hdr2;
   Itemset *ljoin, *ejoin, *mjoin;
   int lsup, esup, msup;

   //join with sequences
   hdr2 = hdr1;
   for (j=i; hdr2; j++, hdr2=hdr2->next()){
      ljoin = prune_decision(hdr1->item(), hdr2->item(), EQ[i]->templ(), LJOIN,1);
      if (hdr2 == hdr1){
         ejoin = mjoin = NULL;
      }
      else{
         ejoin = prune_decision(hdr2->item(), hdr1->item(), EQ[j]->templ2(), EJOIN,1);
         mjoin = prune_decision(hdr2->item(), hdr1->item(), EQ[j]->templ(), MJOIN,1);
         //ejoin = mjoin = (Itemset *)1;
      }
      //cout << "process 1 1 1" << endl;
      lsup = esup = msup = 0;
      if (ljoin || ejoin || mjoin)
         if (hdr1->item()==hdr2->item()) 
             get_tmpEqNewF_intersect(ljoin, lsup, hdr1->item(), iter, 2, 2);
         else 
             get_tmpNewF_intersect(ljoin, ejoin, mjoin, lsup, esup, msup,
                           hdr1->item(), hdr2->item(), iter, 2, 2);
      //cout << "SUPPP " << lsup << " " << esup << " " << msup << endl;
      if (lsup >= MINSUPPORT){
         NumLargeItemset[iter-1]++;
         //cout << "high every"<<endl;
         fill_join(ljoin, hdr1->item(), hdr2->item());
         //cout << "LARGE ";
         if (print_seq) ljoin->print_seq(EQ[i]->templ());
         EQ[i]->append(ljoin);
         //EQ[i]->list()->sortedAscend(ljoin, Itemset::scmp);
         //EQ[i]->list()->sortedAscend(ljoin, Itemset::supportcmp);
      }
      if (esup >= MINSUPPORT){
         NumLargeItemset[iter-1]++;
         fill_join(ejoin, hdr2->item(), hdr1->item());
         //cout << "LARGE ";
         if (print_seq) ejoin->print_seq(EQ[j]->templ2());
         EQ[j]->append2(ejoin);
         //EQ[j]->list2()->sortedAscend(ejoin, Itemset::scmp);
         //EQ[j]->list2()->sortedAscend(ejoin, Itemset::scmp);
      }
      if (msup >= MINSUPPORT){
         NumLargeItemset[iter-1]++;
         fill_join(mjoin, hdr2->item(), hdr1->item());
         //cout << "LARGE ";
         if (print_seq) mjoin->print_seq(EQ[j]->templ());
         EQ[j]->append(mjoin);
         //EQ[j]->list()->sortedAscend(mjoin, Itemset::scmp);
         //EQ[j]->list()->sortedAscend(mjoin, Itemset::supportcmp);
      }
   }
   if ((EQ[i]->list()->size() > 0) || (EQ[i]->list2()->size() > 0)){
      if (recursive){
         //if (use_hash) insert_freqarray(EQ[i]);
         process_cluster1(EQ[i],NULL, iter+1);
         delete EQ[i];
         EQ[i] = NULL;
      }
      else LargeL->append(EQ[i]);
   }
   else{
      //if (use_hash && EQ[i]->list2()->size() == 1)
      //   add_freq(EQ[i]->list2()->head()->item(), EQ[i]->templ2());
      delete EQ[i];
      EQ[i] = NULL;
   }

   //cout << "END cluster 2 : " << MEMUSED << endl;
}

void process_cluster1(Eqclass *cluster, Lists<Eqclass *> *LargeL, int iter)
{
   Eqclass **EQ=NULL;
   ListNodes<Itemset *> *hdr1, *hdr2;
   int i;

   if (cluster->list()->head()){
      EQ = new Eqclass *[cluster->list()->size()];
      if (EQ == NULL){
         perror("memory exceeded");
         exit(errno);
      }
      for (i=0; i < cluster->list()->size(); i++){
         EQ[i] = new Eqclass(iter-1,EQCTYP1);
         if (EQ[i] == NULL){
            perror("memory exceeded");
            exit(errno);
         }
         fill_seq_template(EQ[i], cluster, 1);
      }
   }

   hdr1 = cluster->list()->head();
   for (i=0; hdr1; hdr1=hdr1->next(), i++){
      //if (use_hash && iter > 3) add_freq(hdr1->item(), cluster->templ());
      process_cluster_list2(hdr1, i, EQ, cluster->list(), LargeL, iter,
                            EQCTYP1, cluster);
   }
   if (EQ) delete [] EQ;


   hdr2 = cluster->list2()->head();
   for (; hdr2; hdr2=hdr2->next()){
      //if (use_hash && iter > 3) add_freq(hdr2->item(), cluster->templ2());
      process_cluster_list1(hdr2, cluster->list2(), cluster->list(),
                            LargeL, iter, EQCTYP1, cluster);
   }

   //if (recursive) delete cluster;
   if (maxiter < iter) maxiter = iter;

}


void find_large(Eqclass *cluster, int it)
{
   Lists<Eqclass *> *LargeL, *Candidate;
   ListNodes<Eqclass *> *chd;
   int iter;
   int LargelistSum=0;
   int more;

   more = 1;
   Candidate = new Lists<Eqclass *>;
   Candidate->append(cluster);
   //cout << "MEMFIND " << it << " " << MEMUSED << endl;
   for (iter=3; more; iter++){
      LargeL = new Lists<Eqclass *>;
      chd = Candidate->head();
      for (; chd; chd=chd->next()){
         //cout << "EQCLASS ";
         //chd->item()->print_template();
         //chd->item()->print_list(chd->item()->list());
         //cout << "***\n";
         //chd->item()->print_list(chd->item()->list2());
         //cout << "------------------" << endl;
         process_cluster1(chd->item(), LargeL, iter);
         //cout << "BEF MEMFIND " << it << " " << MEMUSED << endl;
         //reclaim memory for this class immediately
         delete chd->item();
         //cout << "AFT MEMFIND " << it << " " << MEMUSED << endl;
         chd->set_item(NULL);
      }
      Candidate->clear();
      delete Candidate;
      //if (maxiter < iter) maxiter = iter;

      if (use_hash) insert_freqarray(LargeL);
      chd = LargeL->head();
      LargelistSum = 0;
      for (;chd; chd=chd->next()){
         LargelistSum += chd->item()->list()->size();
         if (chd->item()->list2())
            LargelistSum += chd->item()->list2()->size();
      }
      //print_freqary();
      more = (LargelistSum > 0);

      Candidate = LargeL;
      if (memtrace) mout << it << " " << MEMUSED << endl;

      if (!more) {
         LargeL->clear();
         delete LargeL;
      }
      //cout << "AFT DEL " << it << " " << MEMUSED << " " << iter << endl;
   }
   //cout << "MEMLAST " << it << " " << MEMUSED << endl;
}

Eqclass * extract_relevant_items(Eqclass *l2it, Array *cliq)
{
   int i;

   Eqclass *RL2 = new Eqclass(1, EQCTYP1);
   if (RL2 == NULL)
   {
      perror("memory exceeded : ext_class ");
      exit (errno);
   }
   //init seq pattern templates
   RL2->set_templ(1);
   RL2->set_templ2(0);

   ListNodes<Itemset *> *ln= l2it->list()->head();
   for (i=0; ln && i < cliq->size()-1; ){
      //cout << "LN " << (*ln->item())[1] << " " << (*cliq)[i+1] << endl;
      if ((*ln->item())[1] == (*cliq)[i+1]){
         RL2->append(ln->item());
         i++;
      }
      ln = ln->next();
   }
   ln= l2it->list2()->head();
   for (i=0; ln && i < cliq->size()-1; ){
      if ((*ln->item())[1] == (*cliq)[i+1]){
         RL2->append2(ln->item());
         i++;
      }
      ln = ln->next();
   }
   return RL2;
}

void process_class(int it)
{

   interval  = new Array(maxitemsup);
   interval2 = new Array(maxitemsup);

   //form 2-itemsets from ext disk
   Eqclass *large2it = get_ext_eqclass(it);

   delete interval;
   delete interval2;

   if (large2it == NULL) return;

   interv1 = new GArray(N*u2); 
   interv2 = new GArray(N*u2);
   interv3 = new GArray(N*u2);

   inter  = new GArray(N+N/8+1);
   inter2 = new GArray(N+N/8+1);
   inter3 = new GArray(N+N/8+1);

   if (memtrace) mout << it << " " << MEMUSED << endl;
   Eqclass *l2cliq;

   if (use_clique){
      ListNodes<Array *> *clhd = eqgraph[it]->clique()->head();
      //process each clique
      for (;clhd; clhd=clhd->next()){
         //cout << it << " processing clique " << *clhd->item() << endl;
         //construct large k-itemsets, k > 2
         l2cliq = extract_relevant_items(large2it, clhd->item());
         find_large(l2cliq, it);
      }
   }
   else{
      if (recursive){
         process_cluster1(large2it, NULL, 3);
         delete large2it;
      }
      else find_large(large2it, it);
   }
   
   delete interv1;
   delete interv2;
   delete interv3;

   delete inter;
   delete inter2;
   delete inter3;
}

void newSeq()
{
   int i,j;

   if (use_hash)
      FreqArray = (FreqIt **) malloc (FreqArraySz*sizeof(FreqIt*));
   //form large itemsets for each eqclass
   if (use_ascending != -2){
      if (use_ascending == -1){
         for (i=0; i < DBASE_MAXITEM; i++)
            if (eqgraph[i]){
               if (memtrace) mout << i << " " << MEMUSED << endl;
               process_class(i);
               if (memtrace) mout << i << " " << MEMUSED << endl;
            }
      }
      else if (eqgraph[use_ascending])
         process_class(use_ascending);
   }
   else{
      for (i=DBASE_MAXITEM-1; i >= 0; i--){
         if (eqgraph[i]){
            if (memtrace) mout << i << " " << MEMUSED << endl;
            //cout << "PROCESSS ITEM " << i << endl << flush;
            if (use_hash) FreqArrayPos = 0;
            process_class(i);
            if (use_hash){
               if (FreqArrayPos > 0){
                  //cout << "FREQUENT ARRAY3" << endl;
                  FreqIt **fit = new FreqIt *[FreqArrayPos];
                  for (j=0; j < FreqArrayPos; j++){
                     fit[j] = FreqArray[j];
                     //cout << *fit[j];
                  }
                  eqgraph[i]->set_freqarray(fit, FreqArrayPos);
               }
            }
            //cout << " -------- " << endl;
            if (memtrace) mout << i << " " << MEMUSED << endl;
         }
      }
   }
}


void read_files()
{
   int i;
   NumLargeItemset = new int [(int)(DBASE_AVG_TRANS_SZ*30)];

   bzero((char *)NumLargeItemset, sizeof(int)*((int)(DBASE_AVG_TRANS_SZ*30)));
   eqgraph = new EqGrNode *[DBASE_MAXITEM];
   bzero((char *)eqgraph, DBASE_MAXITEM*sizeof(EqGrNode *));

   double t1,t2;
   if (ext_l2_pass){
      seconds(t1);
      NumLargeItemset[0] = make_l1_pass(1);
      seconds(t2);
      EXTL1TIME = t2-t1;
      NumLargeItemset[1] = make_l2_pass();
      //cout << "L2 " <<  NumLargeItemset[1] <<endl;
      seconds(t1);
      EXTL2TIME = t1-t2;
   }
   else{
      seconds(t1);
      NumLargeItemset[0] = make_l1_pass(0);
      seconds(t2);
      EXTL1TIME = t2-t1;
      NumLargeItemset[1] = get_file_l2(it2f, seqf);
      seconds(t1);
      EXTL2TIME = t1-t2;
   }
   //cout << NumLargeItemset[0] << "LARGE 1 ITEMS\n";
   maxitemsup = 0;
   int sup;
   for (i=0; i < DBASE_MAXITEM; i++) {
      sup = partition_get_idxsup(i);
      if (maxitemsup < sup) maxitemsup  = sup;
   }
   //cout << "MAXITEMSUP " << maxitemsup << endl;
   //cout << "MAXEQSZIE " << maxeqsize << " " << t2-ts << endl;
}


//Get the greatest common divisor of two integers
int GCD(int x,int y)
{
   int rem;
   while(y!=0){
      rem=x%y;
      x=y;
      y=rem;
   }
   return x;
}

//Preproccessing: construct permanent auxillary data structures 
//such as GCD Matrix, Support Array, and Reference Array, common factors, etc.
void Build_GCD_Matrix()
{
   int P[8]={2,3,5,7,11,13,17,19};
   int P_Star[256]={0,2,3,5,7,11,13,17,19}; 
   int i,j,l,m,n,o,te,tem,temp,temp1,temp2,k=9;
   int oo, temp11;
   int M; 
   unsigned char FFGreater[9]={0};    
   unsigned char P_Factors[9][256]={{0},{0,1,2,3,4,5,6,7,8}, {0},{0},{0},{0},{0},{0},{0}};
 
   for(i=0;i<=6;i++)
      for (j=i+1; j<=7;j++){
	 temp2=P[i]*P[j];
	 P_Star[k] = temp2;
         P_Factors[1][k]=i+1;
         P_Factors[2][k]=j+1;
	 Sup_Array[k]=2;
	 //Ref_Array[temp2]=k++;
         Ref_Array[temp2%1966]=k;
         P_Rem[1][k]=Ref_Array[(temp2/P[i])%1966];
         P_Rem[2][k]=Ref_Array[(temp2/P[j])%1966];
         k++;
      }	
   for ( i=0; i<=5;i++)
      for( j = i+1; j <= 6; j++) {
	 temp = P[i] * P[j];
	 for (l= j+1; l<=7; l++) {
	    temp2=temp*P[l];
	    P_Star[k]=temp2; 
            P_Factors[1][k]=i+1;
            P_Factors[2][k]=j+1;
            P_Factors[3][k]=l+1;
	    Sup_Array[k]=3;
	    //Ref_Array[temp2]=k++;
            Ref_Array[temp2%1966]=k;
            P_Rem[1][k]=Ref_Array[(temp2/P[i])%1966];
            P_Rem[2][k]=Ref_Array[(temp2/P[j])%1966];
            P_Rem[3][k]=Ref_Array[(temp2/P[l])%1966];
            k++;
	 } 
      }
   for (i=0; i<= 4; i++) 
      for (j=i+1; j <= 5; j++){
	 temp = P[i]*P[j];
	 for (l = j+1; l<=6; l++){
	    temp1 = temp*P[l];
	    for(m=l+1; m<=7; m++) {
		temp2=temp1*P[m];
		P_Star[k]=temp2;
                P_Factors[1][k]=i+1;
                P_Factors[2][k]=j+1;
                P_Factors[3][k]=l+1;
                P_Factors[4][k]=m+1;
		Sup_Array[k]=4; 
		//Ref_Array[temp2/100000 + temp2%100000]=k++;
                Ref_Array[temp2%1966]=k;
                P_Rem[1][k]=Ref_Array[(temp2/P[i])%1966];
                P_Rem[2][k]=Ref_Array[(temp2/P[j])%1966];
                P_Rem[3][k]=Ref_Array[(temp2/P[l])%1966];
                P_Rem[4][k]=Ref_Array[(temp2/P[m])%1966];
                k++;
	    }

	 }
      }
   for (i = 0; i <=3; i++) 
      for (j=i+1; j<=4; j++){
	 tem = P[i]*P[j];
	 for (l= j+1; l<=5; l++){
	    temp = tem*P[l];
	    for (m=l+1; m<=6; m++) {
		temp1 = temp*P[m];
		for(n=m+1; n<=7; n++) {
		   temp2=temp1*P[n];
		   P_Star[k]=temp2; 
                   P_Factors[1][k]=i+1;
                   P_Factors[2][k]=j+1;
                   P_Factors[3][k]=l+1;
                   P_Factors[4][k]=m+1;
                   P_Factors[5][k]=n+1;
		   Sup_Array[k]=5; 
		   //Ref_Array[temp2/100000 + temp2%100000] = k++;
                   Ref_Array[temp2%1966]=k;
                   P_Rem[1][k]=Ref_Array[(temp2/P[i])%1966];
                   P_Rem[2][k]=Ref_Array[(temp2/P[j])%1966];
                   P_Rem[3][k]=Ref_Array[(temp2/P[l])%1966];
                   P_Rem[4][k]=Ref_Array[(temp2/P[m])%1966];
                   P_Rem[5][k]=Ref_Array[(temp2/P[n])%1966];
                   k++;
		}
	    }
	 }
      }

   for (i = 0; i <=2; i++) 
      for (j=i+1; j<=3; j++){
	 te = P[i]*P[j];
	 for (l= j+1; l<=4; l++){
	    tem = te * P[l];
	    for (m=l+1; m<=5; m++){
	       temp = tem * P[m];
	       for (n = m+1; n <=6; n++ ){
	          temp1 = temp*P[n];
		  for(o=n+1; o<=7; o++) {
		     temp2=temp1*P[o]; 
		     P_Star[k]=temp2;
                     P_Factors[1][k]=i+1;
                     P_Factors[2][k]=j+1;
                     P_Factors[3][k]=l+1;
                     P_Factors[4][k]=m+1;
                     P_Factors[5][k]=n+1;
                     P_Factors[6][k]=o+1;
		     Sup_Array[k]=6; 
		     //Ref_Array[temp2/100000 + temp2%100000] = k++;
                     Ref_Array[temp2%1966]=k;
                     P_Rem[1][k]=Ref_Array[(temp2/P[i])%1966];
                     P_Rem[2][k]=Ref_Array[(temp2/P[j])%1966];
                     P_Rem[3][k]=Ref_Array[(temp2/P[l])%1966];
                     P_Rem[4][k]=Ref_Array[(temp2/P[m])%1966];
                     P_Rem[5][k]=Ref_Array[(temp2/P[n])%1966];
                     P_Rem[6][k]=Ref_Array[(temp2/P[o])%1966];
                     k++;
		  }
	       }
	    }
         }
      }


   for (i = 0; i <=1; i++) 
      for (j=i+1; j<=2; j++){
	 te = P[i]*P[j];
	 for (l= j+1; l<=3; l++){
	    tem = te * P[l];
	    for (m=l+1; m<=4; m++){
	       temp = tem * P[m];
	       for (n = m+1; n <=5; n++ ){
	          temp1 = temp*P[n];
		  for(o=n+1; o<=6; o++) {
                   temp11 = temp1*P[o];
                   for(oo=o+1; oo<=7; oo++) {
		     temp2=temp11*P[oo]; 
		     P_Star[k]=temp2;
                     P_Factors[1][k]=i+1;
                     P_Factors[2][k]=j+1;
                     P_Factors[3][k]=l+1;
                     P_Factors[4][k]=m+1;
                     P_Factors[5][k]=n+1;
                     P_Factors[6][k]=o+1;
                     P_Factors[7][k]=oo+1;
		     Sup_Array[k]=7; 
                     Ref_Array[temp2%1966]=k;
                     P_Rem[1][k]=Ref_Array[(temp2/P[i])%1966];
                     P_Rem[2][k]=Ref_Array[(temp2/P[j])%1966];
                     P_Rem[3][k]=Ref_Array[(temp2/P[l])%1966];
                     P_Rem[4][k]=Ref_Array[(temp2/P[m])%1966];
                     P_Rem[5][k]=Ref_Array[(temp2/P[n])%1966];
                     P_Rem[6][k]=Ref_Array[(temp2/P[o])%1966];
                     P_Rem[7][k]=Ref_Array[(temp2/P[oo])%1966];
                     k++;
		  }
	       }
	    }
         }
      }
    }


   temp = P[0]*P[1]*P[2]*P[3]*P[4]*P[5]*P[6]*P[7];
/*
   for (i =0; i<=7; i++)  { 
      temp2 = temp/P[i];
      P_Star[k]=temp2; 
      Ref_Array[temp2%1966]=k;
      for (j=0; j<i; j++){ 
          P_Factors[j+1][k]=j+1; 
          P_Rem[j+1][k]=Ref_Array[(temp2/P[j])%1966]; //for seq mining
      }
      for (j=i+1; j<=7; j++){ 
          P_Factors[j][k]=j+1; 
          P_Rem[j][k]=Ref_Array[(temp2/P[j])%1966];//for seq mining
      }
      Sup_Array[k]=7;
      Ref_Array[temp2%1966] = k++;
   } 
*/
   P_Star[k]=temp;
   for (i=1; i<=8;i++)  P_Factors[i][k]=i;
   for (i=0; i<=7;i++)  P_Rem[i+1][k]=Ref_Array[(temp/P[i])%1966]; 
   Sup_Array[k]=8;
   Ref_Array[temp%1966]=k;

//Get multiplications that are greater than the first factors(for sequence mining) 
   for (i=0; i<= 6; i++){
      M =1; 
      for (j=i+1; j<=7; j++) M *= P[j];
      FFGreater[i+1]= Ref_Array[M%1966];
   }

//Build Greatest common divisor matrix and AspArray
   //for (i=0; i<=256; i++){ 
     // GCD_Matrix[i]= new unsigned char [257]; 
     // (GCD_Matrix[i])[0]=0;
   //}
   for(i=1;i<=255;i++) {
      AspArray[i] = FFGreater[P_Factors[1][i]]; //for seq mining
      for(j=1; j <= 255; j++) { 
	 temp=GCD(P_Star[i], P_Star[j]); 
         //(GCD_Matrix[i])[j] = Ref_Array[temp%1966];
          GCD_Matrix[i][j] = Ref_Array[temp%1966];
      }
   }
////////////////////////////////////////////////////////
   
   //Build common factors matrix
   int f;
   for(i=1; i<=255; i++) {
      //Common_Factors[i]= new  unsigned char *[i+1];
      //for(j=1; j<= i; j++) {  
      //for(j=i; j<= 255; j++) { 
      for(j=1; j<= 255; j++) {       
         n= GCD_Matrix[i][j];
         //Common_Factors[i][j]= new GArray (2*Sup_Array[n]+1); 
         //(Common_Factors[i])[j]= new unsigned char [2*Sup_Array[n]+1]; 
         Common_Factors[i][j]= new unsigned char [2*Sup_Array[n]+1]; 
         f=0; m=0; l=0;
         for (o=1; o<=Sup_Array[n]; o++){ 
            for (f=f+1; f<=Sup_Array[i]; f++) 
               if (P_Factors[o][n] == P_Factors[f][i]) break; 
            for (m=m+1; m<=Sup_Array[j]; m++) 
               if (P_Factors[o][n] == P_Factors[m][j]) break; 
            //(Common_Factors[i][j])->optadd(f);
            //(Common_Factors[i][j])->optadd(m);
            (Common_Factors[i][j])[l++]=f;
            (Common_Factors[i][j])[l++]=m;
         }
      }
   }

/*
   for(i=1;i<=255;i++)
      for(j=1; j <= 255; j++) 
           Temp[i][j]=GCD_Matrix[AspArray[i]][j];
  
*/

/*
for(i=1;i<=255;i++) {
      cout << "i ="<<i<<endl;
      for(j=1; j <= 8; j++) 
         cout << (int)P_Rem[j][i]<<"\t";
      cin.get();
}

for(i=1; i<=8; i++){ 
      for(j=1; j<= 255; j++) {
         //cout << "i ="<<i << " j="<<j<<endl; 
         cout << (int)(P_Factors[i][j])<<" ";
}
cin.get();
}

for(i=1;i<=255;i++) {
      for(j=1; j <= 255; j++) { 
         cout << "i ="<<i << " j="<<j<<endl;
         for (k=0; k< 2*Sup_Array[GCD_Matrix[i][j]];k++)
            cout << (int)(Common_Factors[i][j])[k] <<"\t";
         //cout<<endl;
        
         
}
 cin.get();
}

for(i=1;i<=255;i++) 
      for(j=1; j <= 255; j++){cout << "i ="<<i << " j="<<j<<endl;
 cout<<*(Common_Factors[i][j])<<endl;
}
*/
}

int main(int argc, char **argv)
{
   int i;
   double ts, te;
   double t1,t2;

   seconds(ts);
   Build_GCD_Matrix();
   //cout << "BEGIN MEM " << MEMUSED << endl;
   parse_args(argc, argv);

   partition_alloc(dataf, idxf);
   read_files();
   //cout << "AFTER READFILE " << MEMUSED << endl;
   seconds(t1);
   newSeq();
   seconds(t2);
   double FKtime = t2-t1;
   //print_freqary();
   //cout << "AFTER SEQ " << MEMUSED << endl;
   seconds(te);
   if ((out = fopen("summary.out", "a+")) == NULL){
      perror("can't open summary file");
      exit(errno);
   }
   fprintf (out, "PRIMESPADE ");
   if (use_hash) fprintf (out, "USEHASH ");
   fprintf(out, "%s %f %d %d %f %d (", dataf, MINSUP_PER, MINSUPPORT,
           num_intersect, L2ISECTTIME, L2pruning);
   for (i=0; i < maxiter; i++){
      fprintf(out, "%d ", NumLargeItemset[i]);
      cout << "ITER " <<  i+1 << " " << NumLargeItemset[i] << endl;
   }
   cout << "Total elapsed time " << te-ts 
        << ", NumIntersect " << num_intersect 
        <<  " " << EXTL2TIME << endl;

   fprintf(out, ") %f %f %f %f %d %d\n", EXTL1TIME, EXTL2TIME, FKtime,
           te-ts, u2, N);
   fclose(out);
   partition_dealloc();

   for (i=0; i < DBASE_MAXITEM; i++)
      if (eqgraph[i]) delete eqgraph[i];
   delete [] eqgraph;

   if (memtrace){
      mout <<  MEMUSED << endl;
      mout.close();
   }
   //cout << "LAST MEM " << MEMUSED << endl;
   exit(0);
}


