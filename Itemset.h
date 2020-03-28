#ifndef __ITEMSET_H
#define __ITEMSET_H
#include <iostream>
#include <stdio.h>
#include <errno.h>
#include "GArray.h"
#include "Array.h"
#include "Lists.h"
//#include "Bitvec.h"

using namespace std;

//#define SETBIT(a,b) ((a) |= (1 << (b)))
//#define UNSETBIT(a,b) ((a) &= ~(1 << (b)))
#define SETBIT(a,v,b)  (((v) != 0) ? ((a) | (01 << (b))): ((a) & ~(01 << (b))))
#define GETBIT(a,b) ((a) & (01 << (b)))


class Itemset{
protected:
   Array *theItemset;
   GArray *theIval1; //added by karam
   GArray *theIval2; //added by karam
   //Array *theIval2; //added by karam
   int theSupport;
   int theOnes;
   //int mm;

public:
   Itemset(int it_sz, int ival1_sz, int ival2_sz);
   //Itemset(int it_sz, int ival_sz);
   ~Itemset();

   friend ostream& operator << (ostream& outputStream, Itemset& itemset);
   void intersect_neighbors(Itemset *it1, Itemset *it2);
   int compare(Itemset& ar2, int len);
   int compare(Itemset& ar2);
   int compare(Array& ar2, int len);
   int compare(Itemset& ar2, int len, unsigned int);
   int subsequence(Itemset * ar);
   void print_seq(int itempl);

   GArray *ival1(){return theIval1;}

   GArray *ival2(){return theIval2;}
   //Array *ival2(){return theIval2;}

   unsigned char ival1(int pos){return (*theIval1)[pos];}
   unsigned char ival2(int pos){return (*theIval2)[pos];}
   //int ival2(int pos){return (*theIval2)[pos];}

   int ival1size(){return theIval1->size();}
   int ival2size(){return theIval2->size();}
   int ival1totsize(){return theIval1->totsize();}
   int ival2totsize(){return theIval2->totsize();}
   

   void add_ival1(unsigned char it){theIval1->add(it);}
   void add_ival2(unsigned char it){theIval2->add(it);}
   //void add_ival2(int it){theIval2->add(it);}

   void reallocival1(){theIval1->realloc(ival1size());}
   void reallocival2(){theIval2->realloc(ival2size());}
   
   int operator [] (int pos){return (*theItemset)[pos];};
   
   int ones(){ return theOnes;};
   void set_ones(int sup){theOnes = sup;}
   void increment_ones(){theOnes++;};


   int item(int pos){return (*theItemset)[pos];};  
   void setitem(int pos, int val){theItemset->setitem(pos, val);};  
   void set_itemset (Array *ary){theItemset = ary;}  
   Array *itemset(){return theItemset;};
   void add_item(int val){theItemset->add(val);};   
   int size(){return theItemset->size();};

   int support(){ return theSupport;};
   void set_support(int sup){theSupport = sup;}  
   void increment_support(){theSupport++;};
/*
   void setsz(int val){
   mm=val;
   };
   int getsz(){
   return mm;
   }
*/
   static int intcmp (void *it1, void *it2)
   {
      int i1 = *(int *) it1;
      int i2 = *(int *) it2;
      //printf("cmp %d %d\n", i1->theSupport, 
      if (i1 > i2) return 1;
      else if (i1 < i2) return -1;
      else return 0;
   }
   
   static int supportcmp (void *it1, void *it2)
   {
      Itemset * i1 = (Itemset *)it1;
      Itemset * i2 = (Itemset *)it2;
      //printf("cmp %d %d\n", i1->theSupport, 
      if (i1->theSupport > i2->theSupport) return 1;
      else if (i1->theSupport < i2->theSupport) return -1;
      else return 0;
   }
   

   static int scmp (void *it1, void *it2)
   {
      Itemset * i1 = (Itemset *)it1;
      Itemset * i2 = (Itemset *)it2;
      //printf("cmp %d %d\n", i1->theSupport, 
      if (i1->theIval1->size() > i2->theIval1->size()) return 1;
      else if (i1->theIval1->size() < i2->theIval1->size()) return -1;
      else return 0;
   }

   static int Itemcompare(void * iset1, void *iset2)
   {
      Itemset *it1 = (Itemset *) iset1;
      Itemset *it2 = (Itemset *) iset2;
      return it1->compare(*it2);
   }

   //assume ts1 and ts2 are of same length
   static int compare_seq(void *ts1, void *ts2, int len)
   {
      int *tseq1 = (int *)ts1;
      int *tseq2 = (int *)ts2;
      for (int i=0; i < len; i++){
         if (tseq1[i] > tseq2[i]) return 1;
         else if (tseq1[i] < tseq2[i]) return -1;
      }
      return 0;
   }
   //int find(int , int*);
   //int subsequence(Itemset &);
   //int compare(Itemset &);
};



#endif //__ITEMSET_H

