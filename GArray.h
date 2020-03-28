#ifndef __GARRAY_H
#define __GARRAY_H
#include <stdlib.h>
#include <iostream>
#include <sys/types.h>
#include <malloc.h>

using namespace std;

extern long MEMUSED;

class GArray {
protected:   
   unsigned char *theArray;
   unsigned int theSize;
   unsigned int totSize;
   //unsigned int theIncr;
public:
   
   //GArray (int sz, int incr);
   GArray(int sz);
   ~GArray();
   
   int subsequence(GArray * ar);
   //void add (int, unsigned int);
   void add_ext(unsigned char val, int off, unsigned char *ary)
   {
      ary[off+theSize] = val;
      theSize++;
   }
   
   unsigned char operator [] (unsigned int index)
   {
      return theArray[index];
   };
   
   void setitem(int pos, unsigned char val){
      theArray[pos] = val;
   };
   
   int totsize()
   {
      return totSize;
   }
   void set_totsize(int sz){
      totSize = sz;
   }
   void set_size(int sz){
      theSize = sz;
   }
   void reset()
   {
      theSize = 0;
   }

   unsigned char *array()
   {
      return theArray;
   }
   void set_array(unsigned char *ary){
      theArray = ary;
   }
   //int subsequence(Array&);
   //int compare(Array&);
   friend ostream& operator << (ostream& outputStream, GArray& arr);
   static int Arraycompare(void * iset1, void *iset2)
   {
      GArray *it1 = (GArray *) iset1;
      GArray *it2 = (GArray *) iset2;
      return it1->compare(*it2);
   }
   int compare(GArray& ar2);

   unsigned char item (unsigned int index) 
   {
      return theArray[index];
   }
   
   unsigned int size() 
   {
      return theSize; 
   }
   
   void realloc(int newsz)
   {
      MEMUSED -= totSize*sizeof(unsigned char);
      totSize = newsz;
      //theArray = new unsigned char [totSize];
      theArray = (unsigned char *)::realloc(theArray, totSize*sizeof(unsigned char));
      if (theArray == NULL){
         cout << "MEMORY EXCEEDED\n";
         exit(-1);
      }
      MEMUSED += totSize*sizeof(unsigned char);
      //cout << "high real"<<endl;
   }
   void optadd(unsigned char item)
   {
      theArray[theSize++] = item;
   }
   void add (unsigned char item)
   {
      if (theSize+1 > totSize){
         //totSize = (int) (totSize*1.5);
         //cout << " " << MEMUSED;
         realloc((unsigned char)(totSize*1.5));
         //theArray = (int *)realloc(theArray, totSize*sizeof(int));
         //if (theArray == NULL){
         //   cout << "MEMORY EXCEEDED\n";
         //   exit(-1);
         //}
         //MEMUSED += totSize*sizeof(int);
         //cout << " " << MEMUSED << endl;
      }
      theArray[theSize] = item;
      theSize++;
   }
};
#endif //__ARRAY_H


