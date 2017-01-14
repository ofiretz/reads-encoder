/*

   mapper.c         
   
   Copyright (C) 2016, Ofir Etz-Hadar.
   Website: https://www.cs.bgu.ac.il/~negevcb/about.php
   
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   Date   : Aug 2016
   Author : Ofir Etz-Hadar
   
            This program uses 2BWT library to perform reference-based reads mapping and create an ready-to-compress file
            
   The following command compile the program.
   gcc  mapper.c BWT.o dictionary.o DNACount.o HSP.o HSPstatistic.o iniparser.o inistrlib.o karlin.o MemManager.o MiscUtilities.o QSufSort.o r250.o Socket.o TextConverter.o 2BWT-Interface.o -lm -o mapper
*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "2BWT-Interface.h"
#include <sys/resource.h>


// Read set records
struct Records {
    int _ref_idx;
    int _var_idx;
    int _var_idx2;
    int _flag;
    char _sub[8];
    char _sub2[8];
    int _strand;
    char _seq[300];
    int _readId;
};

// Doubly linked-list 
struct Node {
    struct Records* _records;
    struct Node* next;
    struct Node* prev;
};

struct Node* head; // Global variable - pointer to head node
struct Node* tail; // Global variable - pointer to tail node

// Creates a new Node and returns pointer to it
struct Node* CreateNode(struct Records* records) {
    struct Node* newNode = (struct Node*)malloc(sizeof(struct Node));
    newNode->_records = records;
    newNode->prev = NULL;
    newNode->next = NULL;
    return newNode;
}

// Creates a new Records and returns pointer to it
struct Records* CreateRecords(int ref_idx, int var_idx, int var_idx2, int flag, char sub[], char sub2[], int id, char seq[], int strand) {
    struct Records* newRecord = (struct Records*)malloc(sizeof(struct Records));

    newRecord->_ref_idx = ref_idx-1;
    newRecord->_var_idx = var_idx;
    newRecord->_var_idx2 = var_idx2;
    newRecord->_flag = flag;
    strcpy(newRecord->_sub, sub);
    strcpy(newRecord->_sub2, sub2);
    newRecord->_readId = id;
    strcpy(newRecord->_seq, seq);
    newRecord->_strand = strand;
    return newRecord;
}

// Inserts a Node at head of doubly linked list
void InsertAtHead(struct Records* records) {
    struct Node* newNode = CreateNode(records);
    
    if(head == NULL) {
      head = newNode;
      tail = newNode;
      return;
    }
    
    head->prev = newNode;
    newNode->next = head;
    head = newNode;
}

// Compare Records. First compare by _ref_idx then (if necessary) compare by _var_idx
int compare(const void *r1, const void *r2){
    const struct Records *e1 = *(struct Records **)r1;
    const struct Records *e2 = *(struct Records **)r2;
    if(e1->_ref_idx == e2->_ref_idx){
        if(e1->_var_idx >= e2->_var_idx){
            return 1;
        }else{
            return -1;
        }
    }else if(e1->_ref_idx > e2->_ref_idx){
        return 1;
    }else{
        return -1;
    }
}

// Given a read sequence producing the reverse starnd
void reverseStrand(char * rseq){
    char tmp[150];
    int i;
    int len = strlen(rseq);

    for(i = 0; i < len; i++){
        if(rseq[i] == 'A' || rseq[i] == 'a'){
            rseq[i] = 'T';
        }else if(rseq[i] == 'C' || rseq[i] == 'c'){
            rseq[i] = 'G';
        }else if(rseq[i] == 'G' || rseq[i] == 'g'){
            rseq[i] = 'C';
        }else if(rseq[i] == 'T' || rseq[i] == 't'){
            rseq[i] = 'A';
        }
    }

    strcpy(tmp, rseq);

    for(i = 0; i < 100; i++){
        rseq[i] = tmp[99-i];
    }
}

int exactMap1(char * reference, char * read, int * variationIdx, int * variationIdx2, char * variationSub, char * variationSub2, int j, char * oread){
    int n_miss = 0;
    int i,k;
    k=0;

       while(k < j) // make sure all letters are upper case
     {
        reference[k] = (toupper(reference[k]));
        k++;
     }

    if(strncmp(reference, read, j) == 0){
      return 0; // exaxt match
    }else{
      for(i = 0; i < j && n_miss <= 2; i++){
        
        if(reference[i] != read[i]){
          n_miss++;
          if(n_miss == 1){
            (*variationIdx2) = i;
            variationSub2[0] = tolower(oread[i]);
            variationSub2[1] = tolower(reference[i]);
            variationSub2[2] = '\0'; 
          }else if(n_miss == 2){
            (*variationIdx) = i;
            variationSub[0] = tolower(oread[i]);
            variationSub[1] = tolower(reference[i]);
            variationSub[2] = '\0';
          }
        }
      }
    }

    return n_miss;
}

int exactMap2(char * reference, char * read, int * variationIdx, int * variationIdx2, char * variationSub, char * variationSub2, int j, char * oread){
  int n_miss = 0;
  int i,k;
  k=0;

     while(k < j) // make sure all letters are upper case
   {
      reference[k] = (toupper(reference[k]));
      k++;
   }

  if(strncmp(reference, read, j) == 0){
    return 0; // exaxt match
  }else{
    for(i = 0; i < j && n_miss <= 2; i++){
      
      if(reference[i] != read[i]){
        n_miss++;
        if(n_miss == 1){
          (*variationIdx2) = i;
          variationSub2[0] = tolower(oread[99-i]);
          variationSub2[1] = tolower(reference[i]);
          variationSub2[2] = '\0'; 
        }else if(n_miss == 2){
          (*variationIdx) = i;
          variationSub[0] = tolower(oread[99-i]);
          variationSub[1] = tolower(reference[i]);
          variationSub[2] = '\0';
        }
      }
    }
  }

  return n_miss;
}

int main(int argc, char **argv) {
    int i,j,k,c;

    // Variables for backward and forward search
    unsigned int l,r;
        
    // Variables for result
    unsigned int offset;
    int sequenceId;
    unsigned int idx = -1;
    int foundMatch = 1;
    int exactMatch = 1;
    int variationIdx = -1;
    int variationIdx2 = -1;
    char variationSub[5]; 
    char variationSub2[5];
    int readId = 1;

    // flags, counters and temps
    char tmpvar[5];
    char str[101];
    int match_cnt = 0;
    int patternLength = 0;
    char pattern[300];
    char tmp[300];
    char pattern_cpy[300];
    char buf[300];
    int strand;
    char pattern_reverse[300];
    char tmp_reverse[300];
    int n_miss = 0;
    
    if(argc < 4){
        printf("Wrong: missing arguments\n\n");
        return;
    }

    // Load up the reference index 
    printf("Loading index ... "); 
    fflush(stdout);
    Idx2BWT * idx2BWT = BWTLoad2BWT(argv[1],".sa");
    printf("DONE\n\n"); 

    // Load up the reverse reference index
    printf("Loading reverse index ... "); 
    fflush(stdout);
    Idx2BWT * revidx2BWT = BWTLoad2BWT("chr20_rev.fa.index",".sa"); // TO DO: get the reverse name as command line arg
    printf("DONE\n\n"); 
    
    FILE *reads_in;
    reads_in = fopen(argv[2], "rt");
    
    if(!reads_in)
        return 1;
    
    FILE *reads_out_fn;
    reads_out_fn = fopen("FN.fastq", "wb");
    
    if(!reads_out_fn)
        return 1;

    char reference[63025521]; // TO DO: 
    char rev_reference[63025521];

    FILE *reference_file;
    reference_file = fopen(argv[4], "rt");

    fscanf(reference_file, "%s", reference);
    fclose(reference_file);

    reference_file = fopen("refe.rev", "rt"); // TO DO:
    fscanf(reference_file, "%s", rev_reference);
    fclose(reference_file);
   
    clock_t start = clock();
    clock_t end;
      
    // Get one read from file at each loop
    while(fgets(buf, 300, reads_in) != NULL){
        patternLength = strlen(buf);
        strcpy(pattern, buf);
        strcpy(tmp, buf);
        foundMatch = 1;
        strand = 0;

      if(patternLength > 0 && pattern[patternLength-1] == '\n'){
          pattern[patternLength-1] = '\0';
          tmp[patternLength-1] = '\0';
      }
   
      if(patternLength != atoi(argv[3])){ // Indel error type
          fputs(pattern, reads_out_fn);

      }else{
  
          // Convert the pattern into 2BWT recognised coding scheme
          BWTConvertPattern(idx2BWT, pattern, patternLength, pattern);
  
          // The following performs an up to 2 mismatchs search of the pattern
          // ****

          BWTSARangeInitial(idx2BWT, pattern[patternLength-1], &l, &r); // get initial suffix range

          for(i = patternLength-2; i >= 0; i--) {
              if(l >= r) break;

              BWTSARangeBackward(idx2BWT, pattern[i], &l, &r);
          }

          if(l <= r && i < 0){ // Found without a mismatch
              foundMatch = 0;
              exactMatch = 0;
              BWTRetrievePositionFromSAIndex(idx2BWT, l, &sequenceId, &offset); // get exact offset on the reference

          }else if(l == r){ // suffix range interval == 1 then go for string to string matching
              BWTRetrievePositionFromSAIndex(idx2BWT, l, &sequenceId, &offset);
              strncpy(str, reference+offset-(i+2), 100);
              offset = offset-(i+1); // update offset 
              exactMatch = exactMap1(str, tmp, &variationIdx, &variationIdx2, variationSub, variationSub2, i+1, tmp);
  
              if(exactMatch < 3 ){ // found with up to 2 mismatches
                  foundMatch = 0;
              }

          }

          if(foundMatch == 1 && l > r){ // forward search failed - 5th end strand
              BWTSARangeInitial(revidx2BWT, pattern[0], &l, &r); 

              for(i = 1; i < patternLength; i++){
                  if(l >= r) break;

                  BWTSARangeBackward(revidx2BWT, pattern[i], &l, &r);
              }

              if(l == r){ 
                  BWTRetrievePositionFromSAIndex(revidx2BWT, l, &sequenceId, &offset);
                  strncpy(str, reference+(63025519-(offset-(99-i+2))-99), 100);
                  offset = 63025519-(offset-(99-i+1))-99+2;

                  exactMatch = exactMap1(str, tmp, &variationIdx, &variationIdx2, variationSub, variationSub2, 100, tmp);
                  if(exactMatch < 3 ){
                      foundMatch = 0;
                  }
              }
          }

          if(foundMatch == 1){ // search for the reverse strand
              strand = 1;
              strcpy(pattern_reverse, tmp);
              reverseStrand(pattern_reverse);
              strcpy(tmp_reverse, pattern_reverse);
              // Convert the pattern into 2BWT recognised coding scheme
              BWTConvertPattern(idx2BWT, pattern_reverse, patternLength, pattern_reverse);
          
              // The following performs a exact matching search of the pattern
              // =============================================================
          
              BWTSARangeInitial(idx2BWT, pattern_reverse[patternLength-1], &l, &r);

              for(i = patternLength-2; i >= 0; i--) {
                  if(l >= r) break;

                  BWTSARangeBackward(idx2BWT, pattern_reverse[i], &l, &r);
              }

      if(l <= r && i < 0){ // Found without a mismatch
          foundMatch = 0;
          exactMatch = 0;
          BWTRetrievePositionFromSAIndex(idx2BWT, l, &sequenceId, &offset);
      }else if(l == r){ 
          BWTRetrievePositionFromSAIndex(idx2BWT, l, &sequenceId, &offset);
          strncpy(str, reference+offset-(i+2), 100);
          offset = offset-(i+1);//(i+1)?
          exactMatch = exactMap2(str, tmp_reverse, &variationIdx, &variationIdx2, variationSub, variationSub2, i+1, tmp);
  
          if(exactMatch < 3 ){
              foundMatch = 0;
          }
      }

      if(foundMatch == 1 && l > r){
          BWTSARangeInitial(revidx2BWT, pattern_reverse[0], &l, &r);

          for(i = 1; i < patternLength; i++) {
              if(l >= r) break;

              BWTSARangeBackward(revidx2BWT, pattern_reverse[i], &l, &r);
          }

          if(l == r){ 
              BWTRetrievePositionFromSAIndex(revidx2BWT, l, &sequenceId, &offset);
              strncpy(str, reference+(63025519-(offset-(99-i+2))-99), 100);
              offset = 63025519-(offset-(99-i+1))-99+2;//(i+1)?

              exactMatch = exactMap2(str, tmp_reverse, &variationIdx, &variationIdx2, variationSub, variationSub2, 100, tmp);
              if(exactMatch < 3 ){
                  foundMatch = 0;
              }
          }
      }
  }

  if(foundMatch == 0){
    match_cnt++;

    InsertAtHead(CreateRecords(offset, variationIdx, variationIdx2, exactMatch, variationSub, variationSub2, readId, tmp, strand));
    readId++;   
    foundMatch = 1; 

  }else{

    strcpy(pattern_cpy, buf);
    fputs(pattern_cpy, reads_out_fn); 

    
  }

      
  }
    
    }
    
    fclose(reads_in);
    fclose(reads_out_fn);

    if(match_cnt > 0){
      printf("%d reads found\n", match_cnt);
      struct Records **recordsArray = malloc(match_cnt * sizeof(struct Records *));
      //struct Records *recordsArray[match_cnt]; // Array of reads set records
      struct Records *tmp;

      for(i = 0; i < match_cnt-1; i++){
  tmp = tail->_records;
  recordsArray[i] = tmp;
  tail = tail->prev;
  free(tail->next); 
      }

      recordsArray[i] = tail->_records;
      free(tail);
      
      printf("Start sorting...\n");
      // Sort records in array
      qsort(recordsArray, match_cnt, sizeof(struct Records *), compare);
      printf("Finished sorting\n");

      FILE *positions_file;
      positions_file = fopen("positions.csv", "wb");
      
      if(!positions_file)
  return 1;
      
      int curr_pos = recordsArray[0]->_ref_idx;
      int prev_pos = curr_pos;
      if(recordsArray[0]->_flag == 2){
      fprintf(positions_file, "%d\t%d\t%d\t%s\t[(%d,'%s','%s'), (%d,'%s','%s')]\n", recordsArray[0]->_readId, curr_pos, recordsArray[0]->_strand, recordsArray[0]->_seq, recordsArray[0]->_var_idx2, "S", recordsArray[0]->_sub2,recordsArray[0]->_var_idx, "S", recordsArray[0]->_sub);//, "~");
      }
      else if(recordsArray[0]->_flag == 1){
      fprintf(positions_file, "%d\t%d\t%d\t%s\t[(%d,'%s','%s')]\n", recordsArray[0]->_readId, curr_pos, recordsArray[0]->_strand, recordsArray[0]->_seq, recordsArray[0]->_var_idx2, "S", recordsArray[0]->_sub2);//, "~");
      }else{
        fprintf(positions_file, "%d\t%d\t%d\t%s\t[]\n", recordsArray[0]->_readId, curr_pos, recordsArray[0]->_strand, recordsArray[0]->_seq);
      }

      // Write index differences
      for(i = 1; i < match_cnt; i++){
  curr_pos = recordsArray[i]->_ref_idx;
  if(recordsArray[i]->_flag == 2){
  fprintf(positions_file, "%d\t%d\t%d\t%s\t[(%d,'%s','%s'), (%d,'%s','%s')]\n", recordsArray[i]->_readId, curr_pos, recordsArray[i]->_strand, recordsArray[i]->_seq, recordsArray[i]->_var_idx2, "S", recordsArray[i]->_sub2, recordsArray[i]->_var_idx, "S", recordsArray[i]->_sub);//, "~");
  }
  else if(recordsArray[i]->_flag == 1){
  fprintf(positions_file, "%d\t%d\t%d\t%s\t[(%d,'%s','%s')]\n", recordsArray[i]->_readId, curr_pos, recordsArray[i]->_strand, recordsArray[i]->_seq, recordsArray[i]->_var_idx2, "S", recordsArray[i]->_sub2);//, "~");
  }else{
    fprintf(positions_file, "%d\t%d\t%d\t%s\t[]\n", recordsArray[i]->_readId, curr_pos, recordsArray[i]->_strand, recordsArray[i]->_seq);
  }
      }
      
      fclose(positions_file);
    }
    end = clock();
    printf("Running time (seconds): %f", (double)(end - start)/CLOCKS_PER_SEC);
    // Free up the 2BWT index
    printf("\nFree index ... "); 
    fflush(stdout);
    BWTFree2BWT(idx2BWT);
    printf("DONE\n"); 
    
    return 0;
}