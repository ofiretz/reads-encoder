/*

gcc 2m.c 2m3p.h 2m3p.c BWT.o dictionary.o DNACount.o HSP.o HSPstatistic.o iniparser.o inistrlib.o karlin.o MemManager.o MiscUtilities.o QSufSort.o r250.o Socket.o TextConverter.o 2BWT-Interface.o -lm -o 2m

*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "2BWT-Interface.h"
#include "2m3p.h"
#include "2m4p.h"

struct timeval tv1,tv2,dtv;

struct timezone tz;

void time_start() { 
    gettimeofday(&tv1, &tz); 
}

long time_stop()
{ 
    gettimeofday(&tv2, &tz);
    dtv.tv_sec= tv2.tv_sec -tv1.tv_sec;
    dtv.tv_usec=tv2.tv_usec-tv1.tv_usec;
    if(dtv.tv_usec<0) { 
        dtv.tv_sec--; 
        dtv.tv_usec+=1000000; 
    }
    return dtv.tv_sec*1000+dtv.tv_usec/1000;
}

long sumVisitedVertices = 0;
long saCountSum = 0;
int searchCount = 0;

int firstPart, secondPart, thirdPart;

int main(int argc, char **argv) {
    //Variables for pattern
    char* pattern;
    char* line;
    double times[5];
    double visitedVertices[5];
    //strcpy(pattern, "CTGAACTGAACCAACACTGAACTGAATAAAAAA");
    int searchType = 0;
    int patternLengthToSearch = 30;
    if (argc > 1) {
        patternLengthToSearch = atoi(argv[1]);
    }
    pattern = (char*)malloc(patternLengthToSearch + 1);
    line = (char*)malloc(1024);
    
    fflush(stdout);
    Idx2BWT * idx2BWT = BWTLoad2BWT("chr14.fasta.index",".sa");
    //printf("DONE\n\n"); 
    FILE *outA, *outTimes;
    char prefix[1024];
    strcpy(prefix, "results_2mis_1msearches_final");
    outTimes = fopen(strcat(prefix, "_times.txt"), "w");
    strcpy(prefix, "results_2mis_1msearches_final");
    outA = fopen(strcat(prefix, "_a.txt"), "w");
    char symbol;
    int numberOfSearchesToDo = 1000000;
for(patternLengthToSearch = 12; patternLengthToSearch < 46; patternLengthToSearch += 3) {
fprintf(outA, "\nSearch for length %d\n", patternLengthToSearch);
fprintf(outTimes, "\nSearch for length %d\n", patternLengthToSearch);
for(searchType = -1; searchType < 4; ++searchType) {
    sumVisitedVertices = 0;
    saCountSum = 0;
    searchCount = 0;
    //fprintf(out, "\n    Search type: %d\n", searchType);
    int currentLengthOfPattern = 0;
    time_start();
    FILE *in;
    in = fopen("reads_14chr.fa", "r");
    if (searchType >= 0) {
        switch (searchType) {
            case 0:
                firstPart = patternLengthToSearch / 3;
                secondPart = patternLengthToSearch / 3;
                break;
            case 1:
                switch (patternLengthToSearch) {
                    case 12:
                        firstPart = 4;
                        secondPart = 4;    
                        break;
                    case 15:
                        firstPart = 5;
                        secondPart = 5;    
                        break;
                    case 18:
                        firstPart = 6;
                        secondPart = 6;    
                        break;
                    case 21:
                        firstPart = 8;
                        secondPart = 6;    
                        break;
                    case 24:
                        firstPart = 10;
                        secondPart = 7;    
                        break;
                    case 27:
                        firstPart = 12;
                        secondPart = 7;    
                        break;
                    case 30:
                        firstPart = 12;
                        secondPart = 9;    
                        break;
                    case 33:
                        firstPart = 13;
                        secondPart = 10;    
                        break;
                    case 36:
                        firstPart = 13;
                        secondPart = 11;    
                        break;
                    case 39:
                        firstPart = 13;
                        secondPart = 13;    
                        break;
                    case 42:
                        firstPart = 14;
                        secondPart = 14;    
                        break;
                    case 45:
                        firstPart = 15;
                        secondPart = 15;    
                        break;
                }
                break;
            case 2:
                switch (patternLengthToSearch % 4) {
                    case 0:
                        firstPart = patternLengthToSearch / 4;
                        secondPart = patternLengthToSearch / 4;
                        thirdPart = patternLengthToSearch / 4;
                        break;
                    case 1:
                        firstPart = patternLengthToSearch / 4;
                        secondPart = patternLengthToSearch / 4;
                        thirdPart = patternLengthToSearch / 4;
                        break;
                    case 2:
                        firstPart = patternLengthToSearch / 4 + 1;
                        secondPart = patternLengthToSearch / 4;
                        thirdPart = patternLengthToSearch / 4 + 1;
                        break;
                    case 3:
                        firstPart = patternLengthToSearch / 4 + 1;
                        secondPart = patternLengthToSearch / 4 + 1;
                        thirdPart = patternLengthToSearch / 4 + 1;
                        break;
                }
                break;
            case 3:
                switch (patternLengthToSearch) {
                    case 12:
                        firstPart = 2;
                        secondPart = 5;
                        thirdPart = 1;     
                        break;
                    case 15:
                        firstPart = 3;
                        secondPart = 5;
                        thirdPart = 1;     
                        break;
                    case 18:
                        firstPart = 5;
                        secondPart = 5;
                        thirdPart = 1;     
                        break;
                    case 21:
                        firstPart = 6;
                        secondPart = 4;
                        thirdPart = 3;     
                        break;
                    case 24:
                        firstPart = 7;
                        secondPart = 4;
                        thirdPart = 4;    
                        break;
                    case 27:
                        firstPart = 9;
                        secondPart = 4;
                        thirdPart = 5;    
                        break;
                    case 30:
                        firstPart = 10;
                        secondPart = 4;    
                        thirdPart = 6;
                        break;
                    case 33:
                        firstPart = 11;
                        secondPart = 5;    
                        thirdPart = 6;
                        break;
                    case 36:
                        firstPart = 12;
                        secondPart = 5;
                        thirdPart = 7;    
                        break;
                    case 39:
                        firstPart = 13;
                        secondPart = 6;    
                        thirdPart = 7;
                        break;
                    case 42:
                        firstPart = 14;
                        secondPart = 6;    
                        thirdPart = 8;
                        break;
                    case 45:
                        firstPart = 15;
                        secondPart = 7;
                        thirdPart = 8;    
                        break;
                }
                break;    
        }
    }
    //fprintf(out, "%d %d\n", firstPart, secondPart);
    
    while (fscanf(in, "%c", &symbol) != EOF) {
        int goodLine = 1;
        int numberOfSymbolsInLine = 0;
        do {
            if (symbol == '\n') {
                break;
            } else {
                if (symbol == 'A' 
                    || symbol == 'C' 
                    || symbol == 'G' 
                    || symbol == 'T')
                {
                    line[numberOfSymbolsInLine] = symbol;
                    numberOfSymbolsInLine++;
                }
                else {
                    goodLine = 0;
                }
            }
        } while (fscanf(in, "%c", &symbol) != EOF);
        if (!goodLine) {
            continue;
        }
        currentLengthOfPattern = 0;
        //printf("line: %s...\n", line);
        //printf("%d\n", numberOfSymbolsInLine);
        int i = 0;
        while (i < numberOfSymbolsInLine) {
            if (currentLengthOfPattern == patternLengthToSearch) {
                //printf("Search for: %s...\n", pattern);
                //if (searchCount == 8) {
                //    printf("Search for: %s\n", pattern);
                //}
                //printf("parts length: %d %d\n", firstPart, secondPart);
                if (searchType >= 0) {
                    if (searchType <= 1) {
                        twoMisThreePartSearch(pattern, idx2BWT, firstPart, secondPart, &sumVisitedVertices, &saCountSum);
                    }
                    else {
                        //twoMisFourPartSearch(pattern, idx2BWT, firstPart, secondPart, thirdPart, &sumVisitedVertices, &saCountSum);
                    }
                }
                currentLengthOfPattern = 0;
                searchCount++;                
                
            }
            pattern[currentLengthOfPattern] = line[i];
            currentLengthOfPattern++;
            i++;
            if (searchCount >= numberOfSearchesToDo) {
                break;
            }
            //fprintf(out, "%d %ld %ld\n", searchCount, sumVisitedVertices, saCountSum);
        }
        if (searchCount >= numberOfSearchesToDo) {
            break;
        }
    }
    times[searchType + 1] = time_stop() / 1000.0;
    if (searchType >= 0) {
        times[searchType + 1] -= times[0];
    }
    //fprintf(out, "Time: %lf\n", times[searchType + 1]);
    if (searchType >= 0) {
        times[searchType + 1] -= times[0];
        //fprintf(out, "Time without reading: %lf\n", times[searchType + 1]);
        fprintf(outTimes, "%lf & ", times[searchType + 1]);
    }
    if (searchType >= 1) {
        //fprintf(out, "Time / time for equal parts, percents: %lf\n", 100.0 * times[searchType + 1] / times[1]);
        fprintf(outTimes, "%lf & ", 100.0 * times[searchType + 1] / times[1]);
    }
    visitedVertices[searchType + 1] = ((double)sumVisitedVertices) / searchCount;
    //fprintf(out, "Average visited vertices: %lf\n", visitedVertices[searchType + 1]);
    if (searchType >= 0) {
        fprintf(outA, "%lf & ", visitedVertices[searchType + 1]);
    }
    if (searchType >= 1) {
        //fprintf(out, "Visited vertices / visited vertices for equal parts, percents: %lf\n", 100.0 * visitedVertices[searchType + 1] / visitedVertices[1]);
        fprintf(outA, "%lf & ", 100.0 * visitedVertices[searchType + 1] / visitedVertices[1]);
    }
    //fprintf(out, "Searches count %d\n", searchCount);
    //fprintf(outA, "Sum occurences: %ld\n", saCountSum);
}
}
 //   int i, j, l, r, rev_l, rev_r;
 //   BWTConvertPattern(idx2BWT,pattern,patternLength,pattern);
   
    // Free up the 2BWT index
    //printf("\nFree index ... "); 
    fflush(stdout);
    BWTFree2BWT(idx2BWT);
    //printf("DONE\n"); 
    
    return 0;
}
