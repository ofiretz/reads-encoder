/*

gcc 3m.c 3m4p.h 3m4p.c 3m5p.h 3m5p.c BWT.o dictionary.o DNACount.o HSP.o HSPstatistic.o iniparser.o inistrlib.o karlin.o MemManager.o MiscUtilities.o QSufSort.o r250.o Socket.o TextConverter.o 2BWT-Interface.o -lm -o 3m

*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
//#include "2BWT-Interface.h"
#include "3m4p.h"
#include "3m5p.h"

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

int firstPart, secondPart, thirdPart, forthPart;

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
    Idx2BWT * idx2BWT = BWTLoad2BWT("rand_text.fa.index",".sa");
    //printf("DONE\n\n"); 
    
    FILE *outTimes;
    FILE *outA;
    //out = fopen("results_gen14chr_1m_searches.txt", "w");
    char prefix[1024];
    strcpy(prefix, "results_3mis_1msearches_final");
    outTimes = fopen(strcat(prefix, "_times.txt"), "w");
    strcpy(prefix, "results_3mis_1msearches_final");
    outA = fopen(strcat(prefix, "_a.txt"), "w");
    char symbol;
    int numberOfSearchesToDo = 1000000;
for(patternLengthToSearch = 6; patternLengthToSearch < 40; patternLengthToSearch += 3) {
//fprintf(outTimes, "\nSearch for length %d\n", patternLengthToSearch);
//fprintf(outA, "\nSearch for length %d\n", patternLengthToSearch);
fprintf(outTimes, "\n%d & ", patternLengthToSearch);
fprintf(outA, "\n%d & ", patternLengthToSearch);
for(searchType = -1; searchType < 3; ++searchType) {
    sumVisitedVertices = 0;
    saCountSum = 0;
    searchCount = 0;
    //fprintf(out, "\n    Search type: %d\n", searchType);
    int currentLengthOfPattern = 0;
    time_start();
    FILE *in;
    in = fopen("rand_reads.fa", "r");
    if (searchType >= 0) {
        switch (searchType) {
            case 0:
                switch (patternLengthToSearch % 4) {
                    case 0:
                        firstPart = patternLengthToSearch / 4;
                        secondPart = patternLengthToSearch / 4;
                        thirdPart = patternLengthToSearch / 4;
                        break;
                    case 1:
                        firstPart = patternLengthToSearch / 4 + 1;
                        secondPart = patternLengthToSearch / 4;
                        thirdPart = patternLengthToSearch / 4;
                        break;
                    case 2:
                        firstPart = patternLengthToSearch / 4;
                        secondPart = patternLengthToSearch / 4 + 1;
                        thirdPart = patternLengthToSearch / 4;
                        break;
                    case 3:
                        firstPart = patternLengthToSearch / 4;
                        secondPart = patternLengthToSearch / 4 + 1;
                        thirdPart = patternLengthToSearch / 4 + 1;
                        break;
                }
                break;
            case 1:
                switch (patternLengthToSearch % 5) {
                    case 0:
                        firstPart = patternLengthToSearch / 5;
                        secondPart = patternLengthToSearch / 5;
                        thirdPart = patternLengthToSearch / 5;
                        forthPart = patternLengthToSearch / 5;
                        break;
                    case 1:
                        firstPart = patternLengthToSearch / 5;
                        secondPart = patternLengthToSearch / 5 + 1;
                        thirdPart = patternLengthToSearch / 5;
                        forthPart = patternLengthToSearch / 5;
                        break;
                    case 2:
                        firstPart = patternLengthToSearch / 5;
                        secondPart = patternLengthToSearch / 5 + 1;
                        thirdPart = patternLengthToSearch / 5 + 1;
                        forthPart = patternLengthToSearch / 5;
                        break;
                    case 3:
                        firstPart = patternLengthToSearch / 5;
                        secondPart = patternLengthToSearch / 5 + 1;
                        thirdPart = patternLengthToSearch / 5 + 1;
                        forthPart = patternLengthToSearch / 5;
                        break;
                    case 4:
                        firstPart = patternLengthToSearch / 5 + 1;
                        secondPart = patternLengthToSearch / 5 + 1;
                        thirdPart = patternLengthToSearch / 5 + 1;
                        forthPart = patternLengthToSearch / 5;
                        break;
                }
                break;
            case 2:
                switch (patternLengthToSearch) {
                    case 6:
                        firstPart = 1;
                        secondPart = 1;
                        thirdPart = 2;
                        forthPart = 1;
                        break;
                    case 9:
                        firstPart = 1;
                        secondPart = 2;
                        thirdPart = 3;
                        forthPart = 1;
                        break;
                    case 12:
                        firstPart = 1;
                        secondPart = 2;
                        thirdPart = 6;
                        forthPart = 1;
                        break;
                    case 15:
                        firstPart = 2;
                        secondPart = 2;
                        thirdPart = 6;
                        forthPart = 1;
                        break;
                    case 18:
                        firstPart = 4;
                        secondPart = 2;
                        thirdPart = 6;
                        forthPart = 1;
                        break;
                    case 21:
                        firstPart = 3;
                        secondPart = 6;
                        thirdPart = 4;
                        forthPart = 1;
                        break;
                    case 24:
                        firstPart = 2;
                        secondPart = 9;
                        thirdPart = 3;
                        forthPart = 1;
                        break;
                    case 27:
                        firstPart = 3;
                        secondPart = 10;
                        thirdPart = 3;
                        forthPart = 1;
                        break;
                    case 30:
                        firstPart = 4;
                        secondPart = 10;
                        thirdPart = 4;
                        forthPart = 1;
                        break;
                    case 33:
                        firstPart = 6;
                        secondPart = 9;
                        thirdPart = 6;
                        forthPart = 1;
                        break;
                    case 36:
                        firstPart = 8;
                        secondPart = 8;
                        thirdPart = 8;
                        forthPart = 1;
                        break;
                }
                break;    
            
        }
    }
    //fprintf(out, "%d %d\n", firstPart, secondPart, thirdPart);
    
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
                //printf("parts length: %d %d\n", firstPart, secondPart);
                if (searchType >= 0) {
                    if (searchType == 0) {
                        threeMisFourPartSearch(pattern, idx2BWT, firstPart, secondPart, thirdPart, &sumVisitedVertices, &saCountSum);
                    }
                    else {
                        threeMisFivePartSearch(pattern, idx2BWT, firstPart, secondPart, thirdPart, forthPart, &sumVisitedVertices, &saCountSum);
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
    //fprintf(outTimes, "Sum occurences: %ld\n", saCountSum);
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
