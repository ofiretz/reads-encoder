#include "2m4p.h"

void twoMisFourPartSearch(char *pattern, Idx2BWT *idx2BWT, int firstPart, int secondPart, int thirdPart, long *sumVisitedVertices, long *saCountSum) {

    //FILE *f;
    //f = fopen("out.txt", "w");
    
    int patternLength = strlen(pattern);
    BWTConvertPattern(idx2BWT,pattern,patternLength,pattern);
        
    int i,j,k, temp;
    //Variables for backward and forward search
    unsigned int l,r,rev_l,rev_r;
    int firstSymbol, secondSymbol;
    //Variables for search all sa ranges functions
    unsigned int all_l_first_m[ALPHABET_SIZE];
    unsigned int all_r_first_m[ALPHABET_SIZE];
    unsigned int all_rev_l_first_m[ALPHABET_SIZE];
    unsigned int all_rev_r_first_m[ALPHABET_SIZE];
    unsigned int all_l_second_m[ALPHABET_SIZE];
    unsigned int all_r_second_m[ALPHABET_SIZE];
    unsigned int all_rev_l_second_m[ALPHABET_SIZE];
    unsigned int all_rev_r_second_m[ALPHABET_SIZE];
    
    //Variables for result
    unsigned int offset;
    int sequenceId;
    unsigned int saCount = 0;
    
    int visitedVertices = 0;
    
    //printf("First search\n");
    BWTSARangeInitial(idx2BWT,pattern[0],&l,&r);
    BWTSARangeInitial(idx2BWT,pattern[0],&rev_l,&rev_r);
    for (i = 1; i < firstPart; i++) 
    { 
        visitedVertices++;
        BWTSARangeForward_Bidirection(idx2BWT,pattern[i],&l,&r,&rev_l,&rev_r); 
        if (l > r) {
            break;
        }
    }
    //fprintf(f, "initial l, r: %d %d\n", l, r);
    
    //printf("initialized\n");
    if (l <= r) 
    {
        for (i = firstPart; i < patternLength; i++) 
        {
            visitedVertices += ALPHABET_SIZE;
            BWTAllSARangesForward_Bidirection(idx2BWT,l,r,rev_l,rev_r,all_l_first_m,all_r_first_m,all_rev_l_first_m,all_rev_r_first_m);
            //fprintf(f, "i = %d: ", i);
            for (firstSymbol = 0; firstSymbol < ALPHABET_SIZE; firstSymbol++) 
            {
                if (firstSymbol == pattern[i]) 
                    continue;
                l = all_l_first_m[firstSymbol];
                r = all_r_first_m[firstSymbol];
                rev_l = all_rev_l_first_m[firstSymbol];
                rev_r = all_rev_r_first_m[firstSymbol];
                //fprintf(f, "symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                if (l > r) 
                {
                    continue;
                }
                for(temp = i + 1; temp < firstPart + secondPart + thirdPart; temp++) 
                {   
                    visitedVertices++;
                    if (l > r)
                        break;
                    BWTSARangeForward_Bidirection(idx2BWT,pattern[temp],&l,&r, &rev_l, &rev_r); 
                }
                if (l > r) 
                {
                    continue;
                }         
                int jToStartFrom = firstPart + secondPart + thirdPart;
                if (i + 1 > firstPart + secondPart + thirdPart) {
                    jToStartFrom = i + 1;
                }                
                for(j = jToStartFrom; j < patternLength; ++j) 
                {
                    visitedVertices += ALPHABET_SIZE;
                    BWTAllSARangesForward_Bidirection(idx2BWT,l,r,rev_l,rev_r, all_l_second_m,all_r_second_m,all_rev_l_second_m,all_rev_r_second_m);
                    //fprintf(f, "    j = %d: ", j);
                    for (secondSymbol = 0; secondSymbol < ALPHABET_SIZE; secondSymbol++) 
                    {
                        if (secondSymbol == pattern[j]) 
                            continue;
                        l = all_l_second_m[secondSymbol];
                        r = all_r_second_m[secondSymbol];
                        rev_l = all_rev_l_second_m[secondSymbol];
                        rev_r = all_rev_r_second_m[secondSymbol];
                        //fprintf(f, "    symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                        //printf("%d %d %d %d\n", i, j, l, r);
                        if (l > r) 
                        {
                            continue;
                        }           
                        for (k = j + 1; k < patternLength; k++) 
                        {
                            visitedVertices++;
                            if (l > r) 
                                break;
                            BWTSARangeForward_Bidirection(idx2BWT, pattern[k],
                            &l,&r,&rev_l,&rev_r);
                            //fprintf(f, "        %d %d %d %d %d\n", i, j, k, l, r);
                            //printf("%d %d %d %d %d\n", i, j, k, l, r);
                        }
                        if (l <= r && k >= patternLength) {
                            //printf("found: %d %d %d %d\n", i, j, l, r);
                            //fprintf(f, "FOUND: %d %d %d %d\n", i, j, l, r);                                              
                            saCount += r - l + 1;
                        }
                    }
                    l = all_l_second_m[pattern[j]];
                    r = all_r_second_m[pattern[j]];
                    rev_l = all_rev_l_second_m[pattern[j]];
                    rev_r = all_rev_r_second_m[pattern[j]];
                }
                /*if (l <= r) 
                {
                    printf("found with one mistake: %d %d %d %d\n", i, j, l, r);
                    fprintf(f, "FOUND with one mistake: %d %d %d %d\n", i, j, l, r);                                              
                    saCount += r - l + 1;
                }*/
            }
            l = all_l_first_m[pattern[i]];
            r = all_r_first_m[pattern[i]];
            rev_l = all_rev_l_first_m[pattern[i]];
            rev_r = all_rev_r_first_m[pattern[i]];
        }
        /*if (l <= r) 
        {
            printf("found exact: %d %d %d %d\n", i, j, l, r);
            fprintf(f, "FOUND exact: %d %d %d %d\n", i, j, l, r);                                              
            saCount += r - l + 1;
        }*/
    }
    
    //fprintf(f, "Second search:\n");
    BWTSARangeInitial(idx2BWT,pattern[patternLength - 1],&l,&r);
    BWTSARangeInitial(idx2BWT,pattern[patternLength - 1],&rev_l,&rev_r);
    for (i = patternLength - 2; i >= firstPart + secondPart + thirdPart; i--) 
    { 
        visitedVertices++;
        BWTSARangeBackward(idx2BWT,pattern[i],&l,&r); 
        if (l > r) {
            break;
        }
    }
    //fprintf(f, "initial l, r: %d %d\n", l, r);
    
    if (l <= r) 
    {
        for (i = firstPart + secondPart + thirdPart - 1; i >= 0; i--) 
        {
            //fprintf(f, "i = %d: l = %d, r = %d\n",i, l, r);
            visitedVertices += ALPHABET_SIZE;
            BWTAllSARangesBackward(idx2BWT,l,r,all_l_first_m,all_r_first_m);
            //fprintf(f, "i = %d: ", i);
            for (firstSymbol = 0; firstSymbol < ALPHABET_SIZE; firstSymbol++) 
            {
                if (firstSymbol == pattern[i]) 
                    continue;
                l = all_l_first_m[firstSymbol];
                r = all_r_first_m[firstSymbol];
                //fprintf(f, "    firstSymbol = %d: l = %d, r = %d\n", firstSymbol, l, r);
                //fprintf(f, "symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                if (l > r) 
                {
                    continue;
                }
                for(temp = i - 1; temp >= firstPart + secondPart; temp--) 
                {   
                    visitedVertices++;
                    if (l > r)
                        break;
                    BWTSARangeBackward(idx2BWT,pattern[temp],&l,&r); 
                }
                if (l > r) {
                    continue;
                }
                int jStartFrom = i - 1;
                if (firstPart + secondPart - 1 < jStartFrom) {
                    jStartFrom = firstPart + secondPart - 1;
                }
                for(j = jStartFrom; j >= 0; j--) 
                {
                    //fprintf(f, "        j = %d: l = %d, r = %d\n", j, l, r);                
                    visitedVertices += ALPHABET_SIZE;
                    BWTAllSARangesBackward(idx2BWT,l,r,all_l_second_m,all_r_second_m);
                    //fprintf(f, "    j = %d: ", j);
                    for (secondSymbol = 0; secondSymbol < ALPHABET_SIZE; secondSymbol++) 
                    {
                        if (secondSymbol == pattern[j]) 
                            continue;
                        l = all_l_second_m[secondSymbol];
                        r = all_r_second_m[secondSymbol];
                        //fprintf(f, "        j = %d: l = %d, r = %d\n", j, l, r);                    
                        //fprintf(f, "    symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                        //printf("%d %d %d %d\n", i, j, l, r);
                        if (l > r) 
                        {
                            continue;
                        }           
                        for (k = j - 1; k >= 0; k--) 
                        {
                            visitedVertices++;
                            if (l > r) 
                                break;
                            BWTSARangeBackward(idx2BWT, pattern[k], &l,&r);
                            //fprintf(f, "        %d - %d, %d - %d, %d %d\n", i, firstSymbol, j, secondSymbol, l, r);
                            //printf("%d %d %d %d %d\n", i, j, k, l, r);
                        }
                        if (l <= r && k < 0) {
                            //fprintf(f, "        %d - %d, %d - %d, %d %d\n", i, firstSymbol, j, secondSymbol, l, r);
                            //printf("found: %d %d %d %d\n", i, j, l, r);
                            //fprintf(f, "FOUND: %d %d %d %d\n", i, j, l, r);                                              
                            saCount += r - l + 1;
                        }
                    }
                    l = all_l_second_m[pattern[j]];
                    r = all_r_second_m[pattern[j]];
                }
                /*if (l <= r) 
                {
                    printf("found with one mistake: %d %d %d %d\n", i, j, l, r);
                    fprintf(f, "FOUND with one mistake: %d %d %d %d\n", i, j, l, r);                                              
                    saCount += r - l + 1;
                }*/
            }
            l = all_l_first_m[pattern[i]];
            r = all_r_first_m[pattern[i]];
        }
        /*if (l <= r) 
        {
            printf("found exact: %d %d %d %d\n", i, j, l, r);
            fprintf(f, "FOUND exact: %d %d %d %d\n", i, j, l, r);                                              
            saCount += r - l + 1;
        }*/
    }
    
    //fprintf(f, "Third search:\n");
    BWTSARangeInitial(idx2BWT,pattern[firstPart],&l,&r);
    BWTSARangeInitial(idx2BWT,pattern[firstPart],&rev_l,&rev_r);
    for (i = firstPart + 1; i < firstPart + secondPart + thirdPart; i++) 
    { 
        visitedVertices++;
        BWTSARangeForward_Bidirection(idx2BWT,pattern[i],&l,&r, &rev_l, &rev_r); 
        if (l > r) {
            break;
        }
    }
    //fprintf(f, "initial l, r: %d %d\n", l, r);
    
    if (l <= r) 
    {
        for (i = firstPart + secondPart + thirdPart; i < patternLength; i++) 
        {
            visitedVertices += ALPHABET_SIZE;
            BWTAllSARangesForward_Bidirection(idx2BWT, l, r, rev_l, rev_r, all_l_first_m, all_r_first_m, all_rev_l_first_m, all_rev_r_first_m);
            //fprintf(f, "i = %d: ", i);
            for (firstSymbol = 0; firstSymbol < ALPHABET_SIZE; firstSymbol++) 
            {
                if (firstSymbol == pattern[i]) 
                    continue;
                l = all_l_first_m[firstSymbol];
                r = all_r_first_m[firstSymbol];
                rev_l = all_rev_l_first_m[firstSymbol];
                rev_r = all_rev_r_first_m[firstSymbol];
                //fprintf(f, "symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                if (l > r) 
                {
                    continue;
                }
                for(j = i + 1; j < patternLength; j++) 
                {   
                    visitedVertices++;
                    BWTSARangeForward_Bidirection(idx2BWT,pattern[j],&l,&r, &rev_l, &rev_r); 
                }           
                for(j = firstPart - 1; j >= 0; j--) 
                {
                    visitedVertices += ALPHABET_SIZE;
                    BWTAllSARangesBackward(idx2BWT,l,r,  all_l_second_m, all_r_second_m);
                    //fprintf(f, "    j = %d: ", j);
                    for (secondSymbol = 0; secondSymbol < ALPHABET_SIZE; secondSymbol++) 
                    {
                        if (secondSymbol == pattern[j]) 
                            continue;
                        l = all_l_second_m[secondSymbol];
                        r = all_r_second_m[secondSymbol];
                        //fprintf(f, "    symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                        //printf("%d %d %d %d\n", i, j, l, r);
                        if (l > r) 
                        {
                            continue;
                        }           
                        for (k = j - 1; k >= 0; k--) 
                        {
                            visitedVertices++;
                            if (l > r) 
                                break;
                            BWTSARangeBackward(idx2BWT, pattern[k], &l,&r);
                            //fprintf(f, "        %d %d %d %d %d\n", i, j, k, l, r);
                            //printf("%d %d %d %d %d\n", i, j, k, l, r);
                        }
                        if (l <= r && k < 0) {
                            //printf("found: %d %d %d %d\n", i, j, l, r);
                            //fprintf(f, "FOUND: %d %d %d %d\n", i, j, l, r);                                              
                            saCount += r - l + 1;
                        }
                    }
                    l = all_l_second_m[pattern[j]];
                    r = all_r_second_m[pattern[j]];
                }
                /*if (l <= r) 
                {
                    printf("found with one mistake: %d %d %d %d\n", i, j, l, r);
                    fprintf(f, "FOUND with one mistake: %d %d %d %d\n", i, j, l, r);                                              
                    saCount += r - l + 1;
                }*/
            }
            l = all_l_first_m[pattern[i]];
            r = all_r_first_m[pattern[i]];
            rev_l = all_rev_l_first_m[pattern[i]];
            rev_r = all_rev_r_first_m[pattern[i]];
        }
        /*if (l <= r) 
        {
            printf("found exact: %d %d %d %d\n", i, j, l, r);
            fprintf(f, "FOUND exact: %d %d %d %d\n", i, j, l, r);                                              
            saCount += r - l + 1;
        }*/
    }
    
    //printf("Forth search:\n");
    BWTSARangeInitial(idx2BWT,pattern[0],&l,&r);
    BWTSARangeInitial(idx2BWT,pattern[0],&rev_l,&rev_r);
    for (i = 1; i < firstPart + secondPart; i++) 
    { 
        visitedVertices++;
        BWTSARangeForward_Bidirection(idx2BWT,pattern[i],&l,&r,&rev_l,&rev_r); 
        if (l > r) {
            break;
        }
    }
    //fprintf(f, "initial l, r: %d %d\n", l, r);
    
    if (l <= r) 
    {
        for (i = firstPart + secondPart; i < firstPart + secondPart + thirdPart; i++) 
        {
            visitedVertices += ALPHABET_SIZE;
            BWTAllSARangesForward_Bidirection(idx2BWT,l,r,rev_l,rev_r,all_l_first_m,all_r_first_m,all_rev_l_first_m,all_rev_r_first_m);
            //fprintf(f, "i = %d: ", i);
            for (firstSymbol = 0; firstSymbol < ALPHABET_SIZE; firstSymbol++) 
            {
                if (firstSymbol == pattern[i]) 
                    continue;
                l = all_l_first_m[firstSymbol];
                r = all_r_first_m[firstSymbol];
                rev_l = all_rev_l_first_m[firstSymbol];
                rev_r = all_rev_r_first_m[firstSymbol];
                //fprintf(f, "symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                if (l > r) 
                {
                    continue;
                }
                        
                for(j = i + 1; j < firstPart + secondPart + thirdPart; ++j) 
                {
                    visitedVertices += ALPHABET_SIZE;
                    BWTAllSARangesForward_Bidirection(idx2BWT,l,r,rev_l,rev_r, all_l_second_m,all_r_second_m,all_rev_l_second_m,all_rev_r_second_m);
                    //fprintf(f, "    j = %d: ", j);
                    for (secondSymbol = 0; secondSymbol < ALPHABET_SIZE; secondSymbol++) 
                    {
                        if (secondSymbol == pattern[j]) 
                            continue;
                        l = all_l_second_m[secondSymbol];
                        r = all_r_second_m[secondSymbol];
                        rev_l = all_rev_l_second_m[secondSymbol];
                        rev_r = all_rev_r_second_m[secondSymbol];
                        //fprintf(f, "    symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                        //printf("%d %d %d %d\n", i, j, l, r);
                        if (l > r) 
                        {
                            continue;
                        }           
                        for (k = j + 1; k < patternLength; k++) 
                        {
                            visitedVertices++;
                            if (l > r) 
                                break;
                            BWTSARangeForward_Bidirection(idx2BWT, pattern[k],
                            &l,&r,&rev_l,&rev_r);
                            //fprintf(f, "        %d %d %d %d %d\n", i, j, k, l, r);
                            //printf("%d %d %d %d %d\n", i, j, k, l, r);
                        }
                        if (l <= r && k >= patternLength) {
                            //printf("found: %d %d %d %d\n", i, j, l, r);
                            //fprintf(f, "FOUND: %d %d %d %d\n", i, j, l, r);                                              
                            saCount += r - l + 1;
                        }
                    }
                    l = all_l_second_m[pattern[j]];
                    r = all_r_second_m[pattern[j]];
                    rev_l = all_rev_l_second_m[pattern[j]];
                    rev_r = all_rev_r_second_m[pattern[j]];
                }
                /*if (l <= r) 
                {
                    printf("found with one mistake: %d %d %d %d\n", i, j, l, r);
                    fprintf(f, "FOUND with one mistake: %d %d %d %d\n", i, j, l, r);                                              
                    saCount += r - l + 1;
                }*/
            }
            l = all_l_first_m[pattern[i]];
            r = all_r_first_m[pattern[i]];
            rev_l = all_rev_l_first_m[pattern[i]];
            rev_r = all_rev_r_first_m[pattern[i]];
        }
        /*if (l <= r) 
        {
            printf("found exact: %d %d %d %d\n", i, j, l, r);
            fprintf(f, "FOUND exact: %d %d %d %d\n", i, j, l, r);                                              
            saCount += r - l + 1;
        }*/
    }
    
    //printf("\n%d visited vertices\n", visitedVertices);
    //printf("%u SA-indexes/occurrences were found.\n\n",saCount);
    *sumVisitedVertices += visitedVertices;
    *saCountSum += saCount;
}

