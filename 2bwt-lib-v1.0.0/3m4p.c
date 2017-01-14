#include "3m4p.h"

void threeMisFourPartSearch(char *pattern, Idx2BWT *idx2BWT, int firstPart, int secondPart, int thirdPart, long *sumVisitedVertices, long *saCountSum) {

    //FILE *f;
    //f = fopen("out.txt", "w");
    
    int patternLength = strlen(pattern);
    BWTConvertPattern(idx2BWT,pattern,patternLength,pattern);
        
    int i, j, k, t, temp, i1;
    //Variables for backward and forward search
    unsigned int l, r, rev_l, rev_r;
    int firstSymbol, secondSymbol, thirdSymbol;
    //Variables for search all sa ranges functions
    unsigned int all_l_first_m[ALPHABET_SIZE];
    unsigned int all_r_first_m[ALPHABET_SIZE];
    unsigned int all_rev_l_first_m[ALPHABET_SIZE];
    unsigned int all_rev_r_first_m[ALPHABET_SIZE];
    unsigned int all_l_second_m[ALPHABET_SIZE];
    unsigned int all_r_second_m[ALPHABET_SIZE];
    unsigned int all_rev_l_second_m[ALPHABET_SIZE];
    unsigned int all_rev_r_second_m[ALPHABET_SIZE];
    unsigned int all_l_third_m[ALPHABET_SIZE];
    unsigned int all_r_third_m[ALPHABET_SIZE];
    unsigned int all_rev_l_third_m[ALPHABET_SIZE];
    unsigned int all_rev_r_third_m[ALPHABET_SIZE];
    
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
            if (l > r) {
                break;
            }
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
                for(temp = i + 1; temp < firstPart + secondPart; temp++) 
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
                int jToStartFrom = firstPart + secondPart;
                if (i + 1 > firstPart + secondPart) {
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
                        for(t = j + 1; t < patternLength; ++t) 
                        {
                            visitedVertices += ALPHABET_SIZE;
                            BWTAllSARangesForward_Bidirection(idx2BWT,l,r,rev_l,rev_r, all_l_third_m,all_r_third_m,all_rev_l_third_m,all_rev_r_third_m);
                            //fprintf(f, "    j = %d: ", j);
                            for (thirdSymbol = 0; thirdSymbol < ALPHABET_SIZE; thirdSymbol++) 
                            {
                                if (thirdSymbol == pattern[t]) 
                                    continue;
                                l = all_l_third_m[thirdSymbol];
                                r = all_r_third_m[thirdSymbol];
                                rev_l = all_rev_l_third_m[thirdSymbol];
                                rev_r = all_rev_r_third_m[thirdSymbol];
                                //fprintf(f, "    symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                                //printf("%d %d %d %d\n", i, j, l, r);
                                if (l > r) 
                                {
                                    continue;
                                }           
                                for (k = t + 1; k < patternLength; k++) 
                                {
                                    visitedVertices++;
                                    if (l > r) 
                                        break;
                                    BWTSARangeForward_Bidirection(idx2BWT, pattern[k], &l,&r,&rev_l,&rev_r);
                                    //fprintf(f, "        %d %d %d %d %d\n", i, j, k, l, r);
                                    //printf("%d %d %d %d %d\n", i, j, k, l, r);
                                }
                                if (l <= r && k >= patternLength) {
                                    //printf("found: %d %d %d %d\n", i, j, l, r);
                                    //fprintf(f, "FOUND: %d %d %d %d\n", i, j, l,  r);                                              
                                    saCount += r - l + 1;
                                }
                            }
                            l = all_l_third_m[pattern[t]];
                            r = all_r_third_m[pattern[t]];
                            rev_l = all_rev_l_third_m[pattern[t]];
                            rev_r = all_rev_r_third_m[pattern[t]];
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
    
    //2 search
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
    
    //printf("initialized\n");
    if (l <= r) 
    {
        for (i = firstPart + secondPart + thirdPart - 1; i >= 0; i--) 
        {
            if (l > r) {
                break;
            }
            visitedVertices += ALPHABET_SIZE;
            BWTAllSARangesBackward(idx2BWT,l,r,all_l_first_m,all_r_first_m);
            //fprintf(f, "i = %d: ", i);
            for (firstSymbol = 0; firstSymbol < ALPHABET_SIZE; firstSymbol++) 
            {
                if (firstSymbol == pattern[i]) 
                    continue;
                l = all_l_first_m[firstSymbol];
                r = all_r_first_m[firstSymbol];
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
                if (l > r) 
                {
                    continue;
                }         
                int jToStartFrom = firstPart + secondPart - 1;
                if (i - 1 < jToStartFrom) {
                    jToStartFrom = i - 1;
                }                
                for(j = jToStartFrom; j >= 0; j--) 
                {
                    visitedVertices += ALPHABET_SIZE;
                    BWTAllSARangesBackward(idx2BWT,l,r, all_l_second_m,all_r_second_m);
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
                        for(t = j - 1; t >= 0; t--) 
                        {
                            visitedVertices += ALPHABET_SIZE;
                            BWTAllSARangesBackward(idx2BWT,l,r, all_l_third_m, all_r_third_m);
                            //fprintf(f, "    j = %d: ", j);
                            for (thirdSymbol = 0; thirdSymbol < ALPHABET_SIZE; thirdSymbol++) 
                            {
                                if (thirdSymbol == pattern[t]) 
                                    continue;
                                l = all_l_third_m[thirdSymbol];
                                r = all_r_third_m[thirdSymbol];
                                //fprintf(f, "    symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                                //printf("%d %d %d %d\n", i, j, l, r);
                                if (l > r) 
                                {
                                    continue;
                                }           
                                for (k = t - 1; k >= 0; k--) 
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
                                    //fprintf(f, "FOUND: %d %d %d %d\n", i, j, l,  r);                                              
                                    saCount += r - l + 1;
                                }
                            }
                            l = all_l_third_m[pattern[t]];
                            r = all_r_third_m[pattern[t]];
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

    //3 search
    BWTSARangeInitial(idx2BWT,pattern[firstPart + secondPart - 1],&l,&r);
    BWTSARangeInitial(idx2BWT,pattern[firstPart + secondPart - 1],&rev_l,&rev_r);
    for (i = firstPart + secondPart - 2; i >= firstPart; i--) 
    { 
        visitedVertices++;
        BWTSARangeBackward_Bidirection(idx2BWT,pattern[i],&l,&r,&rev_l,&rev_r); 
        if (l > r) {
            break;
        }
    }
    //fprintf(f, "initial l, r: %d %d\n", l, r);
    
    //printf("initialized\n");
    if (l <= r) 
    {
        for (i = firstPart - 1; i >= 0; i--) 
        {
            if (l > r) {
                break;
            }
            visitedVertices += ALPHABET_SIZE;
            BWTAllSARangesBackward_Bidirection(idx2BWT,l,r,rev_l,rev_r,all_l_first_m,all_r_first_m,all_rev_l_first_m,all_rev_r_first_m);
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
                for(temp = i - 1; temp >= 0; temp--) 
                {   
                    visitedVertices++;
                    if (l > r)
                        break;
                    BWTSARangeBackward_Bidirection(idx2BWT,pattern[temp],&l,&r, &rev_l, &rev_r); 
                }
                if (l > r) 
                {
                    continue;
                }         
                int jToStartFrom = firstPart + secondPart;
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
                        for(t = j + 1; t < patternLength; ++t) 
                        {
                            visitedVertices += ALPHABET_SIZE;
                            BWTAllSARangesForward_Bidirection(idx2BWT,l,r,rev_l,rev_r, all_l_third_m,all_r_third_m,all_rev_l_third_m,all_rev_r_third_m);
                            //fprintf(f, "    j = %d: ", j);
                            for (thirdSymbol = 0; thirdSymbol < ALPHABET_SIZE; thirdSymbol++) 
                            {
                                if (thirdSymbol == pattern[t]) 
                                    continue;
                                l = all_l_third_m[thirdSymbol];
                                r = all_r_third_m[thirdSymbol];
                                rev_l = all_rev_l_third_m[thirdSymbol];
                                rev_r = all_rev_r_third_m[thirdSymbol];
                                //fprintf(f, "    symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                                //printf("%d %d %d %d\n", i, j, l, r);
                                if (l > r) 
                                {
                                    continue;
                                }           
                                for (k = t + 1; k < patternLength; k++) 
                                {
                                    visitedVertices++;
                                    if (l > r) 
                                        break;
                                    BWTSARangeForward_Bidirection(idx2BWT, pattern[k], &l,&r,&rev_l,&rev_r);
                                    //fprintf(f, "        %d %d %d %d %d\n", i, j, k, l, r);
                                    //printf("%d %d %d %d %d\n", i, j, k, l, r);
                                }
                                if (l <= r && k >= patternLength) {
                                    //printf("found: %d %d %d %d\n", i, j, l, r);
                                    //fprintf(f, "FOUND: %d %d %d %d\n", i, j, l,  r);                                              
                                    saCount += r - l + 1;
                                }
                            }
                            l = all_l_third_m[pattern[t]];
                            r = all_r_third_m[pattern[t]];
                            rev_l = all_rev_l_third_m[pattern[t]];
                            rev_r = all_rev_r_third_m[pattern[t]];
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

    //4 search
    BWTSARangeInitial(idx2BWT,pattern[firstPart + secondPart],&l,&r);
    BWTSARangeInitial(idx2BWT,pattern[firstPart + secondPart],&rev_l,&rev_r);
    for (i = firstPart + secondPart + 1; i < firstPart + secondPart + thirdPart; i++) 
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
        for (i = firstPart + secondPart + thirdPart; i < patternLength; i++) 
        {
            if (l > r) {
                break;
            }
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
                for(temp = i + 1; temp < patternLength; temp++) 
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
                int jToStartFrom = firstPart + secondPart - 1;
                for(j = jToStartFrom; j >= 0; j--) 
                {
                    visitedVertices += ALPHABET_SIZE;
                    BWTAllSARangesBackward(idx2BWT,l,r, all_l_second_m,all_r_second_m);
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
                        for(t = j - 1; t >= 0; t--) 
                        {
                            visitedVertices += ALPHABET_SIZE;
                            BWTAllSARangesBackward(idx2BWT, l, r, all_l_third_m, all_r_third_m);
                            //fprintf(f, "    j = %d: ", j);
                            for (thirdSymbol = 0; thirdSymbol < ALPHABET_SIZE; thirdSymbol++) 
                            {
                                if (thirdSymbol == pattern[t]) 
                                    continue;
                                l = all_l_third_m[thirdSymbol];
                                r = all_r_third_m[thirdSymbol];
                                //fprintf(f, "    symbol = %d, l = %d, r = %d\n", firstSymbol, l, r);            
                                //printf("%d %d %d %d\n", i, j, l, r);
                                if (l > r) 
                                {
                                    continue;
                                }           
                                for (k = t - 1; k >= 0; k--) 
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
                                    //fprintf(f, "FOUND: %d %d %d %d\n", i, j, l,  r);                                              
                                    saCount += r - l + 1;
                                }
                            }
                            l = all_l_third_m[pattern[t]];
                            r = all_r_third_m[pattern[t]];
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
    //printf("\n%d visited vertices\n", visitedVertices);
    //printf("%u SA-indexes/occurrences were found.\n\n",saCount);
    *sumVisitedVertices += visitedVertices;
    *saCountSum += saCount;
}

