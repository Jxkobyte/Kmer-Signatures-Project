#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "uthash.h"
#include <chrono>

#include <windows.h>
#include <atomic>
#include <vector>
#include <assert.h>
#include <map>
#include <array>
#include <algorithm>

#define SIGNATURE_LEN 64

struct threadData {
    char buffer[10000];
    int n;
    int ID;
};


struct fileData {
    int ID;
    int doc_sig[SIGNATURE_LEN];
};






const int numThreads = 6;
CRITICAL_SECTION HASH_SECTION;
CRITICAL_SECTION FILE_SECTION;

std::vector<threadData> td;
std::vector<fileData> fd[numThreads];


bool CompareByID(const fileData& a, const fileData& b) {
    return a.ID < b.ID;
}




typedef unsigned char byte;



int DENSITY  = 21;
int PARTITION_SIZE;

int inverse[256];
const char* alphabet = "CSTPAGNDEQHRKMILVFYW";


void seed_random(char* term, int length);
short random_num(short max);
void Init();



int WORDLEN;
FILE *sig_file;



typedef struct
{
    char term[100];
    short sig[SIGNATURE_LEN];
    UT_hash_handle hh;
} hash_term;

hash_term *vocab = NULL;


short* compute_new_term_sig(char* term, short *term_sig)
{
    seed_random(term, WORDLEN);

    int non_zero = SIGNATURE_LEN * DENSITY/100;

    int positive = 0;
    while (positive < non_zero/2)
    {
        short pos = random_num(SIGNATURE_LEN);
        if (term_sig[pos] == 0) 
	{
            term_sig[pos] = 1;
            positive++;
        }
    }

    int negative = 0;
    while (negative < non_zero/2)
    {
        short pos = random_num(SIGNATURE_LEN);
        if (term_sig[pos] == 0) 
	{
            term_sig[pos] = -1;
            negative++;
        }
    }
    return term_sig;
}

short *find_sig(char* term)
{
   
    hash_term *entry;
    
    HASH_FIND(hh, vocab, term, WORDLEN, entry);

    if (entry == NULL)
    {
        EnterCriticalSection(&HASH_SECTION);
        entry = (hash_term*)malloc(sizeof(hash_term));
        strncpy_s(entry->term, sizeof(entry->term), term, WORDLEN);
        memset(entry->sig, 0, sizeof(entry->sig));
        compute_new_term_sig(term, entry->sig);
        HASH_ADD(hh, vocab, term, WORDLEN, entry);
        LeaveCriticalSection(&HASH_SECTION);
    }
    
    return entry->sig;
}


std::atomic<int> doc = 0;

void compute_signature(char* sequence, int length, int ID, int threadId)
{  
    int doc_sig [SIGNATURE_LEN];

    for (int i = 0; i < length - WORDLEN + 1; i++) {
            short* term_sig = find_sig(sequence+i);

            for (int i = 0; i < SIGNATURE_LEN; i++) {
                doc_sig[i] += term_sig[i];
            }           
    }     
    fileData data;
    data.ID = ID;
    memcpy(data.doc_sig, doc_sig, sizeof(doc_sig));
    fd[threadId].push_back(data);
}

#define min(a,b) ((a) < (b) ? (a) : (b))

void partition(char* sequence, int length, int ID, int threadId)
{
    int i=0;
    do
    {
        compute_signature(sequence+i, min(PARTITION_SIZE, length-i), ID, threadId);
        i += PARTITION_SIZE/2;
    }
    while (i+PARTITION_SIZE/2 < length);
    
    doc++;
    
}

int power(int n, int e)
{
    int p = 1;
    for (int j=0; j<e; j++)
        p *= n;
    return p;
}





DWORD WINAPI Foo(LPVOID lpThreadParameter) {

    int threadId = (int)lpThreadParameter;

    // cyclic distribution
    for (int i = threadId; i < td.size(); i += numThreads) {

        char* buffer = td[i].buffer;
        int n = td[i].n;
        int ID = td[i].ID;

        partition(buffer, n, ID, threadId);
    }

    return 0;
}

int main(int argc, char* argv[])
{
    
    //const char* filename = "qut2.fasta";
    const char* filename = "qut3.fasta";
    
    WORDLEN = 3;
    PARTITION_SIZE = 16;
    int WORDS = power(20, WORDLEN);

    for (int i=0; i<strlen(alphabet); i++)
        inverse[alphabet[i]] = i;

    auto start = std::chrono::high_resolution_clock::now();

    FILE* file;
    errno_t OK = fopen_s(&file, filename, "r");

    if (OK != 0)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        return 1;
    }

    char outfile[256];
    sprintf_s(outfile, 256, "%s.part%d_sigs%02d_%d", filename, PARTITION_SIZE, WORDLEN, SIGNATURE_LEN);
    fopen_s(&sig_file, outfile, "w");


    
    
    
    InitializeCriticalSection(&HASH_SECTION);
    
    int ID = 0;


    char buffer[10000];
    while (!feof(file))
    {
        fgets(buffer, 10000, file); // skip meta data line
        fgets(buffer, 10000, file);
        int n = (int)strlen(buffer) - 1;
        buffer[n] = 0;

        threadData data;
        strcpy_s(data.buffer, buffer);
        data.n = n;
        data.ID = ID;
        td.push_back(data);
        ID++;    
    }

    HANDLE* handle = new HANDLE[numThreads];

    for (int i = 0; i < numThreads; i++) {
        handle[i] = CreateThread(NULL, 0, Foo, (LPVOID)i, 0, NULL);
    }


    WaitForMultipleObjects(numThreads, handle, TRUE, INFINITE);
    
    std::vector<fileData> res;

    for (int i = 0; i < numThreads; i++) {
        for (auto& j : fd[i]) {
            res.push_back(j);
        }
    }

    // Sort results to eliminate any race conditions
    std::sort(res.begin(), res.end(), CompareByID);
    
    // loop through sig_results and write to file
    int count = 0;
    
    for (auto& i : res) {

        fwrite(&i.ID, sizeof(int), 1, sig_file);
        auto doc_sig = i.doc_sig;

        for (int i = 0; i < SIGNATURE_LEN; i += 8)
        {
            byte c = 0;
            for (int j = 0; j < 8; j++)
                c |= (doc_sig[i + j] > 0) << (7 - j);
            fwrite(&c, sizeof(byte), 1, sig_file);
        }
        
        count++;
    }
    printf("Printed %d Lines", count);


    fclose(file);

    fclose(sig_file);


    DeleteCriticalSection(&HASH_SECTION);
    

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    printf("%s %f seconds\n", filename, duration.count());

    // qut3.fasta.part16_sigs03_64
    // Open and compare results
    /*
    FILE* sequential;
    FILE* parallel;
    errno_t OK1 = fopen_s(&sequential, "qut3.fasta.part16_sigs03_64_sequential", "r");
    errno_t OK2 = fopen_s(&parallel, "qut3.fasta.part16_sigs03_64", "r");
    assert(OK1 == 0);
    assert(OK2 == 0);

    char buffer1[10000];
    char buffer2[10000];

    while (!feof(sequential)) {
        fgets(buffer1, 10000, sequential);
        fgets(buffer2, 10000, parallel);
        assert(buffer1 == buffer2);
    }
    printf("Files are the same!\n");
    */


    return 0;
}
