#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <math.h>

// Define the k-mer length
#define KMER_LENGTH 8

// Structure to hold a k-mer and its count
typedef struct {
    char kmer[KMER_LENGTH + 1];
    double count; // Change count to double
} KmerCount;

// Function to compare two k-mers (for qsort and bsearch)
int compareKmers(const void *a, const void *b) {
    return strcmp(((KmerCount *)a)->kmer, ((KmerCount *)b)->kmer);
}

// Function to find a k-mer in an array of k-mers
KmerCount *findKmer(KmerCount *kmers, int size, char *kmer) {
    KmerCount key;
    strncpy(key.kmer, kmer, KMER_LENGTH);
    key.kmer[KMER_LENGTH] = '\0';
    return (KmerCount *)bsearch(&key, kmers, size, sizeof(KmerCount), compareKmers);
}

// Function to process a gzipped FASTQ file and count k-mers
void countKmersInFastq(const char *filename, KmerCount *kmers, int kmerCount, int *totalKmers, int *totalReads, int *readLength) {
    gzFile file = gzopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    char line[1024];
    int lineNumber = 0;
    *totalReads = 0;
    *totalKmers = 0;
    *readLength = 0;

    while (gzgets(file, line, sizeof(line))) {
        lineNumber++;
        // Only process sequence lines (every 4th line starting from line 2)
        if (lineNumber % 4 == 2) {
            // Remove newline character
            line[strcspn(line, "\n")] = 0;
            int lineLength = strlen(line);
            // Skip reads with 'N'
            if (strchr(line, 'N') != NULL) {
                continue;
            }
            (*totalReads)++;
            if (*readLength == 0) {
                *readLength = lineLength; // Set read length based on the first read
            }
            // Count all k-mers in the sequence
            for (int i = 0; i <= lineLength - KMER_LENGTH; i++) {
                char kmer[KMER_LENGTH + 1];
                strncpy(kmer, line + i, KMER_LENGTH);
                kmer[KMER_LENGTH] = '\0';
                KmerCount *found = findKmer(kmers, kmerCount, kmer);
                if (found) {
                    found->count++;
                    (*totalKmers)++;
                }
            }
        }
    }

    gzclose(file);
}

// Function to load k-mers from a CSV file
KmerCount *loadKmers(const char *filename, int *kmerCount) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening k-mer file");
        exit(EXIT_FAILURE);
    }

    // Count the number of k-mers
    char line[1024];
    *kmerCount = 0;
    while (fgets(line, sizeof(line), file)) {
        // Remove newline character and trim whitespace and quotes
        char *trimmedLine = line;
        while (*trimmedLine == ' ' || *trimmedLine == '\t' || *trimmedLine == '\n' || *trimmedLine == '"') trimmedLine++;
        char *end = trimmedLine + strlen(trimmedLine) - 1;
        while (end > trimmedLine && (*end == ' ' || *end == '\t' || *end == '\n' || *end == '"')) end--;
        *(end + 1) = '\0';

        if (strlen(trimmedLine) == KMER_LENGTH) {
            (*kmerCount)++;
        } else {
            fprintf(stderr, "Invalid k-mer length in file: '%s' (length: %lu)\n", trimmedLine, strlen(trimmedLine));
        }
    }
    rewind(file);

    // Allocate memory for k-mers
    KmerCount *kmers = malloc(*kmerCount * sizeof(KmerCount));
    if (!kmers) {
        perror("Memory allocation error");
        exit(EXIT_FAILURE);
    }

    // Read k-mers from the file
    int index = 0;
    while (fgets(line, sizeof(line), file)) {
        // Remove newline character and trim whitespace and quotes
        char *trimmedLine = line;
        while (*trimmedLine == ' ' || *trimmedLine == '\t' || *trimmedLine == '\n' || *trimmedLine == '"') trimmedLine++;
        char *end = trimmedLine + strlen(trimmedLine) - 1;
        while (end > trimmedLine && (*end == ' ' || *end == '\t' || *end == '\n' || *end == '"')) end--;
        *(end + 1) = '\0';

        if (strlen(trimmedLine) == KMER_LENGTH) {
            strncpy(kmers[index].kmer, trimmedLine, KMER_LENGTH);
            kmers[index].kmer[KMER_LENGTH] = '\0';
            kmers[index].count = 0.0; // Initialize count as double
            index++;
        } else {
            fprintf(stderr, "Skipping invalid k-mer: '%s' (length: %lu)\n", trimmedLine, strlen(trimmedLine));
        }
    }

    fclose(file);

    // Sort k-mers for binary search
    qsort(kmers, *kmerCount, sizeof(KmerCount), compareKmers);

    return kmers;
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <fastq_file> <kmers_file> <output_file>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // Load k-mers from file
    int kmerCount;
    KmerCount *kmers = loadKmers(argv[2], &kmerCount);

    // Count k-mers in FASTQ file
    int totalKmers, totalReads, readLength;
    countKmersInFastq(argv[1], kmers, kmerCount, &totalKmers, &totalReads, &readLength);

    // Compute total k-mer occurrences
    double totalOccurrences = 0.0;
    for (int i = 0; i < kmerCount; i++) {
        totalOccurrences += kmers[i].count;
    }

    // Normalize counts to sum to 4^k
    double normalizationFactor = pow(4, KMER_LENGTH);
    for (int i = 0; i < kmerCount; i++) {
        kmers[i].count = kmers[i].count*normalizationFactor / totalOccurrences;
    }

    // Open output file
    FILE *outputFile = fopen(argv[3], "w");
    if (!outputFile) {
        perror("Error opening output file");
        free(kmers);
        exit(EXIT_FAILURE);
    }

    // Print normalized k-mer counts to output file and console
    for (int i = 0; i < kmerCount; i++) {
        // Print k-mer:frequency to the output file
        fprintf(outputFile, "%s:%f\n", kmers[i].kmer, kmers[i].count);
    
        // Optional: Print to console for debugging
        printf("%s:%f\n", kmers[i].kmer, kmers[i].count);
    }

fclose(outputFile);
}
