#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* readFile(FILE* fp);

int main() {
    FILE* fp;
    char* file_name = "testing.txt";
    
    fp = fopen(file_name, "r");

    if (!fp) { printf("FASTA file %s not opened\n", file_name); exit(0); }
    printf("FASTA file %s opened\n", file_name);

    while (1) {
        char* a = readFile(fp);
        if (a == NULL) { break; }
        printf("%s\n", a);
        free(a);
    }

    fclose(fp);

}

char* readFile(FILE* fp) {
    char* line = NULL;
    size_t len;

    int i = 0;

    // hard-code for now, change later
    char* read = (char*) malloc(sizeof(char) * 110);

    while (getline(&line, &len, fp) != -1) {
        if (line[0] != '>') {
            while (line[i] != '\0' && line[i] != '\n') {
                read[i] = line[i];
                i++;
            }
        return read;
        }
    }
    free(read);
    return NULL;
}