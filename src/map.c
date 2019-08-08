#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

// direction used for traceback
typedef enum Direction {
    up,
    left,
    diagonal,
    done            // applies to nodes with 0 value
} direction;


// cell structure for use in smith-waterman
typedef struct Cell {
    int S;
    int I;
    int D;              
    int i;              // i index 
    int j;              // j index 
    direction dir;      // traceback direction 
} cell;


// (note: start_index and end_index as described in the assignment are actually map_start and end_start)
typedef struct Node {
    int id;
    int depth;          // how deep we are into the string "banana$". 2 = the first 'n'.
    int label_length;   // how much further we go, >= 1. 
                        // 2 = the second 'a'. Thus label = "na";
    struct Node** children;
    struct Node* parent;
    struct Node* SL;

    // for labels
    int start_index;
    int end_index;

    // for depths
    int map_start;
    int map_end;

} node;

typedef struct Tree {
    struct Node* root;
    struct Node* u;
} tree;


/* Global variables */

// for outputs
int num_nodes = 0;
int num_starts = 0;
int num_reads = 0;


// parameters
int ma;     // max val
int mi;     // min val
int h;      // gap continuation penalty
int g;      // gap opening penalty

int alphabet_length;    // size of alphabet

// x is the minimum concurrent # of characters that must match in a potential read
int x = 25;

// alignment-relevant
int best_optimal_score = 0;
int best_align_start = 0;

// calculations
double percent_identity;
double length_coverage;

// misc
int max_length = 0;
int next_index = 0;
int max = 0;

// tree nodes
int lower_id = 1;   // stores ID of leaves (1 to strlen(S))
int upper_id;       // stores ID of internal nodes


/* Function declarations */

/* Creation */
tree* init();                               
node* createNode();

/* Tree functions*/ 
void sortChildren(node* n, char* S);
node* findPath(tree* t, node* v, char* S, int offset);
int* findLoc(node* root, char* S, char* read, int l);
node* findLocSearch(node* v, char* read, char* S, int read_ptr);
int substitution(char a, char b);

/* Helper functions */ 
int findMaxOfThree(int a, int b, int c);
int findMaxOfFour(int a, int b, int c, int d);
direction findMaxDirection(int S, int D, int I);


/* Tree building*/ 
int* prepareST(tree* t, int n);
void DFS_PrepareST(node* n, int A[]);

/* Primary operations */ 
void mapReads(tree* t, char* S, int* A, char* reads_name);
void align(char* read, int _j, char* S);
/* End function declarations */

/* * * * * * Begin read files * * * * * */

// read our sequence fasta file
char* readFASTA (char* file_name) {
    FILE* fp;
    char* line = NULL;
    size_t len;

    fp = fopen(file_name, "r");

    if (!fp) { printf("FASTA file %s not opened\n", file_name); exit(0); }
    printf("FASTA file %s opened\n", file_name);
    
    int string_length = 0;
    int i = 0;

    while(getline(&line, &len, fp) != -1) {
        if (line[0] != '>') {
            while (line[i] != '\0' && line[i] != '\n') {
                string_length++;
                i++;
            }
            i = 0;
        }
    }
    
    char* S = (char*) malloc(string_length+2);

    rewind(fp);
    
    line = NULL;

    i = 0;
    string_length = 0;

    while(getline(&line, &len, fp) != -1) {
        if (line[0] != '>') {
            while (line[i] != '\0' && line[i] != '\n') {
                S[string_length++] = line[i++];
            }
            i = 0;
        }
    } 
    S[string_length++] = '$';
    S[string_length] = '\0';

    close(fp);
    return S;
}

// read our alphabet file
void readAlphabetFile(char* file){
    FILE* fp;
    char* line = NULL;
    size_t len;

    fp = fopen(file, "r");

    if (!fp) { printf("Alphabet file %s not opened\n", file); exit(0); }
    printf("Alphabet file %s opened\n", file);

    alphabet_length = 0;
    int i = 0;

    while(getline(&line, &len, fp) != -1) {
        if (line[0] != '>') {
            while (line[i] != '\0' && line[i] != '\n') {
                if (line[i] != ' ') {
                    alphabet_length++;
                }
                i++;
            }
            i = 0;
        }
    }
    alphabet_length++;  // for $, per the way I set it up

    fclose(fp);
}

// read our reads file
char* readReadsFile(FILE* fp) {
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
            read[i] = '\0';
            return read;
        }
    }
    free(read);
    return NULL;
}

// read our parameter file
void readParamsFile(char* config_file) {    
    FILE* fp;
    char* line = NULL;
    size_t len;
    
    fp = fopen(config_file, "r");
    if (!fp) { printf("config file not opened\n"); exit(0); }
    printf("config file opened\n");

    // get params from config:
    getline(&line, &len, fp);   
    ma = getParam(line);                 // match
    getline(&line, &len, fp);           
    mi = getParam(line);                 // mismatch
    getline(&line, &len, fp);
    h = getParam(line);                  // opening penalty
    getline(&line, &len, fp);
    g = getParam(line);                  // extension penalty
}

// getParam() is fed a line from the config file and returns the value we care about
int getParam(char* line) {
    char* token  = strtok(line, " ");
    char* last = token;
    for (; token = strtok(NULL, " "); last = token);  // iterate past all tokens
    return atoi(last); // last token evaluated was the int
}

/* * * * * *  End read files  * * * * * */



// substition is the algorithm used by Smith-Waterman to see whether two strings match at a given point
int substitution(char a, char b) {
    if (a == b) {
        return ma;
    }
    else { 
        return mi;
    }
}

// find max of 3
int findMaxOfThree(int a, int b, int c) {
    int max = a;
    if (b > max) { max = b; }
    if (c > max) { max = c; }
    return max;
}

// find max of 4
int findMaxOfFour(int a, int b, int c, int d) {
    int max = a;
    if (b > max) { max = b; }
    if (c > max) { max = c; }
    if (d > max) { max = d; }
    return max;
}

// find max of 3 and give the direction (used in traceback)
direction findMaxDirection(int S, int D, int I) {
    int max = S;
    direction d = diagonal;
    if (D > max) { 
        max = D;
        d = up; 
    }
    if (I > max) { 
        max = I;
        d = left; 
    }
    if (max <= 0) {
        d = done;
    }
    return d;
}

// find time elapsed of function
double FindTime(struct timeval t1, struct timeval t2) {
    return (t2.tv_sec-t1.tv_sec) + ((t2.tv_usec-t1.tv_usec)/100000);
}


/* * * * * * Begin tree files * * * * * */

// creates tree, sets root to null, sets SL to self
tree* init() {
    tree* t = (tree*) malloc(sizeof(struct Tree));

    node* r = createNode();
    r->children = NULL;
    r->parent = NULL;
    r->depth = 0;
    r->id = 0;
    r->SL = r;
    r->parent = r;
    r->start_index = -1;
    r->end_index = -1;

    t->root = r;
    t->u = r;

    return t;
}


// creates a node
node* createNode() {
    node* n = (node*) malloc(sizeof(struct Node));
    n->children = NULL;
    n->depth = 0;
    n->id = 0;
    n->label_length = 0;
    n->start_index = -1;
    n->end_index = -1;
    n->parent = NULL;
    n->SL = NULL;

    num_nodes++;

    return n;
}


// to maintain lexicographical structure
void sortChildren(node* n, char* S) {
    int i,j;

    int len = 0;
    while(n->children[len] != NULL) { len++; }

    if (len < 2) { return; }    // nothing to sort

    // lazy bubble sort
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            if (i == j) { continue; }
            int a = n->children[i]->start_index;
            int b = n->children[j]->start_index;

            if (S[a] < S[b]) {
                node* temp = n->children[i];
                n->children[i] = n->children[j];
                n->children[j] = temp;
            }
        }
    } 
}

// FindPath algorithm. Offset measures where in S we are. 
node* findPath(tree* t, node* v, char* S, int offset) {
    if(v->children != NULL) {
        for (int i = 0; i < alphabet_length; i++) {     // for each child 
            if (v->children[i] != NULL) {                   // if such a child does exist
                node* child = v->children[i];
                if (S[child->start_index] == S[offset]) {           // if the first index of the label matches with s[offset]

                    int length = child->end_index - child->start_index + 1;

                    if(strncmp(S+child->start_index, S+offset, length) == 0) {        // if they are one in the same, we will have to jump to it 
                        return findPath(t,child, S, offset+length);
                    }
                    else {                                  // if they are not one in the same, we must break here along the edge and create two children.

                        // path does not exhaust label -- we must create new nodes
                        int matches = 0;
                        int j = offset;
                        int k = child->start_index;

                        while (S[j++] == S[k++]) {
                            matches++;      // count # matches on the current label
                        }

                        /*              
                                        What we're doing:
                                              parent
                                                | 
                                            new_parent
                                            /       \
                                        child      new_child
                        */

                        node* new_parent = createNode();
                        new_parent->id = upper_id++; // should be upper id
                        new_parent->start_index = child->start_index;
                        new_parent->end_index = child->start_index + matches - 1;
                        new_parent->parent = v;
                        new_parent->SL = NULL;
                        new_parent->depth = v->depth + (new_parent->end_index - new_parent->start_index + 1);

                        node* new_child = createNode();
                        new_child->id = lower_id++;
                        new_child->start_index = offset+matches;
                        new_child->end_index = offset + matches + strlen(S+offset+matches) - 1;
                        new_child->parent = new_parent;
                        new_child->SL = NULL;
                        new_child->depth = new_parent->depth + (new_child->end_index - new_child->start_index + 1);

                        child->parent = new_parent;
                        child->start_index = child->parent->end_index + 1;
                        child->depth = new_parent->depth + (child->end_index - child->start_index + 1);

                        new_parent->children = (node**) malloc(alphabet_length * sizeof(struct Node));
                        new_parent->children[0] = new_child;
                        new_parent->children[1] = child;

                        sortChildren(new_parent, S);

                        t->u = new_parent;

                        v->children[i] = new_parent;
                        sortChildren(v, S);
                        return new_parent;
                    }
                } 
            }
        }
    }


    // if we reach here, node has no adequate children, so we need to make one
    // this node will be a LEAF
    node* temp = createNode();
    temp->id = lower_id++;
    temp->children = NULL;
    temp->start_index = offset;
    temp->end_index = strlen(S) - 1;    // this end index is correct because the node
                                        // is a LEAF, thus label goes to the end of S
    temp->SL = NULL;
    temp->parent = v;
    temp->children = NULL;
    temp->depth = v->depth + (temp->end_index - temp->start_index + 1);

    // insert temp 
    int i = 0;

    if (v->children != NULL) {
        while (v->children[i]) { 
            i++; 
        };
    } 
    else {
        v->children = (node**) malloc(alphabet_length * sizeof(struct Node));
    }

    v->children[i] = temp;
    t->u = v;

    sortChildren(v, S);       
    return v;    // return the parent of the leaf
}


// prepares the suffix tree nodes (note: start_index and end_index as described are actually map_start and end_start)
int* prepareST(tree* t, int n) {
    int* A = (char*) malloc(n * sizeof(int));

    for (int i = 0; i < n; i++) {
        A[i] = -1;  // init content with -1 
    }
    
    DFS_PrepareST(t->root, A);
    return A;
}

// recursive part of tree preparation which starts at root and traverses entire tree
void DFS_PrepareST(node* n, int A[]){
    if (n == NULL) { return; }  
    if (!n->children) {         // if n has no children 
        A[next_index] = n->id;   
        if (n->depth >= x) {
            n->map_start = next_index;
            n->map_end = next_index;
        }
        next_index++;
        return;
    }

    // case: T is an internal node
    int i = 0;
    while (n->children[i] != NULL) {
        DFS_PrepareST(n->children[i], A);
        i++;
    }

    if (n->depth >= x) {
        node* u_left = n->children[0];
        node* u_right = n->children[i-1];

        n->map_start = u_left->map_start;
        n->map_end = u_right->map_end;
    }
}

// mapReads opens the read and output files, goes through all reads, aligns, and outputs
void mapReads(tree* t, char* S, int* A, char* reads_name) {
    FILE* fp, *output;

    // open reads file
    fp = fopen(reads_name, "r");
    if (!fp) { printf("READS file %s not opened\n", reads_name); exit(0); }
    printf("READS file %s opened\n", reads_name);

    // open/create output file for writing
    output = fopen("MappingResults_Peach.txt", "w");

    printf("Begin read alignment.\n");
    struct timeval startA, stopA;
    struct timeval startB, stopB;
    // begin sequencing reads
    for (int r_i = 0; r_i < 2000; r_i++) {
        char* read = readReadsFile(fp);
        if (read == NULL) { break; }

        int l = strlen(read);

        double max_length_coverage = 0.0;
        double max_percent_identity = 0.0;

        gettimeofday(&startA, NULL);
        int* L_i = findLoc(t->root, S, read, l);            // get our start and end indices in an array
        gettimeofday(&stopA, NULL);

        printf("A time=%d\n", stopA.tv_usec-startA.tv_usec);

        num_reads++;

        // read has been made 
        int start, end;
        start = L_i[0];
        end = L_i[1];

        num_starts = end-start+1;

        // if no deepest node found
        if (start <= 0 && end <= 0) {
            free(L_i);
            fprintf(output, "[Read %d] No hit found.\n", r_i+1);
            continue;
        }
        
        // printf("Started alignment:\n");
        for (int i = start; i <= end; i++) {
            gettimeofday(&startB, NULL);
            align(read, A[i]-1, S);  // note: it's A[i]-1 because A is 1-indexed but S is not.
            gettimeofday(&stopB, NULL);

            printf("B time=%d\n", stopB.tv_usec-startB.tv_usec);
            // for this particular read, if it qualifies as being "good"...
            if (percent_identity >= 0.90 && length_coverage >= 0.80) {
                // and if it is "max"...
                if (length_coverage > max_length_coverage) {
                    // set values to these
                    max_percent_identity = percent_identity;
                    max_length_coverage = length_coverage;
                    best_align_start = A[i]-1;
                }
            }
        }

        // if a max length has been set, then obviously we're good
        if (max_length_coverage >= 0.80) {
            fprintf(output, "[Read %d]\t%d\t%d\n", r_i+1, best_align_start, best_align_start+l);
        }
        else {
            fprintf(output, "[Read %d] No hit found. Percent identity = %.2f, length coverage = %.2f\n", r_i+1, percent_identity, length_coverage);
        }
        
        best_align_start = 0;

        percent_identity = 0.0;
        length_coverage = 0.01;

        free(L_i);
        free(read);
    }
    close(fp);
    close(output);
}


// align performs the smith-waterman alignment algorithm and traceback, and sets the values for length coverage and percent identity
void align(char* read, int _j, char* S) {

    int l = strlen(read);

    // S[start, end]
    int start = _j-l;
    int end = _j+l;
    int i, j;

    // fix out of bounds cases 
    if (start < 0) { 
        start = 0; 
    }
    if (end > strlen(S)-1) { 
        end = strlen(S)-1;      // -1 for $ 
    } 

    // create G
    int G_len = end-start+1;
    char* G = (char*) malloc(sizeof(char) * G_len);
    
    
    // populate our G string
    int k = 0;
    for (i = start; i <= end && i < strlen(S)-1; i++) {
        G[k++] = S[i]; 
    }
    G[k] = '\0';

    // m for rows, n for columns
    int m = strlen(read) + 1;
    int n = strlen(G) + 1;
 
    // init table    
    cell** table = (cell**) malloc(sizeof(cell*) * (m));         

    // put cols in table 
    for (i = 0; i < m; i++) { 
        table[i] = (cell*) malloc(sizeof(cell) * (n));            
    }

    // init
    cell* base = &table[0][0];
    base->S = base->I = base->D = base->i = base->j = 0;
    base->dir = done;


    // @ cells (i, 0)
    for (i = 1; i < m; i++) {           // on each col of 0th row
        cell* temp = &table[i][0];
        temp->S = temp->I = temp->D = 0;
        temp->i = i;
        temp->j = 0;
        temp->dir = done;
    }

    // @ cells (0,j) 
    for (j = 1; j < n; j++) {           // on each row of 0th col
        cell* temp = &table[0][j];
        temp->S = temp->D = temp->I = 0;
        temp->i = 0;
        temp->j = j;
        temp->dir = done;
    }

    char* s1 = read;
    char* s2 = G;   

    cell* max_cell = base;

    // forward computation of smith-waterman
    for (i = 1; i < m; i++) { 
        for (j = 1; j < n; j++) {
            cell* curr = &table[i][j];
            curr->i = i;
            curr->j = j;

            cell* temp;

            // set temp to T(i-1,j-1)
            temp = &table[i-1][j-1];
            curr->S = findMaxOfThree(temp->S, temp->D, temp->I) + substitution(s1[i-1], s2[j-1]);
            if (curr->S < 0) { curr->S = 0; }


            // calculate for D by viewing top cell at table[i-1, j]
            temp = &table[i-1][j];
            curr->D = findMaxOfFour(
                                        temp->D + g, 
                                        temp->S + h + g, 
                                        temp->I + h + g,
                                        0
                                    );

            // set temp to above  
            temp = &table[i][j-1];
            curr->I = findMaxOfFour(   
                                        temp->I + g,        
                                        temp->D + h + g, 
                                        temp->S + h + g,
                                        0 
                                    );

            // here we use findMaxOfThree to find the max between S, D, and I
            int max_of_dirs = findMaxOfThree( 
                                        curr->S, 
                                        curr->D, 
                                        curr->I
                                     );

            // if the max of this cell's (S or D or I) is greater than the max cell's 
            if (max_of_dirs > findMaxOfThree(max_cell->D, max_cell->I, max_cell->S)) { 
                max_cell = &table[i][j];        // make this the new max_cell
            }
        }
    }

    int optimal_score = findMaxOfThree(max_cell->D, max_cell->I, max_cell->S);

    // traceback start
    direction d;
    cell* retrace = max_cell;

    // find our first direction just by looking at our neighbors
    if (optimal_score == max_cell->D) {
        d = up;
    }
    else if (optimal_score == max_cell->I) {
        d = left;
    }
    else if (optimal_score == max_cell->S) {
        d = diagonal;
    }
    
    int matches = 0;
    int mismatches = 0;
    int gaps = 0;

    while (1) {

        /* 
        Strategy:
            If max was D:
                Look at cell (i-1,j). 
                For that cell, compute next max using D rules:
                    Max (D(i-1,j) + g, S(i-1,j) + h + g, I(i-1,j) + h + g)
                Update curr_cell = (i-1,j) 
                Set direction = whatever was max (S, I, or D)
                Go onto next iteration 
            If max was I or S:
                Same procedure as above, but with I or S equations
        */

        i = retrace->i;
        j = retrace->j;
        
        if (d == done || i == 0 || j == 0) { 
            break; 
        }

        // D
        if (d == up) {
            retrace = &table[i-1][j];
            gaps++;

            // S, D, I
            d = findMaxDirection(
                                    retrace->S + h + g,
                                    retrace->D + g,
                                    retrace->I + h + g
                                );
        }

        // I
        else if (d == left) {
            retrace = &table[i][j-1];
            gaps++;
            // S, D, I
            d = findMaxDirection(
                                    retrace->S + h + g,
                                    retrace->D + h + g,
                                    retrace->I + g
                                );
        }

        // S
        else if (d == diagonal) {
            retrace = &table[i-1][j-1];
            int sub = substitution(s1[i-1], s2[j-1]);
            if (s1[i-1] == s2[j-1]) {
                matches++;
            }
            else {
                mismatches++;
            }

            // S, D, I
            d = findMaxDirection(
                                    retrace->S + sub,           
                                    retrace->D + sub,
                                    retrace->I + sub
                                );
        }

    }

    // record % identity & length coverage
    int alignlen = matches + mismatches + gaps;
    percent_identity = matches / (double)alignlen;
    length_coverage = alignlen / (double)l;

    for (i = 0; i < m; i++) { 
        free(table[i]);          
    }

    free(table);
    free(G);
}


// findloc returns the start and end node id's for the deepest available node in a subtree, if such a node exists 
int* findLoc(node* root, char* S, char* read, int l) {
    int read_ptr;

    node* deepest_node = NULL;
    for (int i = 0; i < l; i++) {
        read_ptr = 0;
        node* temp = findLocSearch(root, read+i, S, read_ptr);
        if (temp != NULL && max_length >= x) {
            deepest_node = temp;
        }
    }

    // ids = {start, end}
    int* ids = (int*) malloc(sizeof(int) * 2);

    if (deepest_node == NULL) {
        ids[0] = -1;
        ids[1] = -1;
    }
    else { 
        ids[0] = deepest_node->map_start;
        ids[1] = deepest_node->map_end;
    }

    max_length = 0;
    return ids;
}

// recursive read-only findpath used for finding the longest read_ptr of a read compared to the genome. 
node* findLocSearch(node* v, char* read, char* S, int read_ptr) {
    if (read_ptr < strlen(read)) {
        for (int i = 0; i < alphabet_length; i++) {     // for each child 
            if (v->children[i] != NULL) {                   // if such a child does exist    
                node* child = v->children[i];
                if (S[child->start_index] == read[read_ptr]) {           // if the first index of the label matches with s[offset]
                    int length = child->end_index - child->start_index + 1;     // get length of the label 
                    

                    if (strncmp(S+child->start_index, read+read_ptr, length) == 0) {        // if they are one in the same, we will have to jump to it 
                        // if we can move to that child
                        node* temp = findLocSearch(child, read, S, read_ptr+length);     // read_ptr + length of label probably 
                        return temp;
                    }
                    
                    // case B: label breaks between this one and its child, so count the # of matches
                    else { 
                        // let u be the internal node visisted last -- this one!
                        int j = child->start_index;

                        // add to read_ptr the # of matches along the child label 
                        while (S[j] == read[read_ptr]) {
                            read_ptr++;
                            j++;
                        }

                        if (read_ptr > max_length) {
                            max_length = read_ptr;
                            return child;
                        }
                        else {
                            return NULL;
                        }
                    }
                    break;
                } 
            }
        }
    }
    if (read_ptr > max_length) {
        max_length = read_ptr;
        return v;
    }
    else {
        return NULL;
    }
}

//  $ <test executable> <sequnce file> <read file> <alphabet file> 
int main(int argc, char** argv) {

    if (argc < 4) {
        printf("Example functionality: ./a.out <s1.fasta> <reads.fasta> <alphabet.txt>\n");
        exit(0);
    }

    char* sequence_file = argv[1];
    char* reads_file = argv[2];
    char* alphabet_file = argv[3];

    char* S = readFASTA(sequence_file);
    readAlphabetFile(alphabet_file);
    readParamsFile("parameters.config");

    printf("S length + '$' = %d\t alphabet length = %d\n", strlen(S), alphabet_length);
    printf("Config:\n\tma: %d\n\tmi: %d\n\th: %d\n\tg: %d\n", ma, mi, h, g);

    upper_id = strlen(S) + 1;       // upper_id is used for internal nodes -- it'll start at 1 + the number of leaves
    int n = strlen(S);

    // create tree struct
    tree* t = init();

    struct timeval startBuildTree, stopBuildTree, startPrepareST, stopPrepareST, startMapReads, stopMapReads;
     
    printf("Tree is building.\n");
    gettimeofday(&startBuildTree, NULL);
    // loop to insert new nodes into tree
    for (int i = 0; i < n; i++) {
        if (i % 20000 == 0) {
            printf("--- ITERATION #%d, CHARACTER=%c ---\n", i, S[i]);
        }
        findPath(t, t->root, S, i);
    }
    gettimeofday(&stopBuildTree, NULL);
    printf("\n");

    gettimeofday(&startPrepareST, NULL);
    int* A = prepareST(t, n);
    gettimeofday(&stopPrepareST, NULL);

    printf("\n");

    printf("--mapreads:--\n");
    gettimeofday(&startMapReads, NULL);
    mapReads(t, S, A, reads_file);
    gettimeofday(&stopMapReads, NULL);

    double tree_time = FindTime(startBuildTree, stopBuildTree);
    double st_time = FindTime(startPrepareST, stopPrepareST);
    double map_time = FindTime(startMapReads, stopMapReads);

    printf("-- -- -- -- --\n");
    printf("MapReads program finished.\n");
    printf("Statistics:\nx=%d\n", x);
    printf("# internal nodes in tree: %d\n# leaves in tree: %d\n", upper_id-1, lower_id-1);
    printf("Candidate loci per read: %f", num_starts/(double)num_reads);
    printf("# total nodes in tree: %d\n", num_nodes);
    printf("Time elapsed:\n\t Tree build: %.5f\n\tPrepare ST: %.5f\n\tMapReads: %.5f\n", tree_time, st_time, map_time);
    printf("Total starts = %d", num_starts);

    
    /* Sample output results: 
        # internal nodes in tree: 8945575
        # leaves in tree: 5362496
        # total nodes in tree: 8945576
        Candidate loci per read: 1.057386
        Time elapsed:
            Tree build: 349.00000
            Prepare ST: 8.00000
            MapReads: 11489.00000
        Total starts = 528693
    */

    free(t);
    return 0;
}