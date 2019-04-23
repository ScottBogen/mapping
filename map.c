/* 

    What needs to be done:
        Finalize format: (j0, j1) or no read found 
        Write to a file 
        Record times
        Read from parameter file
        Put things into separate files -- ex: .h, multiple .c files

        Clean up functions <-- current

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

// direction used for traceback
typedef enum Direction {
    up,
    left,
    diagonal,
    done            // done only applies to t(0,0)
} direction;

typedef struct Cell {
    int S;
    int I;
    int D;              
    int i;              // i index 
    int j;              // j index 
    direction dir;      // traceback direction 
} cell;

char* reads_file_name = "Peach_simulated_reads.fasta";      // TODO: make this a param


int ma = 1;
int mi = -2;
int h = -5;
int g = -1;

int x = 4;

int duplicate_optimals = 0;
int duplicate_max_cells = 0;


int best_align_start = 0;
int best_optimal_score = 0;


double percent_identity;
double length_coverage;

int max_length = 0;


int next_index = 0;
int last_max = 0;


int lower_id = 1;   // stores ID of leaves (1 to strlen(S))
int upper_id;       // stores ID of internal nodes

int alphabet_length;

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

    // for testing
    int map_start;
    int map_end;
} node;


typedef struct Tree {
    struct Node* root;
    struct Node* u;
} tree;



/* Move to .h files and add all the new functions as well */
tree* init();    // create tree
node* createNode();
void insert(tree* t, char* S, int i);       // insert new node into tree from string S at offset i
                                            // ex: t, "banana$", "3" = insert("nana$")
void sortChildren(node* n, char* S);
node* findPath(tree* t, node* v, char* S, int offset);
node* findLocSearch(node* v, char* read, char* S, int read_ptr);
int* findLoc(node* root, char* S, char* read, int l);


int substitution(char a, char b) {
    if (a == b) {
        return ma;
    }
    else { 
        return mi;
    }
}


int findMax(int a, int b, int c) {
    int max = a;
    if (b > max) { max = b; }
    if (c > max) { max = c; }
    return max;
}

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

int findMaxOfFour(int a, int b, int c, int d) {
    int max = a;
    if (b > max) { max = b; }
    if (c > max) { max = c; }
    if (d > max) { max = d; }
    return max;
}

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
    return n;
}

void writeToFile(node* n, char* S) {
    FILE* fp;
    fp = fopen("output.txt", "w");
    fclose(fp);
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

char* readFile(char* file_name) {
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

int* prepareST(tree* t, int n) {
    int* A = (char*) malloc(n * sizeof(int));

    for (int i = 0; i < n; i++) {
        A[i] = -1;  // init content with -1 
    }
    
    DFS_PrepareST(t->root, A);
    return A;
}

void DFS_PrepareST(node* n, int A[]){
    if (n == NULL) { return; }  // shouldn't ever happen
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
            return read;
        }
    }
    free(read);
    return NULL;
}

void mapReads(tree* t, char* S, int* A, char* reads_name) {
    
    FILE* fp;
    // open reads file
    fp = fopen(reads_name, "r");
    if (!fp) { printf("READS file %s not opened\n", reads_file_name); exit(0); }
    printf("READS file %s opened\n", reads_file_name);

    // begin sequencing reads
    for (int r_i = 0; ; r_i++) {
        
        char* read = readReadsFile(fp);
        if (read == NULL) { break; }


        int l = strlen(read);
        int* L_i = findLoc(t->root, S, read, l);

        // read has been made 
        int start, end;

        start = L_i[0];
        end = L_i[1];

        // if no deepest node found
        if (start <= 0 && end <= 0) {
            free(L_i);
            continue;
        }
        
        // printf("Started alignment:\n");
        for (int i = start; i <= end; i++) {
            align(read, A[i]-1, S, l);  // note: it's A[i]-1 because A is 1-indexed but S is not.
        }

        //printf("Best optimal score: %d\nBest alignment start: %d\n", best_optimal_score, best_align_start);
        //printf("[Read %d] Best start = %d\n", r_i+1, best_align_start);


        if (percent_identity >= .9 && length_coverage >= .8) {

            char* a = "";
            char* b = "";

            if (duplicate_max_cells > 0) { 
                a = "{a}";
            }

            if (duplicate_optimals > 0) {
                b = "{b}";
            }

            printf("[Read %d] Best start = %d, best end = %d\t", r_i+1, best_align_start, best_align_start + (2*l));
            printf("%s\t%s\n", a, b);
        }
        else {
            printf("[Read %d] No hit found\n");
        }
        
        best_optimal_score = 0;
        best_align_start = 0;
        duplicate_optimals = 0;
        duplicate_max_cells = 0;

        free(L_i);
        free(read);
    }

    close(fp);
}


void align(char* read, int _j, char* S, int l) {

    int start = _j-l;    
    int end = _j+l;
    int i, j;

    int max_of_curr, max_of_max;

    // fix out of bounds cases 
    if (start < 0) { 
        start = 0; 
    }
    if (end > strlen(S)-1) { 
        end = strlen(S)-1;      // -1 for $ 
    } 

    int G_len = end-start+1;
    char* G = (char*) malloc(sizeof(char) * G_len);
    
    int k = 0;
    
    // populate G
    for (i = start; i <= end && i < strlen(S)-1; i++) {
        G[k++] = S[i]; 
    }
    
    G[k] = '\0';



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
    int current_max = 0;

    // forward computation of smith-waterman
    for (i = 1; i < m; i++) { 
        for (j = 1; j < n; j++) {
            cell* curr = &table[i][j];
            curr->i = i;
            curr->j = j;

            cell* temp;

            // set temp to T(i-1,j-1)
            temp = &table[i-1][j-1];
            curr->S = findMax(temp->S, temp->D, temp->I) + substitution(s1[i-1], s2[j-1]);
            // S can't be negative on the table, so disallow it
            if (curr->S < 0) {
                curr->S = 0;
            }

            // calculate for D by viewing top cell at table[i-1, j]
            temp = &table[i-1][j];
            curr->D = findMaxOfFour(
                                        temp->D + g, 
                                        temp->S + h + g, 
                                        temp->I + h + g,
                                        0
                                    );

            // set temp to I by viewing left cell at table[i, j-1]  
            temp = &table[i][j-1];
            curr->I = findMaxOfFour(   
                                        temp->I + g,        
                                        temp->D + h + g, 
                                        temp->S + h + g,
                                        0 
                                   );

            // here we use findMax to find the max between S, D, and I
            max_of_curr = findMax(curr->S, curr->D, curr->I);
            max_of_max = findMax(max_cell->D, max_cell->I, max_cell->S);

            // if the max of this cell's (S or D or I) is greater than the max cell's, make this the new max
            if (max_of_curr > max_of_max) { 
                last_max = current_max;
                current_max = max_of_curr;
                max_cell = curr; 
            }
        }
    }   // end forward computation


    // if the last max we calculated was the same as this one, there are multiple max cells in the table, so record that for statistical purposes
    if (last_max == current_max) {
        duplicate_max_cells++;
    }

    // begin retrace

    direction d;
    // set retrace cell to the max cell 
    cell* retrace = max_cell;

    // start the walk back simply by finding where the max came from
    if (max_of_max == max_cell->D) {
        d = up;
    }
    else if (max_of_max == max_cell->I) {
        d = left;
    }
    else if (max_of_max == max_cell->S) {
        d = diagonal;
    }
    
    int matches = 0;
    int mismatches = 0;
    int gaps = 0;

    while (1) {         // breaks upon reaching a cell with value 0
        /* 
        If max was D:
            Look at cell (i-1,j). 
            For that cell, compute next max using D rules:
                Max (D(i-1,j) + g, S(i-1,j) + h + g, I(i-1,j) + h + g)
            Update curr_cell = (i-1,j) 
            Set direction = whatever was max (S, I, or D)
            Go onto next iteration 
        If max was I or S:
            same procedure for above, but using I or S equations
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
    }   // end traceback


    // if the the new max is larger than the old max, make this the new max 
    if (current_max > best_optimal_score) { 
        best_optimal_score = current_max;
        best_align_start = _j;
    }

    // if there's a tie between the old max and new max, we have a duplicate
    else if (current_max == best_optimal_score) {
        best_optimal_score = current_max;
        best_align_start = _j;
        duplicate_optimals++;
    }

    // calculate alignment statistics
    int alignlen = matches + mismatches + gaps;
    percent_identity = (double)matches / (double) alignlen;
    length_coverage = (double) alignlen / l;

    // clear rows of table
    for (i = 0; i < m; i++) { 
        free(table[i]);          
    }

    // clear table and genome subsequence
    free(table);
    free(G);
}

// findloc returns the start and end node id's for the deepest available node in a subtree, if such a node exists 
int* findLoc(node* root, char* S, char* read, int l) {
    int read_ptr;

    // 0: abcde
    // 1:  bcde
    // 2:   cde
    //      ...

    node* deepest_node = NULL;
    for (int i = 0; i < l; i++) {
        read_ptr = 0;
        node* temp = findLocSearch(root, read+i, S, read_ptr);
        if (temp != NULL && max_length >= x) { // edited: used to have && max_length >= x
            //printf("deepest node set at length=%d\n", max_length);
            deepest_node = temp;
        }
    }

    //printf("findLoc says max_len = %d, deepestNode = %d\n", max_length, deepest_node->id);
    //printf("start=%d, end=%d\n", deepest_node->map_start, deepest_node->map_end);

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
    // case A: broke at edge
    if (read_ptr < strlen(read)) {
        for (int i = 0; i < alphabet_length; i++) {     // for each child 
            if (v->children[i] != NULL) {                   // if such a child does exist    
                node* child = v->children[i];
                if (S[child->start_index] == read[read_ptr]) {           // if the first index of the label matches with s[offset]
                    int length = child->end_index - child->start_index + 1;     // get length of the label 
                    
                    if (strncmp(S+child->start_index, read+read_ptr, length) == 0) {        // if they are one in the same, we will have to jump to it 
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
    return NULL;
}

//  $ <test executable> <input file containing sequence s> <input alphabet file> 
int main(int argc, char** argv) {

    if (argc < 4) {
        printf("Example functionality: ./a.out <s1.fasta> <reads.fasta> <alphabet.txt>\n");
        exit(0);
    }

    char* sequence_file = argv[1];
    char* reads_file = argv[2];
    char* alphabet_file = argv[3];

    char* S = readFile(sequence_file);
    readAlphabetFile(alphabet_file);

    printf("S length = %d\t alphabet length = %d\n", strlen(S), alphabet_length);
    
    upper_id = strlen(S) + 1;       // upper_id is used for internal nodes -- it'll start at 1 + the number of leaves
    printf("upper id = %d\n", upper_id);

    int n = strlen(S);

    // create tree struct
    tree* t = init();

    struct timeval start, stop;

    printf("Press ENTER to build tree\n");
    //getchar();
    printf("\n");
     
    gettimeofday(&start, NULL);
    
    // loop to insert new nodes into tree
    for (int i = 0; i < n; i++) {
        // uncomment the two below lines to see if/when the program breaks   
        if (i % 20000 == 0) {
            printf("--- ITERATION #%d, CHARACTER=%c ---\n", i, S[i]);
        }
        findPath(t, t->root, S, i);
    }
    
    gettimeofday(&stop, NULL);
    printf("\n");

    int* A = prepareST(t, n);
     
    printf("\n");
    //writeToFile(t->root, S);

    //searchTree(t->root, S);

    printf("--mapreads:--\n");
    mapReads(t, S, A, reads_file);
    printf("-- -- -- -- --\n");

    printf("McCreight's Algorithm program finished.\n");
    printf("Options used:\n\tx=%d\n", x);
    free(t);
    return 0;
}