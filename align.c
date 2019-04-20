#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

// end results
int matches = 0;
int mismatches = 0;

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


int num_nodes = 0;
int max_depth = 0;
int max_length = 0;

int ma = 1;
int mi = -2;

int x = 0;
int m = 100;

int next_index = 0;

int lower_id = 1;   // stores ID of leaves (1 to strlen(S))
int upper_id;       // stores ID of internal nodes

int alphabet_length;
int line_iterator = 1;  // used to ensure 10 li

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


tree* init();    // create tree
node* createNode();

void insert(tree* t, char* S, int i);       // insert new node into tree from string S at offset i
                                            // ex: t, "banana$", "3" = insert("nana$")
void display(node* u);
void enumerate(node* n, char* S);
void BWT(tree* t);
void sortChildren(node* n, char* S);
node* findPath(tree* t, node* v, char* S, int offset);
void printNodeInfo(node* temp, char* S);
node* nodeHops(node* n, char* S, char* beta, int offset); 
node* nodeHops2(node* n, char* S, char* beta, int offset);
node* findLocSearch(node* v, char* read, char* S, int read_ptr);
int* findLoc(node* root, char* S, char* read);

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

// does ctrl+c, ctrl+v count as code reuse? 
int findMaxLocal(int a, int b, int c, int d) {
    int max = a;
    if (b > max) { max = b; }
    if (c > max) { max = c; }
    if (d > max) { max = d; }
    return max;
}


int main() {


}

void align(char* G, char* read);

    int G_len = strlen(G);
    printf("Glen = %d\n", G_len);
    char* G = (char*) malloc(sizeof(char) * G_len);

    int k = 0;
    for (i = start; i <= end && i < strlen(S)-1; i++) {
        G[k++] = S[i]; 
    }
    G[k] = '\0';

    printf("align: G = %s\n", G);
    printf("align: read = %s\n", read);

    /* params setup:
        0: ma
        1: mi
        2: h
        3: g
    */ 

    int h = -5;
    int g = -1;

    /* 
        create read between read and S[start... end]

        perform SW on it (local alignment)
        calculate matches and alignlen either by forward walk or traceback
    */

    int m = strlen(read) + 1;
    int n = strlen(G) + 1;

    printf("m = %d,   n = %d\n", m, n);

 
    // init table    
    cell** table = (cell**) malloc(sizeof(cell*) * (m));         

    // put cols in table 
    for (i = 0; i < m; i++) { 
        table[i] = (cell*) malloc(sizeof(cell) * (n));            
    }

    printf("\n\n");

    /* 
        keeping i constant and j moving produces a row
        keeping i moving and j constant produces a col

        which means i corresponds to a row#, j corresponds to a col#
    */

    printf("Initialization\n");

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
        //printf("m loop #%d\n", i);
        for (j = 1; j < n; j++) {
            //printf("\tn loop #%d\n", j);
            cell* curr = &table[i][j];
            curr->i = i;
            curr->j = j;

            cell* temp;

            // set temp to T(i-1,j-1)
            temp = &table[i-1][j-1];
            table[i][j].S = findMax(temp->S, temp->D, temp->I) + substitution(s1[i-1], s2[j-1]);

            // calculate for D by viewing top cell at table[i-1, j]
            temp = &table[i-1][j];
            table[i][j].D = findMaxLocal(
                                        temp->D + g, 
                                        temp->S + h + g, 
                                        temp->I + h + g,
                                        0
                                        );

            // set temp to above  
            temp = &table[i][j-1];
            table[i][j].I = findMaxLocal(   
                                        temp->I + g,        
                                        temp->D + h + g, 
                                        temp->S + h + g,
                                        0 
                                        );

            // here we use findMax to find the max between S, D, and I
            int max_of_dirs = findMax( table[i][j].S, 
                                       table[i][j].D, 
                                       table[i][j].I
                                     );
            
            if (i == 4 && j == 3) {
                printf("MAX of dirs = %d\n", max_of_dirs);
                printf(" S = %d \n I = %d \n D = %d \n", table[i][j].S, table[i][j].I, table[i][j].D);
            
                printf("t[i-1][j]: D+g = %d, S+h+g = %d, I+h+g = %d\n", 
                        table[i-1][j].D + g, table[i-1][j].S+h+g, table[i-1][j].I+h+g);
            
            }
            
            if (max_of_dirs < 0) {
                max_of_dirs = 0;
                curr->S = curr->I = curr->D = 0;
            }

            // if the max of this cell's (S or D or I) is greater than the max cell's 
            if (max_of_dirs > findMax(max_cell->D, max_cell->I, max_cell->S)) { 
                max_cell = &table[i][j]; // make this the new max_cell
                printf("max cell changed to table[%d][%d]\n", i, j);
            }

            if (max_of_dirs == 0) {
                table[i][j].dir = done;
            }

            else if (max_of_dirs == table[i][j].D) {
                table[i][j].dir = up;
                printf("set to up\n");
            }
            
            else if (max_of_dirs == table[i][j].I) {
                table[i][j].dir = left;
                printf("set to left\n");
            }

            // none of these cases actually work
            else if (max_of_dirs == table[i][j].S) {
                table[i][j].dir = diagonal;
            }  


        }
    }


    // put cols in table 
    // for (int i = 0; i <= n; i++) { 
    //    table[i] = (cell*) malloc(sizeof(cell) * (m+1));            
    // } 
    

    //printf("Matches: %d\nMismatches: %d\n", matches, mismatches);
    //printf("Identities: %d/%d (%f%)\n", matches, max_length, (float) matches/max_length * 100.0);

    printf("Final table:\n      ");
    for (i = 0; i < strlen(G); i++) {
        printf(" %5c", G[i]);
    }
    printf("\n");

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf(" %5d", findMax(table[i][j].D, table[i][j].I, table[i][j].S));
        }
        printf("\n");
    }
    printf("\n");

    printf("Direction table:\n      ");
    for (i = 0; i < strlen(G); i++) {
        printf(" %5c", G[i]);
    }
    printf("\n");

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            switch(table[i][j].dir) {
                case up:
                    printf("%6s", "U");
                    break;
                case left:
                    printf("%6s", "L");
                    break;
                case diagonal:
                    printf("%6s", "D");
                    break;
                case done:
                    printf("%6s", "0");
                    break;
            }

        }
        printf("\n");
    }
    printf("\n");

    printf("Local optimal score: %d\n", findMax(max_cell->D, max_cell->I, max_cell->S));
    printf("Maximum cell position: Table[%d][%d]\n", max_cell->i, max_cell->j);
    printf("\n");


    // let's try a traceback in a different way -- not using the direction matrix

    // "d is saying look up, but not saying which of the 3 up values to use"
    // "if it came from d i-1,j you just did g, otherwise you did + h + g"


    /* 
        Here's what I think so far:

        Start at max_cell. Pick the max value between S, D, and I, then enter a loop:

        If max was 0, just end the retrace.

        If max was D:
            Look at cell (i-1,j). 
            For that cell, compute next max using D rules:
                Max (D(i-1,j) + g, S(i-1,j) + h + g, I(i-1,j) + h + g)
                Update curr_cell = (i-1,j) 
                Set direction = whatever was max (S, I, or D)
                Go onto next iteration 
        If max was I:
            Look at cell (i, j-1)
            For that cell, compute next max using I rules:
                Max (I(i, j-1) + g, S(i,j-1) + h + g, D(i,j-1) + h + g)
                Update curr_cell = (i, j-1)
                Set direction = whatever was max (S, I, or D)
                Go onto next iteration

        If max was S:
            (same procedure as above)
        
    */

    direction d;
    cell* retrace = max_cell;
    int val = findMax(max_cell->D, max_cell->I, max_cell->S);
    if (val == max_cell->D) {
        d = up;
    }
    else if (val == max_cell->I) {
        d = left;
    }
    else if (val == max_cell->S) {
        d = diagonal;
    }
    
    while (1) {

        /* 
        If max was D:
            Look at cell (i-1,j). 
            For that cell, compute next max using D rules:
                Max (D(i-1,j) + g, S(i-1,j) + h + g, I(i-1,j) + h + g)
            Update curr_cell = (i-1,j) 
            Set direction = whatever was max (S, I, or D)
            Go onto next iteration 
        If max was I:
            Look at cell (i, j-1)
            For that cell, compute next max using I rules:
                Max (I(i, j-1) + g, S(i,j-1) + h + g, D(i,j-1) + h + g)
            Update curr_cell = (i, j-1)
            Set direction = whatever was max (S, I, or D)
            Go onto next iteration

        If max was S:
            (same procedure as above)
        */
        
        i = retrace->i;
        j = retrace->j;
        
        if (d == done) { 
            printf("Done at [%d][%d]\n", i, j);
            break; 
        }

        // D
        if (d == up) {
            retrace = &table[i-1][j];

            // S, D, I
            d = findMaxDirection(
                                    retrace->S + h + g,
                                    retrace->D + g,
                                    retrace->I + h + g
                                );
            printf("Up from [%d][%d]\n", i, j);
        }

        // I
        else if (d == left) {
            retrace = &table[i][j-1];

            // S, D, I
            d = findMaxDirection(
                                    retrace->S + h + g,
                                    retrace->D + h + g,
                                    retrace->I + g
                                );
            printf("Left from [%d][%d]\n", i, j);
        }

        // S
        else if (d == diagonal) {
            retrace = &table[i-1][j-1];
            int sub = substitution(s1[i-2], s2[j-2]);

            // S, D, I
            d = findMaxDirection(
                                    retrace->S + sub,
                                    retrace->D + sub,
                                    retrace->I + sub
                                );
            printf("Diagonal from [%d][%d]\n", i, j);
        }
    }

    free(table);
    free(G);