#include <GraphBLAS.h>
#include <LAGraph.h>
#include <stdio.h>
#include <unistd.h>

#define error(x) fprintf(stderr, "%s\n", x)

int load_dfa(const char *filename, GrB_Matrix *A)
{
    FILE *f = fopen(filename, "r");
    GrB_Info code = LAGraph_MMRead(A, f, NULL);
    fclose(f);
    if (code != GrB_SUCCESS)
    {
        error("failed reading;");
        return code;
    }
    return 0;
}

int load_states(const char *filename, GrB_Vector *V)
{
    GrB_Matrix mtx_data;
    FILE *f = fopen(filename, "r");
    GrB_Info code = LAGraph_MMRead(&mtx_data, f, NULL);
    fclose(f);

    GrB_Index nrows;
    GrB_Matrix_nrows(&nrows, mtx_data);

    code = GrB_Vector_new(V, GrB_BOOL, nrows);

    GrB_Col_extract(*V, NULL, NULL, mtx_data, NULL, nrows, 0, NULL);

    GrB_Matrix_free(&mtx_data);

    return code;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        error("not enough input");
        return -1;
    }

    GrB_init(GrB_NONBLOCKING);
    LAGraph_Init(NULL);

    GrB_Info code;

    // load dfas

    GrB_Matrix first_dfa_label1, first_dfa_label2;
    GrB_Matrix second_dfa_label1, second_dfa_label2;

    if (load_dfa(argv[1], &first_dfa_label1)) return 1;
    if (load_dfa(argv[2], &first_dfa_label2)) return 2;
    if (load_dfa(argv[3], &second_dfa_label1)) return 3;
    if (load_dfa(argv[4], &second_dfa_label2)) return 4;

    GrB_Index first_states, second_states;
    GrB_Matrix_nrows(&first_states, first_dfa_label1);
    GrB_Matrix_nrows(&second_states, second_dfa_label1);
    GrB_Index intersection_states = first_states * second_states;

    GrB_Matrix first_dfa, second_dfa;

    GrB_Matrix_new(&first_dfa, GrB_INT64, first_states, first_states);
    GrB_Matrix_new(&second_dfa, GrB_INT64, second_states, second_states);

    GrB_eWiseAdd(first_dfa, NULL, NULL, GrB_PLUS_INT64, first_dfa_label2, first_dfa_label2, NULL);
    GrB_eWiseAdd(first_dfa, NULL, NULL, GrB_PLUS_INT64, first_dfa, first_dfa_label1, NULL);

    GrB_eWiseAdd(second_dfa, NULL, NULL, GrB_PLUS_INT64, second_dfa_label2, second_dfa_label2, NULL);
    GrB_eWiseAdd(second_dfa, NULL, NULL, GrB_PLUS_INT64, second_dfa, second_dfa_label1, NULL);

    // TODO: load starting and accepting states for non-trivial automaton

    // two boolean Kroneckers for each label and element-wise OR

    GrB_Matrix T;
    GrB_Matrix_new(&T, GrB_BOOL, intersection_states, intersection_states);

    code = GrB_kronecker(T, NULL, NULL, GrB_LAND, first_dfa_label1, second_dfa_label1, NULL);
    if (code != GrB_SUCCESS) { error("kronecker on first label error"); return 1; }
    code = GrB_kronecker(T, NULL, GrB_LOR, GrB_LAND, first_dfa_label2, second_dfa_label2, NULL);
    if (code != GrB_SUCCESS) { error("kronecker on second label error"); return 1; }

    GrB_Index T_nvals;
    GrB_Matrix_nvals(&T_nvals, T);
    //printf("T nvals: %d\n", T_nvals);

    GrB_Matrix_free(&first_dfa_label1);
    GrB_Matrix_free(&first_dfa_label2);
    GrB_Matrix_free(&second_dfa_label1);
    GrB_Matrix_free(&second_dfa_label2);

    // MS BFS over matrix T to get reachable from starting states of automaton

    LAGraph_Graph T_graph;

    code = LAGraph_New(&T_graph, &T, LAGraph_ADJACENCY_DIRECTED, NULL);
    if (code != GrB_SUCCESS) { error("failed to create graph"); return code; }

    GrB_Vector level_f;

    GrB_Vector_new(&level_f, GrB_BOOL, intersection_states);

    //printf("intersection states: %d\n", intersection_states);

    code = LAGr_BreadthFirstSearch(&level_f, NULL, T_graph, 0, NULL);
    if (code != GrB_SUCCESS) { error("failed forward BFS"); return 1; }
    code = GrB_apply(level_f, NULL, NULL, GxB_ONE_BOOL, level_f, NULL);
    if (code != GrB_SUCCESS) { error("failed conversion for forward bfs"); return 1; }

    // is it really needed?

    // BFS over transposed matrix T to get vertices, which reach any accepting states of automaton

    GrB_Vector level_b;

    GrB_Vector_new(&level_b, GrB_BOOL, intersection_states);

    code = GrB_transpose(T_graph->A, NULL, NULL, T_graph->A, NULL);
    if (code != GrB_SUCCESS) { error("transpose failed"); return code; }

    code = LAGr_BreadthFirstSearch(&level_b, NULL, T_graph, 0, NULL);
    if (code != GrB_SUCCESS) { error("failed reversed BFS"); return 1; }
    code = GrB_apply(level_b, NULL, NULL, GxB_ONE_BOOL, level_b, NULL);
    if (code != GrB_SUCCESS) { error("failed conversion for reverse bfs"); return 1; }

    GrB_Matrix transposed_level_b;
    GrB_Matrix_new(&transposed_level_b, GrB_BOOL, 1, intersection_states);

    code = GrB_transpose(transposed_level_b, NULL, NULL, level_b, NULL);
    if (code != GrB_SUCCESS) { error("transposed rev bfs fail"); return code; }

    // get Mask matrix with mxm of two vectors from BFS and BFS over transposed T ---
    // to count kronecker only for non-zero transitions

    GrB_Matrix Mask;
    GrB_Matrix_new(&Mask, GrB_BOOL, intersection_states, intersection_states);

    code = GrB_transpose(T_graph->A, NULL, NULL, T_graph->A, NULL);
    if (code != GrB_SUCCESS) { error("transpose failed"); return code; }

    GrB_Index level_f_nvals, level_b_nvals;
    GrB_Matrix_nvals(&level_f_nvals, level_f);
    GrB_Matrix_nvals(&level_b_nvals, level_b);
    //printf("level_f: %d, level_b: %d\n", level_f_nvals, level_b_nvals);

    code = GrB_kronecker(Mask, T_graph->A, NULL, GrB_LAND, level_f, transposed_level_b, NULL);
    if (code != GrB_SUCCESS) { error("mxm error"); return code; }

    LAGraph_Delete(&T_graph, NULL);

    GrB_Matrix result;
    GrB_Matrix_new(&result, GrB_INT64, intersection_states, intersection_states);

    GrB_kronecker(result, Mask, NULL, GrB_BAND_INT64, first_dfa, second_dfa, NULL);

    GrB_Index mask_nvals;

    GrB_Matrix_nvals(&mask_nvals, Mask);

    GrB_Index res_nvals;

    GrB_Matrix_nvals(&res_nvals, result);

    //printf("mask: %d, res: %d\n", mask_nvals, res_nvals);

    LAGraph_Finalize(NULL);
    GrB_finalize();
    return 0;
}
