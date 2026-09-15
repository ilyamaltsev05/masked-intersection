#include <GraphBLAS.h>
#include <LAGraph.h>
#include <LAGraphX.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>


#define error(x) fprintf(stderr, "%s\n", x)


int64_t second_states_size = 0;


int load_dfa(char **argv, int labels, int automaton_number, GrB_Matrix *A)
{
    FILE *f;
    GrB_Info code = GrB_SUCCESS;
    int offset = labels * (automaton_number - 1) + 2;

    for (int i = 0; i < labels; i++)
    {
        f = fopen(argv[offset + i], "r");
        code = LAGraph_MMRead(A + i, f, NULL);
        if (code != GrB_SUCCESS)
        {
            fprintf(stderr, "failed reading matrix %d of automaton %d, code:%d\n", i + 1, automaton_number, code);
            return code;
        }
        fclose(f);
    }

    return code;
}

// states are represented as a full number_of_sources * 1 vector with values being ids of source vertices
int load_states(char **argv, int offset, GrB_Vector *V)
{
    FILE *f;
    GrB_Info code = GrB_SUCCESS;

    f = fopen(argv[offset], "r");
    code = LAGraph_MMRead(V, f, NULL);
    fclose(f);
    if (code != GrB_SUCCESS)
    {
        error("failed reading states of first automaton");
        return code;
    }

    f = fopen(argv[offset + 1], "r");
    code = LAGraph_MMRead(V + 1, f, NULL);
    fclose(f);
    if (code != GrB_SUCCESS)
    {
        error("failed reading states of second automaton");
        return code;
    }

    return code;
}


void f(void *z, const void *x, const void *y)
{
    int64_t a = *((int64_t *) x);
    int64_t b = *((int64_t *) y);
    *((int64_t *) z) = (a - 1) * second_states_size + (b - 1);
}


void g(void *z, const void *x, const void *y)
{
    *((int64_t *) z) = ((*(int64_t *)x) << 1) + *((int64_t *) y);
}


int main(int argc, char **argv)
{
    // usage: ./a.out alphabeth_cardinality [boolean decomposition of first] [boolean decomposition of second]
    // [vector of start states] [vector of start states]
    if (argc < 6)
    {
        error("not enough input");
        return -1;
    }

    GrB_init(GrB_BLOCKING);
    LAGraph_Init(NULL);

    int decomposition_size = atoi(argv[1]);

    GrB_Info code;

    // load dfas

    GrB_Matrix *first_decomp = malloc(sizeof(GrB_Matrix) * decomposition_size);
    GrB_Matrix *second_decomp = malloc(sizeof(GrB_Matrix) * decomposition_size);

    if (load_dfa(argv, decomposition_size, 1, first_decomp)) { return 1; }
    if (load_dfa(argv, decomposition_size, 2, second_decomp)) { return 2; }

    GrB_Vector *states;
    states = malloc(sizeof(GrB_Vector) * 2);
    load_states(argv, 2 + 2 * decomposition_size, states);

    GrB_Index nrows_states_first, nrows_states_second;
    GrB_Matrix_nrows(&nrows_states_first, states[0]);
    GrB_Matrix_nrows(&nrows_states_second, states[1]);

    GrB_Matrix start_states;
    GrB_Index nrows_start_states = nrows_states_first * nrows_states_second;
    GrB_Matrix_nrows(&second_states_size, second_decomp[0]);
    //second_states_size = nrows_states_second;
    GrB_BinaryOp states_mul = malloc(sizeof(GrB_BinaryOp));
    if (states_mul == NULL) { error("null"); return -1; }
    code = GrB_BinaryOp_new(&states_mul, &f, GrB_INT64, GrB_INT64, GrB_INT64);
    if (code != GrB_SUCCESS) { error("failed to created binary op"); return 1; }
    GrB_Matrix_new(&start_states, GrB_INT64, nrows_start_states, 1);

    code = GrB_kronecker(start_states, NULL, NULL, states_mul, states[0], states[1], NULL);
    if (code != GrB_SUCCESS) { error("failed to get start states of intersection"); return code; }

    GrB_Index first_size, second_size, intersection_states;
    GrB_Matrix_nrows(&first_size, first_decomp[0]);
    GrB_Matrix_nrows(&second_size, second_decomp[0]);
    intersection_states = first_size * second_size;

    GrB_Matrix T;
    GrB_Matrix_new(&T, GrB_BOOL, intersection_states, intersection_states);

    GrB_Matrix first_dfa, second_dfa;

    GrB_Matrix_new(&first_dfa, GrB_INT64, first_size, first_size);
    GrB_Matrix_new(&second_dfa, GrB_INT64, second_size, second_size);

    GrB_BinaryOp shift_and_plus = malloc(sizeof(GrB_BinaryOp));
    code = GrB_BinaryOp_new(&shift_and_plus, &g, GrB_INT64, GrB_INT64, GrB_INT64);

    GrB_Matrix intermediate;

    for (int i = 0; i < decomposition_size; i++)
    {
        GrB_Matrix_new(&intermediate, GrB_BOOL, intersection_states, intersection_states);
        code = GrB_kronecker(intermediate, NULL, NULL, GrB_LAND, first_decomp[i], second_decomp[i], NULL);
        if (code != GrB_SUCCESS) { fprintf(stderr, "kronecker on %d label error", i + 1); return 1; }
        code = GrB_eWiseAdd(T, NULL, NULL, GrB_LOR, T, intermediate, NULL);
        if (code != GrB_SUCCESS) { fprintf(stderr, "ewise add on %d label error", i + 1); return 1; }
        GrB_free(&intermediate);
    }

    GrB_Scalar two;
    GrB_Scalar_new(&two, GrB_INT64);
    GrB_Scalar_setElement_INT64(two, 2);

    for (int i = 0; i < decomposition_size; i++)
    {
        GrB_Matrix_apply_BinaryOp2nd_Scalar(first_dfa, NULL, NULL, GrB_TIMES_INT64, first_dfa, two, NULL);
        code = GrB_eWiseAdd(first_dfa, NULL, NULL, GrB_PLUS_INT64, first_dfa, first_decomp[i], NULL);
        if (code != GrB_SUCCESS) { error("failed to merge another label into first"); return code; }
        GrB_free(first_decomp + i);
        GrB_Matrix_apply_BinaryOp2nd_Scalar(second_dfa, NULL, NULL, GrB_TIMES_INT64, second_dfa, two, NULL);
        code = GrB_eWiseAdd(second_dfa, NULL, NULL, GrB_PLUS_INT64, second_dfa, second_decomp[i], NULL);
        if (code != GrB_SUCCESS) { error("failed to merge another label into second"); return code; }
        GrB_free(second_decomp + i);
    }
    GrB_free(&two);

    free(first_decomp);
    free(second_decomp);

    GrB_Index T_nvals;
    GrB_Matrix_nvals(&T_nvals, T);

    // MS BFS over matrix T to get reachable from starting states of automaton

    LAGraph_Graph T_graph;

    code = LAGraph_New(&T_graph, &T, LAGraph_ADJACENCY_DIRECTED, NULL);
    if (code != GrB_SUCCESS) { error("failed to create graph"); return code; }

    GrB_Matrix level = NULL;
    GrB_Vector reachable;

    GrB_Matrix_new(&level, GrB_INT64, nrows_start_states, intersection_states);
    GrB_Vector_new(&reachable, GrB_BOOL, intersection_states);

    code = LAGraph_MultiSourceBFS(&level, NULL, T_graph, start_states, NULL);
    if (code != GrB_SUCCESS) { error("failed forward BFS"); return code; }
    GrB_Matrix_apply(level, NULL, NULL, GxB_ONE_INT64, level, NULL);
    code = GrB_Matrix_reduce_Monoid(reachable, NULL, NULL, GrB_LOR_MONOID_BOOL, level, GrB_DESC_T0);
    if (code != GrB_SUCCESS) { error("failed to reduce"); return code; }
    GrB_free(&level);

    GrB_Matrix Mask;
    GrB_Matrix_new(&Mask, GrB_BOOL, intersection_states, intersection_states);

    GrB_Matrix diag, temp;
    GrB_Matrix_new(&diag, GrB_BOOL, intersection_states, intersection_states);
    GrB_Matrix_new(&temp, GrB_BOOL, intersection_states, intersection_states);
    GxB_Matrix_diag(diag, reachable, 0, NULL);
    GrB_free(&reachable);

    GrB_mxm(temp, NULL, NULL, GrB_LAND_LOR_SEMIRING_BOOL, diag, T_graph->A, NULL);
    GrB_mxm(Mask, NULL, NULL, GrB_LAND_LOR_SEMIRING_BOOL, temp, diag, NULL);
    GrB_free(&diag); GrB_free(&temp);

    LAGraph_Delete(&T_graph, NULL);

    GrB_Matrix result;
    GrB_Matrix_new(&result, GrB_INT64, intersection_states, intersection_states);

    GrB_kronecker(result, Mask, NULL, GrB_BAND_INT64, first_dfa, second_dfa, GrB_DESC_S);

    FILE *res = fopen(argv[4 + 2 * decomposition_size], "w");
    if (res) { LAGraph_MMWrite(result, res, NULL, NULL); }

    LAGraph_Finalize(NULL);
    GrB_finalize();
    return 0;
}
