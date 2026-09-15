#include <GraphBLAS.h>
#include <LAGraph.h>
#include <stdio.h>
#include <unistd.h>

#define error(x) fprintf(stderr, "%s\n", x)


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
            fprintf(stderr, "failed reading matrix %d of automaton %d\n", i + 1, automaton_number);
            return code;
        }
        fclose(f);
    }

    return code;
}


void g(void *z, const void *x, const void *y)
{
    *((int64_t *) z) = ((*(int64_t *)x) << 1) + *((int64_t *) y);
}


int main(int argc, char **argv)
{
    // usage: ./a.out alphabeth_cardinality [boolean decomposition of first] [boolean decomposition of second]
    if (argc < 2)
    {
        error("not enough input");
        return -1;
    }

    GrB_init(GrB_BLOCKING);
    LAGraph_Init(NULL);

    int decomposition_size = atoi(argv[1]);

    GrB_Info code;

    GrB_Matrix *first_decomp = malloc(sizeof(GrB_Matrix) * decomposition_size);
    GrB_Matrix *second_decomp = malloc(sizeof(GrB_Matrix) * decomposition_size);

    if (load_dfa(argv, decomposition_size, 1, first_decomp)) { return 1; }
    if (load_dfa(argv, decomposition_size, 2, second_decomp)) { return 2; }

    GrB_Index first_states, second_states;
    GrB_Matrix_nrows(&first_states, first_decomp[0]);
    GrB_Matrix_nrows(&second_states, second_decomp[0]);
    GrB_Index intersection_states = first_states * second_states;

    GrB_Matrix first_dfa, second_dfa;

    GrB_Matrix_new(&first_dfa, GrB_INT64, first_states, first_states);
    GrB_Matrix_new(&second_dfa, GrB_INT64, second_states, second_states);

    GrB_BinaryOp shift_and_plus = malloc(sizeof(GrB_BinaryOp));
    code = GrB_BinaryOp_new(&shift_and_plus, &g, GrB_INT64, GrB_INT64, GrB_INT64);

    for (int i = 0; i < decomposition_size; i++)
    {
        code = GrB_eWiseAdd(first_dfa, NULL, NULL, shift_and_plus, first_dfa, first_decomp[i], NULL);
        if (code != GrB_SUCCESS) { error("failed to merge another label into first"); return code; }
        code = GrB_eWiseAdd(second_dfa, NULL, NULL, shift_and_plus, second_dfa, second_decomp[i], NULL);
        if (code != GrB_SUCCESS) { error("failed to merge another label into second"); return code; }
    }

    for (int i = 0; i < decomposition_size; i++)
    {
        GrB_free(first_decomp + i);
        GrB_free(second_decomp + i);
    }

    free(first_decomp);
    free(second_decomp);

    GrB_Matrix result;
    GrB_Matrix_new(&result, GrB_INT64, intersection_states, intersection_states);

    GrB_kronecker(result, NULL, NULL, GrB_BAND_INT64, first_dfa, second_dfa, NULL);

    GrB_free(&first_dfa);
    GrB_free(&second_dfa);
    GrB_free(&result);
    GrB_free(&shift_and_plus);

    LAGraph_Finalize(NULL);
    GrB_finalize();
    return 0;
}
