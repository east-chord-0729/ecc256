#include "arith_lib.h"

void test_add()
{
    FILE* infile_a;
    FILE* infile_b;
    FILE* outfile;
    BN A, B, R;

    infile_a = fopen("testvectors_add_mod/TV_opA.txt", "r");
    infile_b = fopen("testvectors_add_mod/TV_opB.txt", "r");
    outfile  = fopen("testvectors_add_mod/TV_PFADD_TV_res.txt", "w");
    if(infile_a == NULL || infile_b == NULL) goto end;

    while (!feof(infile_a)) {
        // read A and B
        for(int i = NUMWORD - 1; i >= 0; i--) {
            fscanf(infile_a, "%08x", &A.v[i]);
            fscanf(infile_b, "%08x", &B.v[i]);
        }

        // R = A + B mod P
        addp(&R, &A, &B);

        // write R
        for(int i = NUMWORD - 1; i >= 0; i--) fprintf(outfile, "%08X", R.v[i]);
        fwrite("\n\n", sizeof(char), 2, outfile);
    }
end:
    fclose(infile_a);
    fclose(infile_b);
    fclose(outfile);
}
void test_sub()
{
    FILE* infile_a;
    FILE* infile_b;
    FILE* outfile;
    BN A, B, R;

    infile_a = fopen("testvectors_sub_mod/TV_opA.txt", "r");
    infile_b = fopen("testvectors_sub_mod/TV_opB.txt", "r");
    outfile  = fopen("testvectors_sub_mod/TV_PFSUB_TV_res.txt", "w");
    if(infile_a == NULL || infile_b == NULL) goto end;

    while (!feof(infile_a)) {
        // read A and B
        for(int i = NUMWORD - 1; i >= 0; i--) {
            fscanf(infile_a, "%08x", &A.v[i]);
            fscanf(infile_b, "%08x", &B.v[i]);
        }

        // R = A - B mod P
        subp(&R, &A, &B);

        // write R
        for(int i = NUMWORD - 1; i >= 0; i--) fprintf(outfile, "%08X", R.v[i]);
        fwrite("\n\n", sizeof(char), 2, outfile);
    }
end:
    fclose(infile_a);
    fclose(infile_b);
    fclose(outfile);
}
void test_mul_os()
{
    FILE* infile_a;
    FILE* infile_b;
    FILE* outfile;
    BN A, B;
    BN2 R;

    infile_a = fopen("testvectors_mul_os/TV_opA.txt", "r");
    infile_b = fopen("testvectors_mul_os/TV_opB.txt", "r");
    outfile  = fopen("testvectors_mul_os/TV_MUL_TV_res.txt", "w");
    if(infile_a == NULL || infile_b == NULL) goto end;

    while (!feof(infile_a)) {
        // read A and B
        for(int i = NUMWORD - 1; i >= 0; i--) {
            fscanf(infile_a, "%08x", &A.v[i]);
            fscanf(infile_b, "%08x", &B.v[i]);
        }

        // R = A * B
        umul_os(&R, &A, &B);

        // write R
        for(int i = NUMWORD2 - 1; i >= 0; i--) fprintf(outfile, "%08X", R.v[i]);
        fwrite("\n\n", sizeof(char), 2, outfile);
    }
end:
    fclose(infile_a);
    fclose(infile_b);
    fclose(outfile);
}
void test_mul_ps()
{
    FILE* infile_a;
    FILE* infile_b;
    FILE* outfile;
    BN A, B;
    BN2 R;

    infile_a = fopen("testvectors_mul_ps/TV_opA.txt", "r");
    infile_b = fopen("testvectors_mul_ps/TV_opB.txt", "r");
    outfile  = fopen("testvectors_mul_ps/TV_MUL_TV_res.txt", "w");
    if(infile_a == NULL || infile_b == NULL) goto end;

    while (!feof(infile_a)) {
        // read A and B
        for(int i = NUMWORD - 1; i >= 0; i--) {
            fscanf(infile_a, "%08x", &A.v[i]);
            fscanf(infile_b, "%08x", &B.v[i]);
        }

        // R = A * B
        umul_ps(&R, &A, &B);

        // write R
        for(int i = NUMWORD2 - 1; i >= 0; i--) fprintf(outfile, "%08X", R.v[i]);
        fwrite("\n\n", sizeof(char), 2, outfile);
    }
end:
    fclose(infile_a);
    fclose(infile_b);
    fclose(outfile);
}
void test_sqr()
{
    FILE* infile_a;
    FILE* infile_b;
    FILE* outfile;
    BN A, B;
    BN2 R;

    infile_a = fopen("testvectors_sqr/TV_opA.txt", "r");
    infile_b = fopen("testvectors_sqr/TV_opB.txt", "r");
    outfile  = fopen("testvectors_sqr/TV_SQR_TV_res.txt", "w");
    if(infile_a == NULL || infile_b == NULL) goto end;

    while (!feof(infile_a)) {
        // read A and B
        for(int i = NUMWORD - 1; i >= 0; i--) {
            fscanf(infile_a, "%08x", &A.v[i]);
            fscanf(infile_b, "%08x", &B.v[i]);
        }

        // R = A * B
        usqr_ps(&R, &A);

        // write R
        for(int i = NUMWORD2 - 1; i >= 0; i--) fprintf(outfile, "%08X", R.v[i]);
        fwrite("\n\n", sizeof(char), 2, outfile);
    }
end:
    fclose(infile_a);
    fclose(infile_b);
    fclose(outfile);
}
void test_mod()
{
    FILE* infile_a;
    FILE* infile_b;
    FILE* outfile;
    BN A, B, R;
    BN2 AB;

    infile_a = fopen("testvectors_mod/TV_opA.txt", "r");
    infile_b = fopen("testvectors_mod/TV_opB.txt", "r");
    outfile  = fopen("testvectors_mod/TV_PFMUL_TV_res.txt", "w");
    if(infile_a == NULL || infile_b == NULL) goto end;

    while (!feof(infile_a)) {
        // read A and B
        for(int i = NUMWORD - 1; i >= 0; i--) {
            fscanf(infile_a, "%08x", &A.v[i]);
            fscanf(infile_b, "%08x", &B.v[i]);
        }

        // R = A * B
        umul_os(&AB, &A, &B);
        mod_fast(&R, &AB);

        // write R
        for(int i = NUMWORD - 1; i >= 0; i--) fprintf(outfile, "%08X", R.v[i]);
        fwrite("\n\n", sizeof(char), 2, outfile);
    }
end:
    fclose(infile_a);
    fclose(infile_b);
    fclose(outfile);
}
void test_inv_bin()
{
    FILE* infile_a;
    FILE* outfile;
    BN A, R;

    infile_a = fopen("testvectors_inv_bin/TV_opA.txt", "r");
    outfile  = fopen("testvectors_inv_bin/TV_PFINV_TV_res.txt", "w");
    if(infile_a == NULL) goto end;

    while (!feof(infile_a)) {
        // read A and B
        for(int i = NUMWORD - 1; i >= 0; i--) {
            fscanf(infile_a, "%08x", &A.v[i]);
        }

        // R = A^{-1} mod P
        inv_bin(&R, &A);

        // write R
        for(int i = NUMWORD - 1; i >= 0; i--) fprintf(outfile, "%08X", R.v[i]);
        fwrite("\n\n", sizeof(char), 2, outfile);
    }
end:
    fclose(infile_a);
    fclose(outfile);
}

int main(void) {
    test_add();
    //test_sub();
    //test_mul_os();
    //test_mul_ps();
    //test_sqr();
    //test_mod();
    //test_inv_bin();

    return 0;
}