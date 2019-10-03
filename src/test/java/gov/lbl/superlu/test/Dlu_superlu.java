package gov.lbl.superlu.test;


import gov.lbl.superlu.Dlu_supermatrix;
import gov.lbl.superlu.Dlu_supermatrix.NCPformat;
import gov.lbl.superlu.Dlu_supermatrix.NCformat;
import gov.lbl.superlu.Dlu_supermatrix.SCPformat;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import gov.lbl.superlu.Dlu_slu_mt_util.superlu_memusage_t;
import gov.lbl.superlu.Dlu_slu_mt_util.trans_t;

import static gov.lbl.superlu.Dlu_sp_ienv.sp_ienv;

import static gov.lbl.superlu.Dlu_slu_mt_util.trans_t.NOTRANS;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_ABORT;

import static gov.lbl.superlu.Dlu_supermatrix.Stype_t.SLU_NC;
import static gov.lbl.superlu.Dlu_supermatrix.Stype_t.SLU_DN;
import static gov.lbl.superlu.Dlu_supermatrix.Dtype_t.SLU_D;
import static gov.lbl.superlu.Dlu_supermatrix.Mtype_t.SLU_GE;

import static gov.lbl.superlu.Dlu_pdutil.dCreate_CompCol_Matrix;
import static gov.lbl.superlu.Dlu_pdutil.dCreate_Dense_Matrix;
import static gov.lbl.superlu.Dlu_pdutil.dGenXtrue;
import static gov.lbl.superlu.Dlu_pdutil.dFillRHS;
import static gov.lbl.superlu.Dlu_pdutil.dinf_norm_error;

import static gov.lbl.superlu.Dlu_pdmemory.doubleMalloc;
import static gov.lbl.superlu.Dlu_pdmemory.superlu_dQuerySpace;

import static gov.lbl.superlu.Dlu_pmemory.intMalloc;

import static gov.lbl.superlu.Dlu.printf;

import static gov.lbl.superlu.Dlu_get_perm_c.get_perm_c;

import static gov.lbl.superlu.Dlu_pdgssv.pdgssv;

import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

public class Dlu_superlu {

    /**
     * Purpose
     * =======
     *
     * This is the small 5x5 example used in the Sections 2 and 3 of the
     * Users' Guide to illustrate how to call a SuperLU routine, and the
     * matrix data structures used by SuperLU.
     *
     */
    public static final int nprocs = Runtime.getRuntime().availableProcessors();
    public static void main(String[] args) {

        SuperMatrix A = new SuperMatrix();
        SuperMatrix B = new SuperMatrix();
        SuperMatrix L = new SuperMatrix();
        SuperMatrix U = new SuperMatrix();
        double   a[], rhs[];
        double   s, u, p, e, r, l;
        int      asub[], xa[];
        int      perm_r[]; /* row permutations from partial pivoting */
        int      perm_c[]; /* column permutation vector */
        int      nrhs, i, m, n, nnz;
        int[] info = new int[1];
        int      permc_spec;
//        Dlu_superlu_options_t options = new Dlu_superlu_options_t();
//        SuperLUStat_t stat = new SuperLUStat_t();

        /* Initialize matrix A. */
        m = n = 5;
        nnz = 12;
        a = doubleMalloc(nnz);
        asub = intMalloc(nnz);
        xa = intMalloc(n+1);
        s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
        a[0] = s; a[1] = l; a[2] = l; a[3] = u; a[4] = l; a[5] = l;
        a[6] = u; a[7] = p; a[8] = u; a[9] = e; a[10]= u; a[11]= r;
        asub[0] = 0; asub[1] = 1; asub[2] = 4; asub[3] = 1;
        asub[4] = 2; asub[5] = 4; asub[6] = 0; asub[7] = 2;
        asub[8] = 0; asub[9] = 3; asub[10]= 3; asub[11]= 4;
        xa[0] = 0; xa[1] = 3; xa[2] = 6; xa[3] = 8; xa[4] = 10; xa[5] = 12;

        /* Create matrix A in the format expected by SuperLU. */
        dCreate_CompCol_Matrix(A, m, n, nnz, a, asub, xa, Dlu_supermatrix.Stype_t.SLU_NC, Dlu_supermatrix.Dtype_t.SLU_D, Dlu_supermatrix.Mtype_t.SLU_GE);

        /* Create right-hand side matrix B. */
        nrhs = 1;
        rhs = doubleMalloc(m * nrhs);
        for (i = 0; i < m; ++i) rhs[i] = 1.0;
        dCreate_Dense_Matrix(B,m, nrhs, rhs, m, Dlu_supermatrix.Stype_t.SLU_DN, Dlu_supermatrix.Dtype_t.SLU_D, Dlu_supermatrix.Mtype_t.SLU_GE);

        perm_r = intMalloc(m);
        perm_c = intMalloc(n);

        permc_spec = 1;
        get_perm_c(permc_spec, A, perm_c);
//        /* Set the default input options. */
//        set_default_options(options);
//        options.ColPerm = colperm_t.NATURAL;
//
//        /* Initialize the statistics variables. */
//        StatInit(stat);

        /* Solve the linear system. */
        pdgssv(nprocs, A, perm_c, perm_r, L, U, B, info);
    }

}
