//
//  libtest.m
//  NumKit
//
//  Created by Koichi Oshio on 2-21-2013.
//

#import <Foundation/Foundation.h>
#import "NumKit.h"

#import "../../RecKit/RecKit/Src/timer_macros.h"

void    test0();    // machine eps
void    test1();    // unit test
void    test2();    // Num_least_sq
void    test3();    // stat
void    test4();    // t-test, F-test
void    test5();    // complex t-test
void    test6();    // matrix basic
void    test7();    // gaussj
void    test8();    // pinv
void    test9();    // jacobi
void    test10();   // tri-exponential
void    test11();   // Powell
void    test12();   // Runge-Kutta
void    test13();   // radial sample points
void    test14();   // covariance
void    test15();   // Mat <-> RecImage
void    test16();   // phase graph
void    test17();   // eigenvalue
void    test18();   // SVD
void    test19();   // ICA
void    test20();   // ICA
void	test21();	// NumMatrix (solveLinear, solveLinerarNN)
void	test22();	// NumMatrix (low level func)
void	test23();	// NNLS analysis of multi-exponential data

int main(int argc, const char * argv[])
{
    @autoreleasepool {
    //	test0();    //Machine EPS = 1.4013e-45 (iter = 149)
    //	test1();    // brak / brent
    //	test2();	// LS, marq & conj
    //	test3();    // stat
    //	test4();    // t-test / F-test
    //	test5();    // complex t-test
    //	test6();    // matrix basic
    //	test7();	// gaussj
    //	test8();    // pinv
    //	test9();    // jacobi (test case for lapack interface)
    //	test10();   // tri-exponential
    //	test11();   // Powell
	//	test12();	// rk4
	//	test13();   // radial
	//	test14();	// covariance
	//	test15();	// Mat <-> RecImage conversion
	//	test16();	// phase graph
	//	test17();	// eigenvalue
		test18();	// SVD / ICA
	//	test19();   // ICA
	//	test20();   // ICA
	//	test21();	// NumMatrix
	//	test22();	// NumMatrix (low level methods)
	//	test23();	// NNLS analysis of multi-exponential data
   }
    return 0;
}

float testfunc(float a)
{
    float   mn = 1.23456789;
    float   f;

    f = (a - mn);
    return pow(f, 4);
}

void
test0()
{
    float   a, b;
    int     i;

    a = 1.0;
    for (i = 0; i < 200; i++) {
        b = a;
        a /= 2;
        if (b == a || a == 0) {
            printf("Machine EPS = %g (iter = %d)\n", b, i);
            break;
        }
    }
}

void
test1()
{
    float   a, b, c;
    float   fmin, xmin;
	float	pmn, pmx;
    float   tol = 1.0e-7;   // tolerance for brent. 1.0e-8 doesn't work
    float   (^blk)(float);

    a = 0.0;
    b = 2.0;
    c = 10.0;
	pmn = 0;
	pmx = 10;

    blk = ^float(float tt) {
        float   mn = 1.23456789;
        float   f;
        f = (tt - mn);
        return pow(f, 4);
    };

    printf("=== mnbrak ===\n");
    Num_mnbrak(&a, &b, &c, blk);
    printf("a = %f, b = %f, c = %f\n", a, b, c);

    printf("=== brent ===\n");
    fmin = Num_brent(a, b, c, blk, tol, &xmin);
    printf("min = %g at x = %g\n", fmin, xmin);
}

float   tm[10]   = {  30,   50,   80,  100,  200, 400, 1000};
float   data[10] = {2078, 1917, 1843, 1790, 1472, 959,  450};
int     ndata = 7;
int     nparam = 2;

void
test2()
{
	int         i, iter = 0;
	int			mode = 0; // 0:conj gr, 1:marquardt
    float       minval = 0;
    Num_data    *d;
    Num_param   *p;
    Num_range   *r;
    float       xmin, xmax, ymin, ymax;
    float   (^cost)(float *);
    float   (^model)(float x, float *p, float *dydp);
    cost = ^float(float *ptry) {
        int     i;
        float   val, cst = 0;
        for (i = 0; i < ndata; i++) {
            val = model(tm[i], ptry, NULL);
            val -= data[i];
            cst += val*val / ndata;    // MSE
        }
        return cst;    
    };
    model = ^float(float x, float *p, float *dydp) {
		float   e1;
		float   y;
	// eval func val at x
		e1 = exp(-x * p[1]);
		y = p[0] * e1;
	// gradient at x
		if (dydp != NULL) {
			dydp[0] = e1;
			dydp[1] = -p[0] * x * e1;
		}
		return y;
	};

// alloc param struct
    d = Num_alloc_data(ndata);
    p = Num_alloc_param(nparam);
    r = Num_alloc_range(nparam);
// set data
    for (i = 0; i < ndata; i++) {
        d->x[i] = tm[i];
        d->y[i] = data[i];
    }

// normalize
    if (1) {
        Num_normalize_data(d);
    }

// initial estimate (scaled space)
if (1) {
    ymin = d->y[ndata - 1];
    ymax = d->y[0];
    xmin = d->x[0];
    xmax = d->x[ndata - 1];
    p->data[0] = ymax;
	p->min[0] = 0.0;
	p->max[0] = 2.0;
    p->data[1] =  log(ymax / ymin) / (xmax - xmin);
	p->min[1] = 0.1;
	p->max[1] = 10.0;

    printf("initial guess: p0:%8.2g, p1:%8.2g\n", p->data[0], p->data[1]);
} else {
// global search
    printf("global search\n");
    r->min[0] = d->y[0] * 0.8;
    r->inc[0] = 0.1;
    r->steps[0] = 5;
    r->min[1] = d->x[ndata - 1] * 0.2;
    r->inc[1] = 0.2;
    r->steps[1] = 10;

    Num_search_min(d, p, r, model);
    printf("initial guess: p0:%8.2g, p1:%8.2g\n", p->data[0], p->data[1]);
}

// minimization (Numerical Recipe, p321)
    if (mode == 0) {
        iter = Num_least_sq(d, p, model, &minval);   // conj grad
    } else {
     //   iter = Num_marquardt(d, p, model, &minval);  // m & l
    }

    printf("result:    a0 = %f, a1 = %f, min = %f (%%) (%d iterations)\n",
        p->data[0] * d->yscale, p->data[1] / d->xscale * 1000, sqrt(minval) * 100, iter);
    printf("===========================================================\n");
    printf("ref(conj): a0 = 2132.442871, a1 = 1.803731, min = 2.599068 (%%) (3 iterations)\n");
//    printf("ref(ml),  : a0 = 2132.422119, a1 = 1.803796, min = 6.876489 (%%) (3 iterations)\n");
    Num_free_data(d);
    Num_free_param(p);
    Num_free_range(r);
}

//float betacf(float a, float b, float x);
//float betai(float a, float b, float x);
//float tdist(int n, float x);

void
test3()
{
    float	*p1, *p2;
    float   av, vr;
    float   t, f, x, pval;
    int     n1 = 5;
    int     n2 = 5;
    int     i;

//
	if (1) {
		n1 = 100;
		p1 = VECTOR(n1);
		srandom(1);
		for (i = 0; i < n1; i++) {
			p1[i] = Num_nrml(0.0, 1.0);
			printf("%d %14f\n", i, p1[i]);
		}
		av = 0;
		for (i = 0; i < n1; i++) {
			av += p1[i];
		}
		av /= n1;
		printf("av = %f\n", av);

		vr = 0;
		for (i = 0; i < n1; i++) {
			vr += (p1[i] - av) * (p1[i] - av);
		}
		vr /= n1;
		printf("var1 = %f\n", vr);

		vr = 0;
		for (i = 0; i < n1; i++) {
			vr += p1[i] * p1[i];
		}
		vr = vr/n1 - av*av;
		printf("var2 = %f\n", vr);
		exit(0);
	}


// test avevar
	n1 = n2 = 5;
    p1 = VECTOR(n1);
    p2 = VECTOR(n2);
    p1[0] = 6;
    p1[1] = 7;
    p1[2] = 9;
    p1[3] = 15;
    p1[4] = 21;
    
    p2[0] = 20;
    p2[1] = 28;
    p2[2] = 31;
    p2[3] = 38;
    p2[4] = 40;

    Num_avevar(p1, n1, &av, &vr);
    printf("Av/Var = %f / %f\n", av, vr);
    Num_avevar(p2, n2, &av, &vr);
    printf("Av/Var = %f / %f\n", av, vr);

    // regularized incomplete beta function
    // (= cumulative beta distribution)
    if (0) {
        for (i = 0; i <= 200; i++) {
            x = (float)i / 200.0;
            av = betai(2, 5, x);
            printf("%f %f\n", x, av);
        }
    }

    // t-distribution
    if (0) {
        for (i = 1; i < 80; i++) {
            x = (float)(i - 40) / 10.0;
            av = Num_tdist(x, 2);
            printf("%f %f\n", x, av);
        }
    }

    if (0) {
        t = 0;
        // cumulative t-distribution (two-sided)
        for (i = 1; i < 80; i++) {
            float df = n1 * 2 - 2;
            t = (float)(i - 40) / 10.0;
            av = betai(0.5 * df, 0.5, df / (df + t * t));
            printf("%f %f\n", t, av);
        }
        // t-test
        Num_ttest(p1, n1, n2, p2, &t, &pval);
        printf("tval = %f, pval = %f\n", t, pval);
    }

    if (0) {
    // cumulative F-distribution (two-sided)
        for (i = 1; i < 80; i++) {
            float df1 = n1 - 1;
            float df2 = n2 - 1;
        df1 = 20;
        df2 = 20;
            f = (float)i / 10.0;
            pval = 2.0 * betai(0.5 * df2, 0.5 * df1, df2 / (df2 + df1 * f));
            if (pval > 1.0) pval = 2.0 - pval;

            printf("%f %f\n", f, pval);
        }
    // F-test
        Num_ftest(p1, n1, p2, n2, &f, &pval);
        printf("F-val = %f, p-val = %f\n", f, pval);
    }

    free(p1);
    free(p2);
}

void
test4()
{
    float   f, t, p;
    float	*p1, *p2;
    float   m1, m2, v1, v2;
    float   mdif, vdif;
    int     i, k, df, df2;
    int     n = 100;
    int     ntry = 10000;
    float   tval, fval, pval;
    float   pthres = 0.05;
    int     fls, fls_neg;
    BOOL    ttest = NO;

// === test beta
//for (i = 0; i < 100; i++) {
//	t = (float)i * 10 / 100;
//	printf("%f %f\n", t, beta(100, t));
//}
//exit(0);

    printf("test4\n");
    p1 = VECTOR(n);
    p2 = VECTOR(n);

    if (ttest) {
        // check t-table
        t = 2.0;
        df = 10;
        p = Num_t2p(t, df);
        printf("t table\n");
        printf("t(%3.1f, %d), one sided p = %10.6f (should be 0.036694)\n", t, df, p);
    } else {
        f = 2.414;
        df = 7;
        df2 = 10;
        p = Num_f2p(f, df, df2);
        printf("F table\n");
        printf("F(%3.1f, %d, %d), one sided p = %10.6f (should be 0.1)\n", f, df, df2, p);
        f = 1.09295;
        df = 156;
        df2 = 182;
        p = Num_f2p(f, df, df2);
        printf("F(%3.1f, %d, %d), one sided p = %10.6f (should be 0.2811)\n", f, df, df2, p);
    }

        // test t2p
    fls = fls_neg = 0;
    for (k = 0; k < ntry; k++) {
        for (i = 0; i < n; i++) {
            p1[i] = Num_nrml(0.0, 1.0);
            p2[i] = Num_nrml(0.0, 1.0);
        }
        Num_avevar(p1, n, &m1, &v1);
        Num_avevar(p2, n, &m2, &v2);

        mdif = m1 - m2;
        vdif = v1 - v2;

        if (ttest) {
            Num_ttest(p1, n, p2, n, &tval, &pval);        // t-test (one-sided)
            if (pval < pthres && mdif > 0) {
                fls++;   // one-sided
            }
        } else {
            Num_ftest(p1, n, p2, n, &fval, &pval);        // F-test (one-sided)
            if (pval < pthres) {
                if (vdif > 0) {
                    fls++;
                } else {
                    fls_neg++;
                }
            }
        }
    }

    if (ttest) {
        printf("==== t - test ====\n");
        printf("fls = %f (should be %3.3f)\n", (float)fls / ntry, pthres);
    } else {
        printf("==== F - test ====\n");
        printf("fls/flsn = %f/%f (should be %3.3f)\n", (float)fls/ntry, (float)fls_neg/ntry, pthres);
    }

    free(p1);
    free(p2);
}

void
test6()
{
	RecImage	*im;
	float		*p;

    Num_mat     *a, *at, *c;
    Num_vec     *b, *x;
    float       ip;
	Num_mat		*cm;
	Num_vec		*cv;

    printf("test6\n");

    b = Num_new_vec(3);
    x = Num_new_vec(3);
    a = Num_new_mat(3, 4);
    c = Num_new_mat(3, 3);
//    at = Num_new_mat(3, 3);

    x->data[0] = 1;
    x->data[1] = 2;
    x->data[2] = 3;

    a->data[0] = 1; a->data[1] = 2; a->data[2] = 0;
    a->data[3] = 0; a->data[4] = 1; a->data[5] = 1;
    a->data[6] = 0; a->data[7] = 0; a->data[8] = 1;

printf("a\n");
dump_mat(a);
at = Num_trans(a);
printf("a tr\n");
dump_mat(at);
exit(0);


    Num_mvmul(b, a, x);

    dump_vec(b);
    dump_vec(x);
    dump_mat(a);

    ip = Num_dotpr(b, b);
    printf("ip = %f\n", ip);

    printf("\n gaussj \n");
    Num_gaussj(a, NULL);
    dump_mat(a);
    at = Num_trans(a);
    dump_mat(at);

    Num_mmul(c, a, at);
    dump_mat(c);

	cv = Num_new_cvec(3);
	cv->data[0] = 1.0;	cv->data[1] = 0.2;
	cv->data[2] = 2.0;	cv->data[3] = 0.3;
	cv->data[4] = 4.0;	cv->data[5] = 0.4;
	dump_vec(cv);

	cm = Num_new_cmat(2, 2);
	cm->data[0] = 1.0;	cm->data[1] = 0.2;	cm->data[2] = 2.0;	cm->data[3] = 0.2;
	cm->data[4] = 1.0;	cm->data[5] = 0.2;	cm->data[6] = 2.0;	cm->data[7] = 0.2;
	dump_mat(cm);

	im = [RecImage imageOfType:RECIMAGE_REAL xDim:3 yDim:3];
	p = [im data];
	p[0] = 1; p[1] = 2; p[2] = 3;
	p[3] = 4; p[4] = 5; p[5] = 6;
	p[6] = 7; p[7] = 8; p[8] = 9;

	[im saveAsKOImage:@"mimg.img"];

	a = Num_im_to_m(im);
	dump_mat(a);

	im = Num_m_to_im(a);
	[im dumpData];

	Num_free_mat(a);
}

void
test7()
{
    Num_mat     *a;
    Num_vec     *b, *x;

    b = Num_new_vec(3);
    x = Num_new_vec(3);
    a = Num_new_mat(3, 3);
    
    a->data[0] = 1; a->data[1] = 2; a->data[2] = 4;
    a->data[3] = 1; a->data[4] = 3; a->data[5] = 1;
    a->data[6] = 2; a->data[7] = 4; a->data[8] = 5;

    x->data[0] = 1;
    x->data[1] = 2;
    x->data[2] = 3;

    Num_mvmul(b, a, x);
    dump_vec(b);

    Num_gaussj(a, NULL);
    Num_mvmul(x, a, b);
    dump_vec(x);


}

// mat inv
void
test8()
{
    Num_mat     *A, *B, *X;
	RecImage	*a_img, *b_img;
	float		*pa, *pb;
    int         m, n;
	int			mode = 0;	// 0: square, 1: under, 2: over, 3: single b vector, 4: NNLS

	switch (mode) {
	case 0:
	// === square matrix ====
		m = 3, n = 3;
		a_img = [RecImage imageOfType:RECIMAGE_REAL xDim:m yDim:n];
		b_img = [RecImage imageOfType:RECIMAGE_REAL xDim:2 yDim:n];
	 
		pa = [a_img data];
		pb = [b_img data];

	// 3 x 3
		pa[0] = 1;	pa[1] = 2;	pa[2] = 3;
		pa[3] = 5;	pa[4] = 6;	pa[5] = 7;
		pa[6] = 1;	pa[7] = 1;	pa[8] = 3;

	// 3 x 2
		pb[0] = 1;	pb[1] = 2;
		pb[2] = 3;	pb[3] = 1;
		pb[4] = 1;	pb[5] = 1;
		
		A = Num_im_to_m(a_img);
		B = Num_im_to_m(b_img);
		X = Num_new_mat(A->nc, B->nc);

		dump_mat(A);
		dump_mat(B);
		Num_inv(X, A, B);
		dump_mat(X);

		Num_mmul(B, A, X);
		dump_mat(B);
		
		Num_free_mat(A);
		Num_free_mat(B);
		Num_free_mat(X);
		break;
	case 1:	// under-determined ### not correct yet (second colmumn)
		m = 3, n = 4;
		a_img = [RecImage imageOfType:RECIMAGE_REAL xDim:n yDim:m];
		b_img = [RecImage imageOfType:RECIMAGE_REAL xDim:2 yDim:m];
	 
		pa = [a_img data];
		pb = [b_img data];

	// A: 3 x 4
		pa[0] = 1;	pa[1] = 2;	pa[2] = 3;	pa[3] = 4;
		pa[4] = 5;	pa[5] = 6;	pa[6] = 7;	pa[7] = 4;
		pa[8] = 1;	pa[9] = 1;	pa[10] = 2;	pa[11] = 4;

	// B: 3 x 2
		pb[0] = 1;	pb[1] = 2;
		pb[2] = 3;	pb[3] = 1;
		pb[4] = 1;	pb[5] = 1;
		
		A = Num_im_to_m(a_img);
		B = Num_im_to_m(b_img);
		X = Num_new_mat(A->nc, B->nc);

		dump_mat(A);
		dump_mat(B);
		Num_inv(X, A, B);
		dump_mat(X);	// X : 7 x 2

		Num_mmul(B, A, X);
		dump_mat(B);
		break;
	case 2:
		m = 4, n = 3;
		a_img = [RecImage imageOfType:RECIMAGE_REAL xDim:n yDim:m];
		b_img = [RecImage imageOfType:RECIMAGE_REAL xDim:2 yDim:m];
	 
		pa = [a_img data];
		pb = [b_img data];

	// A: 4 x 3
		pa[0] = 1;	pa[1] = 2;	pa[2] = 3;;
		pa[3] = 5;	pa[4] = 6;	pa[5] = 7;
		pa[6] = 1;	pa[7] = 1;	pa[8] = 3;;
		pa[9] = 1;	pa[10] = 3;	pa[11] = 1;;

	// B: 3 x 2
		pb[0] = 1;	pb[1] = 2;
		pb[2] = 3;	pb[3] = 1;
		pb[4] = 1;	pb[5] = 3;
		pb[6] = 1;	pb[7] = 1;
		
		A = Num_im_to_m(a_img);
		B = Num_im_to_m(b_img);
		X = Num_new_mat(A->nc, B->nc);

		dump_mat(A);
		dump_mat(B);
		Num_inv(X, A, B);
		dump_mat(X);	// X : 7 x 2

		Num_mmul(B, A, X);
		dump_mat(B);
		break;
	case 3:
		m = 4, n = 3;
		a_img = [RecImage imageOfType:RECIMAGE_REAL xDim:n yDim:m];
		b_img = [RecImage imageOfType:RECIMAGE_REAL xDim:1 yDim:m];
	 
		pa = [a_img data];
		pb = [b_img data];

	// A: 4 x 3
		pa[0] = 1;	pa[1] = 2;	pa[2] = 3;;
		pa[3] = 5;	pa[4] = 6;	pa[5] = 7;
		pa[6] = 1;	pa[7] = 1;	pa[8] = 3;;
		pa[9] = 1;	pa[10] = 3;	pa[11] = 1;;

	// B: 3 x 1
		pb[0] = 1;
		pb[1] = 3;
		pb[2] = 1;
		pb[3] = 1;
		
		A = Num_im_to_m(a_img);
		B = Num_im_to_m(b_img);
		X = Num_new_mat(A->nc, B->nc);

		dump_mat(A);
		dump_mat(B);
		Num_inv(X, A, B);
		dump_mat(X);	// X : 7 x 2

		Num_mmul(B, A, X);
		dump_mat(B);
		break;
	case 4:
		m = 4, n = 3;
		a_img = [RecImage imageOfType:RECIMAGE_REAL xDim:n yDim:m];
		b_img = [RecImage imageOfType:RECIMAGE_REAL xDim:1 yDim:m];
	 
		pa = [a_img data];
		pb = [b_img data];

	// A: 4 x 3
		pa[0] = 1;	pa[1] = 1;	pa[2] = 1;
		pa[3] = 2;	pa[4] = 4;	pa[5] = 8;
		pa[6] = 3;	pa[7] = 9;	pa[8] = 27;
		pa[9] = 4;	pa[10] = 16;	pa[11] = 64;

	// B: 3 x 1
		pb[0] = 0.73;
		pb[1] = 3.24;
		pb[2] = 8.31;
		pb[3] = 16.72;
		
		A = Num_im_to_m(a_img);
		B = Num_im_to_m(b_img);
		X = Num_new_mat(A->nc, B->nc);

		dump_mat(A);
		dump_mat(B);
		Num_inv(X, A, B);
		dump_mat(X);	// X : 7 x 2

		Num_mmul(B, A, X);
		dump_mat(B);
		break;
	}
}

void
test9()
{
    Num_mat     *a;
    Num_vec     *eval;
	float		*work;
	int			lwork;
    int         n = 3;
	int			info;

    a = Num_new_mat(n, n);	// row-major
	eval = Num_new_vec(n);

// 0.6, 0.4, 1.0; 1.0, -1.0, 0.5;
    a->data[0] =  1; a->data[1] = -1; a->data[2] =  0;
    a->data[3] = -1; a->data[4] =  2; a->data[5] = -1;
    a->data[6] =  0; a->data[7] = -1; a->data[8] =  1;

    printf("initial A\n");
    dump_mat(a);

// test LAPACK

	// convert to col major
	Num_tocm(a);
    printf("converted to col-major\n");
	dump_mat(a);

	// get work size
	lwork = -1;
	work = (float *)malloc(sizeof(float) * 1);
	printf("=== ssysev get work size\n");
	ssyev_("V", "U", &n, a->data, &n, eval->data, work, &lwork, &info);
	lwork = work[0];
	printf("work size = %d\n", lwork);
	free(work);

	work = (float *)malloc(sizeof(float) * lwork);
	printf("=== ssysev \n");
	ssyev_("V", "U", &a->nc, a->data, &a->ld, eval->data, work, &lwork, &info);

	Num_torm(a);
    printf("A after ssyev\n");
    dump_mat(a);
    printf("evec and eval\n");
    dump_vec(eval);

    Num_free_vec(eval);
    Num_free_mat(a);
	free(work);
}

void
test10()
{
    float       b;
    float       d[3];
    Num_data    *dd;
    Num_mat     *A, *B;
    Num_vec     *x, *s, *e;
    int         i, j;
    int         np = 3;
    int         nd = 7;

    printf("Tri-exponential, linear fitting\n");

    d[0] = 0.2;
    d[1] = 1.5;
    d[2] = 4.0;

// alloc param struct
    dd = Num_alloc_data(nd);
// set data
    for (i = 0; i < nd; i++) {
        dd->x[i] = tm[i];
        dd->y[i] = data[i];
    }

// model
// s1 = f1 exp(-d1 b1) + f2 exp(-d2 b1) * f3 exp(-d3 b1)
// s2 = f1 exp(-d1 b2) + f2 exp(-d2 b2) * f3 exp(-d3 b2)
// s3 = f1 exp(-d1 b3) + f2 exp(-d2 b3) * f3 exp(-d3 b3)
// ...
// Si = f1 exp(-d1 bi) + f2 exp(-d2 bi) * f3 exp(-d3 bi);
// -> s = A f : A(exp(-d1 b1) 

    A = Num_new_mat(nd, np);
    B = Num_new_mat(np, nd);
    x = Num_new_vec(np);
    e = Num_new_vec(nd);
    s = Num_new_vec(nd);
    
//  A->data[i][j] = exp(-d[j] * x[i]), s[i] = y[i];
    for (i = 0; i < nd; i++) {
        s->data[i] = dd->y[i];
        for (j = 0; j < np; j++) {
            A->data[i * np + j] = exp(- d[j] * dd->x[i] * 1.0e-3);
        }
    }
    dump_mat(A);
    dump_vec(s);

// linear least sqares
//    Num_pinv(B, A);
//	Num_inv(X, A, B);
    //dump_mat(B);

    Num_mvmul(x, B, s);
    dump_vec(x);

// model
    for (i = 0; i < 100; i++) {
        b = i * 10.0;;
        printf("%f %f\n", b,
            x->data[0] * exp(- d[0] * b * 1.0e-3) +
            x->data[1] * exp(- d[1] * b * 1.0e-3) +
            x->data[2] * exp(- d[2] * b * 1.0e-3));
    }

// error
    Num_mvmul(e, A, x);
    dump_vec(s);
    dump_vec(e);
    b = 0;
    for (i = 0; i < nd; i++) {
        b += (e->data[i] - s->data[i]) * (e->data[i] - s->data[i]);
    }
    b /= s->n;
//    b = sqrt(b);
    printf("min = %f\n", b);
        
    Num_free_mat(A);
    Num_free_vec(x);
    Num_free_vec(s);
    Num_free_data(dd);
}

void
test11()   // Powell
{
    float       (^cost)(float *);
    float       minval;
    Num_data    *d;
    Num_param   *p;
    int         i, iter, nd = 7, np = 2;
    float       ymin, ymax, xmin, xmax;

// alloc param struct
    d  = Num_alloc_data(nd);
    p = Num_alloc_param(np);
// set data
    for (i = 0; i < nd; i++) {
        d->x[i] = tm[i];
        d->y[i] = data[i];
    }

    ymin = d->y[nd - 1];
    ymax = d->y[0];
    xmin = d->x[0];
    xmax = d->x[nd - 1];
    p->data[0] = ymax;
	p->min[0] = 0;
	p->max[0] = 2.0;	// normalized
    p->data[1] =  log(ymax / ymin) / (xmax - xmin);
	p->min[1] = 0;
	p->max[1] = 10;

// normalize
    if (1) {
        Num_normalize_data(d);
    }


    cost = ^float(float *p) {
        float   e1;
        float   y;
        float   cst;
        int     i;

        cst = 0;
        for (i = 0; i < nd; i++) {
            e1 = exp(- d->x[i] * p[1]);
            y = p[0] * e1;
            y -= d->y[i];
            y *= y;
            cst += y;
        }
    
        return cst / nd;
    };

    iter = Num_powell(p, cost, &minval);
    printf("result:    a0 = %f, a1 = %f, min = %f (%%) (%d iterations)\n",
        p->data[0] * d->yscale, p->data[1] / d->xscale * 1000, sqrt(minval) * 100, iter);
    printf("===========================================================\n");
    printf("ref(conj): a0 = 2132.442871, a1 = 1.803731, min = 2.599068 (%%) (3 iterations)\n");

    Num_free_data(d);
    Num_free_param(p);
    
}

void
test12()
{
    int     i;
    float   x, y;
    float   max = 4;
    int     n = 10;
    float   step = max / n;
	float	(^deriv)(float x, float y);

    deriv = ^float(float x, float y) {
        return 3 * x * x;
	};

    x = 0;
    y = 0;
    for (i = 0; i < n; i++) {
        Num_rk4(&x, &y, deriv, step);
        printf("%f %f %f\n", x, y, x*x*x);
    }
}

void
test13()
{
    int     i, j, n = 16;
    float   r, th, cs, sn;
/*
    for (i = 0; i < 10 * n; i++) {
        th = (i + 1) * M_PI / n;
        r = sin(th) / th;
        printf("%f %f\n", th / 5, r);
    }
*/
    for (i = 0; i < n; i++) {
        th = i * M_PI / n;
        cs = cos(th);
        sn = sin(th);
        for (j = 0; j < n; j++) {
            r = (j - n/2 + 0.5);
            printf("%f %f\n", r*cs, r*sn);
        }
    }
}

void
test14()
{
    int         i, n = 1000;
    Num_vec     *x, *y, *z;
    float       var;

    x = Num_new_vec(n);
    y = Num_new_vec(n);
    z = Num_new_vec(n);

    for (i = 0; i < n; i++) {
        x->data[i] = Num_nrml(0.0, 1.0);
        y->data[i] = Num_nrml(0.0, 1.0);
        Num_vadd(x, y, z);
    }

    var = Num_var(x->data, 1, n);
    printf("var = %f\n", var);

    var = Num_var(y->data, 1, n);
    printf("var = %f\n", var);

//    var = Num_covar(x->data, 1, z->data, 1, n);
    printf("covar(xy) = %f\n", var);
}

void
test15()
{
	Num_mat		*a;
	RecImage	*img;
	float		*p;
	int			i, j, ix;

	img = [RecImage imageOfType:RECIMAGE_REAL xDim:5 yDim:6];
	p = [img data];
	for (i = 0; i < 6; i++) {
		for (j = 0; j < 5; j++) {
			ix = i*5 + j;
			p[ix] = ix;
		}
	}
	[img saveAsKOImage:@"IMG_mat"];
	a = Num_im_to_m(img);
	dump_mat(a);

	img = Num_m_to_im(a);
	p = [img data];
	printf("====\n");
	for (i = 0; i < 6; i++) {
		for (j = 0; j < 5; j++) {
			ix = i*5 + j;
			printf("%5.2f ", p[ix]);;
		}
		printf("\n");
	}
}

void
test16()
{
	Num_rf		*rf;
	RecImage	*grph;
	int			i, n = 200;	// 200
	int			eint, necho;
	float		c;
	float		ang;
//	int			nrot;
	int			mode = 3;	// 0:stim, 1:steady-state, 2:RF spoil 3:CPMG
	BOOL		relax = NO;
	BOOL		cpmg;
	float		flop = 120;
	float		cpphs = 90;

	rf = Num_new_rf(n);
	if (relax) {
		rf->t1 = 500;	// msec
		rf->t2 = 100;	// msec
	} else {
		rf->t1 = 0; //1000;	// msec
		rf->t2 = 0; //100;	// msec
	}
	rf->tr = 5;	// msec
	rf->te = 2;	// msec

	switch (mode) {
	case 0 :	// stimuated echo
		rf->rho[0]		= 90;
		rf->rho[n / 5]	= 90;
		rf->rho[n / 2]	= 90;
		rf->theta[0]	= 0;
		rf->theta[n / 5] = 90;
		rf->theta[n / 2] = 0;
		cpmg = NO;
		break;
	case 1 :	// steady-state
		for (i = 0; i < n; i++) {
			rf->rho[i] = 90;
			rf->theta[i] = 180 * (i % 2);
		}
		cpmg = NO;
		break;
	case 2 :		// RF spoiling
		for (i = 0; i < n; i++) {
			c = 52;			// quad phase
			rf->rho[i] = 30;	// deg
			ang = c * i * (i + 1) / 2.0;	// quad
		//	ang = c * i * i;				// real quad
		//	ang = c * i;					// linear
			rf->theta[i] = ang;	// deg
		}
		cpmg = NO;
		break;
	case 3 :	// CPMG
		eint = 2;
		necho = n / eint;
		rf->rho[0] = 90;
		rf->theta[0] = 0;
		for (i = 0; i < necho; i++) {
			rf->rho[eint * (i*2 + 1)] = flop;
			rf->theta[eint * (i*2 + 1)] = cpphs;
		}
		cpmg = YES;
		break;
	}
	grph = Num_phase_graph(rf, cpmg);
	grph = Num_phase_graph_s(rf, cpmg);

	Num_free_rf(rf);
}

// eigen value
void
test17()
{
    Num_mat     *A, *Evec;
	int			n = 3;
    Num_vec     *eval;
	Num_vec		*x, *b;

    A = Num_new_mat(n, n);	// row-major
    Evec = Num_new_mat(n, n);	// row-major
	eval = Num_new_vec(n);

// 0.6, 0.4, 1.0; 1.0, -1.0, 0.5;
    A->data[0] =  1; A->data[3] = -1; A->data[6] =  0;
    A->data[1] = -1; A->data[4] =  2; A->data[7] = -1;
    A->data[2] =  0; A->data[5] = -1; A->data[8] =  1;

    printf("initial A\n");
    dump_mat(A);

	Num_evd_ref(A, Evec, eval);

    printf("A after ssyev\n");
    dump_mat(Evec);
    printf("evec and eval\n");
    dump_vec(eval);

// test
	x = Num_new_vec(n);
	b = Num_new_vec(n);
	Num_col_vec(x, Evec, 1);
	Num_mvmul(b, A, x);
	dump_vec(b);

// free mem
    Num_free_mat(A);
    Num_free_mat(Evec);
    Num_free_vec(eval);
	Num_free_vec(x);
	Num_free_vec(b);
}


// SVD / ICA
// add complex version
void
test18()
{
	RecImage		*img, *img_a, *img_u, *img_e, *tmp;
	NumMatrix		*A;
	int				xDim, yDim, nSlc;
	NSDictionary	*res;

	system("rm IMG_*");

// === low level
	if (0) {
		NumMatrix	*B;
		float		*p;
		int			i;
		A = [NumMatrix matrixOfType:NUM_REAL nRow:3 nCol:3];
		p = [A data];
		for (i = 0; i < [A len]; i++) {
			p[i] = i;
		}
		[A dump];
		B = [A colCenter];
		[B dump];
		B = [A rowCenter];
		[B dump];

		exit(0);
	}
	// non-linearity func
	if (0) {
		int		i;
		float	x, y1, y2, y3, y4, alpha = 2.0;
		for (i = 0; i < 100; i++) {
			x = (float)(i - 50) / 25;
			y1 = tanh(alpha * x);
			y2 = 1 / cosh(alpha * x);
			y3 = x * exp(-x*x/2);
			y4 = (1 - x* x) * exp(-x*x/2);
			printf("%f %f %f\n", x, y1, y2);
		}
		exit(0);
	}


// ==== new test image set

//	img = [RecImage imageWithKOImage:@"../RecKit/toshiba_images/DWI-rinkan-2/b200-proc/b200-1.img"];
//	img = [RecImage imageWithKOImage:@"../RecKit/nciRecon/IMG_bold.sav"];
	img = [RecImage imageWithKOImage:@"/Users/oshio/images/NCI/NIRS/2018-0627/results/2/IMG_mg"];
//	[img magnitude];
	xDim = [img xDim];
	yDim = [img yDim];
	nSlc = [img zDim];

	nSlc = 800;
	[img crop:[img zLoop] to:nSlc startAt:0];
	[img saveAsKOImage:@"IMG_img"];

	printf("x/y/slc = %d/%d/%d\n", xDim, yDim, nSlc);
	img_a = [RecImage imageOfType:RECIMAGE_REAL xDim:xDim * yDim yDim:nSlc];

	[img_a copyImageData:img];
	[img_a saveAsKOImage:@"IMG_A"];

//	A = Num_im_to_m(img_a);
	A = [img_a toMatrix];

// PCA (new interface)
	if (0) {
		printf("new interface\n");
		res = [A svd];
//		img_u = Num_m_to_im(sres->U);
		img_u = [[res objectForKey:@"U"] toRecImage];
		[img_u trans];
		[img_u saveAsKOImage:@"IMG_Ut"];	// Ut

//		tmp = Num_m_to_im(sres->Vt);
		tmp = [[ res objectForKey:@"Vt"] toRecImage];
		[tmp saveAsKOImage:@"IMG_Vt"];
		img_e = [RecImage imageWithImage:img];
		[img_e copyImageData:tmp];
		[img_e saveAsKOImage:@"IMG_PCA"];
	}

// == ICA
	if (1) {
		int nc = 20;
		printf("ICA (new interface)\n");
		res = [A icaForNC:nc];

		img_u = [[res objectForKey:@"U"] toRecImage];
		[img_u trans];
		[img_u saveAsKOImage:@"IMG_U"];	// U
		img_u = [[res objectForKey:@"XW"] toRecImage];
		[img_u trans];
		[img_u saveAsKOImage:@"IMG_XW"];	// XW
		img_u = [[res objectForKey:@"W"] toRecImage];
		[img_u saveAsKOImage:@"IMG_W"];	// W
		img_e = [img copy];
		[img_e crop:[img_e zLoop] to:nc startAt:0];
		[img_e copyImageData:[[res objectForKey:@"Y"] toRecImage]]; // -> make "copyMatrixData" method
		[img_e saveAsKOImage:@"IMG_Y"];
	}
}

// example #1
void
test19()
{
	Num_mat			*S, *A, *X;
	int				i, n = 5000;
	Num_ica_result	*res;
	float			hist[100];

	system("rm IMG_*");

	S = Num_new_mat(n, 2);
	X = Num_new_mat(n, 2);
	for (i = 0; i < n * 2; i++) {
		S->data[i] = Num_unif(0, 1);
	//	printf("%d %f\n", i, S->data[i]);
	}
	saveAsKOImage(S, @"IMG_S");

	A = Num_new_mat(2, 2);
	A->data[0] = 0.7;
	A->data[1] = 0.7;
	A->data[2] = -0.7;
	A->data[3] = 0.7;
	Num_mmul(X, S, A);
	saveAsKOImage(X, @"IMG_X");
	histogram(hist, 100, X->data + n, n, -5, 5);
	for (i = 0; i < 100; i++) {
		printf("%d %f\n", i, hist[i]);
	}

	res = Num_ica(X, 2);
	histogram(hist, 100, res->WX->data + 0, n, -10, 10);
	for (i = 0; i < 100; i++) {
		printf("%d %f\n", i, hist[i]);
	}
	exit(0);
}

// example #2
void
test20()
{
	Num_mat			*S, *A, *X;
	int				i, n = 1000;
	Num_ica_result	*res;
	float			th;

	system("rm IMG_*");

	S = Num_new_mat(n, 2);
	X = Num_new_mat(n, 2);
	for (i = 0; i < n; i++) {
		th = (float)i / 20;
		S->data[i] = sin(th);
	}
	for (i = 0; i < n; i++) {
		th = ((float)(i % 200) - 100) / 200;
		S->data[i + n] = th;
	}
	saveAsKOImage(S, @"IMG_S");

	A = Num_new_mat(2, 2);
	A->data[0] = 0.291;
	A->data[1] = 0.6557;
	A->data[2] = -0.5439;
	A->data[3] = 0.5572;
	Num_mmul(X, S, A);
	saveAsKOImage(X, @"IMG_X");

	res = Num_ica(X, 2);
	saveAsKOImage(res->WX, @"IMG_ans");
}

void
test21()
{
    NumMatrix	*A, *B, *X;
	float		*pa, *pb;
    int         m, n;
	int			mode = 0;	// 0: square, 1: under, 2: over, 3: single b vector, 4: NNLS

	switch (mode) {
	case 0:
	// === square matrix ====
		m = 4, n = 3;
		A = [NumMatrix matrixOfType:NUM_REAL nRow:m nCol:n];
		B = [NumMatrix matrixOfType:NUM_REAL nRow:m nCol:2];
	 
		pa = [A data];
		pb = [B data];

	// 4 x 3
		pa[0] = 74;	pa[4] = 71;	pa[8]  = 52;
		pa[1] = 87;	pa[5] = 74;	pa[9]  = 46;
		pa[2] = 72;	pa[6] = 2;	pa[10] = 7;
		pa[3] = 80;	pa[7] = 89;	pa[11] = 71;

	// 4 x 2
		pb[0] = 49;	pb[4] = 1;
		pb[1] = 67;	pb[5] = 3;
		pb[2] = 68;	pb[6] = 1;
		pb[3] = 20;	pb[7] = 1;

		X = [A solveLinear:B];
		[X dump];
		X = [A solveLinearNN:B];
		[X dump];
	}
}

void
test22()	// low level
{
    NumMatrix	*A, *B, *C;
	NumMatrix	*c;
	float		*pa;

	A = [NumMatrix matrixOfType:NUM_REAL nRow:3 nCol:4];
	pa = [A data];
	// 3 x 3
	pa[0] = 1;	pa[3] = 2;	pa[6] = 3;	pa[9] = 1;
	pa[1] = 5;	pa[4] = 6;	pa[7] = 7;	pa[10] = 2;
	pa[2] = 1;	pa[5] = 1;	pa[8] = 3;	pa[11] = 3;

	[A dump];

	B = [A trans];
	[B dump];
	
	C = [A multWithMat:B];
	[C dump];
	
	c = [C colVect:1];
	A = [A trans];
	[A dump];
	[c dump];
	[[A multWithMat:c] dump];
	
	c = [C rowVect:0];
	[c dump];
}

void
test23()
{
	NumMatrix	*A, *B, *X;		// A: exponential, B: data, X: coefficient (solution)
	RecImage	*img;
	int			i, j;
	float		t, tc1, tc2, fr;
	float		*p;
	int			n_time = 10;
	int			n_tc = 10;

	A = [NumMatrix matrixOfType:NUM_REAL nRow:n_time nCol:n_tc];
	B = [NumMatrix matrixOfType:NUM_REAL nRow:n_time nCol:1];

	p = [A data];
	for (i = 0; i < n_tc; i++) {
	//	tc1 = (float)i * 10 + 10;
		tc1 = 10 * pow(1.8, (float)i);
		printf("tc = %f\n", tc1);
		for (j = 0; j < n_time; j++) {
			t = (float)j * 10;
			p[i * n_time + j] = exp(-t / tc1);
		}
	}
	img = [A toRecImage];
	[img saveAsKOImage:@"img_a.img"];
	p = [B data];
	tc1 = 10;
	tc2 = 31;
	fr = 1.0;
	for (i = 0; i < n_time; i++) {
		t = (float)i * 10;
		p[i] = fr * exp(-t / tc1) + (1 - fr) * exp(-t / tc2);
		printf("%f %f\n", t, p[i]);
	}

	X = [A solveLinear:B];
	[X dump];
	X = [A solveLinearNN:B];
	[X dump];
}
