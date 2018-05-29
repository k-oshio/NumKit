//
// NumLinear.m
//
//  == plans ==
//	use BLAS & LAPACK

#import "NumKit.h"

@implementation NumVector

+ (NumVector *)vectorOfType:(int)type length:(int)n
{
	NumVector	*v;
	v = [NumVector alloc];
	return [v initWithType:type length:n];
}

- (NumVector *)initWithType:(int)tp length:(int)len
{
	int		actLen = len;

	self = [super init];
	if (!self) return nil;

	if (tp == NUM_COMPLEX) {
		len *= 2;
	}
	data = [NSMutableData dataWithLength:sizeof(float) * len];
	length = actLen;
	type = tp;

	return self;
}

- (int)length
{
	return length;
}

- (float *)data
{
	return [data mutableBytes];
}

- (void)clear
{
	int		i, n;
	float	*p = [self data];

	n = length;
	if (type == NUM_COMPLEX) {
		n *= 2;
	}
	for (i = 0; i < n; i++) {
		p[i] = 0;
	}
}

- (void)normal
{
	int		i, n;
	float	*p = [self data];

	n = length;
	if (type == NUM_COMPLEX) {
		n *= 2;
	}
	for (i = 0; i < n; i++) {
		p[i] = Num_nrml(0, 1);
	}
}

- (void)ramp
{
	int		i, n;
	float	*p = [self data];

	n = length;
	if (type == NUM_COMPLEX) {
		n *= 2;
	}
	for (i = 0; i < n; i++) {
		p[i] = (float)i;
	}
}

@end


@implementation NumMatrix

+ (NumMatrix *)matrixOfType:(int)tp nrow:(int)nr ncol:(int)nc
{
	NumMatrix	*m;

	m = [NumMatrix alloc];
	return [m initWithType:tp nrow:nr ncol:nc];
}

+ (NumMatrix *)matrixWithNumMat:(Num_mat *)mat
{
	NumMatrix	*m;

	m = [NumMatrix alloc];
	return [m initWithType:mat->type nrow:mat->nr ncol:mat->nc];
}

- (NumMatrix *)initWithType:(int)tp nrow:(int)nr ncol:(int)nc
{
	int		len;

	self = [super init];
	if (!self) return nil;

	len = nr * nc;
	if (tp == NUM_COMPLEX) {
		len *= 2;
	}
	data = [NSMutableData dataWithLength:sizeof(float) *len];
	type = tp;
	nrow = nr;
	ncol = nc;

	return self;
}

- (float *)data
{
	return [data mutableBytes];
}

- (int)type
{
	return type;
}

- (int)nrow
{
	return nrow;
}

- (int)ncol
{
	return ncol;
}

- (void)clear
{
	int		i, n;
	float	*p = [self data];

	n = nrow * ncol;
	if (type == NUM_COMPLEX) {
		n *= 2;
	}
	for (i = 0; i < n; i++) {
		p[i] = 0;
	}
}

- (void)normal
{
	int		i, n;
	float	*p = [self data];

	n = nrow * ncol;
	if (type == NUM_COMPLEX) {
		n *= 2;
	}
	for (i = 0; i < n; i++) {
		p[i] = Num_nrml(0, 1);
	}
}

- (NumMatrix *)copy
{
	NumMatrix	*m = [NumMatrix matrixOfType:type nrow:nrow ncol:ncol];
	int			i, n;
	float		*p1 = [self data];
	float		*p2 = [m data];

	n = nrow * ncol;
	if (type == NUM_COMPLEX) {
		n *= 2;
	}
	for (i = 0; i < n; i++) {
		p2[i] = p1[i];
	}
	return m;
}

- (NumMatrix *)trans
{
	NumMatrix *m;

	return m;
}

- (NumMatrix *)unitMatrixOfDim:(int)n
{
	NumMatrix *m;

	return m;
}

- (NumMatrix *)multWithVector:(NumVector *)v
{
	NumMatrix *m;

	return m;
}

- (NumMatrix *)multWithMatrix:(NumMatrix *)m
{
	NumMatrix *res;

	return res;
}

- (NumMatrix *)multWithMatrix:(NumMatrix *)m transSelf:(BOOL)ts transMat:(BOOL)tm
{
	NumMatrix *res;

	return res;
}

- (NumMatrix *)diagMatrixWithVector:(NumVector *)v
{
	NumMatrix *res;

	return res;
}


- (NSDictionary *)svd	// returns dict with entry: "U", "s", "Vt"
{
	NSDictionary	*dict;

	return dict;
}

- (NSDictionary *)ica	// returns dict with entry: "WX", "W"
{
	NSDictionary	*dict;

	return dict;
}


@end

// conversion to/from RecImage
// real / complex
Num_mat	*
Num_im_to_m(RecImage *im)
{
	Num_mat	*m;

	m = Num_alloc_m_with_im(im);
	Num_copy_im_to_m(m, im);

	return m;
}

Num_vec *
Num_im_to_v(RecImage *im)
{
	Num_vec	*v;

	v = Num_alloc_v_with_im(im);
	Num_copy_im_to_v(v, im);

	return v;
}

RecImage *
Num_m_to_im(Num_mat *m)
{
	RecImage	*im;

	im = Num_alloc_im_with_m(m);
	Num_copy_m_to_im(im, m);

	return im;
}

RecImage *
Num_v_to_im(Num_vec *v)
{
	RecImage	*im;

	im = Num_alloc_im_with_v(v);
	Num_copy_v_to_im(im, v);

	return im;
}

RecImage *
Num_alloc_im_with_m(Num_mat *m)
{
	RecImage	*im;

	if (m->type == RECIMAGE_REAL) {
		im = [RecImage imageOfType:RECIMAGE_REAL xDim:m->nc yDim:m->nr];
	} else {
		im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:m->nc yDim:m->nr];
	}
	return im;
}

RecImage *
Num_alloc_im_with_v(Num_vec *v)
{
	RecImage	*im;

	if (v->type == RECIMAGE_REAL) {
		im = [RecImage imageOfType:RECIMAGE_REAL xDim:v->n];
	} else {
		im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:v->n];
	}
	return im;
}

Num_mat *
Num_alloc_m_with_im(RecImage *im)
{
	Num_mat	*m;

	if ([im type] == RECIMAGE_REAL) {
		m = Num_new_mat([im yDim], [im xDim]);
	} else {
		m = Num_new_cmat([im yDim], [im xDim]);
	}
	return m;
}

Num_vec *
Num_alloc_v_with_im(RecImage *im)
{
	Num_vec	*v;

	if ([im type] == RECIMAGE_REAL) {
		v = Num_new_vec([im xDim]);
	} else {
		v = Num_new_cvec([im xDim]);
	}
	return v;
}

void
Num_copy_im_to_m(Num_mat *m, RecImage *im)
{
	int		c, r, ixc, ixr;
	int		xDim, yDim;
	float	*p, *q;

	xDim = [im xDim];
	yDim = [im yDim];

	if (([im type] == RECIMAGE_REAL && m->type == NUM_COMPLEX)
		|| ([im type] == RECIMAGE_COMPLEX && m->type == NUM_REAL)) {
		printf("Type mismatch\n");
		return;
	}

	if (m->type == NUM_REAL) {
		p = [im data];
		for (c = 0; c < m->nc; c++) {	// x
			for (r = 0; r < m->nr; r++) {	// y
				ixc = c * m->nr + r;
				ixr = r * m->nc + c;
				m->data[ixc] = p[ixr];
			}
		}
	}
	if (m->type == NUM_COMPLEX) {
		p = [im real];
		q = [im imag];
		ixc = 0;
		for (c = 0; c < m->nc; c++) {
			for (r = 0; r < m->nr; r++) {
				ixc = c * m->nr + r;
				ixr = r * m->nc + c;
				m->data[ixc * 2]     = p[ixr];
				m->data[ixc * 2 + 1] = q[ixr];
			}
		}
	}
}

void
Num_copy_im_to_v(Num_vec *v, RecImage *im)
{
	int		i;
	int		xDim;
	float	*p, *q;

	xDim = [im xDim];

	if (([im type] == RECIMAGE_REAL && v->type == NUM_COMPLEX)
		|| ([im type] == RECIMAGE_COMPLEX && v->type == NUM_REAL)) {
		printf("Type mismatch\n");
		return;
	}

	if (v->type == NUM_REAL) {
		p = [im data];
		for (i = 0; i < v->n; i++) {	// x
			v->data[i] = p[i];
		}
	}
	if (v->type == NUM_COMPLEX) {
		p = [im real];
		q = [im imag];
		for (i = 0; i < v->n; i++) {	// x
			v->data[i * 2]     = p[i];
			v->data[i * 2 + 1] = q[i];
		}
	}
}

void
Num_copy_m_to_im(RecImage *im, Num_mat *m)
{
	int			r, c, ixc, ixr;
	int			xDim, yDim;
	float		*p, *q;

	// copy data
	xDim = [im xDim];
	yDim = [im yDim];

	if (m->type == NUM_REAL) {
		p = [im data];
		for (r = 0; r < yDim; r++) {
			for (c = 0; c < xDim; c++) {
				ixc = c * m->nr + r;
				ixr = r * m->nc + c;
				p[ixr] = m->data[ixc];
			}
		}
	}
	if (m->type == NUM_COMPLEX) {
		p = [im real];
		q = [im imag];
		for (r = 0; r < yDim; r++) {
			for (c = 0; c < xDim; c++) {
				ixc = c * m->nr + r;
				ixr = r * m->nc + c;
				p[ixr] = m->data[ixc * 2];
				q[ixr] = m->data[ixc * 2 + 1];
			}
		}
	}
}

void
Num_copy_v_to_im(RecImage *im, Num_vec *v)
{
	int			i;
	int			xDim;
	float		*p, *q;

	// copy data
	xDim = [im xDim];

	if (v->type == NUM_REAL) {
		p = [im data];
		for (i = 0; i < xDim; i++) {
			p[i] = v->data[i];
		}
	}
	if (v->type == NUM_COMPLEX) {
		p = [im data];
		q = p + [im dataLength];
		for (i = 0; i < xDim; i++) {
			p[i] = v->data[i * 2];
			q[i] = v->data[i * 2 + 1];
		}
	}
}

// ========= real ===========
Num_mat *
Num_new_mat(int nr, int nc)
{

    Num_mat		*m;
    m = (Num_mat *)malloc(sizeof(Num_mat));
	m->type = NUM_REAL;
	m->order = NUM_COL_MAJOR;	// just a reminder...
    m->nr = nr;
    m->nc = nc;
	m->ld = nr;
    m->data = VECTOR(nr * nc);

    return m;
}

Num_vec *
Num_new_vec(int n)
{
    Num_vec *v;

    v = (Num_vec *)malloc(sizeof(Num_vec));
	v->type = NUM_REAL;
    v->n = n;
    v->data = VECTOR(n);

    return v;
}

// ========= complex ===========
Num_mat *
Num_new_cmat(int nr, int nc)
{
    Num_mat		*m;
    m = (Num_mat *)malloc(sizeof(Num_mat));
	m->type = NUM_COMPLEX;
	m->order = NUM_COL_MAJOR;	// just a reminder...
    m->nr = nr;
    m->nc = nc;
	m->ld = nr;
    m->data = VECTOR(nr * nc * 2);

    return m;
}

Num_vec *
Num_new_cvec(int n)
{
    Num_vec		*v;

    v = (Num_vec *)malloc(sizeof(Num_vec));
	v->type = NUM_COMPLEX;
    v->n = n;
    v->data = VECTOR(n * 2);

    return v;
}

void
Num_make_cmat(Num_mat *m)
{
	int			i, n;
	float		*data;

	n = m->nr * m->nc;
	data = (float *)calloc(sizeof(float), n * 2);
	for (i = 0; i < n; i++) {
		data[i * 2] = m->data[i];
	}
	free(m->data);
	m->data = data;
	m->type = NUM_COMPLEX;
}

void
Num_make_cvec(Num_vec *v)
{
	int			i, n;
	float		*data;

	n = v->n;
	data = (float *)calloc(sizeof(float), n * 2);
	for (i = 0; i < n; i++) {
		data[i * 2] = v->data[i];
	}
	free(v->data);
	v->data = data;
	v->type = NUM_COMPLEX;
}

void
Num_free_mat(Num_mat *m)
{
    if (m) {
        if (m->data) {
            free(m->data);
        }
        free(m);
    }
}

void
Num_free_vec(Num_vec *v)
{
    if (v) {
        if (v->data) {
            free(v->data);
        }
        free(v);
    }
}

void
Num_copy_vec(Num_vec *a, Num_vec *b)						// b <- a
{
	if (a->type == NUM_REAL && b->type == NUM_COMPLEX) {
		Num_make_cvec(a);
	}
	if (a->type == NUM_COMPLEX && b->type == NUM_REAL) {
		Num_make_cvec(b);
	}
	switch (a->type) {
	case NUM_REAL :
		cblas_scopy(a->n, a->data, 1, b->data, 1);
		break;
	case NUM_COMPLEX :
		cblas_ccopy(a->n, a->data, 1, b->data, 1);
		break;
	}
}

// ### real only
void
Num_vadd(Num_vec *v1, Num_vec *v2, Num_vec *res)
{
	int		n = v1->n;

	if (res == v1 || res == v2) {	// vDSP doesn't handle in-place op
		Num_vec	*buf = Num_new_vec(n);
		vDSP_vadd(v1->data, 1, v2->data, 1, buf->data, 1, n);
		Num_free_vec(buf);
	} else {
		vDSP_vadd(v1->data, 1, v2->data, 1, res->data, 1, n);
	}
}

// ### real only
void
Num_vsma(Num_vec *x, float a, Num_vec *y)    // y = ax + y
{
	cblas_saxpy(x->n, a, x->data, 1, y->data, 1);
}

// ### real only
float
Num_dotpr(Num_vec *v1, Num_vec *v2)
{
    float   sum;

	if (v1->type == NUM_COMPLEX || v2->type == NUM_COMPLEX) {
		Num_error("dotpr is real only");
		return 0;
	}
    if (v1->n != v2->n) {
		Num_error("dotpr dim mismatch");
		return 0;
	}
    if (v1->type != v2->type) {
		Num_error("dotpr type mismatch");
		return 0;
	}
	sum = cblas_sdot(v1->n, v1->data, 1, v2->data, 1);

    return sum;
}

int dbg = 0;
// ### real only
void
Num_mvmul(Num_vec *b, Num_mat *a, Num_vec *x) // b = Ax
{
    if (x->n != a->nc || b->n != a->nr) {
        Num_error("mvmul dim mismatch");
        return;
    }
	// BLAS matrix-matrix multiplication
	cblas_sgemv(CblasColMajor, CblasNoTrans,
		a->nr, a->nc, 1.0, a->data, a->ld, x->data, 1, 0, b->data, 1);
// dbg
	if (dbg) {
		printf("Num_mvmul: A:%d,%d, x:%d, b:%dn", a->nr, a->nc, x->n, b->n);
	}
}

// ### real only
void
Num_mmul(Num_mat *c, Num_mat *a, Num_mat *b)   // C = AB
{
    if (a->nc != b->nr || b->nc != c->nc || a->nr != c->nr) {
        Num_error("mmul dim mismatch");
    }
	// BLAS matrix-matrix multiplication
	cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		a->nr, b->nc, a->nc, 1.0, a->data, a->ld, b->data, b->ld, 0, c->data, c->ld);
// dbg
	if (dbg) {
		printf("Num_mmul: A:%d,%d, B:%d,%d, C:%d,%d\n", a->nr, a->nc, b->nr, b->nc, c->nr, c->nc);
	}
}

// probably ok ...
void
Num_mtmul(Num_mat *c, Num_mat *a, BOOL ta, Num_mat *b, BOOL tb)	// C = A(t)B(t)
{
	int atr, btr;
	int	na, nb, nk;

	if (ta) {
		atr = CblasTrans;
		na = a->nc;
		nk = a->nr;
	} else {
		atr = CblasNoTrans;
		na = a->nr;
		nk = a->nc;
	}
	if (tb) {
		btr = CblasTrans;
		nb = b->nr;
	} else {
		btr = CblasNoTrans;
		nb = b->nc;
	}
	// BLAS matrix-matrix multiplication
	cblas_sgemm(CblasColMajor, atr, btr,
		na, nb, nk, 1.0, a->data, a->ld, b->data, b->ld, 0, c->data, c->ld);
// dbg
	if (dbg) {
		printf("Num_mmul: A:%d,%d, B:%d,%d, C:%d,%d\n", a->nr, a->nc, b->nr, b->nc, c->nr, c->nc);
	}
}

// ### real only
int
Num_gaussj(Num_mat *a, Num_mat *b)  // A: n x n, B: n x m
{
    int     n, m;
    float   *p, *q;
    int     *indxc, *indxr, *ipiv;
    int     i, icol, irow, j, k, l, ll;
    float   big, dum, pivinv;

    if (a->nc != a->nr) {
        Num_error("Num_gaussj: A not square");
        return -1;
    }
    if (b && a->nr != b->nr) {
        Num_error("Num_gaussj: B dim mismatch");
		return -1;
    }

    n = a->nc;
    p = a->data;

    if (b) {
        m = b->nc;
        q = b->data;
    } else {
        m = 0;
    }

    indxc = (int *)malloc(n * sizeof(int));
    indxr = (int *)malloc(n * sizeof(int));
    ipiv  = (int *)malloc(n * sizeof(int));
    for (i = 0; i < n; i++) {
        indxc[i] = indxr[i] = ipiv[i] = 0;
    }
    icol = irow = 0;
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++) {
            if (ipiv[j] != 1) {
                for (k = 0; k < n; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(p[j*n + k]) > big) {
                            big = fabs(p[j*n + k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        ipiv[icol] += 1;
        if (irow != icol) {
            vDSP_vswap(p + irow * n, 1, p + icol * n, 1, n);
            if (b) {
                vDSP_vswap(q + irow * m, 1, q + icol * m, 1, m);
            }
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (p[icol*n + icol] == 0.0) {
    //    if (fabs(p[icol*n + icol]) < 0.0001) {
    //        printf("gaussj() : Singular matrix\n");
            return -1;
        }
        pivinv = 1.0 / p[icol*n + icol];
        p[icol*n + icol] = 1.0;
        vDSP_vsmul(p + icol*n, 1, &pivinv, p + icol*n, 1, n);
        if (b) {
            vDSP_vsmul(q + icol*m, 1, &pivinv, q + icol*m, 1, m);
        }
        for (ll = 0; ll < n; ll++) {
            if (ll != icol) {
                dum = -p[ll*n + icol];
                p[ll*n + icol] = 0.0;
                for (l = 0; l < n; l++) {
                    p[ll*n + l] += p[icol*n + l] * dum;
                }
                if (b) {
                    for (l = 0; l < m; l++) {
                        q[ll*m + l] += q[icol*m + l] * dum;
                    }
                }
            }
        }
    }   // i loop

    for (l = n-1; l >= 0; l--) {
        if (indxr[l] != indxc[l]) {
            vDSP_vswap(p + indxr[l], n, p + indxc[l], n, n);
        }
    }
    free(ipiv);
    free(indxr);
    free(indxc);

	return 0;
}

void
Num_clear_mat(Num_mat *a)                  // A <- 0
{
	switch (a->type) {
	case NUM_REAL :
		vDSP_vclr(a->data, 1, a->nr * a->nc);
		break;
	case NUM_COMPLEX :
		vDSP_vclr(a->data, 1, a->nr * a->nc * 2);
		break;
	}
}

void
Num_unit_mat(Num_mat *a)                        // A <- I
{
    int     i, n;
    n = a->nr < a->nc ? a->nr : a->nc;  // min(nc, nr)
    Num_clear_mat(a);
    for (i = 0; i < n; i++) {
        a->data[i * a->nr + i] = 1;	// imag part is 0, which is ok
    }
}

Num_mat *
Num_trans(Num_mat *a)            // B = At, A:mxn, B:nxm
{
	Num_mat	*b;

	b = Num_new_mat(a->nc, a->nr);
    vDSP_mtrans(a->data, 1, b->data, 1, a->nr, a->nc);
	return b;
}

void
Num_trans_ip(Num_mat *b, Num_mat *a)
{
    vDSP_mtrans(a->data, 1, b->data, 1, a->nr, a->nc);
}

void
Num_hermit(Num_mat *b, Num_mat *a)						// B = At, A:mxn, B:nxm
{
	int		r, c;

// chk header
	if (a->type != b->type) {
		Num_error("Num_trans: type mismatch");
	}
	if (a->nr * a->nc != b->nr * b->nc) {
		Num_error("Num_trans: size mismatch");
	}
// set header
	b->nc = a->nr;
	b->nr = a->nc;
// copy data
//    vDSP_mtrans(a->data, 1, b->data, 1, a->nr, a->nc);  //in-place doesn't work

	for (r = 0; r < a->nr; r++) {
		for (c = 0; c < a->nc; c++) {	// col major
			b->data[(c * a->nr + r) * 2]     =  a->data[(r * a->nc + c) * 2];
			b->data[(c * a->nr + r) * 2 + 1] = -a->data[(r * a->nc + c) * 2 + 1];
		}
	}
}

void
Num_copy_mat(Num_mat *b, Num_mat *a)         // B <- A
{
	int		i, n;
// chk header
	if (a->type != b->type) {
		Num_error("Num_copy: type mismatch");
	}
	if (a->nr * a->nc != b->nr * b->nc) {
		Num_error("Num_copy: size mismatch");
	}
	switch (a->type) {
	case NUM_REAL :
		n = a->nr * a->nc;
		break;
	case NUM_COMPLEX :
		n = a->nr * a->nc * 2;
		break;
	}
	for (i = 0; i < n; i++) {
		b->data[i] = a->data[i];
	}
}

void
Num_copy_sub_mat(Num_mat *b, Num_mat *a, int r0, int rn, int c0, int cn)
{
	int		i, j, na, nb;
// chk header
	if (a->type != b->type) {
		Num_error("Num_copy: type mismatch");
	}
	switch (a->type) {
	case NUM_REAL :
		for (i = 0; i < rn; i++) {
			for (j = 0; j < cn; j++) {
				b->data[(i + r0) * b->nr + j + c0] = a->data[i * a->nr + j];
			}
		}
		break;
	case NUM_COMPLEX :
		na = a->nr * a->nc;
		nb = b->nr * b->nc;
		for (i = 0; i < rn; i++) {
			for (j = 0; j < cn; j++) {
				b->data[(i + r0) * b->nr + j + c0]      = a->data[i * a->nr + j];
				b->data[(i + r0) * b->nr + j + c0 + nb] = a->data[i * a->nr + j + na];
			}
		}
		break;
	}
}

void
Num_col_vec(Num_vec *x, Num_mat *A, int ix)	// copy i-th col-v of A
{
	int		i;

	if (A->nr != x->n) {
		Num_error("Num_col_vec: size mismatch");
	}
	for (i = 0; i < x->n; i++) {
		x->data[i] = A->data[ix * A->nc + i];
	}
}

void
Num_row_vec(Num_vec *x, Num_mat *A, int ix)	// copy i-th row-v of A
{
	int		i;

	if (A->nc != x->n) {
		Num_error("Num_row_vec: size mismatch");
	}
	for (i = 0; i < x->n; i++) {
		x->data[i] = A->data[i * A->nr + ix];
	}
}

void
Num_pinv(Num_mat *b, Num_mat *a)             // B = A+, A:mxn, B:nxm
{
    int     m, n;
    Num_mat *at, *att;

    m = a->nr, n = a->nc;
    if (m == n) {   // square
        Num_copy_mat(b, a);
        Num_gaussj(b, NULL);
    } else
    if (m > n) {    // more eq than unk -> LS
        at = Num_trans(a);
        att = Num_new_mat(n, n);
        Num_mmul(att, at, a);   // ATT = At A
        Num_gaussj(att, NULL);        // (At A)-1
        Num_mmul(b, att, at);   // B = (At A)-1 At
        Num_free_mat(at);
        Num_free_mat(att);
    } else {        // more unk than eq -> min norm
        at = Num_trans(a);
        att = Num_new_mat(m, m);
        Num_mmul(att, a, at);   // ATT = A * At
        Num_gaussj(att, NULL);        // (A At)-1
        Num_mmul(b, at, att);   // B = At (A At)-1
        Num_free_mat(at);
        Num_free_mat(att);
    }
}

void
Num_tocm(Num_mat *m)	// convert to column major
{
	int		i, len;
	int		c, r;
	int		ixr, ixc;
	float	*p, *buf;

	if (m->order == NUM_COL_MAJOR) {
		return;
	}
	// trans data
	len = m->nc * m->nr;
	p = m->data;
	buf = (float *)malloc(sizeof(float) * len);
	// real part
	for (i = 0; i < len; i++) {
		buf[i] = p[i];
	}
	for (c = 0; c < m->nc; c++) {
		for (r = 0; r < m->nr; r++) {
			ixr = r * m->nc + c;
			ixc = c * m->nr + r;
			p[ixc] = buf[ixr];
		}
	}
	// imag part
	p += len;
	if (m->type == NUM_COMPLEX) {
		for (i = 0; i < len; i++) {
			buf[i] = p[i];
		}
		for (c = 0; c < m->nc; c++) {
			for (r = 0; r < m->nr; r++) {
				ixr = r * m->nc + c;
				ixc = c * m->nr + r;
				p[ixc] = buf[ixr];
			}
		}
	}

	// set mat header
	m->order = NUM_COL_MAJOR;
	m->ld = m->nc;

	free(buf);
}

void
Num_torm(Num_mat *m)
{
	int		i, len;
	int		c, r;
	int		ixr, ixc;
	float	*p, *buf;

	if (m->order == NUM_ROW_MAJOR) {
		return;
	}
	// trans data
	len = m->nc * m->nr;
	p = m->data;
	buf = (float *)malloc(sizeof(float) * len);
	// real part
	for (i = 0; i < len; i++) {
		buf[i] = p[i];
	}
	for (c = 0; c < m->nc; c++) {
		for (r = 0; r < m->nr; r++) {
			ixr = r * m->nc + c;
			ixc = c * m->nr + r;
			p[ixr] = buf[ixc];
		}
	}
	// imag part
	p += len;
	for (i = 0; i < len; i++) {
		buf[i] = p[i];
	}
	for (c = 0; c < m->nc; c++) {
		for (r = 0; r < m->nr; r++) {
			ixr = r * m->nc + c;
			ixc = c * m->nr + r;
			p[ixr] = buf[ixc];
		}
	}
	// set mat header
	m->order = NUM_ROW_MAJOR;
	m->ld = m->nr;

	free(buf);
}

// gramm-schmidt orthogonalization
void
Num_grmsch(Num_mat *a)
{
	// ### not done yet
}

// complex not done yet ###
void
Num_evd_ref(Num_mat *A, Num_mat *Evec, Num_vec *eval)
{
	float		*work;
	int			lwork, info;
	int			n;

	if (A->nr != A->nc) {
		printf("A is not square\n");
		return;
	}
	// copy input (no side effects)
	n = A->nr;
	Evec = Num_new_mat(n, n);
	Num_copy_mat(Evec, A);

	// get work size
	lwork = -1;
	work = (float *)malloc(sizeof(float) * 1);
	ssyev_("V", "U", &Evec->nc, Evec->data, &Evec->ld, eval->data, work, &lwork, &info);
	lwork = work[0];
	free(work);

	work = (float *)malloc(sizeof(float) * lwork);
	ssyev_("V", "U", &Evec->nc, Evec->data, &Evec->ld, eval->data, work, &lwork, &info);
	free(work);
}

// "S" option
// A		: input (not changed)
// U, s, Vt : output
void
Num_svd_ref(Num_mat *A, Num_mat *U, Num_vec *s, Num_mat *Vt)
{
	float		*rwork, *cwork, *adata;
	int			lwork, rworklen, info;
	int			i, nr, nc, n, nrc, ld;

	nr = A->nr;
	nc = A->nc;
	nrc = nr * nc;
	ld = nr;
	n = MIN(nr, nc);

	if (A->type == NUM_REAL) {
	// get work size
		lwork = -1;
		rwork = (float *)malloc(sizeof(float) * 1);
		adata = (float *)malloc(sizeof(float) * nrc);
		sgesvd_("S", "S", &nr, &nc, adata, &ld,
			s->data, U->data, &U->ld, Vt->data, &Vt->ld,
			rwork, &lwork, &info);
		lwork = rwork[0];
		free(rwork);

		//printf("sizeof work area %d\n", lwork);
		// actual proc
		rwork = (float *)malloc(sizeof(float) * lwork);
		// copy A
		for (i = 0; i < nrc; i++) {
			adata[i] = A->data[i];
		}
		sgesvd_("S", "S", &nr, &nc, adata, &ld,
			s->data, U->data, &U->ld, Vt->data, &Vt->ld, rwork,
			&lwork, &info);
		free(rwork);
		free(adata);
	} else {	// NUM_COMPLEX
	// get work size
		rworklen = 2 * MIN(nr, nc) * MAX(nr, nc);
		adata = (float *)malloc(sizeof(float) * nrc * 2);
		lwork = -1;
		rwork = (float *)malloc(sizeof(float) * rworklen);
		cwork = (float *)malloc(sizeof(float) * 2);
		cgesvd_("S", "S", &nr, &nc, (__CLPK_complex *)adata, &ld,
			s->data, (__CLPK_complex *)U->data, &U->ld, (__CLPK_complex *)Vt->data, &Vt->ld,
			(__CLPK_complex *)cwork, &lwork, rwork, &info);
		lwork = cwork[0];
		//printf("lwork[0] = %d\n", lwork);
		free(cwork);

		//printf("sizeof work area %d\n", lwork);
		// actual proc
		cwork = (float *)malloc(sizeof(float) * lwork * 2);
		// copy A
		for (i = 0; i < nrc * 2; i++) {
			adata[i] = A->data[i];
		}
		cgesvd_("S", "S", &nr, &nc, (__CLPK_complex *)A->data, &ld,
			s->data, (__CLPK_complex *)U->data, &U->ld, (__CLPK_complex *)Vt->data, &Vt->ld,
			(__CLPK_complex *)cwork, &lwork, rwork, &info);
		free(rwork);
		free(cwork);
		free(adata);
	}
}

// 
int
Num_inv(Num_mat *A, Num_mat *B)
{
	int		n, nrhs;
	int		lda, ldb, info;
	int		*ipiv;

	n = A->nc;
	nrhs = B->nc;
	lda = A->nr;
	ldb = B->nr;
	ipiv = (int *)malloc(sizeof(int) * n);

	sgesv_(&n, &nrhs, A->data, &lda, ipiv, B->data, &ldb, &info);

	free(ipiv);
	return info;
}

Num_svd_result *
Num_new_svd_result(Num_mat *A)
{
	Num_svd_result	*res;
	int				nc, nr, n;

	// dim
	nc = A->nc;
	nr = A->nr;
	n = MIN(nr, nc);

	// alloc res struct
	res = (Num_svd_result *)malloc(sizeof(Num_svd_result));
	if (A->type == NUM_REAL) { // real mat
		res->U  = Num_new_mat(nr, n);
		res->Vt = Num_new_mat(n, nc);
		res->s  = Num_new_vec(n);
	} else {				// complex mat
		res->U  = Num_new_cmat(nr, n);
		res->Vt = Num_new_cmat(n, nc);
		res->s  = Num_new_cvec(n);
	}

	return res;
}

Num_svd_result *
Num_svd(Num_mat *A)
{
	Num_svd_result	*res;

/*
	int				nc, nr, n;

	// dim
	nc = A->nc;
	nr = A->nr;
	n = MIN(nr, nc);

	// alloc res struct
	res = (Num_svd_result *)malloc(sizeof(Num_svd_result));
	if (A->type == NUM_REAL) { // real mat
		res->U  = Num_new_mat(nr, n);
		res->Vt = Num_new_mat(n, nc);
		res->s  = Num_new_vec(n);
	} else {				// complex mat
		res->U  = Num_new_cmat(nr, n);
		res->Vt = Num_new_cmat(n, nc);
		res->s  = Num_new_cvec(n);
	}
*/
	res = Num_new_svd_result(A);

	// call ref
	Num_svd_ref(A, res->U, res->s, res->Vt);

	// negate U & Vt
	Num_negate_mat(res->U);
	Num_negate_mat(res->Vt);

	return res;
}

Num_evd_result *
Num_evd(Num_mat *A)
{
	Num_evd_result	*res;
	int				n = A->nr;	// symmetric

	// alloc res struct
	res = (Num_evd_result *)malloc(sizeof(Num_evd_result));
	if (A->type == NUM_REAL) {	// resl
		res->Evec = Num_new_mat(n, n);
		res->eval = Num_new_vec(n);
	}
	
	// call ref
	Num_evd_ref(A, res->Evec, res->eval);

	return res;
}

// make auc = 1
void
Num_normalize_vec(Num_vec *v)
{
	float	sum = 0;
	int		i;

	for (i = 0; i < v->n; i++) {
	//	sum += v->data[i] * v->data[i];
		sum += v->data[i];
	}
//	sum = 1.0 / sqrt(sum);
	sum = 1.0 / sum;
	for (i = 0; i < v->n; i++) {
		v->data[i] *= sum;
	}

	//chk
	sum = 0;
	for (i = 0; i < v->n; i++) {
		sum += v->data[i];
	}
	printf("sum = %f\n", sum);
}

// orthogonalize square matrix A, using SVD
// A -> A', A'T A' = A' A'T = E
void
Num_orth_mat(Num_mat *A)
{
	Num_svd_result	*res;
	res = Num_svd(A);
	Num_mmul(A, res->U, res->Vt);
	Num_free_svd_result(res);
}

Num_mat	*
Num_diag(Num_vec *v)
{
	Num_mat *m;
	int		i, n = v->n;

	m = Num_new_mat(n, n);
	for (i = 0; i < n; i++) {
		m->data[i * n + i] = v->data[i];
	}
	return m;
}

void
Num_diag_ip(Num_mat *m, Num_vec *v)
{
	int		i, n = v->n;

	Num_clear_mat(m);
	for (i = 0; i < n; i++) {
		m->data[i * n + i] = v->data[i];
	}
}

void
Num_mmul_const(Num_mat *A, float c)
{
	int		i, n = A->nr * A->nc;

	for (i = 0; i < n; i++) {
		A->data[i] *= c;
	}
}

void
Num_col_norm(Num_mat *a)
{
	int		i, j;
	int		nr = a->nr;
	int		nc = a->nc;
	float	m, v, tmp;

	for (i = 0; i < nc; i++) {
		m = v = 0;
		for (j = 0; j < nr; j++) {
			m += a->data[i * nr + j];
			v += a->data[i * nr + j] * a->data[i * nr + j];
		}
		m /= nr;
		tmp = (v - nr * (m * m)) / (nr - 1);
		v = sqrt(tmp);
		for (j = 0; j < nr; j++) {
			a->data[i * nr + j] = (a->data[i * nr + j] - m) / v;
		}
	}
}

void
Num_col_center(Num_mat *a)
{
	int		i, j;
	int		nr = a->nr;
	int		nc = a->nc;
	float	m;

	for (i = 0; i < nc; i++) {
		m = 0;
		for (j = 0; j < nr; j++) {
			m += a->data[i * nr + j];
		}
		m /= nr;
		for (j = 0; j < nr; j++) {
			a->data[i * nr + j] -= m;
		}
	}
}

void
Num_row_norm(Num_mat *a)
{
	int		i, j;
	int		nr = a->nr;
	int		nc = a->nc;
	float	m, v, tmp;

	for (i = 0; i < nr; i++) {
		m = v = 0;
		for (j = 0; j < nc; j++) {
			m += a->data[j * nr + i];
			v += a->data[j * nr + i] * a->data[j * nr + i];
		}
		m /= nc;
		tmp = (v - nc * (m * m)) / (nc - 1);
		v = sqrt(tmp);
		printf("%d %f %f\n", i, m, v);
		for (j = 0; j < nc; j++) {
			a->data[j * nr + i] = (a->data[j * nr + i] - m) / v;
		}
	}
}

void
Num_row_center(Num_mat *a)
{
	int		i, j;
	int		nr = a->nr;
	int		nc = a->nc;
	float	m;

	for (i = 0; i < nr; i++) {
		m = 0;
		for (j = 0; j < nc; j++) {
			m += a->data[j * nr + i];
		}
		m /= nc;
		for (j = 0; j < nc; j++) {
			a->data[j * nr + i] -= m;
		}
	}
}

void
Num_scale_columns(Num_mat *A, Num_vec *v)
{
	int		i, j;

	if (A->nc != v->n) {
		Num_error("Num_scale_columns: size mismatch");
	}
	for (i = 0; i < A->nr; i++) {		// row
		for (j = 0; j < A->nc; j++) {	// col
			A->data[j * A->nr + i] *= v->data[j];
		}
	}
}

void
Num_scale_rows(Num_mat *A, Num_vec *v)
{
	int		i, j;

	if (A->nr != v->n) {
		Num_error("Num_scale_columns: size mismatch");
	}
	for (i = 0; i < A->nr; i++) {		// row
		for (j = 0; j < A->nc; j++) {	// col
			A->data[j * A->nr + i] *= v->data[i];
		}
	}
}

void
Num_negate_mat(Num_mat *A)					// for SVD result
{
	int		i, n = A->nc * A->nr;

	if (A->type == NUM_COMPLEX) {
		n *= 2;
	}
	for (i = 0; i < n; i++) {
		A->data[i] = -A->data[i];
	}
}

Num_mat *
Num_sort_mat(Num_vec *v)
{
	int			n = v->n;
	Num_mat	*m = Num_new_mat(n, n);
	int			i;
	int			(^compar)(const void *, const void *);
	typedef struct	ix {
		int		i;
		float	val;
	} ix;
	ix			*ixp = (ix *)malloc(sizeof(ix) * n);

	for (i = 0; i < n; i++) {
		ixp[i].i = i;
		ixp[i].val = v->data[i];
	}

	compar = ^int(const void *p1, const void *p2) {
		float pp = ((ix *)p1)->val;
		float qq = ((ix *)p2)->val;
		// descending
		if (pp > qq) {
			return -1;
		} else 
		if (pp == qq) {
			return 0;
		} else {
			return 1;
		}
	};
	qsort_b(ixp, n, sizeof(ix), compar);
	for (i = 0; i < n; i++) {
	//	printf("%d %d\n", i, ixp[i].i);
		m->data[ixp[i].i * n + i] = 1;
	}

	free(ixp);
	return m;
}

// dbg
// just copy image data
void
copy_mat_to_img(Num_mat *m, RecImage *img) {
	RecImage	*tmp_img;
	tmp_img = Num_m_to_im(m);
	[img copyImageData:tmp_img];
}

//======================== preserve until successfull port
Num_svd_result *
Num_ica_ref(Num_mat *X, int ncomp, int maxiter)
{
	Num_svd_result	*ires;	
	Num_svd_result	*sres;	
	int				n, p, ncn, ncnc;
	int				i, j;
	int				iter;
	double			tol;

	Num_mat			*Xt;				// p x n
	Num_mat			*V, *D, *K, *U;		// p x p
	Num_mat			*K1;				// ncomp x p	### chk
	Num_mat			*X1;				// ncomp x n
	Num_mat			*tmp_mat;			// n x ncomp
	Num_mat			*WX, *gWX;			// n x ncomp
	Num_mat			*W;					// ncomp x ncomp
	Num_mat			*V1, *V2;			// ncomp x ncomp
	Num_mat			*W1;				// ncomp x ncomp
	double			tmp;
	
	BOOL			row_norm = NO;
//	float			alpha = 1.0;

	n = X->nr;	// n row (observations)
	p = X->nc;	// p column (pixels)

	saveAsKOImage(X, @"IMG_X");


//	w.init <- matrix(rnorm(n.comp^2, n.comp, n.comp)
	W = Num_new_mat(ncomp, ncomp);
	for (i = 0; i < ncomp * ncomp; i++) {
		W->data[i] = Num_nrml(0.0, 1.0);
	}

// === if (method == "R"_ =====
//	X <- scale(X, scale = FALSE)
	Num_col_center(X);		// X <- scale(X, scale = FALSE)
	saveAsKOImage(X, @"IMG_Xcenter");

//	X <- if (row.norm) t(scale(X, scale = row.snorm)) else t(X)
	Xt = Num_trans(X);
	if (row_norm) {
		Num_row_center(X);
		Num_row_norm(X);
		Num_trans_ip(Xt, X);
	}
//	V <- X %*% t(X) / n
	V = Num_new_mat(p, p);
	Num_mmul(V, Xt, X);
	Num_mmul_const(V, 1.0/n);

// s <- La.svd(V)
	sres = Num_svd(V);

// D <- diag(c(1/sqrt(s$d)))
	D = Num_new_mat(p, p);
	for (i = 0; i < p; i++) {
		D->data[i * p + i] = 1.0 / sqrt(sres->s->data[i]);
	}
// K <- D %*% t(s$u)
	K = Num_new_mat(p, p);
	U = Num_trans(sres->U);
	Num_mmul(K, D, U);
// K <- matrix(K[1:n.comp, ], n.comp, p)
// ncomp != p case is not tested yet
	K1 = Num_new_mat(ncomp, p);
	for (i = 0; i < ncomp; i++) {
		for (j = 0; j < p; j++) {
			K1->data[i * ncomp + j] = K->data[i * p + j];
		}
	}
//saveAsKOImage(K, @"IMG_K");
//saveAsKOImage(K1, @"IMG_K1");

// X1 <- K %x% X
	X1 = Num_new_mat(ncomp, n);
	Num_mmul(X1, K1, Xt);

saveAsKOImage(X1, @"IMG_X1");

// ica.R.par()
// W <- w.init
// sW <- La.svd(W)
// W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W
	Num_orth_mat(W);

// if (fun == "exp")

// fixed point iteration
	WX  = Num_new_mat(ncomp, n);
	gWX = Num_new_mat(ncomp, n);
	V1 = Num_new_mat(ncomp, ncomp);
	V2 = Num_new_mat(ncomp, ncomp);
	W1 = Num_new_mat(ncomp, ncomp);
	tmp_mat = Num_new_mat(ncomp, ncomp);
	ncn = ncomp * n;
	ncnc = ncomp * ncomp;
	maxiter = 100;

// while (lim(it] > tol && it < maxit) {
//	float	alpha = 1.0;
	for (iter = 0; iter < maxiter; iter++) {
		// wx <- W %*% X
		Num_mmul(WX, W, X1);
saveAsKOImage(WX, @"IMG_WX");

		// gwx <- wx * exp(-wx^2)/2)
		for (i = 0; i < ncn; i++) {
			tmp = WX->data[i];
			// exp (faster)
			gWX->data[i] = tmp * exp(-(tmp *tmp)/2);				
			// logcosh
		//	gWX->data[i] = tanh(alpha * tmp);
		}
		// V1 = gWX * Xt / p
		Num_mtmul(V1, gWX, NO, X1, YES);
		for (i = 0; i < ncnc; i++) {
			V1->data[i] /= p;
		}
saveAsKOImage(V1, @"IMG_V1");
		// g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
		for (i = 0; i < ncn; i++) {
			tmp = WX->data[i];
			// exp (faster)
			gWX->data[i] = (1 - tmp*tmp) * exp(-(tmp*tmp)/2);
			// logcosh
		//	gWX->data[i] = alpha * (1 - (tanh(alpha * tmp)) * (tanh(alpha * tmp)));
		}

		// V2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
		Num_clear_mat(tmp_mat);
		for (i = 0; i < ncomp; i++) {	// row
			tmp = 0;
			for (j = 0; j < n; j++) {	// col
				tmp += gWX->data[i + j * ncomp];
			}
			tmp /= n;
			tmp_mat->data[i * ncomp + i] = tmp;
		}
saveAsKOImage(tmp_mat, @"IMG_diag");

		Num_mmul(V2, tmp_mat, W);
saveAsKOImage(V2, @"IMG_V2");
		// W1 <- V1 - V2
		for (i = 0; i < ncnc; i++) {
			W1->data[i] = V1->data[i] - V2->data[i];	// need to preserve V1 ###
		}
saveAsKOImage(V1, @"IMG_V12");
		// make orthogonal
		// sW1 <- La.svd(W1)
		// W1 <- sW1$u %*% Diag*1/sW1$d) %*% t(sW1%u) %*% W1
		Num_orth_mat(W1);	// #### chk
saveAsKOImage(V1, @"IMG_V12o");
				
		// calc tol (V1 * Wt)
		Num_mtmul(tmp_mat, W1, NO, W, YES);
saveAsKOImage(tmp_mat, @"IMG_tol");

		tmp = 0;
		for (i = 0; i < ncomp; i++) {
			if (fabs(1 - fabs(tmp_mat->data[i * ncomp + i])) > tmp) {
				tmp = fabs(1 - fabs(tmp_mat->data[i * ncomp + i]));
			}
		}
		tol = tmp;
		printf("%d %10.3e\n", iter, tol);

		Num_copy_mat(W, W1);

		if (tol < 1.0e-6) break;
	}	// iteration

// Wout
	saveAsKOImage(W, @"IMG_Wout");
// WX
	Num_mmul(WX, W, X1);
//copy_mat_to_img(WX, img);
//[img saveAsKOImage:@"IMG_out"];
saveAsKOImage(WX, @"IMG_out");

// sort/set U, s, Vt
	ires = (Num_svd_result *)malloc(sizeof(Num_ica_result));
	ires = (Num_svd_result *)malloc(sizeof(Num_svd_result));
	ires->U  = Num_new_mat(ncomp, ncomp);	// W
	ires->Vt = Num_new_mat(n, ncomp);		// WX
	ires->s  = Num_new_vec(ncomp);

	Num_copy_mat(ires->U, W);
	Num_trans_ip(ires->Vt, WX);

// free mem (chk ###)
	Num_free_mat(Xt);
	Num_free_mat(W);
	Num_free_mat(W1);
	Num_free_mat(WX);
	Num_free_mat(gWX);
	Num_free_mat(V);
	Num_free_mat(V1);
	Num_free_mat(V2);
	Num_free_mat(tmp_mat);

	return ires;
}

//========================
// 2nd attempt
Num_ica_result *
Num_ica(Num_mat *X, int ncomp)
{
	Num_ica_result	*ires;	
	Num_svd_result	*sres;
	int				n, p, ncn, ncnc;
	int				i, j;
	int				iter;
	double			tol;

	Num_mat			*Xt;				// p x n
	Num_mat			*V, *D, *U;			// p x p
	Num_mat			*K;					// p x p
	Num_mat			*K1;				// ncomp x p
	Num_mat			*X1;				// ncomp x n
	Num_mat			*tmp_mat;			// n x ncomp
	Num_mat			*WX, *gWX;			// n x ncomp
	Num_mat			*W;					// ncomp x ncomp
	Num_mat			*V1, *V2;			// ncomp x ncomp
	Num_mat			*W1;				// ncomp x ncomp
	double			tmp;
	
	int				maxiter = 2000;
	BOOL			row_norm = NO;
	BOOL			covariance = NO;
	float			alpha = 1.0;

	n = X->nr;	// n row (observations)
	p = X->nc;	// p column (pixels)

// === if (method == "R"_ =====
//	X <- scale(X, scale = FALSE)
	Num_col_center(X);		// X <- scale(X, scale = FALSE)

//	X <- if (row.norm) t(scale(X, scale = row.snorm)) else t(X)
	Xt = Num_trans(X);
	if (row_norm) {
	//	Num_row_center(X);
		Num_row_norm(X);
		Num_trans_ip(Xt, X);
	}

	if (covariance) {	// orthogonalize X using cov matrix
	//	V <- X %*% t(X) / n
		V = Num_new_mat(p, p);
		Num_mmul(V, Xt, X);
		Num_mmul_const(V, 1.0/n);

	// s <- La.svd(V)
		sres = Num_svd(V);

	// D <- diag(c(1/sqrt(s$d)))
		D = Num_new_mat(p, p);
		for (i = 0; i < p; i++) {
			D->data[i * p + i] = 1.0 / sqrt(sres->s->data[i]);
		}
	// K <- D %*% t(s$u)
		K = Num_new_mat(p, p);
		U = Num_trans(sres->U);
		Num_mmul(K, D, U);
	// K <- matrix(K[1:n.comp, ], n.comp, p)
	// ncomp != p case is not tested yet
		K1 = Num_new_mat(ncomp, p);
		for (i = 0; i < ncomp; i++) {
			for (j = 0; j < p; j++) {
				K1->data[i * ncomp + j] = K->data[i * p + j];
			}
		}
		Num_free_mat(K);
	// X1 <- K %x% X
		X1 = Num_new_mat(ncomp, n);
		Num_mmul(X1, K1, Xt);
//		saveAsKOImage(X1, @"IMG_X1");
	} else {	// ========== use direct SVD to X ==========
		sres = Num_svd(X);

		K1 = Num_new_mat(ncomp, p);	// ok
	//	K = Num_new_mat(n, p);	// first, chk full rank version
	//	Num_copy_mat(K, sres->Vt);
	// use top ncomp entries
		for (i = 0; i < ncomp; i++) {
			for (j = 0; j < p; j++) {
				K1->data[j * ncomp + i] = sres->Vt->data[j * n + i] / sres->s->data[i];
			}
		}
		X1 = Num_new_mat(ncomp, n);
		Num_mmul(X1, K1, Xt);

// scale X1 (to have variance of 1.0)
		alpha = 1.0;
		tmp = 0;
		for (i = 0; i < ncomp * n; i++) {
			tmp += X1->data[i] * X1->data[i];
		}
		tmp /= (ncomp * n);
		tmp = sqrt(tmp);
		tmp = 50.0 / tmp;	// sd = 10.0 ???
		for (i = 0; i < ncomp * n; i++) {
			X1->data[i] *= tmp;
		}

		saveAsKOImage(X1, @"IMG_X1");
		saveAsKOImage(K1, @"IMG_K1");
	}

// if (fun == "exp")

//	w.init <- matrix(rnorm(n.comp^2, n.comp, n.comp)
	W = Num_new_mat(ncomp, ncomp);
//	for (i = 0; i < ncomp * ncomp; i++) {
//		W->data[i] = Num_nrml(0.0, 1.0);
//	}
	// or
	Num_unit_mat(W);
//	for (i = 0; i < ncomp * ncomp; i++) {
//		W->data[i] += Num_nrml(0.0, 1.0) * 0.2;
//	}

// ica.R.par()
// W <- w.init
// sW <- La.svd(W)
// W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W
	Num_orth_mat(W);


// fixed point iteration
	WX  = Num_new_mat(ncomp, n);
	gWX = Num_new_mat(ncomp, n);
	V1 = Num_new_mat(ncomp, ncomp);
	V2 = Num_new_mat(ncomp, ncomp);
	W1 = Num_new_mat(ncomp, ncomp);
	tmp_mat = Num_new_mat(ncomp, ncomp);
	ncn = ncomp * n;
	ncnc = ncomp * ncomp;

// while (lim(it] > tol && it < maxit) {
	for (iter = 0; iter < maxiter; iter++) {
		// wx <- W %*% X
		Num_mmul(WX, W, X1);

		// gwx <- wx * exp(-wx^2)/2)
		for (i = 0; i < ncn; i++) {
			tmp = WX->data[i];
			// exp (faster)
		//	gWX->data[i] = tmp * exp(-(tmp *tmp)/2);				
			// logcosh
			gWX->data[i] = tanh(alpha * tmp);
		}
		// V1 = gWX * Xt / p
		Num_mtmul(V1, gWX, NO, X1, YES);
		for (i = 0; i < ncnc; i++) {
			V1->data[i] /= p;
		}
		// g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
		for (i = 0; i < ncn; i++) {
			tmp = WX->data[i];
			// exp (faster)
		//	gWX->data[i] = (1 - tmp*tmp) * exp(-(tmp*tmp)/2);
			// logcosh
			tmp = gWX->data[i];
			gWX->data[i] = alpha * (1 - (tmp * tmp));
		}

		// V2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
		Num_clear_mat(tmp_mat);
		for (i = 0; i < ncomp; i++) {	// row
			tmp = 0;
			for (j = 0; j < n; j++) {	// col
				tmp += gWX->data[i + j * ncomp];
			}
			tmp /= n;
			tmp_mat->data[i * ncomp + i] = tmp;
		}
		Num_mmul(V2, tmp_mat, W);
		// W1 <- V1 - V2
		for (i = 0; i < ncnc; i++) {
			W1->data[i] = V1->data[i] - V2->data[i];	// need to preserve V1 ###
		}
		// make orthogonal
		// sW1 <- La.svd(W1)
		// W1 <- sW1$u %*% Diag*1/sW1$d) %*% t(sW1%u) %*% W1
		Num_orth_mat(W1);	// #### chk
				
		// calc tol (V1 * Wt)
		Num_mtmul(tmp_mat, W1, NO, W, YES);

		tmp = 0;
		for (i = 0; i < ncomp; i++) {
			if (fabs(1 - fabs(tmp_mat->data[i * ncomp + i])) > tmp) {
				tmp = fabs(1 - fabs(tmp_mat->data[i * ncomp + i]));
			}
		}
		tol = tmp;
printf("%d %10.3e\n", iter, tol);

		Num_copy_mat(W, W1);

		if (tol < 1.0e-6) break;	// 1.0e-7 doesn't work
	}	// iteration

// make U, s, Vt equivalent
	ires = (Num_ica_result *)malloc(sizeof(Num_ica_result));
	ires->W = Num_new_mat(ncomp, ncomp);
	ires->WX = Num_new_mat(ncomp, n);
	ires->WK = Num_new_mat(ncomp, p);

	// U
	Num_mmul(ires->WX, W, X1);

	// scale K1 (back to Vt)
	for (i = 0; i < ncomp; i++) {
		for (j = 0; j < p; j++) {
			K1->data[j * ncomp + i] *= sres->s->data[i];
		}
	}

	Num_mmul(ires->WK, W, K1);

	// W
	Num_copy_mat(ires->W, W);

	// sort U & Vt


// free mem (chk ###)
	Num_free_mat(Xt);
	Num_free_mat(W);
	Num_free_mat(W1);
	Num_free_mat(WX);
	Num_free_mat(gWX);
//	Num_free_mat(V);
	Num_free_mat(V1);
	Num_free_mat(V2);
	Num_free_mat(tmp_mat);

	return ires;
}

//========================

void
Num_free_svd_result(Num_svd_result *r)
{
	if (r) {
		Num_free_mat(r->U);
		Num_free_vec(r->s);
		Num_free_mat(r->Vt);
		free(r);
	}
}

void
Num_free_ica_result(Num_ica_result *r)
{
	if (r) {
		Num_free_mat(r->W);
		Num_free_mat(r->WK);
		Num_free_mat(r->WX);
		free(r);
	}
}

void
Num_free_evd_result(Num_evd_result *r)
{
	if (r) {
		Num_free_mat(r->Evec);
		Num_free_vec(r->eval);
		free(r);
	}
}

//========= ODE ==========
// single var first (then vector func)
void
Num_rk4(float *x, float *y, float (^deriv)(float x, float y), float step)
{
    float       k1, k2, k3, k4;
    float       step2 = step / 2;

    k1 = deriv(*x, *y);
    k2 = deriv(*x + step2, *y + step2 * k1);
    k3 = deriv(*x + step2, *y + step2 * k2);
    k4 = deriv(*x + step,  *y + step  * k3);

    *x += step;
    *y += step / 6.0 * (k1 + k2 * 2 + k3 * 2 + k4);
}


//======== dbg =============
void
dump_vec(Num_vec *v)
{
    int     i;
	int		n = v->n < 10 ? v->n : 10;

	switch (v->type) {
	case NUM_REAL :
	default :
		printf("=== real vector[%d] ====\n", v->n);
		for (i = 0; i < v->n; i++) {
			printf(" %3.2e\n", v->data[i]);
		}
		break;
	case NUM_COMPLEX :
		printf("=== complex vector[%d] ====\n", v->n);
		for (i = 0; i < n; i++) {
			printf(" (%3.2e %3.2ei)\n", v->data[i*2], v->data[i*2 + 1]);
		}
		break;
	}
}

// dump data as stored in mem
void
dump_mat(Num_mat *m)
{
    int     c, r, ix;
	int		nc = m->nc < 10 ? m->nc : 10;
	int		nr = m->nr < 10 ? m->nr : 10;

	switch (m->type) {
	case NUM_REAL :
	default :
		printf("=== real matrix[%d][%d] ====\n", m->nr, m->nc);
		if (m->order == NUM_COL_MAJOR) {
			printf("(stored in column major order\n");
			for (r = 0; r < nr; r++) {
				for (c = 0; c < nc; c++) {
				//	ix = Num_col_m_index(m, c, r);
					ix = r * m->nc + c;
					printf("%5.2f ", m->data[ix]);
				}
				printf("\n");
			}
		} else {
			printf("(stored in row major order\n");
			printf("(not supported\n");
			ix = 0;
			for (r = 0; r < nr; r++) {
				for (c = 0; c < nc; c++) {
					ix = c * m->nr + r;
					printf("%5.2f ", m->data[ix]);
					ix++;
				}
				printf("\n");
			}
		}
		break;
	case NUM_COMPLEX :
		printf("=== complex matrix[%d][%d] ====\n", m->nr, m->nc);
		for (r = 0; r < nr; r++) {
			for (c = 0; c < nc; c++) {
				ix = (r * m->nr + c) * 2;
				printf("(%5.2f, %5.2fi)", m->data[ix], m->data[ix + 1]);
			}
			printf("\n");
		}
		break;
	}
}

int
diff_mat(Num_mat *a, Num_mat *b)
{
	int		i, d = 0;
	if (a->nr != b->nr) return 1;
	if (a->nc != b->nc) return 1;
	for (i = 0; i < a->nr * a->nc; i++) {
		if (a->data[i] != b->data[i]) {
			d = 1;
			break;
		}
	}
	return d;
}

void
saveAsKOImage(Num_mat *m, NSString *path)
{
	RecImage *im = Num_m_to_im(m);
	[im saveAsKOImage:path];
}
