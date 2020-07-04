//
// NumLinear.m
//
//  == plans ==
//	make object version
//		svd ok, implement ica

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

- (NumMatrix *)diagMatrix
{
	int			i, n = [self length];
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:n nCol:n];
	float		*p = [m data];
	float		*q = [self data];

	for (i = 0; i < n; i++) {
		p[i * n + i] = q[i];
	}

	return m;
}

@end


@implementation NumMatrix

+ (NumMatrix *)matrixOfType:(int)tp nRow:(int)nr nCol:(int)nc
{
	NumMatrix	*m;

	m = [NumMatrix alloc];
	return [m initWithType:tp nRow:nr nCol:nc];
}

+ (NumMatrix *)matrixWithMatrix:(NumMatrix *)m
{
	return [NumMatrix matrixOfType:[m type] nRow:[m nRow] nCol:[m nCol]];
}

- (NumMatrix *)initWithType:(int)tp nRow:(int)nr nCol:(int)nc
{
	self = [super init];
	if (!self) return nil;

	len = nr * nc;
	if (tp == NUM_COMPLEX) {
		len *= 2;
	}
	data = [NSMutableData dataWithLength:sizeof(float) * len];
	type = tp;
	nRow = nr;
	nCol = nc;

	return self;
}

+ (NumMatrix *)matrixWithImage:(RecImage *)img
{
	NumMatrix	*m;
	int			m_type;

	switch ([img type]) {
	case RECIMAGE_REAL :
		m_type = NUM_REAL;
		break;
	case RECIMAGE_COMPLEX :
		m_type = NUM_COMPLEX;
		break;
	}
	m = [NumMatrix matrixOfType:m_type nRow:[img yDim] nCol:[img xDim]];
	[m copyImage:img];
	return m;
}

+ (NumMatrix *)unitMatrixOfDim:(int)n
{
	NumMatrix	*m = [NumMatrix matrixOfType:NUM_REAL nRow:n nCol:n];
	int			i;
	float		*p = [m data];

	for (i = 0; i < n; i++) {
		p[i * n + i] = 1;
	}

	return m;
}

+ (NumMatrix *)matrixWithNumMat:(Num_mat *)mt
{
    NumMatrix   *m;
    float       *p, *q;
    int         i, n;

    m = [NumMatrix matrixOfType:mt->type nRow:mt->nr nCol:mt->nc];
    n = [m len];
    if (mt->type == NUM_COMPLEX) {
        n *= 2;
    }
    p = [m data];
    q = mt->data;
    for (i = 0; i < n; i++) {
        p[i] = q[i];
    }
    return m;
}

- (float *)data
{
	return [data mutableBytes];
}

- (float *)real
{
	return [self data];
}

- (float *)imag
{
	return [self data] + len;
}

- (int)type
{
	return type;
}

- (int)nRow
{
	return nRow;
}

- (int)nCol
{
	return nCol;
}

- (int)ld
{
	return nRow;
}

- (int)len
{
	return len;
}

- (void)clear
{
	int		i;
	float	*p = [self data];

	for (i = 0; i < len; i++) {
		p[i] = 0;
	}
}

- (void)normal
{
	int		i;
	float	*p = [self data];

	for (i = 0; i < len; i++) {
		p[i] = Num_nrml(0, 1);
	}
}

- (NumMatrix *)copy
{
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:nRow nCol:nCol];
	int			i;
	float		*p1 = [self data];
	float		*p2 = [m data];

	for (i = 0; i < len; i++) {
		p2[i] = p1[i];
	}
	return m;
}

// self is not altered
- (NumMatrix *)trans
{
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:nCol nCol:nRow];
	float		*p1 = [self data];
	float		*p2 = [m data];
	int			n1 = nRow * nCol;
	int			n2 = [m nRow] * [m nCol];

    vDSP_mtrans(p1, 1, p2, 1, nRow, nCol);
	if (type == NUM_COMPLEX) {
		p1 += n1;
		p2 += n2;
		vDSP_mtrans(p1, 1, p2, 1, nRow, nCol);
	}
	return m;
}

- (void)copyMatrix:(NumMatrix *)mt		// copy data, size can be different
{
	float	*p, *q;
	int		i, j;
	int		m, n;
	int		src_nr = [mt nRow];

	m = MIN(nRow, [mt nRow]);
	n = MIN(nCol, [mt nCol]);

	p = [mt data];
	q = [self data];
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			q[i * nRow + j] = p[i * src_nr + j];
		}
	}
}

- (NumMatrix *)addMat:(NumMatrix *)m
{
	NumMatrix	*res;
	int			i;
	float		*q, *p1, *p2;

	res = [NumMatrix matrixWithMatrix:self];
	q = [res data];
	p1 = [self data];
	p2 = [m data];
	for (i = 0; i < len; i++) {
		q[i] = p1[i] + p2[i];
	}
	return res;
}

- (NumMatrix *)subMat:(NumMatrix *)m
{
	NumMatrix	*res;
	int			i;
	float		*q, *p1, *p2;

	res = [NumMatrix matrixWithMatrix:self];
	q = [res data];
	p1 = [self data];
	p2 = [m data];
	for (i = 0; i < len; i++) {
		q[i] = p1[i] - p2[i];
	}
	return res;
}

- (NumMatrix *)addConst:(float)a	// sub not necessary (just make a negative)
{
	NumMatrix	*res;
	int			i;
	float		*q, *p1;

	res = [NumMatrix matrixWithMatrix:self];
	q = [res data];
	p1 = [self data];
	for (i = 0; i < len; i++) {
		q[i] = p1[i] + a;
	}

	return res;
}

// C = AB
- (NumMatrix *)multByMat:(NumMatrix *)m
{
	NumMatrix *res;
	if (nCol != [m nRow]) {
        Num_error("mmul dim mismatch");
	}
	res = [NumMatrix matrixOfType:type nRow:nRow nCol:[m nCol]];
	cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
		nRow, [m nCol], nCol, 1.0, [self data], nRow, [m data], [m nRow], 0, [res data], [res nRow]);

	return res;
}

- (NumMatrix *)multByConst:(float)a
{
	int			i;
	float		*p, *q;
	NumMatrix	*res = [NumMatrix matrixWithMatrix:self];

	q = [res data];
	p = [self data];
	for (i = 0; i < len; i++) {
		q[i] = p[i] * a;
	}

	return res;
}

- (NumMatrix *)colVect:(int)ix
{
	NumMatrix	*res;
	int			i;
	float		*p, *q;

	res = [NumMatrix matrixOfType:type nRow:nRow nCol:1];
	q = [res data];
	p = [self data];
	for (i = 0; i < nRow; i++) {
		q[i] = p[ix * nRow + i];
	}

	return res;
}

- (NumMatrix *)rowVect:(int)ix
{
	NumMatrix	*res;
	int			i;
	float		*p, *q;

	res = [NumMatrix matrixOfType:type nRow:1 nCol:nCol];
	q = [res data];
	p = [self data];
	for (i = 0; i < nCol; i++) {
		q[i] = p[i * nRow + ix];
	}

	return res;
}

- (NumMatrix *)rowMean
{
	int			i, j;
	float		mn;
	float		*p, *q;
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:1 nCol:nCol];

	q = [m data];
	p = [self data];
	for (i = 0; i < nCol; i++) {
		mn = 0;
		for (j = 0; j < nRow; j++) {
			mn += p[i=j * nCol + i];
		}
		mn /= nRow;
		q[i] = mn;
	}
	return m;
}

- (NumMatrix *)colMean
{
	int			i, j;
	float		mn;
	float		*p, *q;
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:nRow nCol:1];

	q = [m data];
	p = [self data];
	for (i = 0; i < nRow; i++) {
		mn = 0;
		for (j = 0; j < nCol; j++) {
        //    mn += p[j * nCol + i];
            mn += p[i * nRow + j];
		}
		mn /= nCol;
		q[i] = mn;
	}
	return m;
}

- (NumMatrix *)colCenter
{
	int			i, j;
	float		mn;
	float		*p, *q;
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:nRow nCol:nCol];

	q = [m data];
	p = [self data];
	for (i = 0; i < nCol; i++) {
		mn = 0;
		for (j = 0; j < nRow; j++) {
			mn += p[i * nRow + j];
		}
		mn /= nRow;
		for (j = 0; j < nRow; j++) {
			q[i * nRow + j] = p[i * nRow + j] - mn;
		}
	}
	return m;
}

- (NumMatrix *)rowCenter
{
	int			i, j;
	float		mn;
	float		*p, *q;
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:nRow nCol:nCol];

	q = [m data];
	p = [self data];
	for (i = 0; i < nRow; i++) {
		mn = 0;
		for (j = 0; j < nCol; j++) {
			mn += p[j * nRow + i];
		}
		mn /= nCol;
		for (j = 0; j < nCol; j++) {
			q[j * nRow + i] = p[j * nRow + i] - mn;
		}
	}
	return m;
}

- (NumMatrix *)colSD		// ## vec -> diag matrix
{
	int			i, j;
	float		sd;
	float		*p, *q;
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:nCol nCol:nCol];

	q = [m data];
	p = [self data];
	for (i = 0; i < nCol; i++) {
		sd = 0;
		for (j = 0; j < nRow; j++) {
			sd +=  p[i * nRow + j] * p[i * nRow + j];
		}
		q[i * nRow + i] = sqrt(sd / nRow);
	}
	return m;
	
}

- (NumMatrix *)rowSD		// ## vec -> diag matrix
{
	int			i, j;
	float		sd;
	float		*p, *q;
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:nRow nCol:nRow];

	q = [m data];
	p = [self data];
	for (i = 0; i < nRow; i++) {
		sd = 0;
		for (j = 0; j < nCol; j++) {
			sd += p[j * nRow + i] * p[j * nRow + i];
		}
		q[i * nRow + i] = sqrt(sd / nCol);
	}
	return m;
}

// 0-mean columns are expected
- (NumMatrix *)colNorm
{
	int			i, j;
	float		sd;
	float		*p, *q;
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:nRow nCol:nCol];

	q = [m data];
	p = [self data];
	for (i = 0; i < nCol; i++) {
		sd = 0;
		for (j = 0; j < nRow; j++) {
			sd +=  p[i * nRow + j] * p[i * nRow + j];
		}
		sd = sqrt(sd / nRow);
		for (j = 0; j < nRow; j++) {
			q[i * nRow + j] = p[i * nRow + j] / sd;
		}
	}
	return m;
}

// 0-mean rows are expected
- (NumMatrix *)rowNorm
{
	int			i, j;
	float		sd;
	float		*p, *q;
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:nRow nCol:nCol];

	q = [m data];
	p = [self data];
	for (i = 0; i < nRow; i++) {
		sd = 0;
		for (j = 0; j < nCol; j++) {
			sd += p[j * nRow + i] * p[j * nRow + i];
		}
		sd = sqrt(sd / nCol);
		for (j = 0; j < nCol; j++) {
			q[j * nRow + i] = p[j * nRow + i] / sd;
		}
	}
	return m;
}

- (NumMatrix *)orthog
{
	NSDictionary	*res;
	NumMatrix		*m;
	res = [self svd];
	m = [[res objectForKey:@"U"] multByMat:[res objectForKey:@"Vt"]];
	return m;
}

- (NumMatrix *)diagMatrix			// make diag mat with single col or single row mat
{
	NumMatrix	*m = [NumMatrix matrixOfType:type nRow:len nCol:len];
	int			i;
	float		*p = [m data];
	float		*q = [self data];

	for (i = 0; i < len; i++) {
		p[i * len + i] = q[i];
	}
	return m;
}

// solve Ax = b
- (NumMatrix *)solveLinear:(NumMatrix *)B
{
	int			m, n, nrhs;
	int			lda, ldb, info;
	int			*ipiv;
	float		*work;
	int			lwork;
	NumMatrix	*AA, *BB;	// copy mat first. both AA & BB are overwritten by LAPACK
	NumMatrix	*X;			// result

	m = nRow;
	n = nCol;
	lda = m;
	ldb = MAX(m, n);
	nrhs = [B nCol];

	// output
	X = [NumMatrix matrixOfType:type nRow:n nCol:nrhs];

	// copy
	AA = [self copy];
	BB = [B copy];	// returns copy of object
	if (m == n) {
		ipiv = (int *)malloc(sizeof(int) * n);
		sgesv_(&n, &nrhs, [AA data], &lda, ipiv, [BB data], &ldb, &info);
		// copy output to B
		[X copyMatrix:BB];
		free(ipiv);
	} else {
		work = (float *)malloc(sizeof(float) * 1);
		lwork = -1;
		sgels_("N", &m, &n, &nrhs, [AA data], &lda, [BB data], &ldb, work, &lwork, &info);
		lwork = work[0];
		free(work);
		work = (float *)malloc(sizeof(float) * lwork);
		sgels_("N", &m, &n, &nrhs, [AA data], &lda, [BB data], &ldb, work, &lwork, &info);
		free(work);
		// copy output to B
		[X copyMatrix:BB];
	}

	return X;
}

- (BOOL)empty
{
	int		i, n = nCol * nRow;
	float	*p = [self data];

	for (i = 0; i < n; i++) {
		if (p[i] >= 0) {
			return NO;
		}
	}
	return YES;
}

- (float)maxValAt:(int *)ix
{
	int		i, n = nCol * nRow;
	float	*p = [self data];
	float	mx = 0;

	for (i = 0; i < n; i++) {
		if (p[i] > mx) {
			mx = p[i];
			*ix = i;
		}
	}
	return mx;
}

- (NumMatrix *)solvePartial:(NumMatrix *)b list:(float *)flg
{
	NumMatrix	*AP, *x;
	float		*p, *q;
	int			i, j, ix;
	int			nc;

	nc = 0;
	for (i = 0; i < nCol; i++) {
		if (flg[i] > 0) {
			nc++;
		}
	}
	AP = [NumMatrix matrixOfType:[self type] nRow:nRow nCol:nc];
	q = [AP data];
	p = [self data];
	for (i = ix = 0; i < nCol; i++) {
		if (flg[i] == 0) continue;
		
		for (j = 0; j < nRow; j++) {
			q[ix * nRow + j] = p[i * nRow + j];
		}
		ix++;
	}
	x = [AP solveLinear:b];

	return x;
}

- (float)minVal
{
	int		i;
	float	mn;
	float	*p = [self data];

	for (i = 0; i < len; i++) {
		if (i == 0) {
			mn = p[i];
		} else {
			if (mn > p[i]) {
				mn = p[i];
			}
		}
	}
	return mn;
}

- (NumMatrix *)selectCol:(int)col
{
	NumMatrix	*v;
	float		*p, *q;
	int			i;

	v = [NumMatrix matrixOfType:[self type] nRow:nRow nCol:1];
	q = [v data];
	p = [self data];
	for (i = 0; i < nRow; i++) {
		q[i] = p[col * nRow + i];

	}
	return v;
}

- (void)copyVec:(NumMatrix *)v atCol:(int)col
{
	float		*p, *q;
	int			i;

	q = [v data];
	p = [self data];
	for (i = 0; i < nRow; i++) {
		p[col * nRow + i] = q[i];
	}
}

// mostly ok, exit condition of main loop ###
- (NumMatrix *)solveLinearNN:(NumMatrix *)B	// NNLS
{
	int			m, n, nrhs;
	NumMatrix	*AA;					// copy mat first. both AA & BB are overwritten by LAPACK
	NumMatrix	*X;						// result
	NumMatrix	*w, *b, *x, *s;			// vectors
	NumMatrix	*p;						// passive list
	float		*pp;
	float		tol = 1e-7;
	float		*xx, *ss, alpha, a;
	int			i, k, ix;
	int			outer, inner;
	BOOL		all;

	m = nRow;
	n = nCol;
	nrhs = [B nCol];

	// output
	X = [NumMatrix matrixOfType:type nRow:n nCol:nrhs];

	// copy
	AA = [self copy];

	// alloc variables
	x = [NumMatrix matrixOfType:NUM_REAL nRow:n nCol:1];	// solution
	xx = [x data];
	
	p = [NumMatrix matrixOfType:NUM_REAL nRow:n nCol:1];	// passive list
	pp = [p data];

	// set loop
	for (k = 0; k < nrhs; k++) {
		// init
		for (i = 0; i < n; i++) {
			pp[i] = 0;
		}
		b = [B colVect:k];	// i-th column
		[x clear];
		w = [[AA trans] multByMat:b];	// w = At b

		// main loop
		for (outer = 0; outer < n * 2; outer++) {
			//printf("outer = %d, max w = %f\n", outer, [w maxValAt:&ix]);
			if ([w maxValAt:&ix] <= tol) break;
			all = YES;
			for (i = 0; i < n; i++) {
				if (pp[i] == 0) {
					all = NO;
					break;
				}
			}
			if (all) break;

			// add ix to p
			pp[ix] = 1;
			s = [AA solvePartial:b list:pp];
			ss = [s data];

			for (inner = 0; inner < n * 2; inner++) {
				if ([s minVal] >= 0) break;
				alpha = 1.0;	// maximum possible value
				for (i = ix = 0; i < n; i++) {
					if (pp[i] == 0) continue;
					if (ss[ix] < 0) {
						if (xx[i] == ss[ix]) {
							a = 1.0;
						} else {
							a = xx[i] / (xx[i] - ss[ix]);
						}
						if (alpha > a) {	// ### ??? alpha selection is wrong...
							alpha = a;
						}
					}
					ix++;
				}
				for (i = ix = 0; i < n; i++) {
					if (pp[i] == 0) continue;
					xx[i] += (ss[ix] - xx[i]) * alpha;
					if (fabs(xx[i]) < tol) {
						pp[i] = 0;
					}
					ix++;
				}
				s = [AA solvePartial:b list:pp];
				ss = [s data];
			}
	//	printf("outer = %d, inner = %d\n", outer, inner);
			for (i = ix = 0; i < n; i++) {
				if (pp[i] == 0) {
					xx[i] = 0;
					continue;
				}
				xx[i] = ss[ix];
				ix++;
			}
			w = [[AA trans] multByMat:[b subMat:[AA multByMat:x]]];	// w = At (b - Ax)
		}
		// copy x to X
		[X copyVec:x atCol:k];
	}

	return X;
}

- (NSDictionary *)svd	// returns dict with entry: "U", "S"(diag mat), "Vt"
{
	NumMatrix		*U, *S, *Vt;
	NumVector		*sv;
	float			*rwork, *cwork, *adata;
	float			*p;
	int				lwork, rworklen, info;
	int				i, nr, nc, n, nrc, ld;

	nr = nRow;
	nc = nCol;
	nrc = nr * nc;
	ld = nr;
	n = MIN(nr, nc);

// alloc result matrices
//	res = (Num_svd_result *)malloc(sizeof(Num_svd_result));
	U  = [NumMatrix matrixOfType:type nRow:nr nCol:n];
	Vt = [NumMatrix matrixOfType:type nRow:n nCol:nc];
	sv  = [NumVector vectorOfType:type length:n];

	if (type == NUM_REAL) {
	// get work size
		lwork = -1;
		rwork = (float *)malloc(sizeof(float) * 1);
		adata = (float *)malloc(sizeof(float) * nrc);
		sgesvd_("S", "S", &nr, &nc, adata, &ld,
			[sv data], [U data], &nr, [Vt data], &n,
			rwork, &lwork, &info);
		lwork = rwork[0];
		free(rwork);

		// printf("sizeof work area %d\n", lwork);
		// copy A
		p = [self data];
		for (i = 0; i < nrc; i++) {
			adata[i] = p[i];
		}
		// actual proc
		rwork = (float *)malloc(sizeof(float) * lwork);
		sgesvd_("S", "S", &nr, &nc, adata, &ld,
			[sv data], [U data], &nr, [Vt data], &n, rwork,
			&lwork, &info);
		if (info != 0) {
			printf("sgesvd info: %d\n", info);
		}
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
			[sv data], (__CLPK_complex *)[U data], &nr, (__CLPK_complex *)[Vt data], &n,
			(__CLPK_complex *)cwork, &lwork, rwork, &info);
		lwork = cwork[0];
		//printf("lwork[0] = %d\n", lwork);
		free(cwork);

		//printf("sizeof work area %d\n", lwork);
		// actual proc
		cwork = (float *)malloc(sizeof(float) * lwork * 2);
		// copy A
		p = [self data];
		for (i = 0; i < nrc * 2; i++) {
			adata[i] = p[i];
		}
		cgesvd_("S", "S", &nr, &nc, (__CLPK_complex *)adata, &ld,
			[sv data], (__CLPK_complex *)[U data], &nr, (__CLPK_complex *)[Vt data], &n,
			(__CLPK_complex *)cwork, &lwork, rwork, &info);
		if (info != 0) {
			printf("cgesvd info: %d\n", info);
		}
		free(rwork);
		free(cwork);
		free(adata);
	}
	// make diag mat
	S = [sv diagMatrix];

	// make dict
	return [NSDictionary dictionaryWithObjectsAndKeys:U, @"U", S, @"S", Vt, @"Vt", nil];
}

- (NSDictionary *)icaForNC:(int)nc	// returns dict with entry: "XW", "Y", "U", "Vt"
{
	NSDictionary	*res;
	int				n, p, r;
	int				i, j;

// nr x nc, p:channels, pixels, n:time, observations, r:rank min(p, n)
	NumMatrix		*A;					// [n,  p]
    NumMatrix       *At;                // [p,  n]
	NumMatrix		*U;					// [p,  r]
    NumMatrix       *Vt;
	NumMatrix		*X;					// [nc, n] white input
    NumMatrix       *K;                 // [nc, p] whitening matrix
    NumMatrix       *Sg;                // [r,  r]
    NumMatrix       *Sg1;               // [nc, r]
// === 7/3
	NumMatrix		*Y;					// [nc, p] chk
	NumMatrix		*W, *W1;			// [nc,nc]
	NumMatrix		*WX;				// [nc, n] chk
    NumMatrix       *gWX;               // [nc, n] chk
	float			*px, *pg;
	NumMatrix		*V1, *V2;			// [nc,nc]
     NumMatrix		*Vtnc;				// [nc, p] chk
	NumMatrix		*tmpMat;
    float           *p1, *p2, *q;
	
	int				iter, maxIter = 20; // 1000  -> should be smaller
	float			lim, tmp, tol = 1.0e-5;

	n = nRow;	    // number of samples
	p = nCol;	    // number of pixels
    r = MIN(n, p);  // rank

	A = [self colCenter];
 [A saveAsKOImage:@"IMG_A_ica"];   
    
    At = [A trans];
 [At saveAsKOImage:@"IMG_At_ica"]; 
   	res = [At svd];
    U = [res objectForKey:@"U"];
//printf("U\n"); [U dump];
[U saveAsKOImage:@"IMG_U_ica"];
    Sg = [res objectForKey:@"S"];
[Sg saveAsKOImage:@"IMG_sg"];
//printf("s\n"); [Sg dump];
Vt = [res objectForKey:@"Vt"];
[Vt saveAsKOImage:@"IMG_Vt_ica"];

if (0) {    // calc from U
// U is [p, p] -> [p, min(p, n)] = r(ank)
//    K = [NumMatrix unitMatrixOfDim:nc]; // [nc, nc] -> K[r, p]
    K = [NumMatrix matrixOfType:NUM_REAL nRow:r nCol:p];
    p1 = [Sg data]; // [r, r]
//    p2 = [[U trans] data];  // [r, p]
    p2 = [U data];  // [r, p]
    q = [K data];   // [r, p]
    for (i = 0; i < p; i++) {  // col
        for (j = 0; j < r; j++) {  // row
            if (p1[j*r + j] > 0.001) {
            //    q[j*nc + i] = sqrt(n) * p2[j*p + i] / p1[j*r + j];
                q[j*p + i] = sqrt(n) * p2[j*p + i] / p1[j*r + j];
            }
        }
    } 
     printf("K\n"); [K dump];   
     [K saveAsKOImage:@"IMG_K"];

    //printf("obj K1\n"); [K dump];
        X = [K multByMat:At];
    //printf("obj X1\n"); [X dump];
    [X saveAsKOImage:@"IMG_X"];
} else {    // calc direct form Vt (now ok)
    Sg1 = [NumMatrix matrixOfType:NUM_REAL nRow:nc nCol:r];
    p1 = [Sg data];
    p2 = [Sg1 data];
    for (i = 0; i < nc; i++) {
    //    p2[i * nc + i] = sqrt(p1[i * r + i]);  // not necessary ???
        p2[i * nc + i] = sqrt(n) * sqrt(p1[i * r + i]);
    }
    [Sg1 saveAsKOImage:@"IMG_sg1"];
    X = [Sg1 multByMat:Vt];
    [X saveAsKOImage:@"IMG_X"];
    // make sub-mat by taking top nc rows of X

}

	W = [NumMatrix unitMatrixOfDim:nc];
//	[W normal]; W = [W orthog];

	for (iter = 0; iter < maxIter; iter++) {
		WX = [W multByMat:X];
[WX saveAsKOImage:@"IMG_WX"];    
[W saveAsKOImage:@"IMG_W"];    
        gWX = [WX copy];
        px = [WX data];
        pg = [gWX data];
        for (i = 0; i < [WX len]; i++) {
            tmp = px[i];
            pg[i] = tmp * exp(-(tmp*tmp)/2);
        }
        // V1 <- gxw %*% t(X)/p
        V1 = [gWX multByMat:[X trans]];
     //   V1 = [V1 multByConst:1.0/p];
        V1 = [V1 multByConst:1.0/nc];   // ???
[V1 saveAsKOImage:@"IMG_V1"];

	// V2 <= Diag(apply(g.wx, 1, FUN = mean) %*% W

	//	tmpMat = [[g2 colMean] diagMatrix];
 
        for (i = 0; i < [WX len]; i++) {
            tmp = px[i];
            pg[i] = (1.0 - tmp*tmp) * exp(-(tmp*tmp)/2);
        }
//[gWX dump];
//[gWX saveAsKOImage:@"IMG_gwx"];
        tmpMat = [[gWX colMean] diagMatrix];
[tmpMat saveAsKOImage:@"IMG_tmp"];
        V2 = [tmpMat multByMat:W];  // ### X
[V2 saveAsKOImage:@"IMG_V2"];
//printf("Obj colMean\n"); [[gWX colMean] dump];
//printf("Obj V2\n"); [V2 dump];

    
		W1 = [V1 subMat:V2];
		W1 = [W1 orthog];
[W1 saveAsKOImage:@"IMG_W1"];

	// lim <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
        tmpMat = [W1 multByMat:[W trans]];  // trans does not change self

// ============ 7/4

		px = [tmpMat data];
		lim = 0;
		for (i = 0; i < nc; i++) {
			tmp = fabs(fabs(px[i * nc + i]) - 1);
			if (tmp > lim) lim = tmp;
		}
printf("%d %e\n", iter, lim);
		if (lim < tol) break;
		W = W1;
        
	}   // iteration

	// calc "Y"
	Vtnc = [NumMatrix matrixOfType:NUM_REAL nRow:nc nCol:p];
	[Vtnc copyMatrix:[res objectForKey:@"Vt"]];
	
	Y = [[W trans] multByMat:Vtnc];

	return [NSDictionary dictionaryWithObjectsAndKeys:
				[WX trans], @"WX", Y, @"Y", W, @"W", nil]; // X, W, WX, Y
}

// RecImage
- (RecImage *)toRecImage
{
	RecImage	*img;
	int			i_type;
	int			i, j;
	float		*p, *q;

	switch (type) {
	case NUM_REAL :
		i_type = RECIMAGE_REAL;
		break;
	case NUM_COMPLEX :
		i_type = RECIMAGE_COMPLEX;
		break;
	}
	img = [RecImage imageOfType:i_type xDim:nCol yDim:nRow];
	q = [img data];
	p = [self data];
	for (i = 0; i < nCol; i++) {
		for (j = 0; j < nRow; j++) {
			q[j * nCol + i] = p[i * nRow + j];
		}
	}
	if (type == NUM_COMPLEX) {
		q = [img imag];
		p = [self imag];
		for (i = 0; i < nCol; i++) {
			for (j = 0; j < nRow; j++) {
				q[j * nCol + i] = p[i * nRow + j];
			}
		}
	}
	return img;
}

- (void)saveAsKOImage:(NSString *)path
{
	RecImage	*img = [self toRecImage];
	[img saveAsKOImage:path];
}

- (void)copyImage:(RecImage *)img
{
	float	*p, *q;
	int		i, j;

	q = [img data];
	p = [self data];
	for (i = 0; i < nCol; i++) {
		for (j = 0; j < nRow; j++) {
			p[i * nRow + j] = q[j * nCol + i];
		}
	}
	if (type == NUM_COMPLEX) {
		q = [img imag];
		p = [self imag];
		for (i = 0; i < nCol; i++) {
			for (j = 0; j < nRow; j++) {
				p[i * nRow + j] = q[j * nCol + i];
			}
		}
	}
}

- (void)dump
{
    int     c, r, ix;
	int		nc = nCol < 10 ? nCol : 10;
	int		nr = nRow < 10 ? nRow : 10;

	switch (type) {
	case NUM_REAL :
	default :
		printf("=== real matrix[%d][%d] ====\n", nRow, nCol);
		for (r = 0; r < nr; r++) {
			for (c = 0; c < nc; c++) {
			//	ix = Num_col_m_index(m, c, r);
				ix = c * nRow + r;
				printf("%6.5f ", [self data][ix]);
			}
			printf("\n");
		}
		break;
	case NUM_COMPLEX :
		printf("=== complex matrix[%d][%d] ====\n", nRow, nCol);
		for (r = 0; r < nr; r++) {
			for (c = 0; c < nc; c++) {
				ix = (r * nRow + c) * 2;
				printf("(%5.2f, %5.2fi)", [self data][ix], [self data][ix + 1]);
			}
			printf("\n");
		}
		break;
	}
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
		for (c = 0; c < m->nc; c++) {	// x, col
			for (r = 0; r < m->nr; r++) {	// y, row
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
// ### need separate func for complex version (need to return complex number)
float
Num_dotpr(Num_vec *v1, Num_vec *v2)
{
    float   sum;

	if (v1->type == NUM_COMPLEX || v2->type == NUM_COMPLEX) {
		// ### use cblas_cdotc()
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

// ### add suport for not-equally-sized mat using copy_sub_mat
void
Num_copy_mat(Num_mat *b, Num_mat *a)         // B <- A
{
	int		i, j, n;
	int		nr, nc, na, nb;
// chk header
	if (a->type != b->type) {
		Num_error("Num_copy: type mismatch");
	}
	if (a->nr == b->nr && a->nc == b->nc) {
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
	} else {
		nr = MIN(a->nr, b->nr);
		nc = MIN(a->nc, b->nc);
		na = a->nr * a->nc;
		nb = b->nr * b->nc;
		for (i = 0; i < nc; i++) {
			for (j = 0; j < nr; j++) {
				b->data[i * b->nr + j] = a->data[i * a->nr + j];
			}
		}
		if (a->type == NUM_COMPLEX) {
			for (i = 0; i < nc; i++) {
				for (j = 0; j < nr; j++) {
					b->data[i * b->nr + j + nb] = a->data[i * a->nr + j + na];
				}
			}
		}
	}
}

void
Num_copy_sub_mat(Num_mat *b, Num_mat *a, int r0, int rn, int c0, int cn)	// B <- A
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

		// printf("sizeof work area %d\n", lwork);
		// actual proc
		rwork = (float *)malloc(sizeof(float) * lwork);
		// copy A
		for (i = 0; i < nrc; i++) {
			adata[i] = A->data[i];
		}
		sgesvd_("S", "S", &nr, &nc, adata, &ld,
			s->data, U->data, &U->ld, Vt->data, &Vt->ld, rwork,
			&lwork, &info);
		if (info != 0) {
			printf("sgesvd info: %d\n", info);
		}
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
		if (info != 0) {
			printf("cgesvd info: %d\n", info);
		}
		free(rwork);
		free(cwork);
		free(adata);
	}
}

// solve linear equation (solution is returned in X. A and B are preserved.)
int
Num_inv(Num_mat *X, Num_mat *A, Num_mat *B)
{
	int		m, n, nrhs;
	int		lda, ldb, info;
	int		*ipiv;
	float	*work;
	int		lwork;
	Num_mat	*AA, *BB;	// copy mat first. both AA & BB are overwritten by LAPACK

	m = A->nr;
	n = A->nc;
	lda = m;
	ldb = MAX(m, n);
	nrhs = B->nc;

	// copy
	AA = Num_new_mat(A->nr, A->nc);
	BB = Num_new_mat(ldb, B->nc);
	Num_copy_mat(AA, A);
	Num_copy_mat(BB, B);		// works for mat with not equal sizes
	if (m == n) {
		ipiv = (int *)malloc(sizeof(int) * n);
		sgesv_(&n, &nrhs, AA->data, &lda, ipiv, BB->data, &ldb, &info);
		// copy output to B
		Num_copy_mat(X, BB);
		free(ipiv);
	} else {
		work = (float *)malloc(sizeof(float) * 1);
		lwork = -1;
		sgels_("N", &m, &n, &nrhs, AA->data, &lda, BB->data, &ldb, work, &lwork, &info);
		lwork = work[0];
		free(work);
		work = (float *)malloc(sizeof(float) * lwork);
		sgels_("N", &m, &n, &nrhs, AA->data, &lda, BB->data, &ldb, work, &lwork, &info);
		free(work);
		Num_copy_mat(X, BB);
	}
	Num_free_mat(AA);
	Num_free_mat(BB);

	return info;
}

float
Num_max_val(Num_mat *m, int *ix)
{
	int		i;
	float	mx = 0;

	for (i = 0; i < m->nr * m->nc; i++) {
		if (m->data[i] > mx) {
			mx = m->data[i];
			*ix = i;
		}
	}
	return mx;
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

	res = Num_new_svd_result(A);

	// call ref
	Num_svd_ref(A, res->U, res->s, res->Vt);

	// negate U & Vt
//	Num_negate_mat(res->U);
//	Num_negate_mat(res->Vt);

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

// remove mean for each col
void
Num_col_center(Num_mat *a)
{
	int		i, j;
	int		nr = a->nr;
	int		nc = a->nc;
	float	m;

	for (i = 0; i < nc; i++) {  // row direction
		m = 0;
		for (j = 0; j < nr; j++) {  // col direction
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

//=========== preserve until successfull port -> don't touch
Num_ica_result *
Num_ica_ref(Num_mat *X, int ncomp)
{
	Num_ica_result	*ires;	
	Num_svd_result	*sres;	
	int				n, p, ncn, ncnc;
	int				i, j;
	int				iter;
    int             maxiter = 2000;
	double			tol;

    // ncol x nrow, p:channels, pixels, n:time, observations
	Num_mat			*Xt;				// p x n
	Num_mat			*V, *D, *K, *U;		// p x p
	Num_mat			*K1;				// ncomp x p
	Num_mat			*X1;				// ncomp x n
	Num_mat			*tmp_mat;			// n x ncomp
	Num_mat			*WX, *gWX;			// n x ncomp
	Num_mat			*W;					// ncomp x ncomp
	Num_mat			*V1, *V2;			// ncomp x ncomp
	Num_mat			*W1;				// ncomp x ncomp
	double			tmp;
	
	BOOL			row_norm = NO;
    BOOL            dbg = YES;
//	float			alpha = 1.0;

    // input
	n = X->nr;	// n row (observations)
	p = X->nc;	// p column (pixels)

    if (dbg) {
        saveAsKOImage(X, @"IMG_X");
    }

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
	ires = (Num_ica_result *)malloc(sizeof(Num_ica_result));
//	ires = (Num_svd_result *)malloc(sizeof(Num_svd_result));
	ires->W  = Num_new_mat(ncomp, ncomp);	// unmixing mat
    ires->WX = Num_new_mat(n, ncomp);       // time course
    ires->WK = Num_new_mat(p, ncomp);       // map

	Num_copy_mat(ires->W, W);
	Num_copy_mat(ires->WX, Num_trans(WX));
 // WK not set yet ####   

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
// 2nd attempt (current version ###)
// restarted using Num_ica_ref as reference (6-17-2020)
Num_ica_result *
Num_ica_2(Num_mat *X, int ncomp)
{
    Num_ica_result  *ires;    
    Num_svd_result  *sres;    
    int             n, p, ncn, ncnc;
    int             i, j;
    int             iter;
    int             maxiter = 2000;
    double          tol;

    // ncol x nrow, p:channels, pixels, n:time, observations
    Num_mat         *Xt;                // p x n
    Num_mat         *V, *D, *K, *U;     // p x p (usaed for covar/SVD step only)
    Num_mat         *K1;                // ncomp x p
    Num_mat         *X1;                // ncomp x n
    Num_mat         *tmp_mat;           // n x ncomp
    Num_mat         *WX, *gWX;          // n x ncomp
    Num_mat         *W;                 // ncomp x ncomp
    Num_mat         *V1, *V2;           // ncomp x ncomp
    Num_mat         *W1;                // ncomp x ncomp
    double          tmp;
    
    BOOL            dbg = YES;

    // input
    n = X->nr;    // n row (observations)
    p = X->nc;    // p column (pixels, channels)


    if (dbg) {
        printf("n/p/ncomp = %d/%d/%d\n", n, p, ncomp);
        saveAsKOImage(X, @"IMG_X"); // input
    }

    Num_col_center(X);        // remove mean along col
    Xt = Num_trans(X);  // ## remove trans at some point (everything has to be trans'ed)

    if (0) {
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
printf("U covar\n");
dump_mat(sres->U);
printf("1/sqrt(s)\n");
dump_mat(D);
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
printf("K1 covar\n");
dump_mat(K1);
//        saveAsKOImage(K1, @"IMG_K1");
        Num_free_mat(V);
        Num_free_mat(D);
        Num_free_mat(K);
        Num_free_mat(U);
    } else {    // direct SVD ok (6-19-2020)
    K1 = Num_trans(X);  // Xt
saveAsKOImage(K1, @"IMG_K1");
 //   Num_mmul_const(K1, 1.0/sqrt(n)); // this can be done later to sigma
    sres = Num_svd(K1);
//    printf("U svd\n"); dump_mat(sres->U);
//    printf("s svd\n"); dump_vec(sres->s);
//    for (i = 0; i < 2; i++) {
//        printf("%f\n", 1.0 / (sres->s->data[i]));
//    }
//    dump_mat(sres->Vt);
    //    K1 = Num_new_mat(ncomp, p);
        K1 = Num_new_mat(ncomp, ncomp);
        for (i = 0; i < ncomp; i++) {
         //   for (j = 0; j < p; j++) {
            for (j = 0; j < ncomp; j++) {
                K1->data[j * ncomp + i] = sres->U->data[j * p + i] / (sres->s->data[i]) * sqrt(n); // NOT sqrt(s) ??? 
            }
        }
   printf("K1 svd\n"); dump_mat(K1);
    exit(0);
    //    saveAsKOImage(K1, @"IMG_K1");
    }

// X1 <- K %x% X
    X1 = Num_new_mat(ncomp, n);
    Num_mmul(X1, K1, Xt);
//printf("X1 svd\n"); dump_mat(X1);
saveAsKOImage(X1, @"IMG_X1");


//    w.init <- matrix(rnorm(n.comp^2, n.comp, n.comp)
    W = Num_new_mat(ncomp, ncomp);
    for (i = 0; i < ncomp * ncomp; i++) {
        W->data[i] = Num_nrml(0.0, 1.0);
    }
Num_unit_mat(W);
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
//saveAsKOImage(WX, @"IMG_WX");
//printf("ref WX\n"); dump_mat(WX);

        // gwx <- wx * exp(-wx^2)/2)
        for (i = 0; i < ncn; i++) {
            tmp = WX->data[i];
            // exp (faster)
            gWX->data[i] = tmp * exp(-(tmp *tmp)/2);                
            // logcosh
        //    gWX->data[i] = tanh(alpha * tmp);
        }
        // V1 = gWX * Xt / p
        Num_mtmul(V1, gWX, NO, X1, YES);
//printf("ref gWX\n"); dump_mat(gWX);
        for (i = 0; i < ncnc; i++) {
            V1->data[i] /= p;
        }
//saveAsKOImage(V1, @"IMG_V1");
//printf("ref V1\n"); dump_mat(V1);

        // g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
        for (i = 0; i < ncn; i++) {
            tmp = WX->data[i];
            // exp (faster)
            gWX->data[i] = (1 - tmp*tmp) * exp(-(tmp*tmp)/2);
            // logcosh
        //    gWX->data[i] = alpha * (1 - (tanh(alpha * tmp)) * (tanh(alpha * tmp)));
        }
//dump_mat(gWX);
        // V2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
        Num_clear_mat(tmp_mat);
        for (i = 0; i < ncomp; i++) {    // row
            tmp = 0;
            for (j = 0; j < n; j++) {    // col
                tmp += gWX->data[i + j * ncomp];
            }
            tmp /= n;
            tmp_mat->data[i * ncomp + i] = tmp;
        }
//saveAsKOImage(tmp_mat, @"IMG_diag");
//printf("ref tmp_mat\n"); dump_mat(tmp_mat);

        Num_mmul(V2, tmp_mat, W);
//saveAsKOImage(V2, @"IMG_V2");
//printf("ref V2\n"); dump_mat(V2);
//exit(0);
        // W1 <- V1 - V2
        for (i = 0; i < ncnc; i++) {
            W1->data[i] = V1->data[i] - V2->data[i];    // need to preserve V1 ###
        }
//saveAsKOImage(V1, @"IMG_V12");
        // make orthogonal
        // sW1 <- La.svd(W1)
        // W1 <- sW1$u %*% Diag*1/sW1$d) %*% t(sW1%u) %*% W1
        Num_orth_mat(W1);    // #### chk
//saveAsKOImage(V1, @"IMG_V12o");
                
        // calc tol (V1 * Wt)
        Num_mtmul(tmp_mat, W1, NO, W, YES);
//saveAsKOImage(tmp_mat, @"IMG_tol");

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

//exit(0); // ### checking ...
    }    // iteration

// Wout
    saveAsKOImage(W, @"IMG_Wout");
// WX
    Num_mmul(WX, W, X1);
//copy_mat_to_img(WX, img);
//[img saveAsKOImage:@"IMG_out"];
saveAsKOImage(WX, @"IMG_out");

// sort/set U, s, Vt
    ires = (Num_ica_result *)malloc(sizeof(Num_ica_result));
//    ires = (Num_svd_result *)malloc(sizeof(Num_svd_result));
    ires->W  = Num_new_mat(ncomp, ncomp);    // unmixing mat
    ires->WX = Num_new_mat(n, ncomp);       // time course
    ires->WK = Num_new_mat(p, ncomp);       // map

    Num_copy_mat(ires->W, W);
    Num_copy_mat(ires->WX, Num_trans(WX));
 // WK not set yet ####   

// free mem (chk ###)
    Num_free_mat(Xt);
    Num_free_mat(W);
    Num_free_mat(W1);
    Num_free_mat(WX);
    Num_free_mat(gWX);
    Num_free_mat(V1);
    Num_free_mat(V2);
    Num_free_mat(tmp_mat);

    return ires;
}

//========================

Num_ica_result *
Num_ica(Num_mat *X, int ncomp)
{
//    return Num_ica_ref(X, ncomp);
    return Num_ica_2(X, ncomp);
}

//============================
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

// 
float
Num_rk4_area(float (^deriv)(float x, float y), float st, float ed, int n)    // not tested yet
{
	float	x, y, dx;
	int		i;

	x = st;
	dx = (ed - st) / n;
	for (i = 0; i < n; i++) {
		Num_rk4(&x, &y, deriv, dx);
	}
	return y;
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
			printf(" %5.4e\n", v->data[i]);
		}
		break;
	case NUM_COMPLEX :
		printf("=== complex vector[%d] ====\n", v->n);
		for (i = 0; i < n; i++) {
			printf(" (%5.4e %5.4ei)\n", v->data[i*2], v->data[i*2 + 1]);
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
			printf("(stored in col major order\n");
			for (r = 0; r < nr; r++) {
				for (c = 0; c < nc; c++) {
				//	ix = Num_col_m_index(m, c, r);
					ix = c * m->nr + r;
					printf("%6.5e ", m->data[ix]);
				}
				printf("\n");
			}
		} else {
			printf("(stored in row major order\n");
			printf("(not supported\n");
			ix = 0;
			for (r = 0; r < nr; r++) {
				for (c = 0; c < nc; c++) {
					ix = r * m->nc + c;
					printf("%6.5e ", m->data[ix]);
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
