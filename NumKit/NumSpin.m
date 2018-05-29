//
//  NumSpin
//
//	=== plan ===
//	phase graph lib
//		*Scheffler version first
//			*chk with non-zero phi -> Zy is not zero !!! (contradicts with text)
//		*Oshio version (spherical PO)
//	bloch simulator lib using PO
//		define block-wise constant sequence structure
//		simulator calculates rot/relaxation/daq
//


#import "NumKit.h"

// GAMMA : Hz / mT
#define GAMMA	42570.0

Num_spin_s *
Num_new_spin_s(int order)
{
	Num_spin_s	*sp;
	int			i;

	sp = (Num_spin_s *)malloc(sizeof(Num_spin_s));
	sp->n = order;
	sp->Fx  = (float *)malloc(sizeof(float) * order * 2);
	sp->Fy  = (float *)malloc(sizeof(float) * order * 2);
	sp->Fxp = (float *)malloc(sizeof(float) * order * 2);
	sp->Fyp = (float *)malloc(sizeof(float) * order * 2);
	sp->Zx  = (float *)malloc(sizeof(float) * order);
	sp->Zy  = (float *)malloc(sizeof(float) * order);
	sp->Zxp = (float *)malloc(sizeof(float) * order);
	sp->Zyp = (float *)malloc(sizeof(float) * order);
	for (i = 0; i < order * 2; i++) {
		sp->Fx[i] = sp->Fy[i] = sp->Fxp[i] = sp->Fyp[i] = 0;
	}
	for (i = 0; i < order; i++) {
		sp->Zx[i] = sp->Zy[i] = sp->Zxp[i] = sp->Zyp[i] = 0;
	}

	return sp;
}

Num_spin *
Num_new_spin(int order)
{
	Num_spin	*sp;
	int			i;

	sp = (Num_spin *)malloc(sizeof(Num_spin));
	sp->n = order;
	sp->Ipr = (float *)malloc(sizeof(float) * order * 2);
	sp->Ipi = (float *)malloc(sizeof(float) * order * 2);
	sp->Imr = (float *)malloc(sizeof(float) * order * 2);
	sp->Imi = (float *)malloc(sizeof(float) * order * 2);
	sp->I0r = (float *)malloc(sizeof(float) * order * 2);
	sp->I0i = (float *)malloc(sizeof(float) * order * 2);
	sp->nIpr = (float *)malloc(sizeof(float) * order * 2);
	sp->nIpi = (float *)malloc(sizeof(float) * order * 2);
	sp->nImr = (float *)malloc(sizeof(float) * order * 2);
	sp->nImi = (float *)malloc(sizeof(float) * order * 2);
	sp->nI0r = (float *)malloc(sizeof(float) * order * 2);
	sp->nI0i = (float *)malloc(sizeof(float) * order * 2);
	for (i = 0; i < order * 2; i++) {
		sp->Ipr[i] = sp->Ipi[i] = 0;
		sp->Imr[i] = sp->Imi[i] = 0;
		sp->I0r[i] = sp->I0i[i] = 0;
		sp->nIpr[i] = sp->nIpi[i] = 0;
		sp->nImr[i] = sp->nImi[i] = 0;
		sp->nI0r[i] = sp->nI0i[i] = 0;
	}

	return sp;
}

void
Num_free_spin_s(Num_spin_s *sp)
{
	if (sp) {
		if (sp->Fx)  free(sp->Fx);
		if (sp->Fy)  free(sp->Fy);
		if (sp->Fxp) free(sp->Fxp);
		if (sp->Fyp) free(sp->Fyp);
		if (sp->Zx)  free(sp->Zx);
		if (sp->Zxp) free(sp->Zxp);
		if (sp->Zy)  free(sp->Zy);
		if (sp->Zyp) free(sp->Zyp);
		free(sp);
	}
}

void
Num_free_spin(Num_spin *sp)
{
	if (sp) {
		if (sp->Ipr)  free(sp->Ipr);
		if (sp->Ipi)  free(sp->Ipi);
		if (sp->Imr)  free(sp->Imr);
		if (sp->Imi)  free(sp->Imi);
		if (sp->I0r)  free(sp->I0r);
		if (sp->I0i)  free(sp->I0i);
		if (sp->nIpr)  free(sp->nIpr);
		if (sp->nIpi)  free(sp->nIpi);
		if (sp->nImr)  free(sp->nImr);
		if (sp->nImi)  free(sp->nImi);
		if (sp->nI0r)  free(sp->nI0r);
		if (sp->nI0i)  free(sp->nI0i);
		free(sp);
	}
}

Num_rf *
Num_new_rf(int n)
{
	Num_rf	*rf;
	int		i;

	rf = (Num_rf *)malloc(sizeof(Num_rf));
	rf->n = n;
	rf->rho   = (float *)malloc(sizeof(float) * n);
	rf->theta = (float *)malloc(sizeof(float) * n);
	for (i = 0; i < n; i++) {
		rf->rho[i] = rf->theta[i] = 0;
	}
	return rf;
}

void
Num_free_rf(Num_rf *rf)
{
	if (rf) {
		if (rf->rho)   free(rf->rho);
		if (rf->theta) free(rf->theta);
	}
}

void dump_rf(Num_rf *rf, int ix)
{
	printf("%d %f %f\n", ix, rf->rho[ix], rf->theta[ix]);
}

void
rotate_spin_s(Num_spin_s *sp, int ix, float alpha, float phi) 
{
	int		j, n = sp->n;
	int		ip, im;
	float	a, b, c, d, e, f, g, h, hb, gb, ec, fc;

	a = cos(alpha / 2); a *= a;
	b = sin(alpha / 2); b *= b;
	c = sin(alpha);
	d = 2 * cos(alpha);
	e = sin(phi);
	f = cos(phi);
	g = sin(phi * 2);
	h = cos(phi * 2);
	hb = h * b;
	gb = g * b;
	ec = e * c;
	fc = f * c;

	for (j = 0; j <= ix; j++) {	// -> move loop within func
		ip = n + j;
		im = n - j;
		sp->Fxp[ip] = a  * sp->Fx[ip]
					+ hb * sp->Fx[im]	+ gb * sp->Fy[im]
					+ ec * sp->Zx[j]	+ fc * sp->Zy[j];
		sp->Fyp[ip] = a  * sp->Fy[ip]
					- hb * sp->Fy[im]	+ gb * sp->Fx[im]
					- fc * sp->Zx[j]	+ ec * sp->Zy[j];
		sp->Fxp[im] = hb * sp->Fx[ip]	+ gb * sp->Fy[ip]
					+  a * sp->Fx[im]
					+ ec * sp->Zx[j]	- fc * sp->Zy[j];
		sp->Fyp[im] = gb * sp->Fx[ip]	- hb * sp->Fy[ip] 
					+  a * sp->Fy[im]
					- fc * sp->Zx[j]	- ec * sp->Zy[j];
		sp->Zxp[j]	= -ec * sp->Fx[ip] + fc * sp->Fy[ip]
					- ec * sp->Fx[im] + fc * sp->Fy[im]
					+ d * sp->Zx[j];
		sp->Zxp[j] /= 2;
		sp->Zyp[j]	= -fc * sp->Fx[ip] - ec * sp->Fy[ip]
					+ fc * sp->Fx[im] + ec * sp->Fy[im]
					+ d * sp->Zy[j];
		sp->Zyp[j] /= 2;
	}
}

// oshio version
void
rotate_spin(Num_spin *s, int ix, float alpha, float phi) 
{
	int		j, n = s->n;
	float	ca, sa, cp, sp;	// cos(alpha), sin(alpha), cos(phi), sin(phi)
	float	c2p, s2p;		// cos(2 * phi), sin(2 * phi)
	float	cap1, cam1;		// 0.5 * (cos(alpha) + 1), 0.5 * cos(alpha) -1)
	float	rs2;				// 1 / sqrt(2)

phi = -phi;	// becomes equivalent to Schefler version
	ca = cos(alpha);
	sa = sin(alpha);
	cap1 = 0.5 * (ca + 1);
	cam1 = 0.5 * (ca - 1);
	cp = cos(phi);
	sp = sin(phi);
	c2p = cos(2 * phi);
	s2p = sin(2 * phi);
	rs2 = sa / sqrt(2.0);
	
	for (j = n - ix; j <= n + ix; j++) {
		s->nI0r[j]	= ca * s->I0r[j]
					+ rs2 * ( sp * s->Ipr[j] + cp * s->Ipi[j])
					+ rs2 * (-sp * s->Imr[j] + cp * s->Imi[j]);

		s->nI0i[j]	= ca * s->I0i[j]
					+ rs2 * (-cp * s->Ipr[j] + sp * s->Ipi[j])
					+ rs2 * (-cp * s->Imr[j] - sp * s->Imi[j]);

		s->nIpr[j]	= rs2 * (-sp * s->I0r[j] + cp * s->I0i[j])
					+ cap1 * s->Ipr[j]
					+ cam1 * (c2p * s->Imr[j] + s2p * s->Imi[j]);

		s->nIpi[j]	= rs2 * (-cp * s->I0r[j] - sp * s->I0i[j])
					+ cap1 * s->Ipi[j]
					+ cam1 * (-s2p * s->Imr[j] + c2p * s->Imi[j]);

		s->nImr[j]	= rs2 * ( sp * s->I0r[j] + cp * s->I0i[j])
					+ cam1 * (c2p * s->Ipr[j] - s2p * s->Ipi[j])
					+ cap1 * s->Imr[j];

		s->nImi[j]	= rs2 * (-cp * s->I0r[j] + sp * s->I0i[j])
					+ cam1 * (s2p * s->Ipr[j] + c2p * s->Ipi[j])
					+ cap1 * s->Imi[j];
	}
}

void
shift_spin_s(Num_spin_s *sp, int ix, float e1, float e2, BOOL cpmg)
{
	int		np = sp->n;
	int		j;

	for (j = 0; j < np * 2; j++) {
		sp->Fx[j + 1] = sp->Fxp[j] * e2;
		sp->Fy[j + 1] = sp->Fyp[j] * e2;
	}
	for (j = 0; j < np; j++) {
		sp->Zx[j] = sp->Zxp[j] * e1;
		sp->Zy[j] = sp->Zyp[j] * e1;
		if (j == 0 && !cpmg) {
			sp->Zx[0] += 1.0 - e1;
		}
	}
}

void
shift_spin(Num_spin *s, int ix, float e1, float e2, BOOL cpmg)
{
	int		n = s->n;
	int		j;

// dbg
//e1 = e2 = 1.0;
	for (j = 1; j < n * 2 - 1; j++) {
		s->Ipr[j] = s->nIpr[j-1] * e2;
		s->Ipi[j] = s->nIpi[j-1] * e2;
		s->Imr[j] = s->nImr[j+1] * e2;
		s->Imi[j] = s->nImi[j+1] * e2;
	}
	for (j = 0; j < n * 2; j++) {
		if (j == n && !cpmg) {
		//	s->I0r[j] = 1.0 - (1.0 - s->nI0r[j]) * e1;
			s->I0r[j] = s->nI0r[j] * e1 + 1.0 - e1;	// same as above
		} else {
			s->I0r[j] = s->nI0r[j] * e1;
			s->I0i[j] = s->nI0i[j] * e1;
		}
	}
}

void
update_image_s(RecImage *imgF, RecImage *imgZ, Num_spin_s *sp, int ix)
{
	float	*p, *q;
	int		j, y;
	int		xDim = [imgF xDim];
	int		np = sp->n;

	p = [imgF data];
	q = p + [imgF dataLength];

	for (j = 0; j < np*2; j++) {
		y = np * 2 - j - 1;
		p[y * xDim + ix] = sp->Fx[j];
		q[y * xDim + ix] = sp->Fy[j];
	}
	p = [imgZ data];
	q = p + [imgZ dataLength];
	for (j = 0; j < np; j++) {
		y = np - j - 1;
		p[y * xDim + ix] = sp->Zx[j];
		q[y * xDim + ix] = sp->Zy[j];
	}
}

void
update_image(RecImage *i0, RecImage *ip, RecImage *im, Num_spin *s, int ix)
{
	float	*p, *q;
	int		j, y;
	int		xDim = [i0 xDim];
	int		n = s->n;

// I0
	p = [i0 data];
	q = p + [i0 dataLength];
	for (j = 0; j < n * 2; j++) {
		y = n * 2 - j - 1;
		p[y * xDim + ix] = s->I0r[j];
		q[y * xDim + ix] = s->I0i[j];
	}
// Ip
	p = [ip data];
	q = p + [ip dataLength];
	for (j = 0; j < n * 2; j++) {
		y = n * 2 - j - 1;
		p[y * xDim + ix] = s->Ipr[j];
		q[y * xDim + ix] = s->Ipi[j];
	}
// Im
	p = [im data];
	q = p + [im dataLength];
	for (j = 0; j < n * 2; j++) {
		y = n * 2 - j - 1;
		p[y * xDim + ix] = s->Imr[j];
		q[y * xDim + ix] = s->Imi[j];
	}
}

void
dump_spin_s(Num_spin_s *sp)
{
	int	i, i2, n = sp->n;

	printf("========================\n");
	for (i = 0; i < n; i++) {
		i2 = i + n + 1;
		printf("%d %f %f %f %f\n", i, sp->Fx[i2], sp->Fy[i2], sp->Zx[i], sp->Zy[i]);
	}
}

void
dump_spin(Num_spin *sp)
{
	int	i, n = sp->n;

	printf("========================\n");
	for (i = 0; i < n*2; i++) {
		printf("%d %f %f %f %f %f %f\n", i, sp->Ipr[i], sp->Ipi[i], sp->Imr[i], sp->Imi[i], sp->I0r[i], sp->I0i[i]);
	}
}

void
dump_spin_mag(Num_spin *sp, int ix)
{
	int		i, n = sp->n;
	float	rp, rm, r0;
	float	ip, im, i0;
	float	mp, mm, m0, mg;

	rp = rm = r0 = 0;
	ip = im = i0 = 0;
	for (i = 0; i < n*2; i++) {
		rp += sp->Ipr[i];
		rm += sp->Imr[i];
		r0 += sp->I0r[i];
		ip += sp->Ipi[i];
		im += sp->Imi[i];
		i0 += sp->I0i[i];
	}
	mp = sqrt(rp*rp + ip*ip);
	mm = sqrt(rm*rm + im*im);
	m0 = sqrt(r0*r0 + i0*i0);
	mg = sqrt(mp*mp + mm*mm + m0*m0);
	printf("%d %f %f %f %f\n", ix, mg, mp, mm, m0);
}

RecImage *
Num_phase_graph_s(Num_rf *rf, BOOL cpmg)
{
	Num_spin_s	*sp;
	int			i, n = rf->n;
	float		alpha, phi;
	float		e1, e2;
	RecImage	*imgF, *imgZ;

	if (rf->t1 == 0) {
		e1 = 1.0;
	} else {
		e1 = exp(-rf->tr / rf->t1);	// msec / msec
	}
	if (rf->t2 == 0) {
		e2 = 1.0;
	} else {
		e2 = exp(-rf->te / rf->t2);	// msec / msec
	}

//	printf("e1 : e2 = %f : %f\n", e1, e2);
	sp = Num_new_spin_s(n);
	imgF = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:n yDim:n * 2];
	imgZ = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:n yDim:n * 2];

	sp->Zx[0] = 1.0;

	for (i = 0; i < n; i++) {	// pulse #
		// update states
		update_image_s(imgF, imgZ, sp, i);
		alpha   = rf->rho[i]   * M_PI / 180;		// deg -> rad
		phi = rf->theta[i] * M_PI / 180;	// deg -> rad
		rotate_spin_s(sp, i, alpha, phi);
		shift_spin_s(sp, i, e1, e2, cpmg);

	}
//	dump_spin_s(sp);
[imgF saveAsKOImage:@"IMG_graphF"];
[imgZ saveAsKOImage:@"IMG_graphZ"];

	Num_free_spin_s(sp);
	return imgF;
}

RecImage *
check_zsymm(RecImage *z) {
	RecImage	*ic;
	int			i, j, y1, y2, xDim, yDim;
	float		*p, *q;

	ic = [z copy];
	xDim = [ic xDim];
	yDim = [ic yDim];
	p = [ic data];
	q = p + [ic dataLength];

	for (i = 1; i < yDim/2; i++) {
		y1 = yDim/2 - i - 1; y2 = yDim/2 + i - 1;
		for (j = 0; j < xDim; j++) {
			p[y1 * xDim + j] -= p[y2 * xDim + j];
			q[y1 * xDim + j] += q[y2 * xDim + j];
		}
	}
	return ic;
}

RecImage *
check_isymm(RecImage *ip, RecImage *im) {
	RecImage	*ic;
	int			i, j, y1, y2, xDim, yDim;
	float		*p1, *q1;	// src (im)
	float		*p2, *q2;	// dst (ic)

	ic = [RecImage imageWithImage:im];
	xDim = [ic xDim];
	yDim = [ic yDim];
	p1 = [im data];
	q1 = p1 + [im dataLength];
	p2 = [ic data];
	q2 = p2 + [ic dataLength];
// flip & conj
	for (i = 0; i < yDim-2; i++) {
		y1 = yDim - i - 2; y2 = i;
		for (j = 0; j < xDim; j++) {
			if (i == yDim/2 && j == 10) continue;
			p2[y1 * xDim + j] = -p1[y2 * xDim + j];
			q2[y1 * xDim + j] =  q1[y2 * xDim + j];
		}
	}
	[ic subImage:ip];
	p2[0] = 1;
	p2[1] = -1;

	return ic;
}


// Oshio version
RecImage *
Num_phase_graph(Num_rf *rf, BOOL cpmg)
{
	Num_spin	*sp;
	int			i, n = rf->n;
	float		rho, theta;
	float		e1, e2;
	RecImage	*i0, *ip, *im;
	RecImage	*ic;	// ip conj

	if (rf->t1 == 0) {
		e1 = 1.0;
	} else {
		e1 = exp(-rf->tr / rf->t1);	// msec / msec
	}
	if (rf->t2 == 0) {
		e2 = 1.0;
	} else {
		e2 = exp(-rf->te / rf->t2);	// msec / msec
	}

// printf("e1 : e2 = %f : %f\n", e1, e2);
	sp = Num_new_spin(n);
	i0 = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:n yDim:n * 2];
	ip = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:n yDim:n * 2];
	im = [RecImage imageOfType:RECIMAGE_COMPLEX xDim:n yDim:n * 2];

	sp->I0r[n] = 1.0;

	for (i = 0; i < n; i++) {	// pulse #
		update_image(i0, ip, im, sp, i);
		rho   = rf->rho[i]   * M_PI / 180;		// deg -> rad
		theta = rf->theta[i] * M_PI / 180;	// deg -> rad
		rotate_spin(sp, i, rho, theta);
		shift_spin(sp, i, e1, e2, cpmg);
	}
//	dump_spin(sp);

[i0 saveAsKOImage:@"IMG_graph0"];
[ip saveAsKOImage:@"IMG_graphp"];
[im saveAsKOImage:@"IMG_graphm"];
[im copyLoopsOf:ip];
[i0 copyLoopsOf:ip];

ic = [ip copy];
[ic addImage:i0];
[ic saveAsKOImage:@"IMG_graph0p"];

if (1) {	// test for symmetry
	ic = check_isymm(ip, im);
	[ic saveAsKOImage:@"IMG_dif"];
	ic = check_zsymm(i0);
	[ic saveAsKOImage:@"IMG_dif0"];
}
	Num_free_spin(sp);
	return ip;
}


