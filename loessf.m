***************************************************************
* LOESS   smoothing scattered data in one or more variables   *
*         documentation of Fortran routines                   *
*         Cleveland, Devlin, Grosse, Shyu                     *
***************************************************************

1. The typical program would call lowesd, set tolerances in iv,v if
   desired, then call lowesb and lowese.
2. To save the k-d tree, call lowesd, lowesb and then losave; subsequent
   programs would call lohead, set liv and lv, then call loread and lowese.
3. For statistics, get diagL and then call lowesa or get the full hat
   matrix and call lowesc.  Robustness iterations can take advantage of
   lowesw and lowesp.

lowesd(106,iv,liv,lv,v,d,n,f,tdeg,nvmax,setLf)	setup workspace
lowesf(x,y,w,iv,liv,lv,v,m,z,L,hat,s)		slow smooth at z
lowesb(x,y,w,diagL,infl,iv,liv,lv,v)		build k-d tree
lowesr(y,iv,liv,lv,v)				rebuild with new data values
						(does not change y)
lowese(iv,liv,lv,v,m,z,  s)			evaluate smooth at z
lowesl(iv,liv,lv,v,m,z,  L)			explicit hat matrix,
						which maps from y to z
lofort(iunit,iv,liv,lv,v)			save k-d tree as Fortran
losave(iunit,iv,liv,lv,v)			save k-d tree in file
lohead(iunit,d,vc,nc,nv)			read d,vc,nc,nv from file
	liv = 50+(vc+3)*nc			determine space
	lv = 50+(2*d+1)*nv+nc				requirements
loread(iunit,d,vc,nc,nv,iv,liv,lv,v)		finish reading k-d tree,
							ready for lowese
lowesa(trL,n,d,tau,nsing,  del1,del2)		approximate delta
lowesc(n,L,LL,  trL,del1,del2)			exact delta
lowesp(n,y,yhat,w,rw,  pi,ytilde)		pseudo-values
lowesw(res,n,  rw,pi)				robustness weights

=== arguments ===
d	number of independent variables [integer]  (called "p" elsewhere)
del1,del2	delta1, delta2
diagL	diagonal of hat matrix, only set if infl=.true.    (n)
f	fraction of points to use in local smooth  (called "alpha" elsewhere)
fc	don't refine cells with less than fc*n points;   ordinarily=.05
hat	is hat matrix desired?  [integer]
	0 = none
	1 = diagonal only
	2 = full matrix
infl	is diagonal of hat matrix desired?	[logical]
iunit	Fortran unit number for i/o
iv	workspace  (liv)
L	hat matrix (m,n)   [real]
	in lowesf, only computed if hat nonzero;  if hat=1 only size (n)
LL	workspace (n,n)
liv	50+(2^d+4)*nvmax+2*n
	if setLf, add nf*nvmax
lv	50+(3*d+3)*nvmax+n+(tau0+2)*nf
	if setLf, add (d+1)*nf*nvmax
m	number of points to smooth at;   ordinarily=n
n	number of observations
nf	min(n,floor(n*f))
nsing	if 0, print warning in lowesa when trL<tau;  typically nsing=iv(30)
nvmax	limit on number of vertices for kd-tree; e.g. max(200,n)
pi	workspace (n)  [integer]
res	residual  yhat-y  (n)
rw	robustness weights  (n)
s	smoothed values at z  (m)
setLf	in lowesb, save matrix factorizations  [logical]
	(needed for lowesr and lowesl)
tau	dimension of local model = iv(DIM);
	=d+1 for linear, (d+2)(d+1)/2 for quadratic
		reduced if dropping squares
tau0	=d+1 for linear, (d+2)(d+1)/2 for quadratic
tdeg    polynomials to fit;  0=constants, 1=linear, 2=quadratics
trL	trace L = sum diagL
v	workspace  (lv)
w	weights  (n)    local regression: min sum wi * (f(xi)-yi)^2
x	sample locations  (n,d)
y	observations  (n)
yhat	smoothed y  (n)
ytilde	pseudo y  (n)
z	locations where smooth is desired  (m,d)

If using the double precision version, [real] above should be understood
as Fortran "double precision".

The first argument to lowesd is a version number, updated when calling
sequences change.

If you peek inside the fortran, you will quickly notice that it
was machine generated;  the typeset original (in the language "pine")
is much easier to read.

=== iv indices ===
1	INFO	return code (not currently used)
2	D	number of independent variables
3	N	number of observations
4	VC	2^d  (number of vertices of a cell)
5	NC	number of k-d cells
6	NV	number of k-d vertices
7	A1	starting index in iv of a
8	C1	starting index in iv of c
9	HI1	starting index in iv of hi
10	LO1	starting index in iv of lo
11	V1	starting index in v of vertices
12	XI1	starting index in v of cut values
13	VV1	starting index in v of vertex values
14	NVMAX	maximum allowed value of nv
15	WORK1	starting index in v of workspace
16	WORK2	starting index in v of workspace
17	NCMAX	maximum allowed value of nc
18	WORK3	starting index in v of workspace
19	NF	floor(n*f) (number of points used as neighborhood)
20	KERNEL	1=tricube, 2=unif
21	KIND	1=k-d,cubic blend, (not implemented:2=quadtree,3=triangulation)
22	PI1	starting index in iv of tree permutation
23	VH	starting index in iv of vhit
24	VV2	starting index in v of work vval used in trL computation
25	LQ	starting index in iv of Lq
26	WORK4	starting index in v of workspace
27	PSI1	starting index in iv of workspace permutation
28	SEQ	sequence number, to check if routines called out of order
		takes on values:
		171	after lowesd
		172	after lowesf
		173	after lowesb
29	DIM	dimension of local regression
		1		constant
		d+1		linear   (default)
		(d+2)(d+1)/2	quadratic
		Modified by ehg127 if cdeg<tdeg.
30	SING	number of times singular tolerance was met in l2fit, l2tr
31	PRINT	verbose output?
32	DEG	total degree (of polynomial for local model)
33	NDIST	dd = variables 1:dd enter into distance calculation
34	LF	starting index in v of Lf
35..40		reserved for future use
41..49	CDEG	componentwise degree
iv(A1)	a	coordinate of cut; 0 for leaf  (nc)
iv(C1)	c	pointers to corners (index into vertex array v)  (vc,nc)
iv(HI1)	hi	right subcell  (nc)
iv(LO1)	lo	left subcell  (nc)
		Leaf cell j encloses points x(pi(i),), lo(j)<=i<=hi(j).
		Also, iv(C1),...,iv(PI1-1) is used as workspace (t) by l2fit
------------------------eval only needs workspace up to here
iv(PI1)	pi	permutation of 1:n for listing points in cells
iv(VH)	vhit	cell whose subdivision creates vertex (nv)
		0 if vertex is corner of original bounding box.
iv(LQ)	Lq	active point indices for block of Lf    (nvmax*nf)
iv(PSI1) psi	workspace permutation of 1:n for sorting distances

=== v indices ===
1	F	fraction of n to be used as neighborhood.   See also iv(19).
2	FCELL	no refinement if #points <= fcell * n
		default .05
3	FDIAM	no refinement if diameter is fdiam * overall bounding box
		default 0;    Warning: reset to 0 by ehg142 when nsteps>0.
4	RCOND	reciprocal condition number
... 49		reserved for future use
iv(V1)	v	vertices  (nv,d)
iv(VV1)	vval	vertex values  (0:d,nv)
iv(XI1)	xi	cut values  (nc)
------------------------eval only needs workspace up to here
iv(WORK1)	workspace  (n)  l2fit:dist
iv(WORK2)	workspace  (nf) l2fit:eta
iv(WORK3)	workspace  (dim,nf)   l2fit:X
iv(VV2)	vval2	workspace  ((d+1)*nv)  pseudo-vval for trL computation
iv(LF) 	Lf	hat matrix (data to vertex)   ((d+1)*nvmax*nf)
iv(WORK4)	workspace  (nf) l2fit:w

Internal routine names have been hidden as follows:
ehg106  select q-th smallest by partial sorting
ehg124  rbuild
ehg125  cpvert
ehg126  bbox
ehg127  l2fit,l2tr computational kernel
ehg128  eval
ehg129  spread
ehg131  lowesb after workspace expansion
ehg133  lowese after workspace expansion
ehg134  abort by calling S Recover function
ehg136  l2fit with hat matrix L
ehg137  vleaf
ehg138  descend
ehg139  l2tr
ehg140(w,i,j)   w(i)=j   used when w is declared real, but should store an int
ehg141  delta1,2 from trL
ehg142  robust iteration
ehg144  now called lowesc
ehg152  like ehg142, but for lowesf
ehg167  kernel for losave
ehg168  kernel for loread
ehg169  compute derived k-d tree information
ehg170  generate Fortran
ehg176,ehg177,ehg178,ehg179,ehg180,ehg181  loeval for delta
ehg182  ehgdie(i)
ehg183	warning(message,i,n,inc)
ehg184	warning(message,x,n,inc)
ehg190  now called lowesa, with slight change in calling sequence
ehg191  lowesl after workspace expansion
ehg192  lowesr after workspace expansion
ehg196(tau,d,f,trl)	trL approximation
ehg197    for deg=1,2
m9rwt	now called lowesw
pseudo	now called lowesp
