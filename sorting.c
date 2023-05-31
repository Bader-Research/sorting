#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <limits.h>
#include <strings.h>

#define DEFAULT_ELEMS  25
#define LOOP_CNT       25
#define SORT_CUTOFF   100

/* #define DO_ISORT        0 */

#ifdef GCC
#define INLINE inline
/* #define INLINE */
#else
#define INLINE
#endif

#define DATA_TYPE int
#define ODD(n) ((n)&1)==1
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

struct timeval  tp;
struct timezone tzp;

#define get_seconds()   (gettimeofday(&tp, &tzp), \
                        (double)tp.tv_sec + (double)tp.tv_usec / 1000000.0)


void assert_malloc(void *ptr) {
    if (ptr==NULL) {
	fprintf(stderr,"ERROR: Null pointer\n");
	exit(1);
    }
}

int check_sort(int *list, int n) {

  int i;
  int val = 1;

  for (i=1 ; i<n ; i++) 
    if (list[i] < list[i-1]) {
      fprintf(stderr,"ERROR: list[%d] < list[%d]  (%d < %d)\n",
	      i,i-1, list[i],list[i-1]);
      val = 0;
    }

  return val;
}

static INLINE unsigned bits(unsigned x, int k, int j) {
/* Return the j bits which appear k bits from the right in x */
    return (x>>k) & ~(~0<<j);
}

/* static */ int intcompare(const void *i, const void *j)
{
    return(*(int *)i - *(int *)j);
}

static INLINE void radixsort_3(int *a, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
#if 0
      w,
#endif
	M;

    int	*count,
	bitArr[n],
	b[n];

#if 0
    w = 32;                 /* The number of bits in the key */
#endif
    M = 1 << 11;
    
    count = (int*)malloc(M*sizeof(int));
    assert_malloc(count);

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    /* Digit 0 (10..0) */
    for (j=0 ; j<M ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[bitArr[i] = bits(a[i],0,11)]++;

    for (j=1 ; j<M ; j++) count[j] += count[j-1];
    for (i=n-1 ; i>=0 ; i--) b[--count[bitArr[i]]] = a[i]; 

    /* Digit 1 (21..11) */
    for (j=0 ; j<M ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[bitArr[i] = bits(b[i],11,11)]++;

    for (j=1 ; j<M ; j++) count[j] += count[j-1];
    for (i=n-1 ; i>=0 ; i--) a[--count[bitArr[i]]] = b[i];

    /* Digit 3 (31..22) */
    M = 1 << 10;
    for (j=0 ; j<M ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[bitArr[i] = bits(a[i],22,10)]++;

    for (j=1 ; j<M ; j++) count[j] += count[j-1];
    for (i=n-1 ; i>=0 ; i--) b[--count[bitArr[i]]] = a[i];

    bcopy(b,a,n*sizeof(int));

    free(count);
}

static INLINE void radixsort_3_rev(int *a, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
#if 0
      w,
#endif
	M;

    int	*count,
	bitArr[n],
	b[n];

#if 0
    w = 32;                 /* The number of bits in the key */
#endif
    M = 1 << 11;
    
    count = (int*)malloc(M*sizeof(int));
    assert_malloc(count);

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    /* Digit 0 (10..0) */
    for (j=0 ; j<M ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[bitArr[i] = bits(a[i],0,11)]++;

    for (j=1 ; j<M ; j++) count[j] += count[j-1];
    for (j=M-1 ; j>0 ; j--) count[j] = count[j-1];
    count[0] = 0;
    for (i=0 ; i<n ; i++)  b[count[bitArr[i]]++] = a[i];
    
    /* Digit 1 (21..11) */
    for (j=0 ; j<M ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[bitArr[i] = bits(b[i],11,11)]++;

    for (j=1 ; j<M ; j++) count[j] += count[j-1];
    for (j=M-1 ; j>0 ; j--) count[j] = count[j-1];
    count[0] = 0;
    for (i=0 ; i<n ; i++) a[count[bitArr[i]]++] = b[i];

    /* Digit 3 (31..22) */
    M = 1 << 10;
    for (j=0 ; j<M ; j++) count[j] = 0;
    for (i=0 ; i<n ; i++) count[bitArr[i] = bits(a[i],22,10)]++;

    for (j=1 ; j<M ; j++) count[j] += count[j-1];
    for (j=M-1 ; j>0 ; j--) count[j] = count[j-1];
    count[0] = 0;
    for (i=0 ; i<n ; i++) b[count[bitArr[i]]++] = a[i];

    bcopy(b,a,n*sizeof(int));

    free(count);
}

static INLINE void radixsort_4(int *a, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
	m,
	w,
	M,
	pass;

    int	nextpass,
	*count,
	bitArr[n],
	b[n];

    w = sizeof(int) << 3;   /* The number of bits in the key */
    m = w >> 2;             /* m = number of bits per pass   */
    M = 1 << m;             /* The range of each pass */
    
    count = (int*)malloc(M*sizeof(int));
    assert_malloc(count);

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    for (pass=0 ; pass<(w/m) ; pass+=2) {
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(a[i],pass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) b[--count[bitArr[i]]] = a[i]; 

	nextpass = pass+1;
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(b[i],nextpass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) a[--count[bitArr[i]]] = b[i];
    }
    free(count);
}

static INLINE void radixsort_4p(int *a, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	m,
	w,
	M,      
	pass;

    int	nextpass,
        *c1, *c2,
	*count,
	bitArr[n],
	b[n];

    w = sizeof(int) << 3;   /* The number of bits in the key */
    m = w >> 2;             /* m = number of bits per pass   */
    M = 1 << m;             /* The range of each pass */
    
    count = (int*)malloc(M*sizeof(int));
    assert_malloc(count);

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    c2 = count + M;
    for (pass=0 ; pass<(w/m) ; pass+=2) {
        c1 = count;
	while (c1 < c2)
	  *c1++ = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(a[i],pass*m,m)]++;
	c1 = count;
	while (++c1 < c2)
	  *c1 += *(c1-1);
	for (i=n-1 ; i>=0 ; i--) b[--count[bitArr[i]]] = a[i]; 

	nextpass = pass+1;
	c1 = count;
	while (c1 < c2)
	  *c1++ = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(b[i],nextpass*m,m)]++;
	c1 = count;
	while (++c1 < c2) 
	  *c1 += *(c1-1);
	for (i=n-1 ; i>=0 ; i--) a[--count[bitArr[i]]] = b[i];
    }
    free(count);
}

INLINE static void fastsort_3(int* arr, int nel) {

    if (nel > SORT_CUTOFF)
	radixsort_3(arr,nel); 
    else qsort(arr,nel,sizeof(int),intcompare); 
}

INLINE static void fastsort_4(int* arr, int nel) {

    if (nel > SORT_CUTOFF)
	radixsort_4(arr,nel); 
    else qsort(arr,nel,sizeof(int),intcompare);
}

INLINE static void fastsort_4p(int* arr, int nel) {

    if (nel > SORT_CUTOFF)
	radixsort_4p(arr,nel); 
    else qsort(arr,nel,sizeof(int),intcompare);
}


#if DO_ISORT
INLINE static void insertsort(DATA_TYPE *A, int n) {

  register DATA_TYPE item;
  register int i,j;

  for (i=1 ; i<n ; i++) {
    item = A[i];
    j = i-1;
    while ((j>=0)&&(item < A[j])) {
      A[j+1] = A[j];
      j--;
    }
    A[j+1] = item;
  }
}
#endif

INLINE static void merge(DATA_TYPE *A, int l, int m, int r, DATA_TYPE *B) {
  register int h,i,j,k;

  i = h = l;
  j = m+1;

  while ((h <= m) && (j <= r)) {
    if (A[h] <= A[j]) {
      B[i] = A[h];
      h++;
    }
    else {
      B[i] = A[j];
      j++;
    }
    i++;
  }
  if (h>m) {
    for (k=j ; k <= r ; k++) {
      B[i] = A[k];
      i++;
    }
  }
  else {
    for (k=h ; k<= m ; k++) {
      B[i] = A[k];
      i++;
    }
  }
  bcopy(B+l,A+l,(r-l+1)*sizeof(DATA_TYPE));
}

static void mergesort_r(DATA_TYPE *A, int l, int r, DATA_TYPE *B) {

  register int m;
  
  if (l < r) {
    m = (l+r) / 2;
    mergesort_r(A, l,   m, B);
    mergesort_r(A, m+1, r, B);
    merge(A, l, m, r, B);
  }
}
		 
INLINE static void mergesort(DATA_TYPE *A, int n) {

  DATA_TYPE *B;
  B = (DATA_TYPE *)malloc(n*sizeof(DATA_TYPE));
  assert_malloc(B);

  mergesort_r(A,0,n-1,B);

  free(B);
}

INLINE static DATA_TYPE * mergesort_nr_2(DATA_TYPE *Buffer1,
					 DATA_TYPE *Buffer2,
					 int col_size) {

#if 0
#define MAX_VAL LONG_MAX 
#define MAX_VAL_m1 (MAX_VAL - 1)
#else
  long MAX_VAL = LONG_MAX;
  long MAX_VAL_m1 = (MAX_VAL - 1);
#endif

  int
    i, j, k,
    values, times,
    adjustments, last1, last2,
    tot, tot1, tot2,
    rem;

  DATA_TYPE
    *start_ptr,
    *finish_ptr,
    *temp_ptr,
    *ex_ptr,
    *s1_ptr, *s2_ptr;

  start_ptr = Buffer1;
  finish_ptr = Buffer2;

  times = (int) ceil(log((double) col_size)/log((double) 2)); 

  adjustments = 0;
  values = col_size;

  if (col_size & 1)
    temp_ptr = start_ptr + col_size - 1;
  else  
    temp_ptr = start_ptr + col_size;    
  while (start_ptr < temp_ptr) {
    if ((*start_ptr) > (*(start_ptr+1))) {
      *(finish_ptr++) = *(start_ptr+1);
      *(finish_ptr) = *(start_ptr);
      if ((*(finish_ptr++)) > MAX_VAL_m1) {
	adjustments++;
	(*(finish_ptr - 1)) = MAX_VAL_m1;
      }
      *(finish_ptr++) = MAX_VAL;
      start_ptr+=2;
    }
    else {
      *(finish_ptr++) = *(start_ptr++);
      *(finish_ptr) = *(start_ptr++);
      if ((*(finish_ptr++)) > MAX_VAL_m1) {
	adjustments++;
	(*(finish_ptr-1)) = MAX_VAL_m1;
	if ((*(finish_ptr - 2)) > MAX_VAL_m1) {
	  adjustments++;
	  (*(finish_ptr-2)) = MAX_VAL_m1; 
	}
      }

      *(finish_ptr++) = MAX_VAL;
    }
  } 

  last1 = last2 = 2;

  if (col_size & 1) {
    *(finish_ptr) = *(temp_ptr);
    if ((*(finish_ptr++)) > MAX_VAL_m1) {
      adjustments++;
      (*(finish_ptr-1)) = MAX_VAL_m1;
    }
    *(finish_ptr) = MAX_VAL;
    last2 = 1;
  }

  ex_ptr = Buffer1;
  Buffer1 = Buffer2;
  Buffer2 = ex_ptr;

  tot = 2;
  tot1 = 1;
  tot2 = col_size >> 1;

  for (i=2;i<times;i++) {
    start_ptr = Buffer1;
    finish_ptr = Buffer2;
    tot <<= 1; 
    tot1 <<= 1;
    tot2 >>= 1;
    for (j=1;j<tot2;j++) {
      s2_ptr = (s1_ptr = (start_ptr)) + tot1 + 1;  
      temp_ptr = finish_ptr + tot;
      while (finish_ptr < temp_ptr) {
	if ((*s1_ptr) < (*s2_ptr))
	  *(finish_ptr++) = *(s1_ptr++);
	else 
	  *(finish_ptr++) = *(s2_ptr++);
      }
      *(finish_ptr++) = MAX_VAL;       
      start_ptr += (tot + 2);
    }
  
    rem = (values & (tot - 1));
  
    if (! rem) {
      s2_ptr = (s1_ptr = (start_ptr)) + tot1 + 1;  
      temp_ptr = finish_ptr + tot;
      while (finish_ptr < temp_ptr) {
	if ((*s1_ptr) < (*s2_ptr))
	  *(finish_ptr++) = *(s1_ptr++);
	else 
	  *(finish_ptr++) = *(s2_ptr++);
      }
      *(finish_ptr++) = MAX_VAL;       
      last1 = last2 = tot;
    }
    else {
      if (rem <= tot1) {
	s1_ptr = start_ptr;
	for (k=0;k<=tot1;k++)
	  *(finish_ptr + k) = *(s1_ptr + k);
	finish_ptr += (tot1 + 1); 
	start_ptr += (tot1 + 1);
	s2_ptr = (s1_ptr = (start_ptr)) + last1 + 1;  
	temp_ptr = finish_ptr + last1 + last2;
	while (finish_ptr < temp_ptr) {
	  if ((*s1_ptr) < (*s2_ptr))
	    *(finish_ptr++) = *(s1_ptr++);
	  else 
	    *(finish_ptr++) = *(s2_ptr++);
	}
	*(finish_ptr++) = MAX_VAL;       
	last2 += last1;
	last1 = tot1;
      }
      else {
	s2_ptr = (s1_ptr = (start_ptr)) + tot1 + 1;  
	temp_ptr = finish_ptr + tot;
	while (finish_ptr < temp_ptr) {
	  if ((*s1_ptr) < (*s2_ptr))
	    *(finish_ptr++) = *(s1_ptr++);
	  else 
	    *(finish_ptr++) = *(s2_ptr++);
	} 
	*(finish_ptr++) = MAX_VAL;       
	start_ptr += (tot + 2);
	s2_ptr = (s1_ptr = (start_ptr)) + last1 + 1;  
	temp_ptr = finish_ptr + last1 + last2;
	while (finish_ptr < temp_ptr) {
	  if ((*s1_ptr) < (*s2_ptr))
	    *(finish_ptr++) = *(s1_ptr++);
	  else 
	    *(finish_ptr++) = *(s2_ptr++);
	}
	*(finish_ptr++) = MAX_VAL;       
	last2 += last1;
	last1 = tot;
      }
    }
    ex_ptr = Buffer1;
    Buffer1 = Buffer2;
    Buffer2 = ex_ptr;
  } 

  start_ptr = Buffer1;
  finish_ptr = Buffer2;
  s2_ptr = (s1_ptr = (start_ptr)) + last1 + 1;  
  temp_ptr = finish_ptr + last1 + last2;

  while (finish_ptr < temp_ptr) {
    if ((*s1_ptr) < (*s2_ptr))
      *(finish_ptr++) = *(s1_ptr++);
    else 
      *(finish_ptr++) = *(s2_ptr++);
  }

  finish_ptr = Buffer2 + values;;
  temp_ptr = finish_ptr - adjustments;

  while ((--finish_ptr) >= temp_ptr) 
    *finish_ptr =  MAX_VAL;

  finish_ptr = Buffer2;

  return(finish_ptr);

}

INLINE static void mergesort_nr(DATA_TYPE *A, int elems) {
  DATA_TYPE *finlist, *B;

  B = (DATA_TYPE *)malloc(2*elems*sizeof(DATA_TYPE));
  assert_malloc(B);
  
  finlist = mergesort_nr_2(A, B, elems);
  if (finlist == B)
    bcopy(B,A,elems*sizeof(DATA_TYPE));

  free(B);
  
}


int * radixsort_h_2(int *Buffer1, int *Buffer2, int col_size, int bits) {

  register int j, v,
  radix, radix_m1,
  digits, values,
  log_value_range, rem, 
  shift1, shift2;

  int
    *Count1, *Count2,
    *start_ptr, *finish_ptr,
    *Count, *Next, *t_ptr,
    *c_ptr, *ex_ptr,
    *Offset, *temp_ptr;

  Count1 = (int *) malloc(((1 << bits)+1)*sizeof(int));
  assert_malloc(Count1);
  Count2 = (int *) malloc(((1 << bits)+1)*sizeof(int));
  assert_malloc(Count2);

  radix = (1 << bits);  
  radix_m1 = radix - 1;
  values = col_size;

  log_value_range = 31;

  digits = (int) ceil(((double) log_value_range)/((double) bits));

  rem = log_value_range % bits;

  if (!rem)
    rem = bits;

  start_ptr = Buffer1;
  finish_ptr = Buffer2;

  Count = (Offset = Count1) + 1;
  Next = Count2 + 1;

  if (digits > 1)
    {t_ptr = (c_ptr = Count) + radix;
     while (c_ptr < t_ptr)
       *(c_ptr++) = 0;
     t_ptr = (c_ptr = Next) + radix;
     while (c_ptr < t_ptr)
       *(c_ptr++) = 0;
     t_ptr = (c_ptr = start_ptr) + values;
     while (c_ptr < t_ptr)
       (*(Count + ((*(c_ptr++)) & radix_m1)))++;
     *(Offset) = 0;
     t_ptr = (c_ptr = Offset) + radix;
     while ((++c_ptr) < t_ptr)
       *c_ptr +=  *(c_ptr - 1);
     t_ptr = (c_ptr = start_ptr) + values;
     while (c_ptr < t_ptr) {
       (*(Next + (((v = *(c_ptr++)) >> bits) & radix_m1)))++;
       *(finish_ptr + ((*(Offset + (v & radix_m1)))++)) = v;
     }
   
     ex_ptr = Buffer1;
     Buffer1 = Buffer2;
     Buffer2 = ex_ptr;
   
     start_ptr = Buffer1;
     finish_ptr = Buffer2;

     temp_ptr = Count;
     Count = Next;
     Next = temp_ptr;
     Offset = Count - 1;
   
     for (j = 1 ; j < (digits - 1) ; j++) {
       shift1 = j * bits;
       shift2 = (j + 1)*bits;
       t_ptr = (c_ptr = Next) + radix;
       while (c_ptr < t_ptr)
	 *(c_ptr++) = 0;
       *(Offset) = 0;
       t_ptr = (c_ptr = Offset) + radix;
       while ((++c_ptr) < t_ptr)
	 *c_ptr +=  *(c_ptr - 1);
       t_ptr = (c_ptr = start_ptr) + values;
       while (c_ptr < t_ptr) {
	 (*(Next + (((v = (*(c_ptr ++))) >> shift2) & radix_m1)))++;
	 *(finish_ptr + ((*(Offset + ((v >> shift1) & radix_m1)))++)) = v;
       }
      
       ex_ptr = Buffer1;
       Buffer1 = Buffer2;
       Buffer2 = ex_ptr;
   
       start_ptr = Buffer1;
       finish_ptr = Buffer2;
 
       temp_ptr = Count;
       Count = Next;
       Next = temp_ptr;
       Offset = Count - 1;
     }

     radix_m1 = (radix = (1 << rem)) - 1;
     shift1 = (digits - 1) * bits;
     *(Offset) = 0;
     t_ptr = (c_ptr = Offset) + radix;
     while ((++c_ptr) < t_ptr)
       *c_ptr +=  *(c_ptr - 1);
     t_ptr = (c_ptr = start_ptr) + values;
     while (c_ptr < t_ptr) {
       v = *(c_ptr++);
       *(finish_ptr + ((*(Offset + ((v >> shift1) & radix_m1)))++)) = v;
     }
   } 
  else {
    radix_m1 = (radix = (1 << rem)) - 1;
    t_ptr = (c_ptr = Count) + radix;
    while (c_ptr < t_ptr)
      *(c_ptr++) = 0;
    t_ptr = (c_ptr = start_ptr) + values;
    while (c_ptr < t_ptr)
      (*(Count + ((*(c_ptr++)) & radix_m1)))++;
    *(Offset) = 0;
    t_ptr = (c_ptr = Offset) + radix;
    while ((++c_ptr) < t_ptr)
      *c_ptr +=  *(c_ptr - 1);
    t_ptr = (c_ptr = start_ptr) + values;
    while (c_ptr < t_ptr) {
      v = *(c_ptr++);
      *(finish_ptr + ((*(Offset + (v & radix_m1)))++)) = v;
    }
  }
  
  finish_ptr = Buffer2;

  free(Count1);
  free(Count2);

  return(finish_ptr);
}


INLINE static void radixsort_h3(DATA_TYPE *A, int elems) {
  DATA_TYPE *finlist, *B;

  B = (DATA_TYPE *)malloc(elems*sizeof(DATA_TYPE));
  assert_malloc(B);
  
  finlist = radixsort_h_2(A, B, elems, 11);
  if (finlist == B)
    bcopy(B,A,elems*sizeof(DATA_TYPE));

  free(B);
}


INLINE static void radixsort_h4(DATA_TYPE *A, int elems) {
  DATA_TYPE *finlist, *B;

  B = (DATA_TYPE *)malloc(elems*sizeof(DATA_TYPE));
  assert_malloc(B);
  
  finlist = radixsort_h_2(A, B, elems, 8);
  if (finlist == B)
    bcopy(B,A,elems*sizeof(DATA_TYPE));

  free(B);
}


/*************************************************************************/
#define         INTEGER         1
#define		MTHRESH		6		/* threshold for median */

static  int		qsz;			/* size of each record */
static  int		(*qcmp)();		/* the comparison routine */

static  int		thresh;			/* THRESHold in chars */
static  int		mthresh;		/* MTHRESHold in chars */

/*
 * qsort:
 * First, set up some global parameters for qst to share.  Then, quicksort
 * with qst(), and then a cleanup insertion sort ourselves.  Sound simple?
 * It's not...
 */


#if !INTEGER
#define THRESH 8
#define COMPARE(i1,i2) (*qcmp)(i1,i2)
#define SWAPP(i1,i2) \
  {ii = qsz;do {c = *(i1);*(i1)++ = *(i2);*(i2)++ = c;} while(--ii);}
#else
#define	THRESH 24
#define COMPARE(i1,i2) (*(int *)(i1) - *(int *)(i2))
#define SWAPP(i1,i2) \
  {register int itmp;itmp=*(int *)(i1);*(int *)(i1)=*(int *)(i2);*(int *)(i2)=itmp;}
#endif

void qst();

void gnu_qsort (base, n, size, compar)
     char *base;
     int n;
     int size;
     int (*compar)();
{
  register char *i, *j, *lo, *hi, *min;
  /*  register int c,ii; */
  char *max;

  if (n <= 1)  return;
  qsz = size;
  qcmp = compar;
  thresh = qsz*THRESH;
  mthresh = qsz*MTHRESH;
  max = base + n*qsz;
  if (n >= THRESH)
    {
      qst (base, max);
      hi = base + thresh;
    }
  else
    {
      hi = max;
    }
  /*
   * First put smallest element, which must be in the first THRESH, in
   * the first position as a sentinel.  This is done just by searching
   * the first THRESH elements (or the first n if n < THRESH), finding
   * the min, and swapping it into the first position.
   */
  for (j = lo = base; (lo += qsz) < hi; )
    {
      if (COMPARE(j, lo) > 0)
	j = lo;
    }
  if (j != base)
    {			/* swap j into place */
      for (i = base, hi = base + qsz; i < hi;)
#if !INTEGER
	{
	  c = *j;
	  *j++ = *i;
	  *i++ = c;
	}
#else
      {
	SWAPP(i,j);
	i+=qsz;j+=qsz;
      }
#endif
    }
  /*
   * With our sentinel in place, we now run the following hyper-fast
   * insertion sort.  For each remaining element, min, from [1] to [n-1],
   * set hi to the index of the element AFTER which this one goes.
   * Then, do the standard insertion sort shift on a character at a time
   * basis for each element in the frob.
   */
  for (min = base; (hi = min += qsz) < max;)
    {
      while ( COMPARE(hi -= qsz, min) > 0);

      if ((hi += qsz) != min)
#if !INTEGER
	{
	  for (lo = min + qsz; --lo >= min;)
	    {
	      c = *lo;
	      for (i = j = lo; (j -= qsz) >= hi; i = j)
		*i = *j;
	      *i = c;
	    }
	}
#else
      {
	register int cI;
	register int *loI,*minI=(int *)min,*hiI=(int *)hi,*iI,*jI;
	for (loI = minI + 1; --loI >= minI;)
	    {
	      cI = *loI;
	      for (iI = jI = loI; (jI -= 1) >= hiI; iI = jI)
		*iI = *jI;
	      *iI = cI;
	    }
      }
#endif
    }
}

/*
 * qst:
 * Do a quicksort
 * First, find the median element, and put that one in the first place as the
 * discriminator.  (This "median" is just the median of the first, last and
 * middle elements).  (Using this median instead of the first element is a big
 * win).  Then, the usual partitioning/swapping, followed by moving the
 * discriminator into the right place.  Then, figure out the sizes of the two
 * partions, do the smaller one recursively and the larger one via a repeat of
 * this code.  Stopping when there are less than THRESH elements in a partition
 * and cleaning up with an insertion sort (in our caller) is a huge win.
 * All data swaps are done in-line, which is space-losing but time-saving.
 * (And there are only three places where this is done).
 */

void qst (base, max)
     char *base, *max;
{
  register char *i, *j, *jj, *mid;
  /* register int ii, c; */
  char *tmp;
  int lo, hi;

  lo = max - base;		/* number of elements as chars */
  do
    {
/*
 * At the top here, lo is the number of characters of elements in the
 * current partition.  (Which should be max - base).
 * Find the median of the first, last, and middle element and make that the
 * middle element.  Set j to largest of first and middle.  If max is larger
 * than that guy, then it's that guy, else compare max with loser of first
 * and take larger.  Things are set up to prefer the middle, then the first
 * in case of ties.
 */
      mid = i = base + qsz * ((lo/qsz) >> 1);
      if (lo >= mthresh)
	{
#if !INTEGER
	  j = (COMPARE((jj = base), i) > 0 ? jj : i);
	  if (COMPARE(j, (tmp = max - qsz)) > 0)
#else
	    jj = base;
	  j = (COMPARE(jj, i) > 0 ? jj : i);
	  tmp = max - qsz;
	  if (COMPARE(j, tmp) > 0)
#endif
	    {
	      j = (j == jj ? i : jj);	/* switch to first loser */
	      if (COMPARE(j, tmp) < 0)
		j = tmp;
	    }
	  if (j != i)
	    {
#if !INTEGER
	      ii = qsz;
	      do
		{
		  c = *i;
		  *i++ = *j;
		  *j++ = c;
		}
	      while(  --ii  );
#else
	      SWAPP(i,j);
	      i+=qsz;j+=qsz;
#endif
	    }
	}
      /*
       * Semi-standard quicksort partitioning/swapping
       */
      for (i = base, j = max - qsz; ;)
	{
	  while (i < mid && COMPARE(i, mid) <= 0)
	    i += qsz;
	  while (j > mid)
	    {
	      if (COMPARE(mid, j) <= 0)
		{
		  j -= qsz;
		  continue;
		}
	      tmp = i + qsz;		/* value of i after swap */
	      if (i == mid)
		{	/* j <-> mid, new mid is j */
		  mid = jj = j;
		}
	      else
		{			/* i <-> j */
		  jj = j;
		  j -= qsz;
		}
	      goto  swap;
	    }
	  if (i == mid)
	    {
	      break;
	    }
	  else
	    {				/* i <-> mid, new mid is i */
	      jj = mid;
	      tmp = mid = i;		/* value of i after swap */
	      j -= qsz;
	    }
	swap:
#if !INTEGER	  
	  ii = qsz;
	  do
	    {
	      c = *i;
	      *i++ = *jj;
	      *jj++ = c;
	    }
	  while (--ii);
#else
	  SWAPP(i,jj);
	  i+=qsz;jj+=qsz;
#endif
	  i = tmp;
	}
      /*
       * Look at sizes of the two partitions, do the smaller one first by
       * recursion, then do the larger one by making sure lo is its size,
       * base and max are update correctly, and branching back.
       * But only repeat (recursively or by branching) if the partition is
       * of at least size THRESH.
       */
      i = (j = mid) + qsz;
      if ((lo = j - base) <= (hi = max - i))
	{
	  if (lo >= thresh)
	    qst (base, j);
	  base = i;
	  lo = hi;
	}
      else
	{
	  if (hi >= thresh)
	    qst (i, max);
	  max = j;
	}
    }
  while (lo >= thresh);
}

/*************************************************************************/

INLINE static void quicksort_moret(int N, int* A)
{
  int SMALLFILE = 16;
  /* Sizes between 6 and 25 show less than 10% difference in running times */
  int SMALLFILEminus1 = 15;
  int STACKSIZE = 100;
  int X,temp,top,i,j,l,r,m;
  int *stack;

  /*  stack: array [1..STACKSIZE] of integer; */
  stack = (int *)malloc(STACKSIZE*sizeof(int));

  if (N <= SMALLFILE) {
    r = N;
    goto final; /* just do insertion sort */
  }

  top = 0; /* stack */
  l = 0;
  r = N-1;

 normal:
  /* l, m and r are all distinct */
  /* Put the median element in A[l] so that we can start the hole at
     the left end.  */
  m = (l+r)/2;
  if (A[m] > A[l]) {
    temp = A[l]; A[l] = A[m]; A[m] = temp; 
  }
  if (A[l] > A[r]) {
    temp = A[l]; A[l] = A[r]; A[r] = temp;
  }
  if (A[m] > A[l]) {
    temp = A[l]; A[l] = A[m]; A[m] = temp;
  }
  X = A[l]; /* Make a hole at A[l] */

  /* One of  i  and  j  points to the hole.  The other points to the last
     element of its portion of A[l..r].  Initially  i  points to the hole and
     j  points the last element of the portion that contains elements >= X */
  i = l;
  j = r;

 partition:
  j--;
  while (A[j] > X) j--;
  if (j > i) {
    A[i] = A[j]; /* Reestablish invariant by switching i and j */
  }
  else {
    /* j = i-1  or possibly j = i */
    A[i] = X;  /* Fill in the hole */
    i++;       /* Element in A[i] is correctly placed */
    goto recurse;
  }
  i++;
  while (A[i] < X) i++;
  if (j > i) {
    A[i] = A[j]; /* Reestablish invariant by switching i and j */
    goto partition;
  }
  else {
    /* i = j+1  or possibly j = i */
    A[j] = X;  /* Fill in the hole */
    j--;       /* Element in A[j] is correctly placed */
  }

recurse:
    /* Usually the element X is in neither part, but it is possible in rare
       circumstances that i-j = 1.  In this case no harm is done. */
    /* stack the larger problem and remove tail recursion to solve the
       smaller, giving up entirely if the remaining pieces are small */
    if (r-i > SMALLFILEminus1) {
       if (j-l > SMALLFILEminus1) {
          top = top+2;
          if (j-l > r-i) {
             stack[top-2] = l;
             stack[top-1] = j;
             l = i;
             goto normal;
	  }
          else {
             stack[top-2] = i;
             stack[top-1] = r;
             r = j;
             goto normal;
          }
       }
       else {
          /* left subproblem is small -- ignore it */
          l = i;
          goto normal;
       }
    }
    /* else right subproblem is small -- ignore it */
    else {
       if (j-l > SMALLFILEminus1) {
          r = j;
          goto normal;
       }
       else { /* retrieve an unsolved problem from the stack */
          if (top > 0) {
             l = stack[top-2];
             r = stack[top-1];
             top = top-2;
             goto normal;
          }
       }
    }
    /* The quicksort phase is finished */

    /* Find the smallest element in the file so that we can use sentinel
       during insertion sort.  It must be in the first chunk */
    r = SMALLFILE;

final:
    X = A[0];
    j = 0;
    for (i = 1; i < r; i++) {
      if (A[i] < X) {
         j = i;
	 /* moves := moves+1; */
         X = A[j];
      }
    }
    A[j] = A[0];
    A[0] = X;

    /* By placing smallest item first the basis case is two elements sorted */
    for (r = 2; r <N; r++) {
        X = A[r];
        i = r-1;
        while (A[i] > X) { /* guaranteed sentinel */
            A[i+1] = A[i];
            i--;
	}
	A[i+1] = X;
    }
    
    free(stack);
}
/*************************************************************************/

/* Debian "potato" libc6 (2.1.3-19) stdlib/qsort.c */
/* Copyright (C) 1991, 1992, 1996, 1997 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Written by Douglas C. Schmidt (schmidt@ics.uci.edu).

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with the GNU C Library; see the file COPYING.LIB.  If not,
   write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */


/* Byte-wise swap two items of size SIZE. */
#define DEBSWAP(a, b, size)						      \
  do									      \
    {									      \
      register size_t __size = (size);					      \
      register char *__a = (a), *__b = (b);				      \
      do								      \
	{								      \
	  char __tmp = *__a;						      \
	  *__a++ = *__b;						      \
	  *__b++ = __tmp;						      \
	} while (--__size > 0);						      \
    } while (0)

/* Discontinue quicksort algorithm when partition gets below this size.
   This particular magic number was chosen to work best on a Sun 4/260. */
#define DEB_MAX_THRESH 4

/* Stack node declarations used to store unfulfilled partition obligations. */
typedef struct
  {
    char *lo;
    char *hi;
  } stack_node;

/* The next 4 #defines implement a very fast in-line stack abstraction. */
#define DEB_STACK_SIZE	(8 * sizeof(unsigned long int))
#define DEB_PUSH(low, high)	((void) ((top->lo = (low)), (top->hi = (high)), ++top))
#define	DEB_POP(low, high)	((void) (--top, (low = top->lo), (high = top->hi)))
#define	DEB_STACK_NOT_EMPTY	(stack < top)


/* Order size using quicksort.  This implementation incorporates
   four optimizations discussed in Sedgewick:

   1. Non-recursive, using an explicit stack of pointer that store the
      next array partition to sort.  To save time, this maximum amount
      of space required to store an array of MAX_INT is allocated on the
      stack.  Assuming a 32-bit integer, this needs only 32 *
      sizeof(stack_node) == 136 bits.  Pretty cheap, actually.

   2. Chose the pivot element using a median-of-three decision tree.
      This reduces the probability of selecting a bad pivot value and
      eliminates certain extraneous comparisons.

   3. Only quicksorts TOTAL_ELEMS / DEB_MAX_THRESH partitions, leaving
      insertion sort to order the DEB_MAX_THRESH items within each partition.
      This is a big win, since insertion sort is faster for small, mostly
      sorted array segments.

   4. The larger of the two sub-partitions is always pushed onto the
      stack first, with the algorithm then concentrating on the
      smaller partition.  This *guarantees* no more than log (n)
      stack size is needed (actually O(1) in this case)!  */

void
debian_quicksort (pbase, total_elems, size, cmp)
     void *const pbase;
     size_t total_elems;
     size_t size;
     int (*cmp)();
     /*     __compar_fn_t cmp; */
{
  register char *base_ptr = (char *) pbase;

  /* Allocating SIZE bytes for a pivot buffer facilitates a better
     algorithm below since we can do comparisons directly on the pivot. */
#ifdef GCC
  char *pivot_buffer = (char *) alloca (size);
#else
  char *pivot_buffer = (char *) malloc (size);
#endif
  const size_t max_thresh = DEB_MAX_THRESH * size;

  if (total_elems == 0)
    /* Avoid lossage with unsigned arithmetic below.  */
    return;

  if (total_elems > DEB_MAX_THRESH)
    {
      char *lo = base_ptr;
      char *hi = &lo[size * (total_elems - 1)];
      /* Largest size needed for 32-bit int!!! */
      stack_node stack[DEB_STACK_SIZE];
      stack_node *top = stack + 1;

      while (DEB_STACK_NOT_EMPTY)
        {
          char *left_ptr;
          char *right_ptr;

	  char *pivot = pivot_buffer;

	  /* Select median value from among LO, MID, and HI. Rearrange
	     LO and HI so the three values are sorted. This lowers the
	     probability of picking a pathological pivot value and
	     skips a comparison for both the LEFT_PTR and RIGHT_PTR. */

	  char *mid = lo + size * ((hi - lo) / size >> 1);

	  if ((*cmp) ((void *) mid, (void *) lo) < 0)
	    DEBSWAP (mid, lo, size);
	  if ((*cmp) ((void *) hi, (void *) mid) < 0)
	    DEBSWAP (mid, hi, size);
	  else
	    goto jump_over;
	  if ((*cmp) ((void *) mid, (void *) lo) < 0)
	    DEBSWAP (mid, lo, size);
	jump_over:;
	  memcpy (pivot, mid, size);
	  pivot = pivot_buffer;

	  left_ptr  = lo + size;
	  right_ptr = hi - size;

	  /* Here's the famous ``collapse the walls'' section of quicksort.
	     Gotta like those tight inner loops!  They are the main reason
	     that this algorithm runs much faster than others. */
	  do
	    {
	      while ((*cmp) ((void *) left_ptr, (void *) pivot) < 0)
		left_ptr += size;

	      while ((*cmp) ((void *) pivot, (void *) right_ptr) < 0)
		right_ptr -= size;

	      if (left_ptr < right_ptr)
		{
		  DEBSWAP (left_ptr, right_ptr, size);
		  left_ptr += size;
		  right_ptr -= size;
		}
	      else if (left_ptr == right_ptr)
		{
		  left_ptr += size;
		  right_ptr -= size;
		  break;
		}
	    }
	  while (left_ptr <= right_ptr);

          /* Set up pointers for next iteration.  First determine whether
             left and right partitions are below the threshold size.  If so,
             ignore one or both.  Otherwise, push the larger partition's
             bounds on the stack and continue sorting the smaller one. */

          if ((size_t) (right_ptr - lo) <= max_thresh)
            {
              if ((size_t) (hi - left_ptr) <= max_thresh)
		/* Ignore both small partitions. */
                DEB_POP (lo, hi);
              else
		/* Ignore small left partition. */
                lo = left_ptr;
            }
          else if ((size_t) (hi - left_ptr) <= max_thresh)
	    /* Ignore small right partition. */
            hi = right_ptr;
          else if ((right_ptr - lo) > (hi - left_ptr))
            {
	      /* Push larger left partition indices. */
              DEB_PUSH (lo, right_ptr);
              lo = left_ptr;
            }
          else
            {
	      /* Push larger right partition indices. */
              DEB_PUSH (left_ptr, hi);
              hi = right_ptr;
            }
        }
    }

  /* Once the BASE_PTR array is partially sorted by quicksort the rest
     is completely sorted using insertion sort, since this is efficient
     for partitions below DEB_MAX_THRESH size. BASE_PTR points to the beginning
     of the array to sort, and END_PTR points at the very last element in
     the array (*not* one beyond it!). */

#define deb_min(x, y) ((x) < (y) ? (x) : (y))

  {
    char *const end_ptr = &base_ptr[size * (total_elems - 1)];
    char *tmp_ptr = base_ptr;
    char *thresh = deb_min(end_ptr, base_ptr + max_thresh);
    register char *run_ptr;

    /* Find smallest element in first threshold and place it at the
       array's beginning.  This is the smallest array element,
       and the operation speeds up insertion sort's inner loop. */

    for (run_ptr = tmp_ptr + size; run_ptr <= thresh; run_ptr += size)
      if ((*cmp) ((void *) run_ptr, (void *) tmp_ptr) < 0)
        tmp_ptr = run_ptr;

    if (tmp_ptr != base_ptr)
      DEBSWAP (tmp_ptr, base_ptr, size);

    /* Insertion sort, running from left-hand-side up to right-hand-side.  */

    run_ptr = base_ptr + size;
    while ((run_ptr += size) <= end_ptr)
      {
	tmp_ptr = run_ptr - size;
	while ((*cmp) ((void *) run_ptr, (void *) tmp_ptr) < 0)
	  tmp_ptr -= size;

	tmp_ptr += size;
        if (tmp_ptr != run_ptr)
          {
            char *trav;

	    trav = run_ptr + size;
	    while (--trav >= run_ptr)
              {
                char c = *trav;
                char *hi, *lo;

                for (hi = lo = trav; (lo -= size) >= tmp_ptr; hi = lo)
                  *hi = *lo;
                *hi = c;
              }
          }
      }
  }
#ifndef GCC
  free(pivot_buffer);
#endif
}

/*************************************************************************/

INLINE static void bubbleSort(int N, int* A)
{
  register int i, j, temp;

  for (i = 0; i < N - 1; i++) {
    for (j = 0; j < N - i - 1; j++) {
      if (A[j] > A[j + 1]) {
	temp = A[j];
	A[j] = A[j + 1];
	A[j + 1] = temp;
      }
    }
  }
}


/*************************************************************************/

INLINE static void selectionSort(int N, int* A)
{
  register int i, j, minIndex, temp;
  
  for (i = 0; i < N - 1; i++) {
    minIndex = i;
    for (j = i + 1; j < N; j++) {
      if (A[j] < A[minIndex])
	minIndex = j;
    }
    temp = A[minIndex];
    A[minIndex] = A[i];
    A[i] = temp;
  }
}

/*************************************************************************/

INLINE static void insertionSort2(int N, int* A)
{
  register int i, j, key;
  
  for (i = 1; i < N; i++) {
    key = A[i];
    j = i - 1;
    while (j >= 0 && A[j] > key) {
      A[j + 1] = A[j];
      j = j - 1;
    }
    A[j + 1] = key;
  }
}

 

/*************************************************************************/

void swap_3(int* a, int* b) {
  int temp = *a;
  *a = *b;
  *b = temp;
}

int partition(int* A, int low, int high) {
  register int i, j;
  register int pivot = A[high];

  i = (low - 1);

  for (j = low; j <= high - 1; j++) {
    if (A[j] <= pivot) {
      i++;
      swap_3(&A[i], &A[j]);
    }
  }
  swap_3(&A[i + 1], &A[high]);
  return (i + 1);
}


void quickSort_3(int* A, int low, int high) {
  register int pi;
  
  if (low < high) {
    pi = partition(A, low, high);
    quickSort_3(A, low, pi - 1);
    quickSort_3(A, pi + 1, high);
  }
}

INLINE static void run_quickSort_3(int N, int* A) {
  register int pi;
  
  pi = partition(A, 0, N-1);
  quickSort_3(A, 0, pi - 1);
  quickSort_3(A, pi + 1, N-1);
}


/*************************************************************************/

void create_input(DATA_TYPE *list, int elems) {
  int i;

  for (i=0 ; i<elems ; i++) {
#if 1
    list[i] = (DATA_TYPE)random();
#endif
#if 0
    list[i] = (DATA_TYPE)i;
#endif
#if 0
    list[i] = (DATA_TYPE)(i % DEFAULT_R);
#endif
#if 0
    list[i] = (DATA_TYPE)1;
#endif
  }

#if 0
  for (i=0 ; i<elems; i++)
    fprintf(stdout,"List[%5d]: %12d\n",i,list[i]);
#endif

}


void test_code() {
#if 0
  printf("(double)((int) 2*3500) f %f\n",(double)((int) 2*3500));
  printf("(double)( 2*3500) f %f\n",(double)( 2*3500));
  printf("(int) 2*3500 f %f\n",(int) 2*3500);
  printf("(int) 2*3500 d %d\n",(int) 2*3500);
  printf("2*3500 f %f\n",2*3500);
  printf("2*3500 d %d\n",2*3500);

  srandom(time(0));
  printf("(double)((int) 2*random()) f %f\n",(double)((int) 2*random()));
  printf("(double)( 2*random()) f %f\n",(double)( 2*random()));
  printf("(int) 2*random() f %f\n",(int) 2*random());
  printf("(int) 2*random() d %d\n",(int) 2*random());
  printf("2*random() f %f\n",2*random());
  printf("2*random() d %d\n",2*random());
#endif
}

int
main(int argc, char **argv) {
  DATA_TYPE 
    *originalList,
    *list;
  int 
    elems,
    loop;
  double 
    total_time,
    over_time;

  int err;

  if (argc <= 1)
    elems = DEFAULT_ELEMS;
  else
    elems = atoi(argv[1]);

  test_code();
  
  fprintf(stdout,"Sorting [%12d]\n",elems);
  
#if 0
  list = (DATA_TYPE*)malloc(elems*sizeof(DATA_TYPE));
  assert_malloc(list);

  originalList = (DATA_TYPE*)malloc(elems*sizeof(DATA_TYPE));
  assert_malloc(originalList);
#else
  
  list = (DATA_TYPE*)malloc(2*elems*sizeof(DATA_TYPE));
  assert_malloc(list);

  originalList = (DATA_TYPE*)malloc(2*elems*sizeof(DATA_TYPE));
  assert_malloc(originalList);

#endif
  
  srandom(time(0));

  create_input(originalList, elems);

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    radixsort_3(list,elems);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with radixsort_3\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t  radixsort3: %f\n",
	  elems,total_time);

/******************************************************************/

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    radixsort_3_rev(list,elems);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with radixsort_3_rev\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t  radixs3rev: %f\n",
	  elems,total_time);

/******************************************************************/

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    radixsort_4(list,elems);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with radixsort_4\n");


  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t  radixsort4: %f\n",
	  elems,total_time);

/******************************************************************/

#if 0
/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    radixsort_4p(list,elems);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with radixsort_4p\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t  radixsrt4p: %f\n",
	  elems,total_time);

/******************************************************************/
#endif

#if DO_ISORT
/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    insertsort(list,elems);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with insertsort\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t  insertsort: %f\n",
	  elems,total_time);

/******************************************************************/
#endif
 
/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    qsort(list,elems,sizeof(int),intcompare);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with qsort()\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t       qsort: %f\n",
	  elems,total_time);

/******************************************************************/

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    mergesort(list,elems);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with mergesort\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t   mergesort: %f\n",
	  elems,total_time);

/******************************************************************/

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    mergesort_nr(list,elems);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with mergesort_nr\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t    msort_nr: %f\n",
	  elems,total_time);

/******************************************************************/

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    radixsort_h3(list,elems);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with radixsort_h3\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t   rsort_h_3: %f\n",
	  elems,total_time);

/******************************************************************/


/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    radixsort_h4(list,elems);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with radixsort_h4\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t   rsort_h_4: %f\n",
	  elems,total_time);

/******************************************************************/

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    gnu_qsort(list,elems,sizeof(int),intcompare);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with gnu_qsort()\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t   gnu_qsort: %f\n",
	  elems,total_time);

/******************************************************************/
  
/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    quicksort_moret(elems, list);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with quicksort_moret\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t   qsort-mor: %f\n",
	  elems,total_time);

/******************************************************************/
  
/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    debian_quicksort(list,elems,sizeof(int),intcompare);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with debian_quicksort\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t   deb-qsort: %f\n",
	  elems,total_time);

/******************************************************************/

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    bubbleSort(elems,list);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with bubbleSort\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t  bubblesort: %f\n",
	  elems,total_time);

/******************************************************************/

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    selectionSort(elems,list);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with selectionSort\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t  selectsort: %f\n",
	  elems,total_time);

/******************************************************************/

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    insertionSort2(elems,list);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with insertionSort2\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t insertsort2: %f\n",
	  elems,total_time);

/******************************************************************/

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
    run_quickSort_3(elems,list);
  }
  total_time = get_seconds() - total_time;
  err = check_sort(list,elems);
  if (!err) fprintf(stderr,"ERROR with quickSort3\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    bcopy(originalList,list,elems*sizeof(DATA_TYPE));
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," elems: %12d \t  quicksort3: %f\n",
	  elems,total_time);

/******************************************************************/

  free(originalList);
  free(list);
  return(0);
}


