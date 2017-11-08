/*
*
* This software is being provided to you, the LICENSEE, by the Massachusetts
* Institute of Technology (M.I.T.) under the following license.  By
* obtaining, using and/or copying this software, you agree that you have
* read, understood, and will comply with these terms and conditions:  
*
* Permission to use, copy, modify and distribute this software and its
* documentation for any purpose and without fee or royalty is hereby granted,
* provided that you agree to comply with the following copyright notice and
* statements, including the disclaimer, and that the same appear on ALL
* copies of the software and documentation, including modifications that you
* make for internal use or for distribution:
*
* Copyright 1992 by the Massachusetts Institute of Technology.  All rights
* reserved.  
*
* THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. MAKES NO REPRESENTATIONS OR
* WARRANTIES, EXPRESS OR IMPLIED.  By way of example, but not limitation,
* M.I.T. MAKES NO REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS
* FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE OR
* DOCUMENTATION WILL NOT INFRINGE ANY THIRD PARTY PATENTS, COPYRIGHTS,
* TRADEMARKS OR OTHER RIGHTS.   
*
* The name of the Massachusetts Institute of Technology or M.I.T. may NOT
* be used in advertising or publicity pertaining to distribution of the
* software.  Title to copyright in this software and any associated
* documentation shall at all times remain with M.I.T., and USER agrees to
* preserve same.  
*
* Written by: K. Nabors, T. Korsmeyer, and J. White
*
*/

/* # ***** sort to /src/misc
   # ***** */
/* 
 memory allocator for fastcap/fastlap
 - almost identical to Kernigan & Ritchie sec 8.7
 - differs in that there is no free list since fastcap never frees any memory
 - also the amount of memory sbrk()'d is exactly equal (within a unit) to 
   the requested amount if its more than NALLOC units
   this cuts down on unused blocks usually added to the free list
   while still keeping the number of sbrk()'s down 
 - the regular allocation functions are still availible for stdio to use
   although free() is rewritten to do nothing - stdio buffers are passed
   over when new memory is set up
 - machine dependancy should be confined to the type `ALIGN' and
   `MORECORE()' (i.e. non-UNIX machines won't have brk(), sbrk())
 - no attempt is made to make allocation efficient in terms of virtual pages
*/
#include <stdio.h>

#define NALLOC 8192		/* >= sizeof(HEADER)*NALLOC bytes sbrk()'d */
#define MAGICN 0xaaaaaaaaL	/* used to check fidelity of allocated blks */
#define UGDEBG 0		/* 0=> most efficient, no verification
				   1=> ualloc_verify() checks magicno, size
				   2=> ualloc_verify() prints lots to stdout
				   - can reduce efficiency due to big HEADER */

                                /* cftk: This is the old align definition, but 
                                   alignment on Snakes was always a problem.
				   The definition below should cure that, and
                                   one wonders if it will cure problems with
                                   free-ing memory on repeated calls for 
                                   some applications.
typedef int *ALIGN;		   all allocated mem starts on an ALIGN bdry */
typedef double ALIGN;		/* all allocated mem starts on an ALIGN bdry */

union header {			/* allocated block header - all mem allocated
				   in units of sizeof(HEADER) bytes */
  struct {
#if UGDEBG == 1 || UGDEBG == 2
    long magicno;		/* to check for overwrites */
    union header *ptr;		/* next allocated block */
#endif
#if UGDEBG == 2
    union header *enlarged;	/* limit of newly enlarged block */
    union header *enlgdfrom;	/* base pointer of newly added block */
    unsigned int osize;		/* size before enlargement/truncation */
    unsigned int request;      	/* nbytes requested from allocator */
#endif
    unsigned int size;		/* block size (memory + 1 or 2 header units) */
  } s;
  ALIGN x;			/* for alignment only, never referenced */
};

typedef union header HEADER;

/* 
  uses calloc to get a new block of memory
  - any header space put in by calloc is wasted
  - zeros memory
  - an alternative to mocore() but should only be used if sbrk() doesnt zero
*/
#define MORECORE(SIZE) (HEADER *)calloc(1, SIZE*sizeof(HEADER))
char *calloc();
char *malloc();

static HEADER *base = NULL;    	/* base of allocated block list */
static HEADER *allocp = NULL;	/* last allocated block */
static unsigned int sizeofHDR = sizeof(HEADER);

/*
  Restarts the storage allocator. 
*/
void restartStore()
{
extern long memQ2M;
extern long memQ2L;
extern long memQ2P;
extern long memL2L;
extern long memM2M;
extern long memM2L;
extern long memM2P;
extern long memL2P;
extern long memMSC;

  brk(base);
  base = NULL;
  allocp = NULL;

memQ2M = 0L;
memQ2L = 0L;
memQ2P = 0L;
memL2L = 0L;
memM2M = 0L;
memM2L = 0L;
memM2P = 0L;
memL2P = 0L;
memMSC = 0L;
}

/*
  asks operating system for more memory which is added to the top block
  - memory not zeroed out
*/
static HEADER *mocore(nu)
unsigned int nu;
{
  char *sbrk();
  char *cp;

  cp = sbrk(nu*sizeofHDR);
  if((int)cp == -1) return(NULL);
  return((HEADER *)cp);
}

/* 
  ugly storage allocator 
  - no frees
  - allocates in blocks at least NALLOC*sizeof(HEADER) bytes long
  - ultimately uses mocore(), since no frees are done (sbrk() zeros added
    memory) this allocator performs like calloc() w/no explicit assigns to 0
*/
char *ualloc(nbytes)
unsigned int nbytes;
{
  HEADER *mocore();
  HEADER *p, *q;
  int nunits;			/* size in number of sizeof(HEADER)'s */
  int brkunits;			/* number to add to heap */
#if UGDEBG == 2
  HEADER *dummy = 0x0;
  HEADER *dummyret = 0x0;
#endif
#if UGDEBG == 0
  static unsigned int topblksize;  	/* size of current top block */
#endif

  if(nbytes == 0) return(NULL);	/* should probably return something else */

#if UGDEBG == 1 || UGDEBG == 2
  nunits = 3+(nbytes-1)/sizeofHDR; /* rm for 2 hdrs too */
#else
  nunits = 2+(nbytes-1)/sizeofHDR; /* rm for 1 hdr too */
#endif

  if((q = allocp) == NULL) {	/* no allocation yet */
    if((q = allocp = base = mocore(1)) == NULL) return(NULL);
    topblksize = 1;
#if UGDEBG == 1 || UGDEBG == 2
    base->s.size = 1;
    base->s.ptr = base;
    base->s.magicno = MAGICN;
#endif
  }

  /* 
    check previously allocated block for room
    - if it's big enough, use head (not tail) part of it
    - if it's not, break out more memory and discard the previous block
  */
  if(topblksize >= nunits) {	/* if it's big enough */
#if UGDEBG == 1 || UGDEBG == 2
    p = q;			/* copy old top block pointer */
    q += (nunits - 1);	/* make q new top block pointer */
    q->s.size = p->s.size - nunits + 1; /* upper block size */
#else
    allocp = q + nunits - 1;
#endif
    topblksize -= (nunits - 1);
#if UGDEBG == 2
    p->s.osize = p->s.size;	/* save size before new block carved out */
    p->s.request = nbytes;	/* actual #bytes requested */
    p->s.enlarged = dummy;
    p->s.enlgdfrom = dummyret;
#endif
#if UGDEBG == 1 || UGDEBG == 2
    q->s.ptr = base;	/* top block connected back to base */
    p->s.ptr = q;		/* old top block pointer to new top block */
    if(allocp != p) allocp->s.ptr = p;	/* link previous old block */
    q->s.magicno = MAGICN;
    p->s.size = nunits;	/* set up returned block size */
    allocp = q;		/* save the top block pointer */
    return((char *) (p+1));	/* return pntr to first unit beyond header */
#else
    return((char *) q);	/* use old header location */
#endif
  }
  else {			/* get more memory, add to top block beyond 
				   any stdio buffers made since last ualloc */
    brkunits = (nunits > NALLOC) ? nunits : NALLOC;
    if((q = mocore(brkunits)) == NULL) return(NULL);
#if UGDEBG == 1 || UGDEBG == 2
    q->s.size = brkunits;	/* setup size, etc. of new top block */
    q->s.magicno = MAGICN;
#endif
    topblksize = brkunits;
#if UGDEBG == 2
    dummy = q + q->s.size - 1; /* limit pointer of new block */
    dummyret = q;		/* base pointer of new block */
#endif
    /* repeat above code to aviod having a loop */
#if UGDEBG == 1 || UGDEBG == 2
    p = q;
    q += (nunits - 1);	/* make q new top block pointer */
    q->s.size = p->s.size - nunits + 1; /* upper block size */
#else
    allocp = q + nunits - 1;
#endif
    topblksize -= (nunits - 1);
#if UGDEBG == 2
    p->s.osize = p->s.size;	/* save size before new block carved out */
    p->s.request = nbytes;	/* actual #bytes requested */
    p->s.enlarged = dummy;
    p->s.enlgdfrom = dummyret;
#endif
#if UGDEBG == 1 || UGDEBG == 2
    q->s.ptr = base;	/* top block connected back to base */
    p->s.ptr = q;		/* old top block pointer to new top block */
    if(allocp != p) allocp->s.ptr = p;	/* link previous old block */
    q->s.magicno = MAGICN;
    p->s.size = nunits;	/* set up returned block size */
    allocp = q;		/* save the top block pointer */
    return((char *) (p+1));	/* return pntr to first unit beyond header */
#else
    return((char *) q);	/* use old header location */
#endif
  }
}

/*
  checks for overwrites at end of blocks
  - checks magic numbers and sizes (UGDEBG == 1 or 2)
  - checks if length corresponds to pointers (UGDEBG == 1 or 2)
  - prints information about list of allocated blocks (UGDEBG == 2)
*/
void ualloc_verify()
{
  HEADER *p;
  int cnt = 1;

#if UGDEBG == 1 || UGDEBG == 2
  for(p = base;; p = p->s.ptr, cnt++) {
    if(p == base && cnt > 1) break;
#endif
#if UGDEBG == 2
    fprintf(stdout, 
	   "%d 0x%x 0x%x %u %u bytes (osize %u enlarged from 0x%x to 0x%x)\n", 
	    cnt, p, p->s.ptr, p->s.size, p->s.request,
	    p->s.osize, p->s.enlgdfrom, p->s.enlarged);
#endif
#if UGDEBG == 1 || UGDEBG == 2
    if(p->s.magicno != MAGICN) {
      fflush(stdout);
      fprintf(stderr, "ualloc_verify: bad block %d magic number\n", cnt);
      fflush(stderr);
      abort();
    }
    if(p->s.ptr - p + 1 < p->s.size) {
      fflush(stdout);
      fprintf(stderr, "ualloc_verify: bad block size ", cnt);
      fprintf(stderr, "(from 0x%x to 0x%x size %u)\n", p, p->s.ptr, p->s.size);
      fflush(stderr);
      abort();
    }
  }
#endif
}

/*
  compares the total address range broken out to amount requested
  - efficiency is 100*memcount/(sbrk(0)-base) = 100*requested/(ttl broken out)
  - UGDEBG == 2 prints waste values: 
      ualloc waste = 100*(mem discarded)/(ttl broken out)
      stdio waste = 100*(mem discarded to apease stdio usage)/(ttl broken out)
       (stdio waste is often zero and ualloc waste does not include the
        memory lost in each header struct (a problem with many small things)
  - if base == NULL (not using ugly allocator), final break value is printed
*/
void uallocEfcy(memcount)
long memcount;
{
#if UGDEBG == 2
  HEADER *p;
  unsigned int waste = 0;
  unsigned int stdiowaste = 0;
  int first = 1;
#endif
  int total;
  char *sbrk();

  total = (int)(sbrk(0) - (char *)base);

  if(base == NULL) fprintf(stdout, "(top of memory = 0x%x", sbrk(0));
  else fprintf(stdout, "(%.3g%% efficiency",
	       100*((double)memcount)/((double)total));

#if UGDEBG == 2
  /* loop through blocks, get total btyes wasted (by allocator) */
  for(p = base;; p = p->s.ptr) {
    if(p == base && first == 0) break;
    if(p == NULL) break;	/* compatability when ualloc not used */
    if(p->s.request == 0) {
      if((p->s.ptr)->s.size >= p->s.size) { /* if block dropped bec. too sm. */
	if(p != allocp) waste += p->s.size - 1; 
	else waste += p->s.size;
      }
      else {
	if(p != allocp) stdiowaste += p->s.size - 1; 
      }	
    }
    first = 0;
  }
  if(p != NULL) 
      fprintf(stdout, 
	      ", waste: %.2g%% (ualloc), %.2g%% (stdio)", 
	      100*(double)(waste*sizeof(HEADER))/((double)total),
	      100*(double)(stdiowaste*sizeof(HEADER))/((double)total));
#endif
  fprintf(stdout, ")\n");

}

/*
  so that stdio free list won't get in the way

void free(ptr)
char *ptr;
{}
*/

char *whatsit(num, size)
unsigned num, size;
{
char *foo;
  foo = calloc(num,size);
  return(foo);
}
