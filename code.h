
#define MCQLI24 1  // used by ICE
//#define RJ1 1
//#define RJ2 1
//#define BJ24 1
//#define MJ 1
//#define MCQLI24 1 // k=24 r=1/2 code for ICE
//#define JQLIODP48 1
//#define MCQLI48 1
//#define OT24 1
//#define QR24 1
//#define LL 1
//#define J60 1
//#define J50 1

//#define MCQLI32 1 // NASA standard
//#define JP24 1
//#define JSODP47 1
//#define BLLF47 1

// Convolutional coding polynomials. All are rate 1/2

#if	defined(MCQLI32)
//"NASA standard" code by Massey & Costello
// Shortened from MCQLI-48:
// Nonsystematic, quick look-in, dmin=11, dfree=23
// used on Pioneer 10-12, Helios A,B
#define	POLY1	0xbbef6bb7      // 111 011 011 101 011 011 110 111 110 111 010 = 73353367672
#define	POLY2	0xbbef6bb5      // 101 011 011 101 011 011 110 111 110 111 010 = 53353367672
#define K       32
#define G1FLIP  0
#define G2FLIP  0

#elif	defined(MJ)
// Massey-Johannesson code
// Nonsystematic, quick look-in, dmin=13, dfree>=23
// Purported to be more computationally efficient than Massey-Costello
// Appears related to  1JQLIODP-24
#define	POLY1	0xb840a20f      // 111 100 000 100 010 100 000 010 000 111 010 = 74042402072
#define POLY2	0xb840a20d      // 101 100 000 100 010 100 000 010 000 111 010 = 54042402072
#define K       32
#define G1FLIP  0
#define G2FLIP  0

#elif	defined(LL)
// Layland-Lushbaugh code
// Nonsystematic, non-quick look-in, dmin=?, dfree=?
#define	POLY1	0xf2d05351      // 100 010 101 100 101 000 001 011 010 011 110 = 42545013236
#define	POLY2	0xe4613c47      // 111 000 100 011 110 010 000 110 001 001 110 = 70436206116
#define K       32
#define G1FLIP  0
#define G2FLIP  0


#elif defined(MCQLI24)
// MCQLI-24
// k=24 r=1/2 Massey QLI code for ISEE-3/International Comet Explorer
// The ref gives the polynomials as 073353367 and 053353367 (octal)
// Shortened from MCQLI-48
#define POLY1   073665667
#define POLY2   073665665
#define K       24
#define G1FLIP  0
#define G2FLIP  1           // Invert the second symbol

// 1JQLIODP-24
// Johannassen QLI ODP k=24
#elif defined(RJ1)
#define POLY1   074121017
#define POLY2   074121015
#define K       24
#define G1FLIP  0
#define G2FLIP  0

// 2JQLIODP-24
// Johannassen QLI ODP k=24
#elif defined(RJ2)
#define POLY1   073541017
#define POLY2   073541015
#define K       24
#define G1FLIP  0
#define G2FLIP  0

// BJ-24
// Bahl-Jelinek k=24 complementary code
#elif defined(BJ24)
#define POLY1   054220245
#define POLY2   063557533
#define K       24
#define G1FLIP  0
#define G2FLIP  0

// QR-24 code
// Massey, Costello, Justesen Quadratic k=24 residue code
#elif defined(QR24)
#define POLY1   026241177
#define POLY2   037620515
#define K       24
#define G1FLIP  0
#define G2FLIP  0

// OT-24
// Massey k=24 optimally truncatable code
#elif defined(OT24)
#define POLY1   062650457
#define POLY2   062650455
#define K       24
#define G1FLIP  0
#define G2FLIP  0

// MCQLI-48
// Massey-Costelli QLI code k=48
#elif defined(MCQLI48)
#define POLY1   06556767373665667LL
#define POLY2   06556767373665665LL
#define K       48
#define G1FLIP  0
#define G2FLIP  0

// JQLIODP-48
// Johannesson QLI ODP code with k=48
#elif defined(JQLIODP48)
#define POLY1    05634247020121017LL
#define POLY2    05634247020121015LL
#define K        48
#define G1FLIP   0
#define G2FLIP   0

//BLLF-47
// Bussgang, Lin, Lynn, Forney k=47
// Systematic
#elif defined(BLLF47)
#define POLY1   1
#define POLY2   0531746407671547LL
#define K       45
#define G1FLIP  0
#define G2FLIP  0

// JSODP-47
// Johannesson's ODP k=47
// Systematic
#elif defined(JSODP47)
#define POLY1    1
#define POLY2    03331355751514473LL
#define K        47
#define G1FLIP   0
#define G2FLIP   0

// JP-24
// Mentioned in addendum to Massey
// Non-QLI, non-systematic
#elif defined (JP24)
#define POLY1    052431655
#define POLY2    061411757
#define K        24
#define G1FLIP   0
#define G2FLIP   0

// Monster k=60 Johannesson systematic
#elif defined (J60)
#define POLY1    1
#define POLY2    073607331355751514473LL
#define K        60
#define G1FLIP   0
#define G2FLIP   0

// Monster k=50 Johannesson QLI
#elif defined (J50)
#define POLY1    075634247020121017
#define POLY2    075634247020121015
#define K        50
#define G1FLIP   0
#define G2FLIP   0


#endif


int encode(
   unsigned char *symbols,	// Output buffer, 2*8*nbytes
   const unsigned char *data,	// Input buffer, nbytes
   unsigned int nbytes);	// Number of bytes in data
