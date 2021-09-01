#include <stdio.h>
#include <inttypes.h>
#include "io_util.h"

/**
 * \brief Subroutine that reads an (optionally) segmented block of data written by a Fortran program
 *
 * Please note that the file must be opened in read mode and enough memory be allocated for buf.
 * This implies you have to know beforehand the size of your data (the same happens in Fortran).
 *
 * \param f	The file object.
 * \param buf0	Pointer to where the data will be written.
 * \param ts	Size of the data type to be read (this is needed for byte-swapping).
 * \param len	Number of elements that will be read.
 * \param skip	Number of elements to skip since the beginning.
 * \param swap	Indicates if swapping should be performed.
 *
 * \return Number of elements actually read.
 *
 * See http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/fortran/lin/compiler_f/bldaps_for/common/bldaps_rectypes.htm for a descriotion of Variable-Length Records
 *
 * Francisco Martinez-Serrano, 2008, 2010
 */

#define MAX(A,B)        ((A)>(B)?(A):(B))
#define MIN(A,B)        ((A)<(B)?(A):(B))
#define ABS(A)        ((A)< 0?-(A):(A))

extern int io_util_readfortran(FILE * const f, void * const buf0, size_t const ts, size_t const len, size_t const skip, int const swap)
{
	int32_t s1,s1p,s2,s2p; //Fortran heading and trailing record sizes
	size_t i,read=0,readt=0,remaining=len*ts,toRead,toSkip,remainingSkip=skip*ts; //counter
	unsigned char * buf= (unsigned char *) buf0;

	do
	{
		if (fread(&s1,sizeof(int32_t),1,f) != 1) return readt;
		//		printf("record size: %d ",s1);
		s1p=ABS(s1);
		if(swap) io_util_sexchange(&s1,sizeof(int32_t));
		toRead = s1p;
		toSkip = MIN(remainingSkip,toRead);
		//		printf("will skip: %ld ",toSkip);
		fseek(f,toSkip,SEEK_CUR);
		remainingSkip -= toSkip;
		toRead -= toSkip;
		toRead = MIN(toRead,remaining);
		//		printf("will read: %ld ",toRead);
		if(toRead > 0){
			if (fread(buf,1,toRead,f) != toRead) return readt;
			buf += toRead;
			read += toRead;
			readt=read/ts;
			remaining -= toRead;
		}
		//		printf("will skip %ld until the end. ",(s1p-(toSkip+toRead)));
		fseek(f,(s1p-(toSkip+toRead)),SEEK_CUR); //Skip to the end of the current slice
		if(fread(&s2,sizeof(int32_t),1,f)!=1) return readt;
		s2p=ABS(s2);
		//		printf("record size: %d ",s2);
		if(swap) io_util_sexchange(&s2,sizeof(int32_t));
		if(s1p != s2p){
			printf("Something went wrong: |%d| != |%d|",s1,s2);
			return readt; //Sanity check
		}

	} while (s1 < 0 /*&& toRead > 0*/);

	if(swap){
		buf=buf0;
		for(i=0;i<len;i++) 
		{
			io_util_sexchange(buf,ts);
			buf += ts;
		}
	}

	if(readt!=len)
		printf("read != len, %ld != %ld\n",readt,len);
	return readt;

}
