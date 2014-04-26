#ifndef IO_UTIL_H
#define IO_UTIL_H

/* $Id: io_util.h,v 1.5 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file util.h
 *
 * Provides utility function for input/output operations.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdbool.h>

#include "io_file.h"


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/**
 * \brief Reads a possibly byteswapped float from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readfloat(f,n,b) \
  io_util_readbytes((f), (void *)(n), sizeof(float), (b))

/**
 * \brief Reads a possibly byteswapped double from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readdouble(f,n,b) \
  io_util_readbytes((f), (void *)(n), sizeof(double), (b))

/**
 * \brief Reads a possibly byteswapped int from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readint(f,n,b) \
  io_util_readbytes((f), (void *)(n), sizeof(int), (b))

/**
 * \brief Reads a possibly byteswapped long from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readlong(f,n,b) \
  io_util_readbytes((f), (void *)(n), sizeof(long), (b))

/**
 * \brief Reads a possibly byteswapped unsigned long from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readulong(f,n,b) \
  io_util_readbytes((f), (void *)(n), sizeof(unsigned long), (b))

/**
 * \brief Reads a possibly byteswapped long long from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readlonglong(f,n,b) \
  io_util_readbytes((f), (void*)(n), sizeof(long long), (b))

/**
 * \brief Reads a possibly byteswapped int8 from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readint8(f,n,b) \
  io_util_readbytes((f), (void*)(n), sizeof(int8_t), (b))

/**
 * \brief Reads a possibly byteswapped int16 from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readint16(f,n,b) \
  io_util_readbytes((f), (void*)(n), sizeof(int16_t), (b))

/**
 * \brief Reads a possibly byteswapped int32 from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readint32(f,n,b) \
  io_util_readbytes((f), (void*)(n), sizeof(int32_t), (b))

/**
 * \brief Reads a possibly byteswapped int64 from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readint64(f,n,b) \
  io_util_readbytes((f), (void*)(n), sizeof(int64_t), (b))

/**
 * \brief Reads a possibly byteswapped uint8 from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readuint8(f,n,b) \
  io_util_readbytes((f), (void*)(n), sizeof(uint8_t), (b))

/**
 * \brief Reads a possibly byteswapped uint16 from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readuint16(f,n,b) \
  io_util_readbytes((f), (void*)(n), sizeof(uint16_t), (b))

/**
 * \brief Reads a possibly byteswapped uint32 from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readuint32(f,n,b) \
  io_util_readbytes((f), (void*)(n), sizeof(uint32_t), (b))

/**
 * \brief Reads a possibly byteswapped uint64 from a stream
 *
 * \param f  The stream from which to read.
 * \param n  A pointer to the variabel which is supposed to hold the
 *           result. 
 * \param b  A toggle to choose between byteswapped or
 *           non-byteswapped data.
 *
 * \return Returns 1 if the value could be read, 0 otherwise.
 */
#define io_util_readuint64(f,n,b) \
  io_util_readbytes((f), (void*)(n), sizeof(uint64_t), (b))


/**********************************************************************\
 *    Prototypes of global functions                                  * 
\**********************************************************************/

/**
 * \brief BRIEF DESCRIPTION
 *
 * LONG DESCRIPTION
 *
 * \param *p  The pointer to the variable that should be byteswapped.
 * \param s   The number of bytes of the input variable.
 *
 * \return Nothing.
 */
extern void
io_util_sexchange(void *p, size_t s);

/**
 * \brief Reads in a possibly byteswapped numerical value
 *
 * \param *f    A pointer in the stream where the value starts
 * \param *n    A pointer to a variable supposed to hold the read value
 * \param len   The number of bytes to be read
 * \param swap  Toggles between swapped reading and unswapped 
 *              reading
 *
 * \return If the value could be read, 1 is returned, 0 otherwise.
 */
extern int
io_util_readbytes(FILE *f, void *n, size_t len, io_file_swap_t swap);

/**
 * \brief Reads in a fixed length string.
 *
 * \param *f  The stream to read from.
 * \param *s  The external string that'll hold the result.
 * \param n   The number of characters to read.
 *
 * \note The String s points to must have enough space to accomodate the
 *       n bytes that'll be read plus the trailing '\0', hence the
 *       string must be at least n+1 bytes long.
 *
 * \return  Returns the number of bytes read.
 */
extern size_t
io_util_readstring(FILE *f, char *s, size_t n);

/**
 * \brief Reads in a fixed length string.
 * 
 * Like readstring, but will also stop at linebreaks.
 *
 * \param *f  The stream to read from.
 * \param *s  The external string that'll hold the result.
 * \param n   The number of characters to read.
 *
 * \note The String s points to must have enough space to accomodate the
 *       n bytes that'll be read plus the trailing '\0', hence the
 *       string must be at least n+1 bytes long.
 *
 * \return  Returns the number of bytes read.
 */
extern size_t
io_util_readline(FILE *f, char *s, size_t n);

/**
 * \brief Function to replicate the behaviour of strdup which is
 *        not a C99 function.
 *
 * \param *str  The string to duplicate.
 *
 * \return Duplicated (including final \0) string.
 */
extern char *
io_util_strdup(const char *str);

/**
 * \brief Splits a string into a path and a filename part.
 *
 * \param *str    The string to split. Will not be touched.
 * \param **path  Pointer to the variable supposed to hold the path.
 * \param **stem  Pointer to the variable supposed to hold the stem.
 *
 * \return Returns a pointer to the path or NULL in case of errors.
 */
extern char *
io_util_split_pathfname(const char *str, char **path, char **stem);

/**
 * \brief Finds all files in a given directory matching a given stem by
 *        adding numbers.
 *
 * \param *path      Where to look.
 * \param *stem      What to look for.
 * \param *format    How the numbers are formatted, e.g. '%04i'.
 * \param *suffic    What to append to get a correct filename.
 * \param ***fnames  Will hold the resulting filename array.
 *
 * \return Returns the number of files found.
 */
extern int32_t
io_util_findfiles(const char *path,
                  const char *stem,
                  const char *format,
                  const char *suffix,
                  char ***fnames);


#endif /* IO_UTIL_H */
