//////////////////////////////////////////////////////////////////
/*
   getopt.c -- 
      getopt is used to break up options in command lines for easy parsing
      by shell procedures and to check for legal options.  optstring is a
      string of recognized option letters (see getopt(3C)).  If a letter is
      followed by a colon, the option is expected to have an argument which
      may or may not be separated from it by white space.

      The positional parameters ($1 $2 ...) of the shell are reset so that
      each option is preceded by a - and is in its own positional parameter;
      each option argument is also parsed into its own positional parameter.

      getopt recognizes two hyphens (--) to delimit the end of the options.
      If absent, getopt places -- at the end of the options.

      The most common use of getopt is in the shell's set command (see the
      example below) where getopt converts the command line to a more easily
      parsed form.  getopt writes the modified command line to the standard
      output.
*/
//
//   
//
//  Version:    
//
//  History:    Downloaded Sept 14, 1998
//
//  Compilers:  gcc
//
//  Source: Christopher Fraser and David R. Hanson,
//          A Retargetable C Compiler -- Design and Implementation
//          ISBN: 0-8053-1670-1
//          Addison-Wesley, 1995.
//          http://www.cs.princeton.edu/software/lcc/
// 
////////////////////////////////////////////////////////////
#if defined _WIN32 || defined _WIN64

#include        <stdio.h>
#include        <string.h>

#define EPR                 fprintf(stderr,
#define ERR(str, chr)       if(opterr){EPR "%s%c\n", str, chr);}
/* [<][>][^][v][top][bottom][index][help] */
int     opterr = 1;
int     optind = 1;
int     optopt;
char    *optarg;
// char    *strchr();

int
getopt (int argc, char *const argv[], const char *opts)
/* [<][>][^][v][top][bottom][index][help] */
{
        static int sp = 1;
        register int c;
        register char *cp;
        
        if (sp == 1) {
                if (optind >= argc ||
                   argv[optind][0] != '-' || argv[optind][1] == '\0')
                        return -1;
                else if (strcmp(argv[optind], "--") == 0) {
                        optind++;
                        return -1;
                }
        }
        optopt = c = argv[optind][sp];
        if (c == ':' || (cp=(char *)strchr(opts, c)) == 0) {
                ERR (": illegal option -- ", c);
                if (argv[optind][++sp] == '\0') {
                        optind++;
                        sp = 1;
                }
                return '?';
        }
        if (*++cp == ':') {
                if (argv[optind][sp+1] != '\0')
                        optarg = &argv[optind++][sp+1];
                else if (++optind >= argc) {
                        ERR (": option requires an argument -- ", c);
                        sp = 1;
                        return '?';
                } else
                        optarg = argv[optind++];
                sp = 1;
        } else {
                if (argv[optind][++sp] == '\0') {
                        sp = 1;
                        optind++;
                }
                optarg = 0;
        }
        return c;
}

#endif