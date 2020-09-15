#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// gcc -c -fpic my_ctype.c 
// gcc -shared -o lib_my_ctype.so my_ctype.o 
//
// import cytpes
// my_ctyp = ctypes.CDLL("./lib_my_ctype.so")
// my_ctyp.my_fput("this is a test\n","test.dat")

void my_fput(char *str, char *filename)
{
    FILE *fp = fopen(filename, "w");
    fputs(str, fp);
    fclose(fp);
}

void my_update(int *a, int *b)
{
	*a += 5;
	*b += 6;
}
