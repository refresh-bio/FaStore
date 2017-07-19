#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>

#define LINE_MAX 10000

// Input arguments:
// 1: Original FASTQ file
// 2: File with new quals
// 3: Output file

int main(int argc, char *argv[])
{
  FILE  *f1, *f2, *fout;
  char line[LINE_MAX];

  f1 = fopen(argv[1],"r");
  f2 = fopen(argv[2],"r");
  fout = fopen(argv[3], "w");

  while (fgets(line,LINE_MAX, f1)!=NULL){
      fprintf(fout,"%s",line);
      
      fgets(line,LINE_MAX,f1);
      fprintf(fout,"%s",line);

      fgets(line,LINE_MAX,f1);
      fprintf(fout,"%s",line);
      
      fgets(line,LINE_MAX,f1);

      fgets(line,LINE_MAX,f2);
      fprintf(fout,"%s",line);
  }
  fclose(f1);
  fclose(f2);
  fclose(fout);
  return 0;
}
