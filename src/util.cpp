#include <ctype.h>
#include <string.h>

#include "util.h"

int bode::splits(char *str,char **dest,int max) {
  int index = 0;
  char *start = str;
  if (str == NULL) {
    return 0;
  }
  while (*str != '\0' && index < max) {
    if (isspace(*str)) {
      dest[index++] = start;
      start = str+1;
      if (index < max) {
        *str = '\0';
      }
    }
    str++;
  }
  if (index < max) {
    dest[index++] = start;
  }
  return index;
}
