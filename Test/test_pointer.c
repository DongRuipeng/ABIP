#include <stdio.h>

void main()
{
  int *pn = 0;
  if (pn)
  {
    printf("the value of null pointer is not zero !\n");
  }
  else
  {
    printf("the value of null pointer is zero.\n");
  }
}