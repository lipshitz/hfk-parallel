#include <stdio.h>
#include "half-integer.h"

int main( int argc, char **argv ) {
  HalfInteger a = 1;
  HalfInteger b(1.5);
  HalfInteger c = a + b;
  HalfInteger d = a - b;
  HalfInteger e = d - a;
  a++;
  printf("%f\n", a.print());
  printf("%f\n", b.print());
  printf("%f\n", c.print());
  printf("%f\n", d.print());
  printf("%f\n", e.print());
  printf("%d\n", a == (a+b)-b);
}
