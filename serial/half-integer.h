#ifndef HALF_INTEGER_H
#define HALF_INTEGER_H
/*
  A Simple class to handle arithmetic of half integers.  It should interoperate relatively seamlessly with integers and doubles; but note that operations between a half-integer and a double will return a half-integer, so there may be (a lot of) precision loss if that isn't what you want.
 */

const double TOL = 0.01;

class HalfInteger {
 public:
  HalfInteger() { doubleValue = 0; }
  HalfInteger( int value ) { doubleValue = 2*value; }
  HalfInteger( double value ) { 
    if( value >= 0 )
      doubleValue = (int) (2*value+TOL);
    else
      doubleValue = (int) (2*value-TOL);
  }
 
  HalfInteger( int dv1, int dv2 ) { doubleValue = dv1+dv2; } // to help with addition

  HalfInteger operator+(const HalfInteger &h2) const { return HalfInteger(doubleValue, h2.doubleValue); }
  HalfInteger operator-(const HalfInteger &h2) const { return HalfInteger(doubleValue, -1*h2.doubleValue); }
  HalfInteger operator-() const { return HalfInteger(0,-1*doubleValue); }
  bool operator==(const HalfInteger &h2) const { return (doubleValue == h2.doubleValue ); }
  bool operator!=(const HalfInteger &h2) const { return (doubleValue != h2.doubleValue ); }
  bool operator<(const HalfInteger &h2 ) const { return (doubleValue < h2.doubleValue ); }
  bool operator>(const HalfInteger &h2 ) const { return (doubleValue > h2.doubleValue ); }
  bool operator<=(const HalfInteger &h2 ) const { return (doubleValue <= h2.doubleValue ); }
  bool operator>=(const HalfInteger &h2 ) const { return (doubleValue >= h2.doubleValue ); }
  void operator+=(const HalfInteger &h2) { doubleValue += h2.doubleValue; }
  void operator-=(const HalfInteger &h2) { doubleValue -= h2.doubleValue; }
  void operator++(int) { doubleValue+=2; }
  void operator++() { doubleValue+=2; }
  void operator--(int) { doubleValue-=2; }
  void operator--() { doubleValue-=2; }

  // Don't do arithmetic on the output of this, that's why it is called print
  double print() const {
    return doubleValue / 2.;
  }
 private:
  int doubleValue;
};

#endif
