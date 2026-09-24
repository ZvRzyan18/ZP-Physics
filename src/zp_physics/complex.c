#include "zp_physics/complex.h"

#ifdef zp_has_complex_extension
zp_complex zp_complex_identity = 1.0f + 0.0fi;
#else
zp_complex zp_complex_identity = {.x = 1.0f, .y = 0.0f};
#endif


/*
 src : https://en.wikipedia.org/wiki/Euler%27s_formula
 
 the idea is based on the Euler's formula 

  e ^ i*θ = cos(θ) + i sin(θ)
 
  e = euler's number (≈2.71828)
  i = imaginary (√-1 or (complex){0, 1})
 
 
 when you use logarithm, you can express it as
 inverse trigonometry

  ln(e ^ i*radians) where you can extract the radians as imaginary part


 Multiplication :
 
  multiplication of a unit complex numbers means adding the
  radians angle.
  
  similar to exponent identity :
  (a ^ b) * (a ^ c) = a ^ (b + c)
  
  for complex :
  (e ^ iA) * (e ^ iB) = e ^ (iA + iB)
  the angle of both complex numbers (iA and iB imaginary angles) accumulate.


 example C++:

  double theta = M_PI_2;
  std::complex<double> i(0, 1);
  std::complex<double> unit_vector = std::exp(i * theta);
  std::complex<double> inverse = std::log(unit_vector);

  std::cout << unit_vector << '\n';
  std::cout << inverse << '\n';

*/


/*
 Minimal Version of Euler's Formula works
 by avoiding the expensive exponential and logarithm
 computation and exploit the fact that e ^ 0 == 0 
 and ln(1) == 0. both of them can be implicit and cancels
 out each other.
*/


/* 
 complex rotate uses a simplifcation that eliminates the 
 exponential.

 complex exponential :
  input (a, b)
  output : (e ^ a * cos(b), e ^ a * sin(b))
 
 idea :
  e ^ i*θ = cos(θ) + i sin(θ)
 
 simplification : 
  since (e ^ 0 = 1.0) you can completely eliminate the exponential
              leaving only the trigonometry part of imaginary
*/
zp_complex zp_crotate(const float radians) {
 zp_vec2 out;
 zp_sincos(radians, &out.y, &out.x);
 return zp_as_complex2(out);
}

/*
 complex angle uses simplification (same tricks, to eliminate the logarithm)
 complex logaritm:
 
 
  input (input_complex)(a, b)
  
            ln(dot(input_complex, input_complex))
  output : (_____________________________________,  atan2(b, a))
                             2.0

  output : (ln(a * a + b * b) * 0.5,  atan2(b, a))
  
  idea :
  
   ln(e ^ i*θ) = iθ
 
 simplification : 
  earlier, in this formula (e ^ i*θ), we input the angle_theta as imaginary (i*θ)
  so we expect to return it as imaginary. because of that, we can ignore the
  logarithm part completely.
*/
float zp_cangle(const zp_complex c) {
 return zp_atan2(zp_cimag(c), zp_creal(c));
}









zp_quaternion zp_quaternion_identity = {.x = 0.0f, .y = 0.0f, .z = 0.0f, .w = 1.0f};


zp_quaternion zp_qfromaxisangle(const float a, const float x, const float y, const float z) {
 zp_quaternion out;
 float s;
 zp_sincos(a * 0.5f, &s, &out.w);
	out.x = x * s;
 out.y = y * s;
	out.z = z * s;
	return out;
}

zp_vec4 zp_qtoaxisangle(const zp_quaternion a) {
 zp_vec4 out;
 out.w = zp_asin(a.w);
 out.w += out.w;
 const float b = zp_rsqrt(1.0f - a.w * a.w);
 out.x = a.x * b;
 out.y = a.y * b;
 out.z = a.z * b;
 return out;
}




