#ifndef MLN_COLORS_LAB_HPP
# define MLN_COLORS_LAB_HPP

# include <mln/core/vec_base.hpp>
# include <mln/core/colors.hpp>
# include <mln/colors/xyz.hpp>
# include <cmath>
// FIXME: optimize this out (slow because of floats and saturations)

namespace mln
{

  struct lab_tag {};

  template <typename T>
  using lab = internal::vec_base<T, 3, lab_tag>;

  typedef lab<uint8> lab8;

  namespace internal
  {
    template <>
    struct vec_base_traits<lab_tag>
    {
      static const bool is_additive = true;
      static const bool is_additive_ext = true;
      static const bool is_multiplicative = false;
      static const bool is_multiplicative_ext = true;
      static const bool is_less_than_comparable = true;
      static const bool is_equality_comparable = true;
    };

  }

  template <typename T>
  lab<float> rgb2lab(const rgb<T>& v);

  template <typename T>
  rgb<T> rgb2lab_normalize(const lab<T>& v);

  template <typename T>
  lab<float> rgb2lab_minh(const rgb<T>& v);

  template <typename T>
  rgb<T> lab2rgb(const lab<T>& v);


  /*********************/
  /*** Implementation **/
  /*********************/

  template <typename T>
  inline
  lab<float>
  rgb2lab(const rgb<T>& v)
  {
    auto xyz = rgb2xyz(v);

    constexpr float xn = (0.4125 + 0.3576 + 0.1804) * value_traits<T>::max();
    constexpr float yn = (0.2127 + 0.7152 + 0.0722) * value_traits<T>::max();
    constexpr float zn = (0.0193 + 0.1192 + 0.9502) * value_traits<T>::max();

    auto f = [] (float t) {
      return (t > .008856f) ? std::cbrt(t) : (t * 7.787f + 16.0f/116.0f);
    };

    float xr = f(xyz[0] / xn);
    float yr = f(xyz[1] / yn);
    float zr = f(xyz[2] / zn);

    float l = 116 * yr - 16;
    float a = 500 * (xr - yr);
    float b = 200 * (yr - zr);

    return {l, a, b};
  }


  template <typename T>
  inline
  lab<float>
  rgb2lab_minh(const rgb<T>& v)
  {

      float
      R = float(v[0])   / 255.0,
      G = float(v[1]) / 255.0,
      B = float(v[2])  / 255.0;

      float X, Y,Z;
      float r1,g1,b1;
      float epsilon = 0.008856;	//actual CIE standard
      float kappa   = 903.3;		//actual CIE standard
      float Xr = 0.950456;	//reference white
      float Yr = 1.0;		//reference white
      float Zr = 1.088754;	//reference white
      double xr,yr,zr;
      double fx, fy, fz;

      if(R <= 0.04045)	r1 = R/12.92;
      else				r1 = pow((R+0.055)/1.055,2.4);
      if(G <= 0.04045)	g1 = G/12.92;
      else				g1 = pow((G+0.055)/1.055,2.4);
      if(B <= 0.04045)	b1 = B/12.92;
      else				b1 = pow((B+0.055)/1.055,2.4);


      X = r1*0.4124564 + g1*0.3575761 + b1*0.1804375;
      Y = r1*0.2126729 + g1*0.7151522 + b1*0.0721750;
      Z = r1*0.0193339 + g1*0.1191920 + b1*0.9503041;

      //------------------------
      // XYZ to LAB conversion
      //------------------------
      xr = X/Xr;
      yr = Y/Yr;
      zr = Z/Zr;
      if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
      else				fx = (kappa*xr + 16.0)/116.0;
      if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
      else				fy = (kappa*yr + 16.0)/116.0;
      if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
      else				fz = (kappa*zr + 16.0)/116.0;


      float l = 116.0*fy-16.0;
      float a = 500.0*(fx-fy);
      float b = 200.0*(fy-fz);



      return {l, a, b};
  }
  
  
    template <typename T>
  inline
  lab<float>
  rgb2lab_normalize(const rgb<T>& v)
  {
      float
      R = float(v[0])   / 255.0,
      G = float(v[1]) / 255.0,
      B = float(v[2])  / 255.0;

      float X, Y,Z;
      float r1,g1,b1;
      float epsilon = 0.008856;	//actual CIE standard
      float kappa   = 903.3;		//actual CIE standard
      float Xr = 0.950456;	//reference white
      float Yr = 1.0;		//reference white
      float Zr = 1.088754;	//reference white
      double xr,yr,zr;
      double fx, fy, fz;

      if(R <= 0.04045)	r1 = R/12.92;
      else				r1 = pow((R+0.055)/1.055,2.4);
      if(G <= 0.04045)	g1 = G/12.92;
      else				g1 = pow((G+0.055)/1.055,2.4);
      if(B <= 0.04045)	b1 = B/12.92;
      else				b1 = pow((B+0.055)/1.055,2.4);


      X = r1*0.4124564 + g1*0.3575761 + b1*0.1804375;
      Y = r1*0.2126729 + g1*0.7151522 + b1*0.0721750;
      Z = r1*0.0193339 + g1*0.1191920 + b1*0.9503041;

      //------------------------
      // XYZ to LAB conversion
      //------------------------
      xr = X/Xr;
      yr = Y/Yr;
      zr = Z/Zr;
      if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
      else				fx = (kappa*xr + 16.0)/116.0;
      if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
      else				fy = (kappa*yr + 16.0)/116.0;
      if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
      else				fz = (kappa*zr + 16.0)/116.0;


      float l = 116.0*fy-16.0;
      l = (l/100*255 );

      float a = 500.0*(fx-fy);
      a = (a+86.1735)*255/184.4197;

      float b = 200.0*(fy-fz);
      b = (b + 107.874) *255 / 202.353;

      return {l, a, b};
  }


  template <typename T>
  inline
  rgb<T>
  lab2rgb(const lab<T>& v)
  {
    constexpr float xn = (0.4125 + 0.3576 + 0.1804) * value_traits<T>::max();
    constexpr float yn = (0.2127 + 0.7152 + 0.0722) * value_traits<T>::max();
    constexpr float zn = (0.0193 + 0.1192 + 0.9502) * value_traits<T>::max();

    auto f = [] (float t) {
      return (t > .2068965f) ? std::pow(t,3) : .12842f * (t - 4.0f /29.0f);
    };

    float K = .00862069f * (v[0] + 16.0f);
    float Y = yn * f(K);
    float X = xn * f(K + 0.002 * v[1]);
    float Z = zn * f(K - 0.005 * v[2]);

    return xyz2rgb(xyz<float>{X, Y, Z});
  }

}

#endif // ! MLN_COLORS_LAB_HPP
