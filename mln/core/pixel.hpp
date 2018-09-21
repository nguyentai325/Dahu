#ifndef MLN_CORE_PIXEL_HPP
# define MLN_CORE_PIXEL_HPP

namespace mln {

  // template <typename Pixel>
  // struct pixel_iterator
  // {
  //   typedef typename Pixel::point_type       point_type;
  //   typedef typename Pixel::value_type       value_type;

  // protected:
  //   void set_point(const point_type& p) { pixel_.p_ = p; }
  //   void set_value(value_type* v) { pixel_.v_ = v; }


  //   Pixel pixel_;
  // };


  template <typename Point, typename Value>
  struct pixel
  {
    typedef Point       point_type;
    typedef Value       value_type;

    Point point() const { return p_; }
    Value& val() const { return *v_; }


    //friend class pixel_iterator< pixel<Point, Value> >;

  private:
    Point p_;
    Value* v_;
  };

} // end of namespace mln

#endif //!MLN_CORE_PIXEL_HPP
