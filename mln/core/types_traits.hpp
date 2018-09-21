#ifndef MLN_CORE_TYPES_TRAITS_HPP
# define MLN_CORE_TYPES_TRAITS_HPP

# define REGISTER_TYPE(T, tname, ID)		\
  template <>					\
  type_info<T> get_info_from_type<T>()		\
  {						\
    static type_info<T> info = {tname, ID};	\
    return info					\
  }						\
						\
  template <>					\
  type_info<T> get_info_from_id<ID>()		\
  {						\
    static type_info<T> info = {tname, ID};	\
    return info;				\
  }

namespace mln
{
  struct type_info_base
  {
    virtual const char* name() const = 0;
    virtual int id() const = 0;
  };

  template <typename T>
  struct type_info : type_info_base
  {
    typedef T	type;

    virtual const char* name() const { return name_; }
    virtual int id() const { return id_; }

    const char* name_;
    int		id_;
  };

  template <typename T>
  type_info<T> get_info_from_type()
  {
    static type_info<void> info;
    return info;
  }


  template <int>
  type_info<void> get_info_from_id()
  {
    static type_info<void> info;
    return info;
  };

  REGISTER_TYPE(void, "unkwown", 0)

}

#endif
