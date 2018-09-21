#ifndef MLN_CORE_ITERATOR_ITERATOR_BASE_HPP
# define MLN_CORE_ITERATOR_ITERATOR_BASE_HPP

# include <mln/core/concept/iterator.hpp>
# include <mln/core/iterator/iterator_traits.hpp>
# include <memory>
# include <type_traits>

namespace mln
{

  struct use_default {};

  namespace internal
  {

    template <typename Reference>
    struct pointer_wrapper
    {
      explicit pointer_wrapper(Reference x) : x_ (x) {}
      Reference* operator->() { return std::addressof(x_); }
      operator Reference* () { return std::addressof(x_); }
    private:
      Reference x_;
    };


    template <typename T>
    struct make_pointer
    {
      typedef pointer_wrapper<T> type;
      static pointer_wrapper<T> foo(T x) { return pointer_wrapper<T> (x); }
    };

    template <typename T>
    struct make_pointer<T&>
    {
      typedef T* type;
      static T* foo(T& x) { return std::addressof (x); }
    };

  };


  ///
  /// Helper class for iterators
  ///
  template <typename Derived, typename Value,
            typename Reference = Value&>
  struct iterator_base : Iterator<Derived>
  {
    typedef Derived iterator;
    typedef Derived const_iterator;
    typedef std::false_type has_NL;

    typedef typename std::remove_const<Value>::type value_type;
    typedef Reference reference;
    typedef typename std::conditional< std::is_reference<Reference>::value,
				       typename std::remove_reference<Reference>::type*,
				       internal::pointer_wrapper<Reference> >::type pointer;

    Derived& iter()
    {
      return *(this->derived());
    }

    Derived iter() const
    {
      return *(this->derived());
    }


    reference
    operator* () const
    {
      // It may not be necessary to start iteration to have valid
      // member (e.g. the image pointer in pixels), thus we should
      // assert this:
      // mln_precondition(not this->derived()->finished());
      return this->derived()->dereference();
    }

    pointer
    operator-> () const
    {
      return internal::make_pointer<reference>::foo(this->derived()->dereference());
    }

    void set_dejavu_(bool)
    {
    }

  private:
    const Derived* derived() const
    {
      return reinterpret_cast<const Derived*>(this);
    }

    Derived* derived()
    {
      return reinterpret_cast<Derived*>(this);
    }

  };

} // end of namespace mln

#endif //!MLN_CORE_ITERATOR_ITERATOR_BASE_HPP
