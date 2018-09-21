#ifndef MLN_CORE_IREF_HPP
# define MLN_CORE_IREF_HPP

# include <memory>
# include <mln/core/assert.hpp>

namespace mln
{

  /// \brief Implement an image reference wrapper.
  ///
  /// iref is:
  /// * default constructible
  /// * copy constructible
  /// * copy assignable
  template <class I>
  struct iref;


  /// \brief Specialization for a reference to temporary image
  /// The temporary image is copied and destroyed when iref is destroyed.
  template <class I>
  struct iref<I&&>;

  /********************************/
  /**   Implementation          ***/
  /********************************/

  template <class I>
  struct iref<I&>
  {
    iref() : m_ptr(nullptr) {}
    iref(I& x) : m_ptr(&x) {}
    iref(I&& x) = delete;

    I& get() {
      mln_precondition(m_ptr != nullptr);
      return *m_ptr;
    }

    const I& get() const {
      mln_precondition(m_ptr != nullptr);
      return *m_ptr;
    }

  private:
    I* m_ptr;
  };

  template <class I>
  struct iref<I&&>
  {
    iref() = default;

    template <class U>
    iref(U&& x) : m_x(std::forward<U>(x)) {}

    I& get() { return m_x; }
    const I& get() const { return m_x; }

  private:
    I m_x;
  };

}

#endif // !MLN_CORE_IREF_HPP
