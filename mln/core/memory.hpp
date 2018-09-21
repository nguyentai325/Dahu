#ifndef MLN_MEMORY_HH
# define MLN_MEMORY_HH

# include <cstdlib>
# include <mln/core/assert.hpp>

namespace mln
{


	/*!
	 * \brief
	 * Allocate \p n bytes of data such that the \p offset th byte of data is 128-bits aligned.
	 * 
	 * \param n Size of the buffer to allocate
	 * \param offset The \p offset th byte to be 128-bits aligned 
	 * 
	 * \returns Pointer to the beginning of the buffer
	 *
	 * \see aligned_free
	 */
	inline
	void*
	aligned_malloc(size_t n, size_t offset)
	{
          char* buffer = (char*)std::malloc(n + 16);
          if (!buffer)
            return buffer;
          char* ptr = (char*)(((size_t)(buffer + offset) & ~((size_t)15)) + 16) - offset;
          char* adr = (char*) ptr - 1;
          *adr = (ptr - buffer);
          mln_assertion(adr >= buffer);
          mln_assertion(ptr >= buffer);
          mln_assertion((ptr + n) <= (buffer + n + 16));
          return ptr;
	}

	inline
	void
	aligned_free(void* ptr)
	{
	  if (ptr) {
	    char* p = (char*) ptr - 1;
            std::free((char*)ptr - *p);
	  }
	}

}


#endif
