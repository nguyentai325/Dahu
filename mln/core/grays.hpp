#ifndef MLN_CORE_GRAYS_HH
# define MLN_CORE_GRAYS_HH

# include <mln/core/value/int.hpp>


/**
   \file This file provides typedefs for commons buitins types
*/

namespace mln
{

  enum MLN_C_GRAYS {
    MLN_UINT8 = 0x01,	///< Constant for an unsigned 8-bits integer.
    MLN_INT8 = 0x02,	///< Constant for an signed 8-bits integer.
    MLN_UINT16 = 0x03,	///< Constant for an unsigned 16-bits integer.
    MLN_INT16 = 0x04, 	///< Constant for an signed 16-bits integer.
    MLN_UINT32 = 0x05,	///< Constant for an unsigned 32-bits integer.
    MLN_INT32 = 0x06,	///< Constant for an signed 32-bits integer.
    MLN_UINT64 = 0x07,
    MLN_INT64 = 0x08,
    MLN_BOOL = 0x09,
    MLN_FLOAT = 0x10,
    MLN_DOUBLE = 0x11
  };

}

#endif /* !MLN_CORE_GRAYS_HH */
