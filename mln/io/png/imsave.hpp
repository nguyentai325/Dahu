#ifndef MLN_IO_PNG_IMSAVE_HPP
# define MLN_IO_PNG_IMSAVE_HPP

# include <FreeImage.h>
# include <mln/core/grays.hpp>
# include <boost/static_assert.hpp>


namespace mln
{

  namespace io
  {

    namespace png
    {

    namespace internal
    {

      template <typename Image>
      void imsave(const Image& ima, const char* path, raw_image_tag)
      {
        BOOST_STATIC_ASSERT((Image::ndim == 2));
        typedef mln_value(Image) V;

        FreeImage_Initialise();
        FREE_IMAGE_TYPE fit;

	if (std::is_same<V, bool>::value or std::is_same<V, uint8>::value or std::is_same<V, rgb8>::value)
	  fit = FIT_BITMAP;
	else if (std::is_same<V, uint16>::value)
	  fit = FIT_UINT16;
	else
	  throw MLNIOException("Unhandled image type.");


        int h = ima.domain().shape()[0], w = ima.domain().shape()[1];
        FIBITMAP* dib;
	if (std::is_same<V, rgb8>::value)
	  dib = FreeImage_AllocateT(fit, w, h, sizeof(V) * 8, 0xFF000000, 0x00FF00, 0xFF000000);
	else if (std::is_same<V, bool>::value) {
	  dib = FreeImage_AllocateT(fit, w, h, 1);
	  RGBQUAD *pal = FreeImage_GetPalette(dib);
	  mln_assertion(pal != NULL);
	  pal[1].rgbRed = 0xFF;
	  pal[1].rgbGreen = 0XFF;
	  pal[1].rgbBlue = 0xFF;
	  mln_assertion(FreeImage_GetColorType(dib) == FIC_MINISBLACK);
	} else
	  dib = FreeImage_AllocateT(fit, w, h, sizeof(V) * 8);


        const size_t* strides = ima.strides();
        const char* ptr = (const char*) &ima(ima.domain().pmin);
	unsigned bpp = sizeof(V);
        for(int y = h-1; y >= 0; --y, ptr += strides[0]) {
          BYTE *bits = FreeImage_GetScanLine(dib, y);
	  if (std::is_same<V, rgb8>::value)
	    {
	      for (int j = 0; j < w; ++j) {
		*(bits + j * bpp + 0) = *(ptr + j * bpp + 2);
		*(bits + j * bpp + 1) = *(ptr + j * bpp + 1);
		*(bits + j * bpp + 2) = *(ptr + j * bpp + 0);
	      }
	    }
	  else if (std::is_same<V, bool>::value)
	    {
	      for (int j = 0, z = 0; j < w; ++z)
		for (int b = 7; b >= 0; --b, ++j)
		  bits[z] |= ptr[j] << b;
	    }
	  else
	    std::memcpy(bits, ptr, w * sizeof(V));
        }

        //FreeImage_DeInitialise();
        bool res = FreeImage_Save(FIF_PNG, dib, path);
	if (!res)
	  throw MLNIOException("Unable to save the image.");
        FreeImage_Unload(dib);
      }
    }

      template <typename Image>
      void imsave(const Image& ima, const char* path)
      {
	internal::imsave(ima, path, typename image_traits<Image>::category ());
      }
    }


  } // end of namespace mln::io
} // end of namespace mln



#endif //!MLN_IO_PNG_IMSAVE_HPP
