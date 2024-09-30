/*
  * Copyright Ali Ghaffaari.
  * SPDX-License-Identifier: MIT
  *
  * Licensed under the MIT License. See the LICENSE accompanying this file
  * for the specific language governing permissions and limitations under
  * the License.
  */


#ifndef  KSEQPP_SEQIO_HPP__
#define  KSEQPP_SEQIO_HPP__

#include <zlib.h>

#include "kseq++.hpp"

namespace klibpp {
  class SeqStreamIn
    : public KStreamIn< gzFile, int(*)(gzFile_s*, void*, unsigned int) > {
    public:
      /* Typedefs */
      typedef KStreamIn< gzFile, int(*)(gzFile_s*, void*, unsigned int) > base_type;
      /* Lifecycle */
      SeqStreamIn( const char* filename )
        : base_type( gzopen( filename, "r" ), gzread, gzclose )
      { }

      SeqStreamIn( int fd )
        : base_type( gzdopen( fd, "r" ), gzread, gzclose )
      { }
  };

  class SeqStreamOut
    : public KStreamOut< gzFile, int(*)(gzFile_s*, const void*, unsigned int) > {
    public:
      /* Typedefs */
      typedef KStreamOut< gzFile, int(*)(gzFile_s*, const void*, unsigned int) > base_type;
      /* Lifecycle */
      SeqStreamOut( const char* filename, bool compressed=false,
                    format::Format fmt=base_type::DEFAULT_FORMAT )
        : base_type( gzopen( filename, ( compressed ? "w" : "wT" ) ), gzwrite, fmt, gzclose )
      { }

      SeqStreamOut( int fd, bool compressed=false,
                    format::Format fmt=base_type::DEFAULT_FORMAT )
        : base_type( gzdopen( fd, ( compressed ? "w" : "wT" ) ), gzwrite, fmt, gzclose )
      { }

      SeqStreamOut( const char* filename, format::Format fmt )
        : base_type( gzopen( filename, "wT" ), gzwrite, fmt, gzclose )
      { }

      SeqStreamOut( int fd, format::Format fmt )
        : base_type( gzdopen( fd, "wT" ), gzwrite, fmt, gzclose )
      { }
  };
}  /* -----  end of namespace klibpp  ----- */
#endif  /* ----- #ifndef KSEQPP_SEQIO_HPP__  ----- */
