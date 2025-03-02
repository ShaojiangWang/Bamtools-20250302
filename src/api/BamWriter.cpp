// ***************************************************************************
// BamWriter.cpp (c) 2009 Michael Str�mberg, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#include "api/BamAlignment.h"
#include "api/BamWriter.h"
using namespace BamTools;
using namespace std;

/*! \class BamTools::BamWriter
    \brief Provides write access for generating BAM files.
*/
/*! \enum BamTools::BamWriter::CompressionMode
    \brief This enum describes the compression behaviors for output BAM files.
*/
/*! \var BamWriter::CompressionMode BamWriter::Compressed
    \brief Use normal BAM compression
*/
/*! \var BamWriter::CompressionMode BamWriter::Uncompressed
    \brief Disable BAM compression

    Useful in situations where the BAM data is streamed (e.g. piping).
    It would be wasteful to compress, and then immediately decompress
    the data.
*/

BamWriter::~BamWriter(void) {
   d->Close();
    delete d;
    d = nullptr;
}

/*! \fn BamWriter::Close(void)
    \brief Closes the current BAM file.
    \sa Open()
*/
void BamWriter::Close(void) {
    d->Close();
}

/*! \fn std::string BamWriter::GetErrorString(void) const
    \brief Returns a human-readable description of the last error that occurred

    This method allows elimination of STDERR pollution. Developers of client code
    may choose how the messages are displayed to the user, if at all.

    \return error description
*/
std::string BamWriter::GetErrorString(void) const {
    return d->GetErrorString();
}

/*! \fn bool BamWriter::IsOpen(void) const
    \brief Returns \c true if BAM file is open for writing.
    \sa Open()
*/
bool BamWriter::IsOpen(void) const {
    return d->IsOpen();
}

/*! \fn void BamWriter::SetCompressionMode(const BamWriter::CompressionMode& compressionMode)
    \brief Sets the output compression mode.

    Default mode is BamWriter::Compressed.

    \note Changing the compression mode is disabled on open files (i.e. the request will
    be ignored). Be sure to call this function before opening the BAM file.

    \code
        BamWriter writer;
        writer.SetCompressionMode(BamWriter::Uncompressed);
        writer.Open( ... );
        // ...
    \endcode

    \param[in] compressionMode desired output compression behavior
    \sa IsOpen(), Open()
*/
void BamWriter::SetCompressionMode(const BamWriter::CompressionMode& compressionMode) {
    d->SetWriteCompressed( compressionMode == BamWriter::Compressed );
}
