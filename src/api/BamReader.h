// ***************************************************************************
// BamReader.h (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 18 November 2012 (DB)
// ---------------------------------------------------------------------------
// Provides read access to BAM files.
// ***************************************************************************

#ifndef BAMREADER_H
#define BAMREADER_H

#include "api/api_global.h"
#include "api/BamAlignment.h"
#include "api/BamIndex.h"
#include "api/SamHeader.h"
#include <string>
#include "api/internal/bam/BamReader_p.h"

using namespace std;

namespace BamTools {
//using namespace BamTools;

namespace Internal {
    class BamReaderPrivate;
} // namespace Internal
  
using namespace BamTools::Internal;
//using namespace BamTools::Internal;

/**
 * Got design problems for this class.
 * TODO: Should get rid of this class in new design.
 *       Use more C++ design patterns.
 *
 * BamRegion is defined in BamAux.h
 * BamAux contains several data structures shared by
 * the API implementations and interfaces.
 */
class API_EXPORT BamReader {
    // constructor / destructor
    public:
        /**
        * Default constructor.
        * There is no other constructors declared but could be
        * implemented by the compiler implicitly.
        * This method creates a new pointer.
        * So copy constructor and assignment operators
        * should all be disabled otherwise will get a double
        * delete probem.
        *
        * Note: We cannot inline this function because of the
        * cross mentioning of each other.
        */
        BamReader();

        /** 
         *  Destructor
         */
        ~BamReader(void);

    // public interface
    public:

        // ----------------------
        // BAM file operations
        // ----------------------

        /** 
         *  Closes the current BAM file.
         *
         *  Also clears out all header and reference data.
         *
         *  @return true if file closed OK
         *  @see IsOpen(), Open()
        */
        bool Close(void);
         bool close(void) {
            return d->Close();
         }
        /** 
         *  Returns name of current BAM file.
         *
         *  Retrieved filename will contain whatever was passed via Open().
         *  If you need full directory paths here, be sure to include them
         *  when you open the BAM file.
         *
         *  @return name of opened BAM file. If no file is open, returns an empty string.
         *  @see IsOpen()
         *
         *  TODO: should remove the const return type
        */
        const std::string& GetFilename(void) const;
        const std::string& getFilename(void) const {
           return GetFilename();
        }
        /** 
         *  @return true if a BAM file is open for reading.
         *  An BamReader may be in its unopen state because
         *  it closed the connection to the file or has not
         *  started connecting with a file.
        */
        bool IsOpen(void) const;
        bool isOpen(void) const {
           return IsOpen();
        }

        // returns true if a BAM file is open for reading
        /** 
         *  Performs a random-access jump within BAM file.
         *
         *  This is a convenience method, equivalent to calling SetRegion()
         *  with only a left boundary specified.
         *
         *  @param[in] refID    left-bound reference ID
         *  @param[in] position left-bound position
         *
         *  @return true if jump was successful
         *  @see HasIndex()
         */
        bool Jump(int refID, int position=0);
        bool jump(int refID, int position) {
            return d->SetRegion(BamRegion(refID, position));
        }
        /** 
         *  Opens a BAM file.
         *
         *  If BamReader is already opened on another file, this function closes
         *  that file, then attempts to open requested \a filename.
         *
         *  @param[in] filename name of BAM file to open
         *
         *  @return true if BAM file was opened successfully
         *  @see Close(), IsOpen(), OpenIndex()
        */
        bool Open(const std::string& filename);
        bool open(const std::string& filename) {
            return d->Open(filename);
        }

        /** 
         *  Returns the internal file pointer to the first alignment record.
         *
         *  Useful for performing multiple sequential passes through a BAM file.
         *  Calling this function clears any prior region that may have been set.
         *
         *  @note This function sets the file pointer to first alignment record
         *  in the BAM file, NOT the beginning of the file.
         *
         *  @returns \c true if rewind operation was successful
         *  @see Jump(), SetRegion()
        */
        bool Rewind(void);
        bool rewind(void) {
            return d->Rewind();
        }
        /** 
         *  Sets a target region of interest
         *
         *  Region in command option: chr1:20..999
         *  Requires that index data be available. Attempts a random-access
         *  jump in the BAM file, near \a region left boundary position.
         *
         *  Subsequent calls to GetNextAlignment() or GetNextAlignmentCore()
         *  will only return \c true when alignments can be found that overlap
         *  this \a region.
         *
         *  A \a region with no right boundary is considered open-ended, meaning
         *  that all alignments that lie downstream of the left boundary are
         *  considered valid, continuing to the end of the BAM file.
         *
         *  @warning BamRegion now represents a zero-based, HALF-OPEN interval [L, R).
         *  In previous versions of BamTools (0.x & 1.x) all intervals were treated
         *  as zero-based, CLOSED.
         *
         *  @param[in] region desired region-of-interest to activate
         *
         *  @returns true if reader was able to jump successfully to the region's left boundary
         *  @see HasIndex(), Jump()
        */
        bool SetRegion(const BamRegion& region);
        bool setRegion(const BamRegion& region) {
            return d->SetRegion(region);
        }
        /** 
         *  Sets a target region of interest.
         * 
         *  This is an overloaded function.
         * 
         *  @warning This function expects a zero-based, HALF-OPEN interval.
         *  In previous versions of BamTools (0.x & 1.x) all intervals were treated
         *  as zero-based, CLOSED.
         * 
         *  @param[in] leftRefID     referenceID of region's left boundary
         *  @param[in] leftPosition  position of region's left boundary
         *  @param[in] rightRefID    reference ID of region's right boundary
         *  @param[in] rightPosition position of region's right boundary
         * 
         *  @returns true if reader was able to jump successfully to the region's left boundary
         * 
         *  @see HasIndex(), Jump()
         *  For single chromosome such as chr4 refid=15, not sure what is the internal
         *    representation 15:-1, 15:-1?  
         */
        bool SetRegion(const int& leftRefID, const int& leftPos,
                       const int& rightRefID, const int& rightPos);
        bool setRegion(const int& leftRefID, const int& leftPos,
                       const int& rightRefID, const int& rightPos)
        {
            return d->SetRegion(BamRegion(leftRefID, leftPos, rightRefID, rightPos));
        }
        /**
         * This version is for human interface and use 1-based coordinates.
         * [leftPos, rightPos) is 1-based half open interval.
         */
        bool setRegion(const string& leftName, const int& leftPos,
                       const string& rightName, const int& rightPos)
        {
            return d->SetRegion(BamRegion(getReferenceID(leftName), leftPos-1, getReferenceID(rightName), rightPos-1));
        }

        // ----------------------
        // access alignment data
        // ----------------------
        /** Retrieves next available alignment.
         *
         *  Attempts to read the next alignment record from BAM file, and
         *  checks to see if it overlaps the current region. If no region is
         *  currently set, then the next alignment available is always
         *  considered valid.
         *
         *  If a region has been set, via Jump() or SetRegion(), an alignment
         *  is only considered valid if it overlaps the region. If the actual
         *  'next' alignment record in the BAM file does not overlap this
         *  region, then this function will read sequentially through the file
         *  until the next alignment that overlaps this region is found.  Once
         *  the region has been exhausted (i.e. the next alignment loaded is
         *  beyond the region), the function aborts and returns \c false. In
         *  this case, there is no point to continue reading, assuming properly
         *  sorted alignments.
         *
         *  This function fully populates all of the alignment's available data
         *  fields, including the string data fields (read name, bases,
         *  qualities, tags, filename).  If only positional data (refID,
         *  position, CIGAR ops, alignment flags, etc.) are required, consider
         *  using GetNextAlignmentCore() for a significant performance boost.
         *
         *  This method act as an iterator.
         *
         *  @param[out] alignment destination for alignment record data
         *  @returns \c true if a valid alignment was found
         *  @see SetRegion()
         */
        bool GetNextAlignment(BamAlignment& alignment);
        bool getNextAlignment(BamAlignment& alignment) {
        // uses Internal::BamReaderPrivate to do the work
            return d->GetNextAlignment(alignment);
        }
        bool nextAlignment(BamAlignment& alignment) {
            return d->GetNextAlignment(alignment);
        }

        /**
         * TODO: implement other typs of BamAlighment such as
         * a derived calss for MultimappedAlign, or ChimeraAlign.
         * This method could be a factory method.
         *
         * @return next BamAlignment pointer if exists
         *  otherwise return nullptr
         */
        BamAlignment* next();
        BamAlignment* nextCore();
        /** 
         * TODO: Should create a base class to read core only.
         * Retrieves next available alignment, without populating the
         * alignment's string data fields.
         *
         * Equivalent to GetNextAlignment() with respect to what is a valid
         * overlapping alignment. However, this method does NOT populate the
         * alignment's string data fields (read name, bases, qualities, tags,
         * filename). This provides a boost in speed when these fields are not
         * required for every alignment.  These fields, excluding filename, can
         * be populated 'lazily' (as needed) by calling
         * BamAlignment::BuildCharData() later.

         * @param[out] alignment destination for alignment record data
         * @returns true if a valid alignment was found
         * @see SetRegion()
         */
        bool GetNextAlignmentCore(BamAlignment& alignment);
        bool getNextAlignmentCore(BamAlignment& alignment) {
            return d->GetNextAlignmentCore(alignment);
        }
        bool nextAlignmentCore(BamAlignment& alignment) {
            return d->GetNextAlignmentCore(alignment);
        }

        // ----------------------
        // access header data
        // ----------------------

        /**
         *  Returns const reference to SAM header data.
         *
         *  Allows for read-only queries of SAM header data.
         *
         *  If you do not need to modify the SAM header, use this method to avoid the
         *  potentially expensive copy used by GetHeader().
         *
         *  @return const reference to header data object
         *  @see GetHeader(), GetHeaderText()
        */
        const SamHeader& GetConstSamHeader(void) const {
            return d->GetConstSamHeader();
        }
        const SamHeader& getSamHeader(void) const {
            return d->getSamHeader();
        }
        SamHeader& getSamHeader(void) {
            return d->getSamHeader();
        }
        /** 
         *  Returns SAM header data stand-alone object.
         *
         *  Header data is wrapped in a SamHeader object that can be
         *  conveniently queried and/or modified.  If you only need read
         *  access, consider using GetConstSamHeader() instead.
         *
         *  Note: Modifying the retrieved SamHeader object does NOT affect the
         *  current BAM file. This file has been opened in a read-only mode.
         *  However, your modified SamHeader object can be used in conjunction
         *  with BamWriter to generate a new BAM file with the appropriate
         *  header information.
         *
         *  @return an editable copy of SAM header data
         *  @see GetConstSamHeader(), GetHeaderText()
        */
        const SamHeader& GetHeader(void) const {
            return d->getSamHeader();
        }
        const SamHeader& getHeader(void) const {
            return d->getSamHeader();
        }
        SamHeader& getHeader(void) {
            return d->getSamHeader();
        }
        /** 
         *  Returns SAM header data, as SAM-formatted text.
         *
         *  Note: Modifying the retrieved text does NOT affect the current
         *  BAM file. This file has been opened in a read-only mode. However,
         *  your modified header text can be used in conjunction with BamWriter
         *  to generate a new BAM file with the appropriate header information.
         *
         *  @return SAM-formatted header text
         *  @see GetHeader()
        */
        std::string GetHeaderText(void) const;

        // ----------------------
        // access reference data
        // ----------------------

        /** 
         *  @return number of reference sequences.
         *  Mentioned in the header or actually used in the bamfile?
        */
        int GetReferenceCount(void) const;
        int getReferenceCount(void) const {
           return GetReferenceCount();
        }

        /**
         * @see RefData RefData is a simple typedef in BamAux.h
         *    typedef std::vector<RefData> RefVector;
         *    RefData: has two fields 
         *     { RefName (string), RefLength(int32_t) }
         *    
         * @returns all reference meta data
         *
         * The index in the vector (RefVector) has a pattern
         * chrM 0, chr1 1, ..., chr22 22, chrX 23, chrY 24.
         * Not sure how others are named, so you have to
         * rely on the index to get the string to be safe.
         * The number after 24 is unreliable.
         */
        const RefVector& GetReferenceData(void) const {
           return d->GetReferenceData();
        }
        /**
         * @return a const reference to vector<RefData> stored in 
         *  the BamReaderPrivate d pointer.
         */
        const RefVector& getReferenceData(void) const {
           return d->getReferenceData();
        }
        /**
         * You should never modify the refvector. This maybe dangerous.
         * Only library writer need this function.
         */
        RefVector& getReferenceData(void) {
           return d->getReferenceData();
        }

        /**
         * @return the reference meta data: ref_name, ref_length
         *   without actual sequence information. This is
         *   extracted from the header section of the Bam file.
         *
         *  the refid in the bamfile can directly index into
         *  the vector returned. refid of zero is the mitochrondria,
         *  ref 1 is the first human chromosome.
         *  ref 24 is the Y chromosome.
         *  BUT the above is not all the case, for a different
         *  refgenome. 0 is (chr)1, 1 is (chr)2,  ....
         *  So could not assume anything.
         *
         *  This method will be useful for interacting with the
         *  stdandard library of C++.
         *  Can only return a new object, there is no reference
         *  in this object.
         * would be more portable and more readable.
         * convert the d->GetReferenceData() into RefVector;
         */
        vector<pair<string,int>> getReferenceMetaData() const;
        /** 
         *  @return the ID of the reference with this name.
         *
         *  If \a refName is not found, returns -1.
         *
         *  @param[in] refName name of reference to look up
        */
        int GetReferenceID(const std::string& refName) const;
        /**
         * More efficient implementation
         */
        int getReferenceID(const std::string& refName) const {
            return d->getReferenceID(refName);
        }
        int getReferenceId(const std::string& refName) const {
            return d->getReferenceID(refName);
        }
        /**
         * @return reference name given refid
         */
        const string& getReferenceName(int refid) const;
        /**
         * @return random_refid => main_refid (chr2:1) for GRCh38decoy or empty
         *   for reference genome GRCh37decoy
         */
        map<int,int> getRefidMatch() const {
           return d->getRefidMatch();
        }
        // ----------------------
        // BAM index operations
        // ----------------------

        /** 
         *  Creates an index file for current BAM file.
         *
         *  @param[in] type file format to create, see BamIndex::IndexType for available formats
         *  @return true if index created OK
         *  @see LocateIndex(), OpenIndex()
        */
        bool CreateIndex(const BamIndex::IndexType& type=BamIndex::STANDARD);
        bool createIndex(const BamIndex::IndexType& type=BamIndex::STANDARD) {
            return d->CreateIndex(type);
        }
        /** 
         *  @return true if index data is available.
        */
        bool HasIndex(void) const;
         bool hasIndex(void) const {
            return d->HasIndex();
         }
        /** 
         *  Looks in BAM file's directory for a matching index file.
         *
         *  Use this function when you need an index file, and perhaps have a
         *  preferred index format, but do not depend heavily on which format
         *  actually gets loaded at runtime.
         *
         *  This function will defer to your \a preferredType whenever possible.
         *  However, if an index file of \a preferredType can not be found, then
         *  it will look for any other index file that corresponds to this BAM file.
         *
         *  If you want precise control over which index file is loaded, use OpenIndex()
         *  with the desired index filename. If that function returns false, you can use
         *  CreateIndex() to then build an index of the exact requested format.
         *
         *  @param[in] preferredType desired index file format, see BamIndex::IndexType for available formats
         *
         *  @return true if (any) index file could be found
        */
        bool LocateIndex(const BamIndex::IndexType& preferredType = BamIndex::STANDARD);
        bool locateIndex(const BamIndex::IndexType& preferredType = BamIndex::STANDARD) {
            return d->LocateIndex(preferredType);
        }
        /** 
         *  Opens a BAM index file.
         *
         *  @param[in] indexFilename name of BAM index file to open
         *
         *  @returns true if BAM index file was opened & data loaded successfully
         *  @see LocateIndex(), Open(), SetIndex()
        */
        bool OpenIndex(const std::string& indexFilename);
        bool openIndex(const std::string& indexFilename) {
            return d->OpenIndex(indexFilename);
        }
        /** 
         *  Sets a custom BamIndex on this reader.
         *
         *  Only necessary for custom BamIndex subclasses. Most clients should
         *  never have to use this function.
         *
         *  Example:
         *  \code
         *      BamReader reader;
         *      reader.SetIndex(new MyCustomBamIndex);
         *  \endcode
         *
         *  \note BamReader takes ownership of \a index - i.e. the BamReader will
         *  take care of deleting it when the reader is destructed, when the current
         *  BAM file is closed, or when a new index is requested.
         *
         *  \param[in] index custom BamIndex subclass created by client
         *  \sa CreateIndex(), LocateIndex(), OpenIndex()
        */
        void SetIndex(BamIndex* index);
        void setIndex(BamIndex* index) {
            d->SetIndex(index);
        }

        // ----------------------
        // error handling
        // ----------------------

        /* 
         *  Returns a human-readable description of the last error that occurred
         *
         *  This method allows elimination of STDERR pollution. Developers of client code
         *  may choose how the messages are displayed to the user, if at all.
         *
         *  @return error description
        */
        //std::string GetErrorString(void) const;
        
    // private implementation
    private:
        Internal::BamReaderPrivate* d;
};

} // namespace BamTools

#endif // BAMREADER_H
