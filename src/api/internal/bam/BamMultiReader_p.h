// ***************************************************************************
// BamMultiReader_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 14 January 2013 (DB)
// ---------------------------------------------------------------------------
// Functionality for simultaneously reading multiple BAM files
// *************************************************************************

#ifndef BAMMULTIREADER_P_H
#define BAMMULTIREADER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/SamHeader.h"
#include "api/BamMultiReader.h"
#include "api/internal/bam/BamMultiMerger_p.h"
#include <string>
#include <vector>
#include <thread> // swang

namespace BamTools {
namespace Internal {

static const int MAX_NUM_THREADS = 10; // swang , should be configuable,and used later
class BamMultiReaderPrivate {

    // typedefs
    public:
        typedef std::pair<BamReader*, BamAlignment*> ReaderAlignment;

    // constructor / destructor
    public:
        BamMultiReaderPrivate(void);
        ~BamMultiReaderPrivate(void);

    // public interface
    public:

        // file operations
        bool Close(void);
        // swang not needed ? bool CloseFile(const std::string& filename);
        const std::vector<std::string> Filenames(void) const;
        bool Jump(int refID, int position = 0);
        bool Open(const std::vector<std::string>& filenames);

        // swang begin
        bool CreateReaderThreads(const std::vector<std::string>& filenames);
        std::vector<std::thread> m_threads;
        // swang end

        //bool OpenFile(const std::string& filename); // swang not used ??
        bool Rewind(void);
        bool SetRegion(const BamRegion& region);

        // access alignment data
        BamMultiReader::MergeOrder GetMergeOrder(void) const;
        bool GetNextAlignment(BamAlignment& al);
        bool GetNextAlignmentCore(BamAlignment& al);
        bool HasOpenReaders(void);
        bool SetExplicitMergeOrder(BamMultiReader::MergeOrder order);

        // access auxiliary data
        SamHeader GetHeader(void) const;
        std::string GetHeaderText(void) const;
        int GetReferenceCount(void) const;
        const BamTools::RefVector GetReferenceData(void) const;
        int GetReferenceID(const std::string& refName) const;

        // BAM index operations
        bool CreateIndexes(const BamIndex::IndexType& type = BamIndex::STANDARD);
        bool HasIndexes(void) const;
        bool LocateIndexes(const BamIndex::IndexType& preferredType = BamIndex::STANDARD);
        bool OpenIndexes(const std::vector<std::string>& indexFilenames);

        // error handling
        std::string GetErrorString(void) const;

    // 'internal' methods
    // swang public:
    private: // swang
        bool CloseFiles(const std::vector<std::string>& filenames); // swang will let the thread to do it
        IMultiMerger* CreateAlignmentCache(void);
        bool PopNextCachedAlignment(BamAlignment& al, const bool needCharData);
        bool RewindReaders(void);
        void SaveNextAlignment(BamReader* reader, BamAlignment* alignment);
        void SetErrorString(const std::string& where, const std::string& what) const; //
        bool UpdateAlignmentCache(void);
        bool ValidateReaders(void) const;

    // data members
    private: // swang
        std::vector<MergeItem> m_readers; // swang the vector should be used to create threads
    public:
        IMultiMerger* m_alignmentCache;

        bool m_hasUserMergeOrder;
        BamMultiReader::MergeOrder m_mergeOrder;

        mutable std::string m_errorString;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMMULTIREADER_P_H
