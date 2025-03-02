// ***************************************************************************
// BamMultiReader.h (c) 2010 Erik Garrison, Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 14 January 2013 (DB)
// ---------------------------------------------------------------------------
// Convenience class for reading multiple BAM files.
// ***************************************************************************

#ifndef BAMMULTIREADER_H
#define BAMMULTIREADER_H

#include "api/api_global.h"
#include "api/BamReader.h"
#include <map>
#include <sstream>
#include <string>
#include <utility>

namespace BamTools {

namespace Internal {
    class BamMultiReaderPrivate;
} // namespace Internal

/** 
 * Convenience class for reading multiple BAM files.
 * Each alignment remembers its file name.
 * If member files sorted, the reading from them will
 * also be in sorted order.
 */
class API_EXPORT BamMultiReader {

    // enums
    public:
        // possible merge order strategies
        enum MergeOrder { RoundRobinMerge=0, MergeByCoordinate, MergeByName };

    // constructor / destructor
    public:
        /**
         * Default constructor.
         */
        BamMultiReader(void);
        ~BamMultiReader(void);

    // public interface
    public:

        // ----------------------
        // BAM file operations
        // ----------------------

        // closes all open BAM files
        bool Close(void);
        // close only the requested BAM file
        // swang not used ? bool CloseFile(const std::string& filename);
        // returns list of filenames for all open BAM files
        const std::vector<std::string> Filenames(void) const;
        // returns curent merge order strategy
        BamMultiReader::MergeOrder GetMergeOrder(void) const;
        // returns true if multireader has any open BAM files
        bool HasOpenReaders(void) const;
        // performs random-access jump within current BAM files
        bool Jump(int refID, int position = 0);
        /** opens BAM files.
         * @param filenames a list of file names stored in a vector.
         */
        bool Open(const std::vector<std::string>& filenames);

        // swang begin
        bool CreateThreads(const std::vector<std::string>& filenames);  // swang addition to replace Open() above
        /***
        after all threads finished we need to call follow:
        // creates index files for current BAM files
        bool CreateIndexes(const BamIndex::IndexType& type = BamIndex::STANDARD);
        // returns true if all BAM files have index data available
        bool HasIndexes(void) const;
        // looks for index files that match current BAM files
        bool LocateIndexes(const BamIndex::IndexType& preferredType = BamIndex::STANDARD);
        // opens index files for current BAM files.
        bool OpenIndexes(const std::vector<std::string>& indexFilenames);
         *
         */
        // swang end

        // opens a single BAM file, adding to any other current BAM files
        // swang not used ? bool OpenFile(const std::string& filename);
        // returns file pointers to beginning of alignments
        bool Rewind(void);
        // sets an explicit merge order, regardless of the BAM files' SO header tag
        bool SetExplicitMergeOrder(BamMultiReader::MergeOrder order);
        // sets the target region of interest
        bool SetRegion(const BamRegion& region);
        // sets the target region of interest
        bool SetRegion(const int& leftRefID,
                       const int& leftPosition,
                       const int& rightRefID,
                       const int& rightPosition);

        // ----------------------
        // access alignment data
        // ----------------------

        /** 
         *  Retrieves next available alignment.
         *
         *  Equivalent to BamReader::GetNextAlignment() with respect to what is a valid
         *  overlapping alignment and what data gets populated.
         *
         *  This method takes care of determining which alignment actually is 'next'
         *  across multiple files, depending on their sort order.
         *
         *  @param[out] alignment destination for alignment record data
         *  @returns \c true if a valid alignment was found
         *  @see GetNextAlignmentCore(), SetExplicitMergeOrder(), SetRegion(), BamReader::GetNextAlignment()
         *
         *  This method can be used to perform merge sort
         *  in essence.
        */
        bool GetNextAlignment(BamAlignment& alignment);
        // retrieves next available alignment (without populating the alignment's string data fields)
        bool GetNextAlignmentCore(BamAlignment& alignment);

        // ----------------------
        // access auxiliary data
        // ----------------------

        // returns unified SAM header for all files
        SamHeader GetHeader(void) const;
        // returns unified SAM header text for all files
        std::string GetHeaderText(void) const;
        // returns number of reference sequences
        int GetReferenceCount(void) const;
        // returns all reference sequence entries.
        const BamTools::RefVector GetReferenceData(void) const;
        // returns the ID of the reference with this name.
        int GetReferenceID(const std::string& refName) const;
        int getReferenceId(const std::string& refName) const {
           return GetReferenceID(refName);
        }

        // ----------------------
        // BAM index operations
        // ----------------------

        // creates index files for current BAM files
        bool CreateIndexes(const BamIndex::IndexType& type = BamIndex::STANDARD);
        // returns true if all BAM files have index data available
        bool HasIndexes(void) const;
        // looks for index files that match current BAM files
        bool LocateIndexes(const BamIndex::IndexType& preferredType = BamIndex::STANDARD);
        // opens index files for current BAM files.
        bool OpenIndexes(const std::vector<std::string>& indexFilenames);

        // ----------------------
        // error handling
        // ----------------------

        // returns a human-readable description of the last error that occurred
        std::string GetErrorString(void) const;

    // private implementation
    private:
        /** private implementation.
         * The sorrounding class use the private class
         * as a delegator.  This is a design pattern commonly
         * used in Java.  It looks not very good in C++.
         */
        Internal::BamMultiReaderPrivate* d;
        //Internal::BamReaderPrivate* m_p_bamReaderPrivate;  // swang
};

} // namespace BamTools

#endif // BAMMULTIREADER_H
