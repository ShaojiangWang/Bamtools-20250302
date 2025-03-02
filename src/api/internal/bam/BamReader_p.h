// ***************************************************************************
// BamReader_p.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 18 November 2012 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for reading BAM files
// ***************************************************************************

#ifndef BAMREADER_P_H
#define BAMREADER_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#include "api/BamAlignment.h"
#include "api/BamIndex.h"
//#include "api/BamReader.h" // include this will cause recursive definition!
#include "api/SamHeader.h"
#include "api/internal/bam/BamHeader_p.h"
#include "api/internal/bam/BamRandomAccessController_p.h"
#include "api/internal/io/BgzfStream_p.h"
#include <string>

using namespace std;
//using namespace BamTools;
namespace BamTools {
class BamReader;
namespace Internal {

/**
 * Should have a hiarchy of different readers
 * Base: Sequential readers
 * Random readers for different situations.
 */
class BamReaderPrivate {
    // ctor & dtor
    public:
       /**
        * Hold a pointer to the BamReader parent
        * So that BamReader could have multiple 
        * objects such as random controller in addtion to this one?
        */
        //BamReaderPrivate(BamTools::BamReader* parent); // swang: why do we need m_parent which is not used in the class ?
        BamReaderPrivate(); // swang:
        ~BamReaderPrivate(void);

    // BamReader interface
    public:
        // file operations
        bool Close(void);
        /**
         * retrieve the file name
         */
        const std::string& Filename(void) const;
        std::string& Filename(void) {
           return m_filename;
        }
        bool IsOpen(void) const;
        /**
         * Open the stream only without opening the index.
         */
        bool Open(const std::string& filename);
        bool Rewind(void);
        bool SetRegion(const BamRegion& region);
        /**
         * Take alignment data from file stream and 
         * fill up the Alignment object.
         * @return true if obtained the next valid alignment
         *    false otherwise.
         */
        bool GetNextAlignment(BamAlignment& alignment);
        /**
         * Collect only the core part of the alignment.
         */
        bool GetNextAlignmentCore(BamAlignment& alignment);

        // access auxiliary data
        std::string GetHeaderText(void) const;
        /**
         * @return a const reference to the internal SamHeader object.
         */
        const SamHeader& GetConstSamHeader(void) const;
        /**
         * @return a const reference to the internal SamHeader stored in
         *   m_header.
         */
        const SamHeader& getSamHeader() const {
           return m_header.getSamHeader();
        }
        /**
         * @return a reference to the SamHeader stored in m_header.
         *   The caller can modify the return object which if using
         *   reference will change this object.
         */
        SamHeader& getSamHeader() {
           return m_header.getSamHeader();
        }
        /**
         * Avoid using this object if you only want to read.
         */
        SamHeader GetSamHeader(void) const;
        int GetReferenceCount(void) const;
        /**
         * @return a const reference to m_references
         */
        const RefVector& GetReferenceData(void) const {
            return m_references;
        }
        const RefVector& getReferenceData(void) const {
            return m_references;
        }
        RefVector& getReferenceData(void) {
            return m_references;
        }
        int GetReferenceID(const std::string& refName) const;
        /**
         * Better implementation, more efficient.
         */
        int getReferenceID(const std::string& refName) const;

        // index operations
        /**
         * Create an index of type and attach to RandomAccessController.
         */
        bool CreateIndex(const BamIndex::IndexType& type);
        /**
         * Delegate to RandomAccessController that has
         * information about index or mot.
         */
        bool HasIndex(void) const;
        /**
         * delegate to RandomAccessControler
         */
        bool LocateIndex(const BamIndex::IndexType& preferredType);
        /**
         * Open a given index file.
         * RandomAcessController does the opening
         * and store the index information.
         * @return true if successful.
         */
        bool OpenIndex(const std::string& indexFilename);
        void SetIndex(BamIndex* index);
        /**
         * Return a look up table for refid => master_refid
         * As an example:
         * chr2_KI270715v1_random (refid:34) => chr2 (refid:1)
         * For GRCh38decoy 
         * chr1, ..., chr21, chr22, chrX, chrY, chrM, 
         * chr1_KI270706v1_random,  chr1_KI270707v1_random, chr1_KI270708v1_random,
         * , ... ,
         * chrUn_JTFH01001992v1_decoy can be ignored
         * chrEBV  AC:AJ507799.2  gi:86261677  LN:171823  rl:decoy  M5:6743bd63b3ff2b5b8985d8933c53290a  SP:Human_herpesvirus_4  tp:circular
         * ========================
         * For human_g1k_v37_decoy.fasta the format is different 
         * and there is no such naming method of random chromosomes
         * so can ignore.
         * @return chr_random => chr refid mapping for GRCh38decopy only
         *    For human_g1k_v37_decoy will return empty container.
         */
        map<int,int> getRefidMatch() const;

    // internal methods, but available as a BamReaderPrivate 'interface'
    //
    // these methods should only be used by BamTools::Internal classes
    // (currently only used by the BamIndex subclasses)
    public:
        // retrieves header text from BAM file
        void LoadHeaderData(void);
        /**
         * retrieves BAM alignment under file pointer
         * (does no overlap checking or character data parsing)
         * Parse the Cigar data into CigarData field.
         *
         * @return true if success false if failure.
         */
        bool LoadNextAlignment(BamAlignment& alignment);
        // builds reference data structure from BAM file
        bool LoadReferenceData(void);
        // seek reader to file position
        bool Seek(const int64_t& position);
        // return reader's file position
        int64_t Tell(void) const;

    // data members

    public:
        // general BAM file data
        int64_t     m_alignmentsBeginOffset;
        /**
         * The input bamf file name or path
         */
        std::string m_filename;
        /**
         * vector of [refname, reflen] index by refid
         * RefVector is a typedef of vector<RefData> in BamAux.h
         * RefData { string, int32_t }
         */
        RefVector m_references;
        /** 
         * system data
         * TODO: Should auto detect.
         */
        bool m_isBigEndian;
        /**
         * parent BamReader
         * TODO: poor design should remove.
         */
        //BamReader* m_parent; // swang: why do we need m_parent which is not used in the class ?

        // BamReaderPrivate components
        /**
         * Holds the BamFile header
         */
        BamHeader m_header;
        /**
         * Random access controller
         */
        BamRandomAccessController m_randomAccessController;
        /**
         * Bgzfile stream
         */
        BgzfStream m_stream;
};

} // namespace Internal
} // namespace BamTools

#endif // BAMREADER_P_H
