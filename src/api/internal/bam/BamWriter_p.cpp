// ***************************************************************************
// BamWriter_p.cpp (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 18 November 2012 (DB)
// ---------------------------------------------------------------------------
// Provides the basic functionality for producing BAM files
// ***************************************************************************

#include "api/BamAlignment.h"
#include "api/BamConstants.h"
#include "api/IBamIODevice.h"
#include "api/internal/bam/BamWriter_p.h"
#include "api/internal/utils/BamException_p.h"
using namespace BamTools;
using namespace BamTools::Internal;

#include <cstdlib>
#include <cstring>

using namespace std;


// ctor
BamWriterPrivate::BamWriterPrivate(void)
    : m_isBigEndian( BamTools::SystemIsBigEndian() )
{ }

// dtor
BamWriterPrivate::~BamWriterPrivate(void) {
    Close();
}

// calculates minimum bin for a BAM alignment interval [begin, end)
uint32_t BamWriterPrivate::CalculateMinimumBin(const int begin, int end) const {
    --end;
    if ( (begin >> 14) == (end >> 14) ) return 4681 + (begin >> 14);
    if ( (begin >> 17) == (end >> 17) ) return  585 + (begin >> 17);
    if ( (begin >> 20) == (end >> 20) ) return   73 + (begin >> 20);
    if ( (begin >> 23) == (end >> 23) ) return    9 + (begin >> 23);
    if ( (begin >> 26) == (end >> 26) ) return    1 + (begin >> 26);
    return 0;
}

// closes the alignment archive
void BamWriterPrivate::Close(void) {
    // skip if file not open
    if ( !IsOpen() ) return;

    // close output stream
    try {
        m_stream.Close();
    } catch ( BamException& e ) {
        m_errorString = e.what();
    }
}

// creates a cigar string from the supplied alignment
void BamWriterPrivate::CreatePackedCigar(const vector<CigarOp>& cigarOperations, string& packedCigar) {
    // initialize
    const size_t numCigarOperations = cigarOperations.size();
    packedCigar.resize(numCigarOperations * Constants::BAM_SIZEOF_INT);
    // pack the cigar data into the string
    unsigned int* pPackedCigar = (unsigned int*)packedCigar.data();
    // iterate over cigar operations
    vector<CigarOp>::const_iterator coIter = cigarOperations.begin();
    vector<CigarOp>::const_iterator coEnd  = cigarOperations.end();
    for ( ; coIter != coEnd; ++coIter ) {
        // store op in packedCigar
        uint8_t cigarOp;
        switch ( coIter->Type ) {
            case (Constants::BAM_CIGAR_MATCH_CHAR)    : cigarOp = Constants::BAM_CIGAR_MATCH;    break;
            case (Constants::BAM_CIGAR_INS_CHAR)      : cigarOp = Constants::BAM_CIGAR_INS;      break;
            case (Constants::BAM_CIGAR_DEL_CHAR)      : cigarOp = Constants::BAM_CIGAR_DEL;      break;
            case (Constants::BAM_CIGAR_REFSKIP_CHAR)  : cigarOp = Constants::BAM_CIGAR_REFSKIP;  break;
            case (Constants::BAM_CIGAR_SOFTCLIP_CHAR) : cigarOp = Constants::BAM_CIGAR_SOFTCLIP; break;
            case (Constants::BAM_CIGAR_HARDCLIP_CHAR) : cigarOp = Constants::BAM_CIGAR_HARDCLIP; break;
            case (Constants::BAM_CIGAR_PAD_CHAR)      : cigarOp = Constants::BAM_CIGAR_PAD;      break;
            case (Constants::BAM_CIGAR_SEQMATCH_CHAR) : cigarOp = Constants::BAM_CIGAR_SEQMATCH; break;
            case (Constants::BAM_CIGAR_MISMATCH_CHAR) : cigarOp = Constants::BAM_CIGAR_MISMATCH; break;
            default:
                const string message = string("invalid CIGAR operation type") + coIter->Type;
                throw BamException("BamWriter::CreatePackedCigar", message);
        }

        *pPackedCigar = coIter->Length << Constants::BAM_CIGAR_SHIFT | cigarOp;
        pPackedCigar++;
    }
}

// encodes the supplied query sequence into 4-bit notation
void BamWriterPrivate::EncodeQuerySequence(const string& query, string& encodedQuery) {
    // prepare the encoded query string
    const size_t queryLength = query.size();
    const size_t encodedQueryLength = static_cast<size_t>((queryLength+1)/2);
    encodedQuery.resize(encodedQueryLength);
    char* pEncodedQuery = (char*)encodedQuery.data();
    const char* pQuery = (const char*)query.data();
    // walk through original query sequence, encoding its bases
    unsigned char nucleotideCode;
    bool useHighWord = true;
    while ( *pQuery ) {
        switch ( *pQuery ) {
            case (Constants::BAM_DNA_EQUAL) : nucleotideCode = Constants::BAM_BASECODE_EQUAL; break;
            case (Constants::BAM_DNA_A)     : nucleotideCode = Constants::BAM_BASECODE_A;     break;
            case (Constants::BAM_DNA_C)     : nucleotideCode = Constants::BAM_BASECODE_C;     break;
            case (Constants::BAM_DNA_G)     : nucleotideCode = Constants::BAM_BASECODE_G;     break;
            case (Constants::BAM_DNA_T)     : nucleotideCode = Constants::BAM_BASECODE_T;     break;
            case (Constants::BAM_DNA_N)     : nucleotideCode = Constants::BAM_BASECODE_N;     break;
            case (Constants::BAM_DNA_M)     : nucleotideCode = Constants::BAM_BASECODE_M;     break;
            case (Constants::BAM_DNA_R)     : nucleotideCode = Constants::BAM_BASECODE_R;     break;
            case (Constants::BAM_DNA_S)     : nucleotideCode = Constants::BAM_BASECODE_S;     break;
            case (Constants::BAM_DNA_V)     : nucleotideCode = Constants::BAM_BASECODE_V;     break;
            case (Constants::BAM_DNA_W)     : nucleotideCode = Constants::BAM_BASECODE_W;     break;
            case (Constants::BAM_DNA_Y)     : nucleotideCode = Constants::BAM_BASECODE_Y;     break;
            case (Constants::BAM_DNA_H)     : nucleotideCode = Constants::BAM_BASECODE_H;     break;
            case (Constants::BAM_DNA_K)     : nucleotideCode = Constants::BAM_BASECODE_K;     break;
            case (Constants::BAM_DNA_D)     : nucleotideCode = Constants::BAM_BASECODE_D;     break;
            case (Constants::BAM_DNA_B)     : nucleotideCode = Constants::BAM_BASECODE_B;     break;
            default:
                const string message = string("invalid base: ") + *pQuery;
                throw BamException("BamWriter::EncodeQuerySequence", message);
        }
        // pack the nucleotide code
        if ( useHighWord ) {
            *pEncodedQuery = nucleotideCode << 4;
            useHighWord = false;
        } else {
            *pEncodedQuery |= nucleotideCode;
            ++pEncodedQuery;
            useHighWord = true;
        }
        // increment the query position
        ++pQuery;
    }
}

// returns a description of the last error that occurred
std::string BamWriterPrivate::GetErrorString(void) const {
    return m_errorString;
}

// returns whether BAM file is open for writing or not
bool BamWriterPrivate::IsOpen(void) const {
    return m_stream.IsOpen();
}

// opens the alignment archive
bool BamWriterPrivate::Open(const string& filename, const string& samHeaderText,
                            const RefVector& referenceSequences)
{
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " filename :" << filename << ", referenceSequences.size() = " << referenceSequences.size() << endl; // 7

    try {
        // open the BGZF file for writing
        m_stream.Open(filename, IBamIODevice::WriteOnly);
        // write BAM file 'metadata' components
        WriteMagicNumber();
        WriteSamHeaderText(samHeaderText);
        WriteReferences(referenceSequences);
        // return success
        return true;

    } catch ( BamException& e ) {
        m_errorString = e.what();
        return false;
    }
}

// saves the alignment to the alignment archive
bool BamWriterPrivate::SaveAlignment(const BamAlignment& al) {
    //if (!al.validCigar()) {
    //  cerr << __FILE__ << ":" << __LINE__ << ": invalid cigar when saving\n";
    //  throw logic_error("cigar is inconsistent with query or reference");
    //}
    try {
        // if BamAlignment contains only the core data and a raw char data buffer
        // (as a result of BamReader::GetNextAlignmentCore())
        if ( al.SupportData.HasCoreOnly )
            WriteCoreAlignment(al);
        // otherwise, BamAlignment should contain character in the standard fields: Name, QueryBases, etc
        // (resulting from BamReader::GetNextAlignment() *OR* being generated directly by client code)
        else WriteAlignment(al);
        // if we get here, everything OK
        return true;
    } catch ( BamException& e ) {
        m_errorString = e.what();
        return false;
    }
    catch (exception& err) {
       cerr << endl << string(50, '!') << endl << al << endl;
       cerr << err.what() << endl;
       cerr << __FILE__ << ":" << __LINE__ << ":" << __func__
          << " Failed to write object inconsistent value.\n";
       throw;
    }
}

void BamWriterPrivate::SetWriteCompressed(bool ok) {
    // modifying compression is not allowed if BAM file is open
    if ( !IsOpen() )
        m_stream.SetWriteCompressed(ok);
}

void BamWriterPrivate::WriteAlignment(const BamAlignment& al) {
    // calculate char lengths
    const unsigned int nameLength         = al.Name.size() + 1;
    const unsigned int numCigarOperations = al.CigarData.size();
    const unsigned int queryLength        = ( (al.QueryBases == "*") ? 0 : al.QueryBases.size() );
    const unsigned int tagDataLength      = al.TagData.size();
    // no way to tell if alignment's bin is already defined (there is no default, invalid value)
    // so we'll go ahead calculate its bin ID before storing
    const uint32_t alignmentBin = CalculateMinimumBin(al.Position, al.GetEndPosition());
    // create our packed cigar string
    string packedCigar;
    CreatePackedCigar(al.CigarData, packedCigar);
    const unsigned int packedCigarLength = packedCigar.size();
    // encode the query
    unsigned int encodedQueryLength = 0;
    string encodedQuery;
    if ( queryLength > 0 ) {
        EncodeQuerySequence(al.QueryBases, encodedQuery);
        encodedQueryLength = encodedQuery.size();
    }
    // write the block size
    const unsigned int dataBlockSize = nameLength + packedCigarLength +
         encodedQueryLength +
         queryLength +         // here referring to quality length
         tagDataLength;
    unsigned int blockSize = Constants::BAM_CORE_SIZE + dataBlockSize;
    if (m_isBigEndian) BamTools::SwapEndian_32(blockSize);
    m_stream.Write((char*)&blockSize, Constants::BAM_SIZEOF_INT);

    // assign the BAM core data
    uint32_t buffer[Constants::BAM_CORE_BUFFER_SIZE];
    buffer[0] = al.RefID;
    buffer[1] = al.Position;
    buffer[2] = (alignmentBin << 16) | (al.MapQuality << 8) | nameLength;
    buffer[3] = (al.AlignmentFlag << 16) | numCigarOperations;
    buffer[4] = queryLength;
    buffer[5] = al.MateRefID;
    buffer[6] = al.MatePosition;
    buffer[7] = al.InsertSize;
    // swap BAM core endian-ness, if necessary
    if ( m_isBigEndian ) {
        for ( int i = 0; i < 8; ++i )
            BamTools::SwapEndian_32(buffer[i]);
    }
    // write the BAM core
    m_stream.Write((char*)&buffer, Constants::BAM_CORE_SIZE);
    // write the query name
    m_stream.Write(al.Name.c_str(), nameLength);
    // write the packed cigar
    if ( m_isBigEndian ) {
        char* cigarData = new char[packedCigarLength]();
        memcpy(cigarData, packedCigar.data(), packedCigarLength);
        if ( m_isBigEndian ) {
                for ( size_t i = 0; i < packedCigarLength; ++i )
                BamTools::SwapEndian_32p(&cigarData[i]);
        }
        m_stream.Write(cigarData, packedCigarLength);
        delete[] cigarData; // TODO: cleanup on Write exception thrown?
    }
    else
        m_stream.Write(packedCigar.data(), packedCigarLength);
    if ( queryLength > 0 ) {
        // write the encoded query sequence
        //cerr << "encodedQuery.len=" << encodedQuery.size() 
        //    << " encodedQueryLength=" << encodedQueryLength << endl;
        m_stream.Write(encodedQuery.data(), encodedQueryLength);
        //cerr << "written to m_stream\n";
        // write the base qualities as phred score (+33)
        char* pBaseQualities = new char[queryLength]();
        if (al.Qualities.empty() || (al.Qualities.size() == 1 && al.Qualities[0] == '*') || al.Qualities[0] == (char)0xFF)
            memset(pBaseQualities, 0xFF, queryLength); // if missing or '*', fill with invalid qual
        else {
            //cerr << " queryLength=" << queryLength << " al.Qualitities.len="
            //    << al.Qualities.size() << endl;
#ifdef DEBUG
            if (queryLength != al.Qualities.size()) {
                throw logic_error("queryLength " + to_string(queryLength) +
                        " different from quality length " + to_string(al.Qualities.size()));
            }
#endif
            for ( size_t i = 0; i < queryLength; ++i ) {
                //pBaseQualities[i] = al.Qualities.at(i) - 33; // FASTQ ASCII -> phred score conversion
                pBaseQualities[i] = al.Qualities[i] - 33; // FASTQ ASCII -> phred score conversion
            }
        }
        m_stream.Write(pBaseQualities, queryLength);
        //std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__
        //   << " pBaseQualities: " << pBaseQualities << endl;
        delete[] pBaseQualities;
    }
    //cerr << __LINE__ << ": writting tag data\n";
    // write the tag data
    try {
    if (m_isBigEndian) {
        char* tagData = new char[tagDataLength]();
        memcpy(tagData, al.TagData.data(), tagDataLength);
        size_t i = 0;
        while ( i < tagDataLength ) { // this loop swaps numeric type if is BigEndian
            i += Constants::BAM_TAG_TAGSIZE;  // skip tag chars (e.g. "RG", "NM", etc.)
            const char type = tagData[i];     // get tag type at position i
            ++i;
            switch ( type ) {
                case(Constants::BAM_TAG_TYPE_ASCII) :
                case(Constants::BAM_TAG_TYPE_INT8)  :
                case(Constants::BAM_TAG_TYPE_UINT8) :
                    ++i;
                    break;
                case(Constants::BAM_TAG_TYPE_INT16)  :
                case(Constants::BAM_TAG_TYPE_UINT16) :
                    BamTools::SwapEndian_16p(&tagData[i]);
                    i += sizeof(uint16_t);
                    break;
                case(Constants::BAM_TAG_TYPE_FLOAT)  :
                case(Constants::BAM_TAG_TYPE_INT32)  :
                case(Constants::BAM_TAG_TYPE_UINT32) :
                    BamTools::SwapEndian_32p(&tagData[i]);
                    i += sizeof(uint32_t);
                    break;
                case(Constants::BAM_TAG_TYPE_HEX) :
                case(Constants::BAM_TAG_TYPE_STRING) :
                    // no endian swapping necessary for hex-string/string data
                    while (tagData[i]) ++i;
                    // increment one more for null terminator
                    ++i;
                    break;
                case(Constants::BAM_TAG_TYPE_ARRAY) :
                {
                    // read array type
                    const char arrayType = tagData[i];
                    ++i;
                    // swap endian-ness of number of elements in place, then retrieve for loop
                    BamTools::SwapEndian_32p(&tagData[i]);
                    int32_t numElements;
                    memcpy(&numElements, &tagData[i], sizeof(uint32_t));
                    i += sizeof(uint32_t);
                    // swap endian-ness of array elements
                    for ( int j = 0; j < numElements; ++j ) {
                        switch (arrayType) {
                            case (Constants::BAM_TAG_TYPE_INT8)  :
                            case (Constants::BAM_TAG_TYPE_UINT8) :
                                // no endian-swapping necessary
                                ++i;
                                break;
                            case (Constants::BAM_TAG_TYPE_INT16)  :
                            case (Constants::BAM_TAG_TYPE_UINT16) :
                                BamTools::SwapEndian_16p(&tagData[i]);
                                i += sizeof(uint16_t);
                                break;
                            case (Constants::BAM_TAG_TYPE_FLOAT)  :
                            case (Constants::BAM_TAG_TYPE_INT32)  :
                            case (Constants::BAM_TAG_TYPE_UINT32) :
                                BamTools::SwapEndian_32p(&tagData[i]);
                                i += sizeof(uint32_t);
                                break;
                            default:
                                delete[] tagData;
                                const string message = string("invalid binary array type: ") + arrayType;
                                throw BamException("BamWriter::SaveAlignment", message);
                        }
                    }

                    break;
                }
                default :
                    delete[] tagData;
                    const string message = string("invalid tag type: ") + type;
                    throw BamException("BamWriter::SaveAlignment", message);
            }
        }
        m_stream.Write(tagData, tagDataLength);
        delete[] tagData; // TODO: cleanup on Write exception thrown?
    }
    else
        m_stream.Write(al.TagData.data(), tagDataLength);
    }
    catch(const exception& err) {
       cerr << __LINE__ << ": faied to write tag data\n" << err.what() << endl;
       throw;
    }
}

void BamWriterPrivate::WriteCoreAlignment(const BamAlignment& al) {

    // write the block size
    unsigned int blockSize = al.SupportData.BlockLength;
    if ( m_isBigEndian ) BamTools::SwapEndian_32(blockSize);
    m_stream.Write((char*)&blockSize, Constants::BAM_SIZEOF_INT);

    // re-calculate bin (in case BamAlignment's position has been previously modified)
    const uint32_t alignmentBin = CalculateMinimumBin(al.Position, al.GetEndPosition());

    // assign the BAM core data
    uint32_t buffer[Constants::BAM_CORE_BUFFER_SIZE];
    buffer[0] = al.RefID;
    buffer[1] = al.Position;
    buffer[2] = (alignmentBin << 16) | (al.MapQuality << 8) | al.SupportData.QueryNameLength;
    buffer[3] = (al.AlignmentFlag << 16) | al.SupportData.NumCigarOperations;
    buffer[4] = al.SupportData.QuerySequenceLength;
    buffer[5] = al.MateRefID;
    buffer[6] = al.MatePosition;
    buffer[7] = al.InsertSize;

    // swap BAM core endian-ness, if necessary
    if ( m_isBigEndian ) {
        for ( int i = 0; i < 8; ++i )
            BamTools::SwapEndian_32(buffer[i]);
    }

    // write the BAM core
    m_stream.Write((char*)&buffer, Constants::BAM_CORE_SIZE);

    // write the raw char data
    m_stream.Write((char*)al.SupportData.AllCharData.data(),
                   al.SupportData.BlockLength-Constants::BAM_CORE_SIZE);
}

void BamWriterPrivate::WriteMagicNumber(void) {
    std::unique_lock<std::mutex> lock(m_mutex_write); // swang
    // write BAM file 'magic number'
    m_stream.Write(Constants::BAM_HEADER_MAGIC, Constants::BAM_HEADER_MAGIC_LENGTH);
}

void BamWriterPrivate::WriteReferences(const BamTools::RefVector& referenceSequences) {

    // write the number of reference sequences
    uint32_t numReferenceSequences = referenceSequences.size();
    cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " numReferenceSequences " << numReferenceSequences << endl; // swang
    if ( m_isBigEndian ) BamTools::SwapEndian_32(numReferenceSequences);
    m_stream.Write((char*)&numReferenceSequences, Constants::BAM_SIZEOF_INT);

    // foreach reference sequence
    RefVector::const_iterator rsIter = referenceSequences.begin();
    RefVector::const_iterator rsEnd  = referenceSequences.end();
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " referenceSequences " << referenceSequences.size() << endl; // swang
    for ( ; rsIter != rsEnd; ++rsIter ) {

        // write the reference sequence name length (+1 for terminator)
        const uint32_t actualNameLen = rsIter->RefName.size() + 1;
        uint32_t maybeSwappedNameLen = actualNameLen;
        if ( m_isBigEndian ) BamTools::SwapEndian_32(maybeSwappedNameLen);
        m_stream.Write((char*)&maybeSwappedNameLen, Constants::BAM_SIZEOF_INT);

        // write the reference sequence name
        m_stream.Write(rsIter->RefName.c_str(), actualNameLen);

        // write the reference sequence length
        int32_t referenceLength = rsIter->RefLength;
        if ( m_isBigEndian ) BamTools::SwapEndian_32(referenceLength);
        m_stream.Write((char*)&referenceLength, Constants::BAM_SIZEOF_INT);
    }
}

// swang begin

bool BamWriterPrivate::MultiThreadsWrap(const std::string& filename,
            const std::string& samHeaderText,
            const BamTools::RefData& referenceSequence, int index)
{
    std::unique_lock<std::mutex> lock(mtx);
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " index = " << index << " next_index = " << next_index << endl; // swang
    cv.wait(lock, [this, index]() { return index == next_index; });  // Ensure correct order
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " wake up  ===> " << endl; // swang
    m_stream.Open(filename, IBamIODevice::WriteOnly);
    // write BAM file 'metadata' components
    WriteMagicNumber();
    WriteSamHeaderText(samHeaderText);
    WriteReference(std::ref(referenceSequence)); // index);

    next_index++; // Move to the next sequence
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " after wake up next_index " << next_index << endl;
    cv.notify_all(); // Wake up the next waiting thread
    return true;
}

// opens the alignment archive
bool BamWriterPrivate::OpenMultiThreads(const string& filename, const string& samHeaderText,
                            const RefVector& referenceSequences)
{
    try {
        // open the BGZF file for writing
        int sz = referenceSequences.size();
        //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " sz " << sz << endl; // swang

        std::vector<std::thread> threads;

        /***/
        if (1)
        {
            //m_stream.Open(filename, IBamIODevice::WriteOnly);
            for (int ix = 0; ix < sz; ix++)
            {
                int n  = ix;
                threads.emplace_back([this, n, &filename, &samHeaderText, &referenceSequences]() {
                    MultiThreadsWrap(filename,samHeaderText,referenceSequences[n],n);
                    });
            }
            // Join all threads
            for (auto& th : threads) {
                if (th.joinable()) {
                    th.join();
                }
            }
        }
        else  // debug code
        {
            for (int n = 0; n < sz; n++)
            {
                cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " n  referenceSequences " << n << endl;
                if(0)
                {
                    m_stream.Open(filename, IBamIODevice::WriteOnly);
                    // write BAM file 'metadata' components
                    WriteMagicNumber();
                    WriteSamHeaderText(samHeaderText);
                    WriteReference(referenceSequences[n]);
                }
                else
                {
                    MultiThreadsWrap(filename, samHeaderText,referenceSequences[n], n);
                }
            }
        }

        // return success
        return true;

    } catch ( BamException& e ) {
        m_errorString = e.what();
        return false;
    }
}

void BamWriterPrivate::WriteReference(const BamTools::RefData& referenceSequence)
{
    // write the number of reference sequences
    uint32_t numReferenceSequences = 1;// referenceSequences.size();
    if ( m_isBigEndian ) BamTools::SwapEndian_32(numReferenceSequences);
    m_stream.Write((char*)&numReferenceSequences, Constants::BAM_SIZEOF_INT);

    // write the reference sequence name length (+1 for terminator)
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << "numReferenceSequences = " << numReferenceSequences << endl;
    string refName =  referenceSequence.RefName;
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << "refName = " << refName << endl;
    const uint32_t actualNameLen = referenceSequence.RefName.size() + 1; //->RefName.size() + 1;
    uint32_t maybeSwappedNameLen = actualNameLen;
    if ( m_isBigEndian ) BamTools::SwapEndian_32(maybeSwappedNameLen);
    m_stream.Write((char*)&maybeSwappedNameLen, Constants::BAM_SIZEOF_INT);

    // write the reference sequence name
    //m_stream.Write(rsIter->RefName.c_str(), actualNameLen);
    m_stream.Write(referenceSequence.RefName.c_str(), actualNameLen);
    // write the reference sequence length
    int32_t referenceLength = referenceSequence.RefLength; // ->RefLength;
    if ( m_isBigEndian ) BamTools::SwapEndian_32(referenceLength);
    m_stream.Write((char*)&referenceLength, Constants::BAM_SIZEOF_INT);
}

// swang end

void BamWriterPrivate::WriteSamHeaderText(const std::string& samHeaderText) {

    std::unique_lock<std::mutex> lock(m_mutex_write); // swang
    // write the SAM header  text length
    const uint32_t actualHeaderLen = samHeaderText.size();
    uint32_t maybeSwappedHeaderLen = samHeaderText.size();
    if ( m_isBigEndian ) BamTools::SwapEndian_32(maybeSwappedHeaderLen);
    m_stream.Write((char*)&maybeSwappedHeaderLen, Constants::BAM_SIZEOF_INT);

    // write the SAM header text
    if ( actualHeaderLen > 0 )
        m_stream.Write(samHeaderText.data(), actualHeaderLen);

}
