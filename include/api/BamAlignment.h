// ***************************************************************************
// BamAlignment.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 July 2013 (DB)
// ---------------------------------------------------------------------------
// Provides the BamAlignment data structure
// ***************************************************************************

#ifndef BAMALIGNMENT_H
#define BAMALIGNMENT_H

#include "api/api_global.h"
#include "api/BamAux.h"
#include "api/BamConstants.h"
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <typeinfo>
#include <map>
#include <typeinfo>
//#include <mutex>
#include <mutex>
#include <limits>

// to turn off debug mode uncomment the following line
//#define NDEBUG 
//above will disable assert()
namespace BamTools {

//! \cond
// forward declaration of BamAlignment's "friends"
namespace Internal {
    class BamReaderPrivate;
    class BamWriterPrivate;
} // namespace Internal
//! \endcond

class BamAlignmentException : public logic_error {
   public:
      BamAlignmentException(const string& msg) 
      : logic_error(msg) 
      { }
};

class BamNotagException : public BamAlignmentException {
   public:
      BamNotagException(const string& msg)
         : BamAlignmentException(msg)
      { }
};

class BamTypeException : public logic_error {
   public:
      BamTypeException(const string& msg)
         : logic_error(msg)
      { }
};

/**
 *  BamAlignment data structure to hold a single
 *  Alignment data with summary information to it partner 
 *  (mate or read2) alignment data.
 */
class API_EXPORT BamAlignment {
   // API_EXPORT are constructed used for Window DDL
    /////// constructors & destructor //////////
    public:
        /** 
         *   Default constructor of empty objects.
         */
        BamAlignment()
            : Name(), QueryBases(), AlignedBases(), Qualities(), 
            TagData(), RefID(-1), Position(-1), Bin(0), MapQuality(0),
            AlignmentFlag(0), CigarData(), MateRefID(-1), MatePosition(-1), 
            InsertSize(0), SupportData()           
         { }

        /**
         * Convenient constructor for testing
         * No Cigar string thus no actual alignment information.
         * Also Filename is not given.
         */
        BamAlignment(const std::string& qname, int32_t refid, int32_t refpos, uint32_t alnflag, 
              int32_t mrefid, int32_t mrefpos, const std::string& queryseq, const std::string& qstring)
           : Name(qname), QueryBases(queryseq), AlignedBases(), Qualities(qstring),
             TagData(), RefID(refid), Position(refpos), Bin(0), MapQuality(0),
             AlignmentFlag(alnflag), CigarData(),
            MateRefID(mrefid), MatePosition(mrefpos), InsertSize(0),
            SupportData()
        {
           SupportData.QuerySequenceLength=QueryBases.size(); // is this duplicate of QueryLength?
        }
        BamAlignment(const std::string& qname, int32_t refid, int32_t refpos, uint32_t alnflag, 
              int32_t mrefid, int32_t mrefpos, std::string&& queryseq, std::string&& qstring)
           : Name(qname), QueryBases(std::move(queryseq)), AlignedBases(), Qualities(std::move(qstring)),
             TagData(), RefID(refid), Position(refpos), Bin(0), MapQuality(0),
             AlignmentFlag(alnflag), CigarData(),
            MateRefID(mrefid), MatePosition(mrefpos), InsertSize(0),
            SupportData()
        {
           SupportData.QuerySequenceLength=QueryBases.size(); // is this duplicate of QueryLength?
        }
        /**
         * Convenient constructor with Cigar input for testing.
         * If single then mate_refid is -1.
         */
        BamAlignment(const std::string& qname, int32_t refid, int32_t refpos, uint32_t alnflag, 
              int32_t mrefid, int32_t mrefpos, const std::string& queryseq, const std::string& qstring,
              const string& cigarstr)
           : Name(qname), 
             QueryBases(queryseq), AlignedBases(), Qualities(qstring),
             TagData(), RefID(refid), Position(refpos), Bin(0), MapQuality(0),
             AlignmentFlag(alnflag), CigarData(),
            MateRefID(mrefid), MatePosition(mrefpos), InsertSize(0),
            SupportData()
        { 
           setCigar(cigarstr); 
           SupportData.QuerySequenceLength=QueryBases.size(); // is this duplicate of QueryLength?
        }
        /**
         * Using move as much as possible.
         */
        BamAlignment(const std::string& qname, int32_t refid, int32_t refpos, uint32_t alnflag, 
              int32_t mrefid, int32_t mrefpos, std::string&& queryseq, std::string&& qstring,
              const string& cigarstr)
           : Name(qname), 
             QueryBases(std::move(queryseq)), AlignedBases(), 
             Qualities(std::move(qstring)),
             TagData(), RefID(refid), Position(refpos), Bin(0), MapQuality(0),
             AlignmentFlag(alnflag), CigarData(),
            MateRefID(mrefid), MatePosition(mrefpos), InsertSize(0),
            SupportData()
        { 
           setCigar(cigarstr); 
           SupportData.QuerySequenceLength=QueryBases.size(); // is this duplicate of QueryLength?
        }
        /** 
         *  Copy constructor
         */
        BamAlignment(const BamAlignment& other);
        /**
         * Move constructor
         */
        BamAlignment(BamAlignment&& other);
        /**
         * This class does not have dynamicly allocated data.
         * Destructor plannng to derive from this class:
         * MultiAln, ChimearAln
         */
        virtual ~BamAlignment() { };
        /**
         * Test version for looking at the critical information about the
         * alignment.  Should output a human readable tab-dlimited forms with
         * most of the essential information about an alignment.
         */
        friend std::ostream& operator<<(std::ostream& ous, const BamAlignment& ba);
        /**
         * Assignment operator
         */
        BamAlignment& operator=(const BamAlignment& other);
        /**
         * move assignment operator
         */
        BamAlignment& operator=(BamAlignment&& other);
        /**
         * First compare by [begin,end], then first mate < second mate
         * Mainly used for sorting out alignments on the reference.
         */
        bool operator<(const BamAlignment& other) const;
        bool operator>(const BamAlignment& other) const;
        /**
         * Two objects are considered identical if
         * same name, same mate. Mainly using this to sort out
         * same query mapped to multiple locations.
         */
        bool operator==(const BamAlignment& other) const;
        bool operator<=(const BamAlignment& other) const {
           return operator<(other) || sameLocation(other);
        }
        bool operator>=(const BamAlignment& other) const {
           return operator>(other) || sameLocation(other);
        }
        bool sameLocation(const BamAlignment& o) const;
        /**
         * interface function should never be used by this class
         * Should be implemented by derived classes.
         */
        virtual bool operator()(BamAlignment* ba) { 
           if (ba != nullptr) {
              throw runtime_error("base class cannot add true pointer");
           }
            return true;
        }

        //// informational methods //////
        /**
         * Should be derived by subclasses
         */
        virtual bool isFull() const {
           return true;
        }
        /**
         * Should be derived by subclasses
         */
        virtual void findBestLocation() {
        }
        /**
         * This object has not been populated with real data.
         */
        bool empty() const {
           return getLength() == 0;
        }
        bool isEmpty() const { return getLength() == 0; }
        /**
         * Not sure the condition is sufficient, need some testing.
         * The caller needs to make sure the query name is the same or 
         * not. This function only cares about the reference location.
         * @return true if two aligns have the same refid, strand, Postion, endPosition
         */
        bool sameHit(const BamAlignment& other) const {
           return getReferenceId() == other.getReferenceId()
              && getStrand() == other.getStrand()
              && getPosition() == other.getPosition()
              && getEndPosition() == other.getEndPosition();
        }
        uint32_t getAlignmentFlag() const {
            return AlignmentFlag;
        }
         /** 
          * This cannot be relied on.
          * @return true if this read is a PCR duplicate
         */
        bool IsDuplicate(void) const;
        bool IsFailedQC(void) const;          // returns true if this read failed quality control
        /**
         * Mate one is determined by the sequencer operation. It 
         * read the same read twice. The second mate is reading from
         * the other strand.
         * @return true if alignment is first mate on paired-end read
         * First mate is also first read: Two words meaning the same thing.
         * First read usually have better quality than second read.
         */
        bool IsFirstMate(void) const;         
        /**
         * @return true if is the First Read of a paired end reads.
         * Alias for IsFirstMate.
         * Use this version for carmel casing.
         * @see isFirstRead
         */
        bool isFirstMate(void) const { 
           return (AlignmentFlag & READ_1) == READ_1; 
        }
        /**
         * Alias for isFirstMate()
         * First read usually has slightly better quality than the second read
         * in the read pair.
         */
        bool isFirstRead(void) const { 
           return (AlignmentFlag & READ_1) == READ_1; 
        }
        /** 
         * @returns true if alignment is second mate on paired-end read
         */
        bool IsSecondMate(void) const;        
        /**
         * @return true if it is the second read (mate)
         */
        bool isSecondMate(void) const { 
           return (AlignmentFlag & READ_2) == READ_2; 
        }
        /**  in C++ true is 1 false is 0
         * Alias for isSecondMate()
         */
        bool isSecondRead(void) const { 
           return (AlignmentFlag & READ_2) == READ_2; 
        }
        /**
         * @return 1 for first mate, 2 for second mate, 
         *     and 0 for unknown mate or not pair-end read.
         */
        int getMate() const { 
           if (isFirstMate()) return 1; 
           else if (isSecondMate()) return 2;
           else return 0;
        }
        /**
         * @return 0 if is unpaired or one of the mate is unmapped
         *    otherwise return 1 for read (mate) 1, 2 for read (mate) 2
         */
        int getMappingStatus() const {
           if (!isPaired() || isMateUnmapped() || isUnmapped()
                 || !mateOnSameReference()) return 0;
           return getMate();
        }
        /**
         * @return true if query is mapped to references
         * @see isMapped is an alias conforming to convension.
         */
        bool IsMapped(void) const;            
        /** 
         * @return true if alignment is mapped
         */
        bool isMapped(void) const {
           return !((AlignmentFlag & UNMAPPED) == UNMAPPED);            
        }
        /** 
         * @return true if alignment is not mapped
         */
        bool isUnmapped() const {
           return (AlignmentFlag & UNMAPPED) == UNMAPPED;
        }
        /** 
         * Is a flag check
         * @return true if alignment's mate is mapped
         */
        bool IsMateMapped(void) const;        
        bool isMateMapped(void) const {
            return !((AlignmentFlag & MATE_UNMAPPED) == MATE_UNMAPPED);
        }
        /**
         * @return true if Mate is unmapped
         */
        bool isMateUnmapped() const {
           return (AlignmentFlag & MATE_UNMAPPED) == MATE_UNMAPPED;
        }
        /** 
         *  If the read is mapped to the reference then there is
         *  only two choices +/-.  - is the reverse strand.
         *
         *  @return true if alignment mapped to reverse strand of refseq.
         */
        bool IsReverseStrand(void) const;     
        /**
         * @return true if alignment mapped to reverse strand of refseq.
         */
        bool isReverseStrand(void) const {
           return (AlignmentFlag & REVERSE_STRAND) == REVERSE_STRAND; 
        }
        /**
         * The query sequence is in the original forward direction.
         */
        bool isForwardStrand() const {
           return !((AlignmentFlag & REVERSE_STRAND) == REVERSE_STRAND); 
        }
        /**
         * The query sequence is the reverse complement of the original
         * sequence.
         */
        void setReverseStrand() {
           AlignmentFlag |= REVERSE_STRAND;
        }
        void setForwardStrand() {
           AlignmentFlag &= ~(REVERSE_STRAND);
        }
        /**
         * There are only two states; cannot be zero.
         * @return -1 reverse strand, +1 for forward strand; 
         */
        int getStrand() const {
           if (isReverseStrand()) return -1;
           return 1;
           //if (isForwardStrand()) return 1;
           //return 0; // never reach this line
        }
        char getStrandChar() const {
           if (isReverseStrand()) return '-';
           return '+';
           //if (isForwardStrand()) return '+';
           //return '?'; // nover reach this line
        }
        /** 
         * @return true if alignment's mate mapped to reverse strand
         */
        bool IsMateReverseStrand(void) const; 
        bool isMateReverseStrand() const {
            return (AlignmentFlag & MATE_REVERSE_STRAND) == MATE_REVERSE_STRAND;
        }
        bool isMateForwardStrand() const {
            return !((AlignmentFlag & MATE_REVERSE_STRAND) == MATE_REVERSE_STRAND);
        }
        void setMateForwardStrand() {
           AlignmentFlag &= ~MATE_REVERSE_STRAND;
        }
        void setMateReverseStrand() {
           AlignmentFlag |= MATE_REVERSE_STRAND;
        }
        bool isMateOppositeStrand() const {
           return  (isForwardStrand() && isMateReverseStrand())
              || (isReverseStrand() && isMateForwardStrand()); 
        }
        bool isMateSameStrand() const {
           return  (isForwardStrand() && isMateForwardStrand())
              || (isReverseStrand() && isMateReverseStrand()); 
        }
        /**
         * Get a value to represent the strand.
         * For simple alignment of single sequence it is either
         * +1 or -1.  For Merged sequence which is reprented as
         * having the XO tag this value could be [0, 1)
         * a fractional number between 0 and 1.
         * If the overlap is 100% of the read length, then it is zero.
         * The higher the value approaching 1 the lower the overlap.
         * We use a single return value to record
         * this multiple situations.
         */
        double getFractionStrand() const;
        /** @returns true if alignment part of paired-end read
         */
        bool IsPaired(void) const; 
        /**
         * If is not paired then READ1 or READ2 is meaningless
         * @return this alignment is paired or single.
         */
        bool isPaired() const { 
            return (AlignmentFlag & PAIRED) == PAIRED;
        }
        bool isUnpaired() const { 
            return (AlignmentFlag & PAIRED) != PAIRED;
        }
        /** 
         * @return true if reported position is primary alignment
         */
        bool IsPrimaryAlignment(void) const;  
        /**
         * Test SECONDARY flag is not set.
         */
         bool isPrimaryAlignment(void) const  {
            return !((AlignmentFlag & SECONDARY) == SECONDARY);
         }
        /**
         * @return true if is secondary alignment
         *    This indicates that the query was mapped
         *    multiple times in the genome. The choice of
         *    primary/secondary is usually arbitrary!
         */
        bool isSecondaryAlignment() const { 
           return (AlignmentFlag & SECONDARY) == SECONDARY; }
        /**
         * Alignment is part of a chimera. The choice
         * of representative/supplementary is arbitrary.
         */
        bool isSupplementaryAlignment() const { 
           return (AlignmentFlag & SUPPLEMENTARY) == SUPPLEMENTARY; 
        }
        /**
         * Alias for isSupplementaryAlignment
         */
        bool isSupplementary() const { 
           return (AlignmentFlag & SUPPLEMENTARY) == SUPPLEMENTARY; 
        }
        void unsetSupplementary() {
           AlignmentFlag &= ~SUPPLEMENTARY; 
        }
        /** 
         * @return true if alignment is part of read that satisfied paired-end resolution
         */
        bool IsProperPair(void) const;        
        /**
         * Test against PROPER_PAIR flag. This is usually set by an aligner.
         * For a negation test the result will include:
         *    1. Singleton
         *    2. pairs located to multi-references
         *    3. Inert size way out of range compared to the
         *       mean std (3 or 4 std away from mean?) This is
         *       all determined by the aligner.
         */
        bool isProperPair() const {
            return (AlignmentFlag & PROPER_PAIR) == PROPER_PAIR;
        }
        bool isNotProperPair() const {
            return !((AlignmentFlag & PROPER_PAIR) == PROPER_PAIR);
        }
        bool isImproperPair() const {
            return !((AlignmentFlag & PROPER_PAIR) == PROPER_PAIR);
        }

    ////// manipulate alignment flags ///////
    public:        
        void setAlignmentFlag(uint32_t flag) {
           AlignmentFlag = flag;
        }
        /** sets value of "PCR duplicate" flag
         */
        void SetIsDuplicate(bool ok);         
        void SetIsFailedQC(bool ok);          // sets value of "failed quality control" flag
        void SetIsFirstMate(bool ok);         // sets value of "alignment is first mate" flag
        void SetIsMapped(bool ok);            // sets value of "alignment is mapped" flag
        void SetIsMateMapped(bool ok);        // sets value of "alignment's mate is mapped" flag
        /**
         * Only update the flag. For full operation you need to call
         * makeUnmapped()
         */
        void setUnmapped() {
            AlignmentFlag |= UNMAPPED;
        }
        void setMapped() {
            AlignmentFlag &= ~UNMAPPED;
        }
        void setMateUnmapped() {
           AlignmentFlag |=  MATE_UNMAPPED;
        }
        void setMateMapped() {
           AlignmentFlag &= ~MATE_UNMAPPED;
        }
        /** 
         * sets value of "alignment mapped to reverse strand" flag
         */
        void SetIsReverseStrand(bool ok);     
        /**  
         * sets value of "alignment's mate mapped to reverse strand" flag
         */
        void SetIsMateReverseStrand(bool ok); 
        /** sets value of "alignment part of paired-end read" flag
         */
        void SetIsPaired(bool ok);            
        void setUnpaired() {
           AlignmentFlag &= (~PAIRED);
        }
        void SetIsPrimaryAlignment(bool ok);  // sets value of "position is primary alignment" flag
        /** sets value of "alignment is part of read that satisfied paired-end 
         *  resolution" flag
         */
        void SetIsProperPair(bool ok);        
        void setProperPair() {
           AlignmentFlag |=  PROPER_PAIR;
        }
        void setImproperPair() {
           AlignmentFlag &= ~PROPER_PAIR;
        }
        void SetIsSecondMate(bool ok);        // sets value of "alignment is second mate on read" flag

        // convenient constants octal number, first digit is 0
         static const uint32_t PAIRED              = 0x0001;
         static const uint32_t PROPER_PAIR         = 0x0002;
         static const uint32_t UNMAPPED            = 0x0004;
         static const uint32_t MATE_UNMAPPED       = 0x0008;
         static const uint32_t REVERSE_STRAND      = 0x0010;
         static const uint32_t MATE_REVERSE_STRAND = 0x0020;
         static const uint32_t READ_1              = 0x0040;
         static const uint32_t READ_2              = 0x0080;
         static const uint32_t SECONDARY           = 0x0100;
         static const uint32_t QC_FAILED           = 0x0200;
         static const uint32_t DUPLICATE           = 0x0400;
         static const uint32_t SUPPLEMENTARY       = 0x0800; 

    //////// tag data access methods //////////
    public:
        unsigned int getTagDataSize() const {
           return TagData.size();
        }
        bool isTagDataEmpty() const {
           return TagData.empty();
        }
        static bool isValidTagName(const string& tag) {
           return tag.size() == Constants::BAM_TAG_TAGSIZE && isalpha(tag[0]) && isalpha(tag[1]);
        }
        /**
         * @return true if p points to a valid tag followed by a type char
         *   and if array followed by a correct atomic type.
         */
        bool isValidTag(const char* p) const;
        /**
         * A valid array tag is one that after moving
         * tag_width to the right, p will land on the 
         * first char of another tag or on the '\0' terminator.
         * [TAG][B][T][L]{ data }
         * |
         * p
         * @param p is at start of the TAG.
         * @return true if the array tag pointed to by p
         *   is valid.  
         */
        bool isValidArrayTag(const char* p) const;
        bool isValidArrayTag(const string& tag) const;
        /** 
         * Adds a field to the BAM tags. Deprecated.
         * Does NOT modify an existing tag - use \link BamAlignment::EditTag() 
         *    \endlink instead.
         * @param[in] tag   2-character tag name
         * @param[in] type  1-character tag type such as Z, H, or i
         * @param[in] value data to store. When data is atomic type the end of 
         *    tag data will be the last byte. When the value is string then
         *    the end will be the null char. Actually the end will have 2 null
         *    chars, one for marking the end of string; the other is the internal
         *    representation of C++ string (since version 11).
         * @return \c true if the \b new tag was added successfully
         * @see \samSpecURL for more details on reserved tag names, 
         *     supported tag types, etc.
         */
        template<typename T> bool AddTag(const std::string& tag, const std::string& type, const T& value);
        /**
         * Version not returning bool.
         * Will throw exception.
         */
        template<typename T> void addTag(const std::string& tag, const std::string& type, const T& value);
        template<typename T> void addTag(const std::string& tag, const char type, const T& value);
        /**
         * Adds a numeric array field to the BAM tags. Deprecated.
         *
         *  Does NOT modify an existing tag - use \link BamAlignment::EditTag() \endlink instead.
         *
         *  @param[in] tag    2-character tag name
         *  @param[in] values vector of data values (type T) to store
         *  @return true if the new tag was added successfully
         *  @sa samSpecURL for more details on reserved tag names, supported tag types, etc.
         */
        template<typename T> bool AddTag(const std::string& tag, const std::vector<T>& values);
        /**
         * Version use exception without return value.
         */
        template<typename T> void addTag(const std::string& tag, const std::vector<T>& values);
        template<typename T> void addOrReplaceTag(const std::string& tag, const std::vector<T>& values) {
           if (hasTag(tag)) {
              removeTag(tag);
           }
           addTag(tag, values);
        }

        /** 
         *  Edits a BAM tag field. Deprecated.
         *
         *  If \a tag does not exist, a new entry is created.
         *  If tag exists then it will be removed first. So the same tag
         *  can store a different type of data.
         *
         *  @param tag[in]   2-character tag name
         *  @param type[in]  1-character tag type (must be "Z", "i", or "H")
         *     Z for string. H for byte array in Hex format.
         *  @param value[in] new data value
         *
         *  @return \c true if the tag was modified/created successfully
         *
         *  @see BamAlignment::RemoveTag()
         *  @see \samSpecURL for more details on reserved tag names, supported tag types, etc.
        */
        template<typename T> bool EditTag(const std::string& tag, const std::string& type, const T& value);
        /**
         * Morden version throwing exception.
         */
        template<typename T> void editTag(const std::string& tag, const std::string& type, const T& value);
        template<typename T> void editTag(const std::string& tag, const char type, const T& value);
        /**
         *  Edits a BAM tag field containing a numeric array.
         *
         *  If \a tag does not exist, a new entry is created.
         *
         *  \param tag[in]   2-character tag name
         *  \param value[in] vector of data values
         *
         *  \return \c true if the tag was modified/created successfully
         *  \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
         */
        template<typename T> bool EditTag(const std::string& tag, const std::vector<T>& values);
        /**
         * Not need to check. Always return. If not will throw exception.
         */
        template<typename T> void editTag(const std::string& tag, const std::vector<T>& values);

        // retrieves tag data
        /** 
         *  Retrieves the value associated with a BAM tag.
         *  There is no excetion, when the function return false
         *  it means no such tag.
         *
         *  @param tag[in]          2-character tag name
         *  @param destination[out] retrieved value, will clear previous values.
         *  @return true if found.
         *
        */
        template<typename T> bool GetTag(const std::string& tag, T& destination) const;
        /**
         * This version should work with both string and atomic types
         * but should not be used for vector type. By checking the 
         * second field of the returned value, there is not need to 
         * call hasTag() method. This will speed up the algorithm.
         * The parameter T can be any type that is wider than than
         * values stores. int32_t can get stored int8_t. So it is
         * save to use T::int32_t to get all integer types of different
         * length.
         * @throw exceptions if anything goes wrong: BamNotagException.
         * @return [value, hastag] tag value store under tag and whether 
         *    the tag is present or not. If no such tag then hastag field 
         *    will be false.
         */
        template<typename T> pair<T,bool> getTag(const std::string& tag) const;
        /**
         * String specialization.
         * @param destination will clear previous values and fill with tag value.
         * @return true if found
         */
        template<typename T> bool GetTag(const std::string& tag, std::vector<T>& destination) const;
        /**
         * get array tag.
         * New implementation should be faster.
         */
        template<typename T> bool getTag(const std::string& tag, std::vector<T>& destination) const;
        /**
         * @return the underlying tag data as a vector of type T
         *    If no such tag then return an empty container.
         */
        template<typename T> vector<T> getArrayTag(const std::string& tag) const;
        /**
         * More convenient version return value directly. If not found
         * for string then return empty string. For integer types
         * will return the smallest integer.
         * TODO: Need testing
         * More convenient method than the reference version
         * @return the tag string value if exists otherwise return empty string
         */
        string getStringTag(const std::string& tag) const;
        /** 
         * retrieves all current tag names
         *
         * When paired with GetTagType() and GetTag(), this method allows you
         * to iterate over an alignment's tag data without knowing the names (or types)
         * beforehand.

         * @return \c vector containing all tag names found (empty if none available)
         * \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
         */
        std::vector<std::string> GetTagNames(void) const;
        /** 
         * retrieves the SAM/BAM type-code for requested tag name
         * @param[in]  tag  2-character tag name
         * @param[out] type retrieved (1-character) type-code
         * @return \c true if found
         * \sa \samSpecURL for more details on reserved tag names, supported tag types, etc.
         */
        bool GetTagType(const std::string& tag, char& type) const;
        /**
         * retrieves the SAM/BAM type-code for the data elements in an array tag
         */
        bool GetArrayTagType(const std::string& tag, char& type) const;
        /** 
         * @returns true if alignment has a record for this tag name
         */
        bool HasTag(const std::string& tag) const;
        /**
         * New implementation
         * @return true if has this tag.
         */
        bool hasTag(const std::string& tag) const;
        /** 
         * Removes a tag if exists otherwise do nothing. 
         * If tag does not exist it is fine.
         */
        void removeTag(const std::string& tag);
        /**
         * Alias for removeTag
         * remove tag if exists; otherwise do nothing
         */
        void RemoveTag(const std::string& tag) {
            removeTag(tag);
        }

    ////// additional methods /////////
    public:
        /**
         * populates alignment string fields
         */
        bool BuildCharData(void);
        /** 
         *  Calculates alignment end position, based on its starting 
         *  position and CIGAR data.
         * 
         *  @warning The position returned now represents a zero-based,
         *      HALF-OPEN interval.  In previous versions of BamTools 
         *      (0.x & 1.x) all intervals were treated as zero-based, CLOSED.
         *  @param[in] usePadded Allow inserted bases to affect the reported
         *      position. Default is false, so that reported position
         *      stays synced with reference coordinates.
         *  @param[in] closedInterval Setting this to true will return a
         *      0-based end coordinate. Default is false, so that his value
         *      represents a standard, half-open interval.
         *  @return alignment end position on the reference.
         *      If unmapped, then return -1.
         */
        int GetEndPosition(bool usePadded = false, bool closedInterval = false) const;
        /**
         * @return end of the mapping index (0-based) on reference.
         *    This is the last mapped base index on reference. 
         *    Closed end of the interval.
         */
        int getEndPosition() const { return GetEndPosition(false, true); }
        /**
         * @return the end position of the read pair.
         *   if calling from this which is mapped to the forward
         *   strand, mate will be mapped to the reverse strand
         *   then the end of the end of mate mapping.
         *  if this read is the reverse strand, then the end is 
         *  this mapping's end.
         *  If the mate is mapped to a different chromosome
         *  then just the end position of this mapping for this read.
         */
        int getPairedEndPosition() const;
        /**
         * @return true if this alignment contains a closed end
         *   range [b,e] on the genomic reference coordinate.
         */
        bool contains(int b, int e) const {
            return getPosition() <= b && getEndPosition() >= e;
        }
        /**
         * @return true if p in in [b,e] on the refgerence coordinate
         */
        bool contains(int p) const {
            return getPosition() <= p && getEndPosition() >= p;
        }
        bool contain(int p) const {
            return getPosition() <= p && getEndPosition() >= p;
        }
        bool contain(int b, int e) const {
            return getPosition() <= b && getEndPosition() >= e;
        }

        /**
         * @return the [start, end] range of the mapping 
         * of reads on the reference in 0-based index.
         * start < end
         */
        std::pair<int,int> getRange() const { 
           return std::pair<int,int>(getPosition(), GetEndPosition(false, true)); 
        }
        /**
         * @return the 0-based close interval on the Query sequence.
         *    --===--
         */
        std::pair<int,int> getQInterval() const { 
           if (isForwardStrand()) {
            return std::pair<int,int>(getFirstSoftclipLength(), getQueryLength()-getLastSoftclipLength()-1); 
           }
           else {
            return std::pair<int,int>(getLastSoftclipLength(), getQueryLength()-getFirstSoftclipLength()-1); 
           }
        }

        /**
         * The first is always smaller than the second position.
         * @return the [begin,end] of the closed range for the mapping of this sequence
         *    on the genomic DNA. 0-based index. 
         *    If unmapped, then return [-1,-1].
         */
        std::pair<int,int> getInterval() const;
        bool sameInterval(const BamAlignment& ba) const;
        /**
         * @return the interval including soft-clipped region if they
         *   exist. A more generic version of getInterval() but 
         *   more costly. So only in special situation you need
         *   to use this version.
         */
        std::pair<int,int> getSoftInterval() const;
        /**
         * The distance coverted by the alignment on the reference.
         * Should be the same as the sum of cigar M+D length.
         * The validCigar() function check the above fact.
         * @return the match length on reference.
         */
        int getReferenceWidth() const {
           if (isUnmapped()) return 0;
           return GetEndPosition(false,false) - getPosition();
        }
        /**
         * For properly mapped pairs: +/- or -/+ configuration
         * <pre>
         * Get the range if it is paired on the same reference
         *    b1      e1     b2      e2
         * 1) ==this==> ---- <==mate==  return [b1, b2] we don't know e2
         *    b2      e2     b1      e1
         * 2) ==mate==> ---- <==this==   return [b2, e1] we know the whole range
         * Otherwise just the range of this alignment.
         * </pre>
         * b could be larger than e
         * Can get the mate alignment length from the MC tag
         */
        std::pair<int,int> getPairedRange() const;
        /**
         * This is the start to end regions coverted by the read.
         * Similar to getPairedRange but wihout the informaiton about the
         * direction of the two reads.
         *  |==R1==  ==R2==|
         *  |<------------>|
         *  b              e
         * Should always be [small, large]
         * |==R1====>|
         * |<=R2=    | R2 might have been trimmed due to low quality
         * @return [b,e] pair.
         */
        std::pair<int,int> getPairedInterval() const;

        /**
         * retrieves the size, read locations and reference locations of soft-clip operations
         */
        bool GetSoftClips(std::vector<int>& clipSizes,
                          std::vector<int>& readPositions,
                          std::vector<int>& genomePositions,
                          bool usePadded = false) const;

    //////// getter and informational methods //////
    // adding these getter methods for security
    // using public member is a security issue and may not pass FDA
    // quality
        /**
         * get the name of the query squence
         */
        const std::string& getQueryName() const { return Name; }
        /**
         * @return the name of the query.
         */
        const std::string& getName() const { return Name; }
        /**
         * Alias for shorter method getLength()
         * getter method for the length of the query sequence (read)
         * @return the query sequence length
         */
        int32_t getQueryLength() const { 
           return getLength(); 
        }
        /**
         * @return the length of the matched part of the
         *    the reference.  This is the sum of cigar
         *    string M 
         * The D segment is not counted.
         *  @see getReferenceWidth that count D
         */
        int getMatchedReferenceLength() const;
        /**
         * @return true if More match than softclip
         */
        bool matchDominate() const {
           return getMatchedReferenceLength() > getSoftclipLength();
        }
        /**
         * 'original' sequence (contained in BAM file)
         */
        const std::string& getQueryBases() const { return QueryBases; }
        std::string& getQueryBases() { return QueryBases; }
        /**
         * @return the query sequence same as getQueryBases.
         * If there is Hard-clip, the only the alignmed part will be 
         * returned. So should tell BWA never use hardclip.
         */
        const std::string& getQuerySequence() const { return QueryBases; }
        std::string& accessSequence() {
           return QueryBases;
        }
        std::string getRevcompQueryBases() const {
           return getRevcompQuerySequence();
        }
        std::string getRevcompQuerySequence() const;
        bool sameQuerySequence(const BamAlignment& ba) const {
           return getQueryBases() == ba.getQueryBases();
        }
        /**
         * aligned sequence (QueryBases plus deletion, padding, clipping chars)
         * Not good for sequence analysis, need a separate method to 
         * extract the aligned without indel symbols.
         * @return the formated aligned query sequence.
         *   if not present, could return empty string.
         * see: getMatchedQuerySequence()
         * TODO: Need a method to build the aligned bases when it is empty.
         */
        const std::string& getAlignedQueryBases() const { 
           //if (AlignedBases.empty()) {
           //   cerr << "empty aligned bases!\n";
           //}
           return AlignedBases; 
        }
        /**
         * Will not depend on whether the AlignedBases has been formated or not.
         * @return the portion of query sequences that matched the reference
         *    without any deletion characters. QueryBases excluding the 
         *    softclip sequences if they exist.
         */
        string getMatchedQuerySequence() const;
        /**
         * Query length excluding S,H clip. Only the M + I length.
         * alignment length - (D H or S). 
         * @return the sum length from M,I segments
         * The naming of this method may cause confusion when used
         * with getMatchedQueryLength(). The mirror image is
         * getReferenceWidth() that is the reference length
         * M+D. This method is M+I.
         */
        int getMatchedQueryLength() const;
        /**
         * Query bases in alignment excluding S, H, I, and D segments.
         */
        int numberBaseAligned() const;
        /**
         * Should call this function after update to CigarData
         */
        void clearAlignedBases() {
           AlignedBases.clear();
        }
        /**
         * Update AlignedBases.
         * This function needs to be called if the alignment
         * part of the object is changed.
         * For simplicity, the alignment can simply be cleared.
         * Most operations will not need this part of the computed
         * result. The result can be computed from the Cigar string.
         * The designer should not have included it as a member of the 
         * object in the first place.
         */
        void setAlignedBases(const std::string &alnedseq) { AlignedBases = alnedseq; }
        /** 
         * @return the FASTQ qualities (ASCII characters, not numeric values)
         * Values are ASCII 33-93
         */
        const std::string& getQuality() const { return Qualities; }
        std::string& getQuality() { return Qualities; }
        /**
         * @return the quality in reverse order.
         */
        string getReverseQuality() const;
        /**
         * @return the fastq quality as integer value from 0 to 93
         * This is the Phred score after int(ASCII) - 33
         */
        vector<int> getQualityScore() const;
        int getAverageQualityScore() const;
        /**
         * Make sure all score are not less than '!' which is 33 and 
         * not more than 126 ~
         */
        bool validQScore() const;
        /**
         * @return ID number for reference sequence
         * use this id and RefVector to get reference name.
         * The RefVector is only available in BamReader's
         * header section.
         * RefVector is vector<RefData>, RefData is { RefName, RefLength }
         * unnecessary typedef, more confusing than help.
         * RefID -> index of the RefVector
         */
        int32_t getReferenceId() const { return RefID; }
        /**
         * save some typing.
         */
        bool sameReferenceId(const BamAlignment& ba) const {
           return getReferenceId() == ba.getReferenceId();
        }
        bool sameReference(const BamAlignment& ba) const {
           return getReferenceId() == ba.getReferenceId();
        }
        /**
         * get the fist mapping position in 0-based index on the reference
         * sequence.
         * Soft clips at the beginning does not count, the 
         * position is that of the reference.
         * @return the 0-based index on the reference.
         *    In case of unmapped, -1 is returned.
         */
        int32_t getPosition() const { return Position; }
        /**
         * @return mapping quality.
         */
        int16_t getMapQuality() const { return MapQuality; }
        /**
         * @return the reference if of the mate
         */
        int32_t getMateReferenceId() const { return MateRefID; }
        /**
         * @return mate mapped position (left mapping point)
         */
        int32_t getMatePosition() const { return MatePosition; }
        int getMateEndPosition() const;
        /**
         * @return the strand as a single char + or -
         */
        char getMateStrandChar() const {
           if (IsMateReverseStrand()) return '-';
           else return '+';
        }
        /**
         * @return true if mate is mapped to the same reference.
         *   If no mate or mate is not mapped will return false.
         */
        bool mateOnSameReference() const { return MateRefID == RefID; }
        bool mateOnDifferentReference() const { return MateRefID != RefID; }
        /**
         * If read is unpaired, the insert size is also zero.
         * @return length of the template positive for
         *    first segment (left most) mate and negative for last (right most)mate.
         *    If the two reads map to different references then
         *    the insert size is set to zero.
         * Note: BamWriter use the InsertSize directly for output.
         */
        int32_t getInsertSize() const { 
           return InsertSize; 
        }
        /**
         * The returned value is >= 0. not sure better
         * to use unsigned int or not.
         *
         * @return the length of the InsertSize if it is not zero
         *   otherwise return the length of the projected length
         *   of the read on the reference. result is non-negative.
         *
         *  This function used getInsertSize()
         */
        int getTemplateLength() const;
        /**
         * @return a const reference to the CIGAR operations for this alignment
         * CigarData is a typedef of vector<CigarOp>
         * CigarOp: struc { char Type; uint32_t Length; }
         */
        const std::vector<CigarOp>& getCigar() const { return CigarData; } 
        std::vector<CigarOp>& getCigar() { return CigarData; } 
        /**
         * Provide a more user-friendly interface
         * for working with other applications.
         */
        vector<pair<char,int> > getCigarOperation() const;
        /**
         * @return a string version of Cigar.
         */
        string getCigarString() const;
        uint32_t getCigarHash() const;
        /**
         * @return true if the CigarData of this object is the same as 
         *    the argument cigar
         */
        bool sameCigar(const vector<pair<char,int>>& cigar) const;
        bool sameCigar(const BamAlignment& ba) const {
           return CigarData == ba.CigarData;
        }
        /**
         * Some alignment's cigar entry is *
         * this is the same as no cigar. This function test this situation
         * @return true if no cigar
         */
        bool lackCigar() const { return CigarData.empty(); }
        bool hasCigar() const { return !CigarData.empty(); }
        /**
         * There is no Deletion segment in the Cigar
         */
        bool lackDCigar() const;
        bool hasDCigar() const;
        /**
         * There is not Insertion segment in the Cigar
         */
        bool lackICigar() const;
        bool hasICigar() const;
        /**
         * @return true of the reference length (M+D) and query length (M+I) 
         *   computed from cigar are the same as computed from query length and
         *   END-BEGIN+1. 
         */
        bool validCigar() const;
        /**
         * @return the  Cigar operation type as a single char
         *    at index i.
         */
        char getCigarType(unsigned int i) const {
           return CigarData[i].getType();
        }
        /**
         * @return the length of the Cigar segment at 0-based index i
         */
        int getCigarLength(unsigned int i) const {
           return CigarData[i].getLength();
        }
        /**
         * Same as getCigarSize()
         */
        unsigned int getCigarOperationCount() const {
           return CigarData.size();
        }
        const CigarOp& getCigarOp(unsigned int i) const {
           return CigarData[i];
        }
        /**
         * @return Number of cigar operations
         */
        int getCigarSize() const {
           return CigarData.size();
        }
        int numberOfCigar() const {
           return CigarData.size();
        }
        /**
         * @return the number of indel + softclip
         */
        int numberOfIndelsoft() const;
        /**
         * @return the number of indels
         */
        int numberOfIndel() const;
        /**
         * Has indel near < 22 nt from the end
         */
        bool hasEndIndel() const;
        bool hasAmbiguousBase() const;
        /**
         * To fix certain aligner's tendency to put two gap
         * when a small region of the sequence has more mismatches
         * @deprecated use fix1M.
         */
        void fixStaggerGap();
        /**
         * Fix cigar with pattern like 57M1I1M1I69M 1M flank I or D
         * After this operation MD tag will be invalid. Need the
         * reference sequence to recalculate this tag. NM value is
         * also approximate. If want to be accurate, also need to
         * recalculate.
         * In the case of 15S2M2D59M, 57M2D1M23S, 15S2M2I59M, 95M1I2M15S,
         * the estimated MD tag may be close.
         * In the case of ...57M1I1M1I69M... editing MD tag is almost
         * impossible. So has to be recomputed.
         * NOTE: need to call recalMD() or recalMDsubseq() if you want the
         * object to be valid.
         */
        bool fix1M();
        /**
         * will convert 12I135M into 12S135M
         * Change ending I to S
         */
        void fixCigarError();

        /**
         * @return true if start with softclip
         */
        bool startWithSoftclip() const {
            return !CigarData.empty() && getCigar().front().Type == 'S';
        }
        /**
         * @return true if ends with softclip
         */
        bool endWithSoftclip() const {
            return !CigarData.empty() && getCigar().back().Type == 'S';
        }
        bool bothEndSoft() const {
            return !CigarData.empty() && getCigar().front().Type == 'S' && getCigar().back().Type == 'S';
        }
        bool bothEndsSoft() const {
            return !CigarData.empty() && getCigar().front().Type == 'S' && getCigar().back().Type == 'S';
        }
        /**
         * @param startR is the starting value for reference index.
         *    default is 0.  Can be any value. But the alternative
         *    meanful value is the absolute value of the genome
         *    reference such as the postion value for the alignment.
         * @return pair[cigar_index, bool]. The second value will be true
         *    if reference deisredR index has an insertion (attached to the
         *    last base of the previous M segment). The first value is the
         *    cigar segment index (0-based). If is not insertion, then, the
         *    cigar index is the one whose start index of the segment is
         *    greater than desiedR.
         */
        pair<int,bool> isInsertionAtRefloc(int desiredR, int startR=0) const;
        /**
         * @return pair[cigar_idx, bool] when desiredR is a deletion then
         *   the second value is true. The desiredR is the index of the 
         *   first base of the D segment.
         */
        pair<int,bool> isDeletionAtRefloc(int desiredR, int startR=0) const;

        /**
         * Alignment has soft clip on either start
         * or end. Not checking middle, which may not make sense.
         * @return true of has soft clip otherwise, including no cigar,
         *    then return false.
         */
        bool hasSoftclip() const;
        bool noSoftclip() const {
            if (CigarData.empty() || CigarData.size() == 1) return true;
            return CigarData.front().Type != 'S' && CigarData.back().Type != 'S';
        }
        /**
         * Calculate identity over the aligned part excluding
         * insert/delete and soft/hard clips.
         * @return ungapped identity.
         */
        float getNGIdentity() const;
        /**
         * This function requires the NM tag being present
         * NM is the mismath+del total base count.
         * @return number of mismatch, alngnment length excluding gaps
         */
        pair<int,int> getMismatchCount() const;
        /**
         * @return the number of exact matches for this alignment.
         */
        int getIdentical() const {
           return getAlignLength() - getNMValue();
        }
        /**
         * @return NM/alnlen as a fraction number that
         *   represents the local alignment identity.
         *   Soft clipped regions are not counted.
         */
        float getIdentity() const;
        /**
         * @return a normalized score comparable between different alignments.
         */
        float getScore() const {
           return getIdentity() * getAlignLength();
        }
        /**
         * Same as horizontal coverage of the query.
         * @return the portion of aligned query sequence.
         */
        float getFractionAligned() const {
            return static_cast<float>(getMatchedQueryLength())/getLength();
        }
        /**
         * @return horizontal coverage of the query between 0 and 1
         */
        float getQCoverage() const {
            return static_cast<float>(getMatchedQueryLength())/getLength();
        }
        pair<int,int> getMatchBound() const;
        int getQueryMatchBegin() const {
           if (CigarData.front().getType() == 'S' ||
               CigarData.front().getType() == 'H')
              return CigarData.front().getLength();
           return 0;
        }
        int getQueryMatchEnd() const {
           if (CigarData.back().getType() == 'S' ||
               CigarData.back().getType() == 'H')
              return getLength() - CigarData.front().getLength() - 1;
           return getLength()-1;
        }
        /**
         * @return the query sequence for the first soft clip.
         *    If there is no soft clip then an empty string is returned.
         */
        string getFirstSoftclip() const;
        string getFirstSoftquality() const;
        /**
         * @return the length of the first softclip.
         *   if not start with softclip, then return 0.
         */
        int getFirstSoftclipLength() const;
        /**
         * @return the last soft clip in query sequence.
         */
        string getLastSoftclip() const;
        string getLastSoftquality() const;
        /**
         * @return the length of the last softclip
         *   if no softclip return 0.
         */
        int getLastSoftclipLength() const;
        /**
         * @return sum of softclip length if both are present.
         */
        int getSoftclipLength() const;
        /**
         * @return the longer of the left and right softclip lengths.
         *   If no softlcip, the return 0
         */
        int getMaxSoftclipLength() const;
        /**
         * Remove the first soft clip so that the alignment
         * appears to be better. The query sequence will also 
         * be changed.
         */
        void chopFirstSoftclip();
        /**
         * Remove the last soft clip
         */
        void chopLastSoftclip();
        /**
         * Chop first and last softclip if they exist.
         */
        void chopSoftclip();

        /**
         * Same as chopFirstSoftclip() except for checking 
         * that first soft is off the start of the reference.
         */
        void chopDangleFrontSoft();
        /**
         * Chop soft clip that is off the end of the reference.
         *  ====== ref ref could be circular or assembly not long enough
         *      ==-- read
         */
        void chopDangleBackSoft();

        ///// setter methods ////
        /**
         * change the name of the query
         */
        void setQueryName(const std::string &qname) { Name = qname; }
        void setQuerySequenceLength(int32_t qlen) {  
           SupportData.QuerySequenceLength = qlen; 
        }
        void setQueryLength(int32_t qlen) {  
           SupportData.QuerySequenceLength = qlen; 
        }
        /**
         * Once query bases is changed the length will also
         * change. There is no need to call setQueryLength()
         * after calling this method.
         * The Cigar string may also need to be update.
         * If paired, the mate needs to set MC tag.
         */
        void setQueryBases(const std::string &qseq) { 
           QueryBases = qseq; 
           setQueryLength(QueryBases.size());
        }
        void setQueryBases(std::string &&qseq) { 
           QueryBases = std::move(qseq); 
           setQueryLength(QueryBases.size());
        }
        void setQuerySequence(const std::string &qseq) { 
           QueryBases = qseq; 
           setQueryLength(QueryBases.size());
        }
        void setQuerySequence(std::string &&qseq) { 
           QueryBases = std::move(qseq); 
           setQueryLength(QueryBases.size());
        }
        void appendQueryBases(const string& tail) {
            QueryBases.append(tail);
            SupportData.QuerySequenceLength += tail.size();
        }
        /**
         * At this point I have only implemented the operation on
         * unmapped object. For mapped, the Cigar, position will 
         * also need to be changed. Not sure this is meaningful or not.
         * After this operation, the object is useless.
         */
        void revcomp();
        /** 
         * set quality from string data. Quality and Query sequence
         * modification should always be done in sync.
         * */
        void setQuality(const std::string &qual) { Qualities = qual; }
        void setQuality(std::string &&qual) { Qualities = std::move(qual); }
        void appendQuality(const string& tail) {
            Qualities.append(tail);
        }
        /**
         * If the Qualities member is shorter than QueryBases
         * then will pad to the tail part and fill with the
         * score+33 ==> char
         */
        void fillQuality(int score);
        /** 
         * integer version.
         * Since the internal representation is with char, you cannot
         * use r_reference.
         * @param qual is a vector of Phred scores from 0-93. 
         *   This function will char(+33) 33-126 all are visible characters.
         */
        void setQuality(const std::vector<int> &qual);
        void setRefID(int32_t refid) { RefID = refid; }
        void setReferenceId(int32_t refid) { RefID = refid; }
        /**
         * Update the refgenome mapping start position
	      * If Cigar string needs to be updated then also do that.
         * Will also clear AlignedBases member.
         */
        void changePosition(int32_t alnstart); 
        /**
         * Simply update the Position member without changing anything
         * else.
         */
        void setPosition(int32_t alnstart) {
            Position = alnstart; 
        }
        /** alias for setPosition */
        void setStart(int32_t alnstart) { Position = alnstart; }
        void setBin(uint16_t indexbin) { Bin = indexbin; }
        void setMapQuality(uint16_t mqual) { MapQuality = mqual; }
        /**
         * update CigarData with a new value.
         * TODO: SupportData.NumCigarOperatons is redundant with
         *       CigarData.size() need to remove
         */
        void setCigarData(const std::vector<CigarOp> &cd) { 
           CigarData = cd; 
           SupportData.NumCigarOperations=cd.size();
        } 
        /** 
         * set CigarData from a vector a pair {char, int}
         *  for easy communication with external world.
         *  @param cd CigarData in more universal std data type
         * */
        void setCigarOperation(const std::vector<pair<char,int> > &cd); 
        void setCigar(const string& cstr);
        void setCigar(const std::vector<pair<char,int> > &cd) {
           setCigarOperation(cd);
        }
        void setCigar(const std::vector<CigarOp> &cd) {
           setCigarData(cd);
           SupportData.NumCigarOperations=cd.size();
        }
        /**
         * @return the sum of M,D,I segment length that is the
         *    total aligned length.
         */
        int getAlignLength() const;
        /**
         * @return the suma of M,I length of query aligned.
         */
        int getQueryAlignLength() const;
        /**
         * @return the suma of M,D length of reference aligned.
         */
        int getReferenceAlignLength() const;
        /**
         * @param materefid Mate reference id, set to -1 if mate unmapped
         */
        void setMateRefID(int32_t materefid)  {  MateRefID = materefid; } 
        void setMateRefid(int32_t materefid)  {  MateRefID = materefid; } 
        void setMateReferenceId(int32_t materefid)  {  MateRefID = materefid; } 
        /**
         * update mate position
         */
        void setMatePosition(int32_t matepos) { MatePosition = matepos; } 
        /**
         * Sets the insert size which is the 
         * length of the template. For left-most alignment, it is positive;
         * for right-most it is negative. This function will automaticall
         * adjust the sign of the value if the user did it wrong.
         * Only properly paired alignment will have non-zero values.
         * Thus the left most is also the forwardstrand, and the 
         * right most is the reverse strand.
         */
        void setInsertSize(int32_t insize); 
        /////// end of setter methods ///////
        //
        // mutation functions
        /**
         * helper function to iterate over the cigar object
         * @param i is the index on the reference. The index
         *    can be based on any value: 0-based on the first
         *    base of the subreference sequence, or the absolute
         *    index for the full chromosome. This function will
         *    make relative move.
         * @param j index on the query. Usually starting from 0
         *   for the first base of the query.
         * @param ci is the index for the cigar segment. Start
         *    from zero for the first cigar element.
         */
        //void nextCigar(int& i, int& j, unsigned int& ci) const;
        void nextCigar(int& i, int& j, int& ci) const;
        /**
         * @return index from reference to that of query.
         */
        int indexRef2Query(int ri) const;
        /**
         * Pick the subsequence based on the 0-based index
         * of query sequence. If the [b,e] fall in a 
         * Softclip region then the returned alignment may not
         * make sense. This function does not check for this
         * special case.
         */
        BamAlignment subsequence(int b, int e) const;
        /**
         * Use the reference (genomic) coordinate to pick subsequence
         * of the query alignment.
         * @return a new BamAlignment from b to e on the reference
         *   [b,e] is an closed end interval.
         */
        BamAlignment subsequenceByRef(int b, int e) const;
        /**
         * If [b,e] contains a deletion then the returned query substring
         * will be shorter than length([b,e]) 
         * @return the substring of the query sequence according
         * to closed range [b,e]
         */
        std::string substringByRef(int b, int e) const;
        /**
         * @param ri is the reference index [pos, endpos]
         * @return the char of query at index b according to the reference
         *    sequence index coordinate system.
         */
        char charAtByRef(int ri) const;
        /**
         * @param len length of the deletion from query.
         * @return true of index is a deletion of query sequence for
         *    length of len.
         */
        bool isDeletionAt(int ri, int len) const;
        /**
         * The position is followed by a insertion of particular
         *   sequence.
         * @param ri is the reference index.
         * @param seq is the sequence of last Base of M
         *   plus the entire bases of I.
         */
        bool isInsertionAt(int ri, const string& seq) const;
        /**
         * @return true if no syntax error in MD tag.
         */
        bool validMD() const;
        /** 
         * even position number, odd position letter last 3 digits
         * for bases 00 A, 01 C, 10 G, 11 T, 100 N, first bit del
         * MD:Z:20^A127 20 match, A in ref deleted, 127 match
         * MD:Z:108^TTCTAAGGCCAGCTCCTGCACC39 =>108 identical deletion of TT...ACC match of 39
         * MD tag ^AC to represent deletion of ref AC in query
         * MD: ^A0T means A in ref deleted, followed by ref base T mismatch
         *      0 is used to separate to consecutive position.
         *     0G154, 0 match, followed by G ref mismatch, then 154 match
         *     154G0, 154 match, ref G mismatch, 0 match
         * But insertion in query is not recorded.
         * Difference in base is represented by Base
         * identical residues are represented by number.
         * If the alignment does not have MD tag then
         * the return vectors will be empty. The caller must
         * test that.
         * Usually the first vector is one element more than the second one.
         * @return the mathed and difference sequence in reference
         *    in two vectors. 
        */ 
        pair<vector<int>, vector<string>> getMDArray() const;
        /**
         * @return the length of the reference width computed from MD tag
         * This number should be identical to that of computed from the
         * Cigar string getReferenceWidth()
         */
        int getMDWidth() const;
        /**
         * validCigar() method checks whether cigar agree with refwidth
         * This method checks MD agree with refwidth or not.
         */
        bool refwidthAgreeWithMD() const;
        string getSAString() const {
           return getReferenceName() + "," +
               to_string(getPosition()) + "," +
               string(1, getStrandChar()) + "," +
               getCigarString() + "," +
               to_string(getMapQuality()) + "," +
               to_string(getNMValue());
        }
        /**
         * @return bwa XA tag (chr, +/-pos, CIGAR, NM-value) 
         * for additional locations of this alignment.
         */
        string getXAString() const {
           return getReferenceName() + ","
              + string(1, getStrandChar()) + to_string(getPosition()) 
              + "," + getCigarString() + "," + to_string(getNMValue());
        }
        /**
         * QC function.
         */
        bool valid() const;
        /**
         * Update the MD tag with the value mdvec
         */
        void updateMDTag(const pair<vector<int>, vector<string>>& mdvec);
        /**
         * Given the refernce sequence for this object. The MD tag will 
         * be refreshed with the recomputed value.
         * @param refsq is the whole chromosome sequence.
         */
        void recalMD(const string& refsq);
        /**
         * @param refsq is a subsequence from [pos to endpos]
         */
        void recalMDSubseq(const string& refsq);
        /**
         * @return the difference betwee reference and the query
         *   in chopped region (before idx).
         *   NM value should be reduce by this amount plus 
         *   the sum of insertion length.
         */
        int chopMDBefore(int idx);
        /**
         * Helper method used by chopAfter()
         */
        int chopMDAfter(int idx);
        /**
         * Helper function.
         * @param mtag is either XM or XW. Actaully it can be any
         *   array tag. TODO: rename to a more generic one.
         */
         void chopMethyTagBefore(const string& mtag, int idx);
         void chopMethyTagAfter(const string& mtag, int idx);
        /**  
         * If align is + strand then insertSize also needs to be updated.
         * This method cannot change that, needs the mate information.
         * Usually done when both mates are available by the caller.
         * @param idx is the 0-based chromosome index and is retained
         *    in the resulting object.
         * The position will be changed. The mate needs to update mate position.
         */
        void chopBefore(int idx);
        /**
         * Remove alignment after idx, idx will be the last aligned base.
         *      idx
         * =====|----
         * =====|
         * If - strand insert size will be changed.
         * @param idx is the 0-based index on the reference.
         */
        void chopAfter(int idx);
        /**
         *      idx
         *    ==|===  -> { ==, |=== }
         *      + goes to the second object
         * @return two alignments by breaking this object into two.
         */
        pair<BamAlignment,BamAlignment> cut(int idx) const;
        /**
         * Helper method not tested
         */
        int countFrontMismatch(int len) const;
        int countBackMismatch(int len) const;
        /**
         * Helper method used by trimFront()
         */
        void chopFront(size_t len, int numMismatch);
        /**
         * Helper used by trimBack()
         * @param numMismatch is the number of mismatch in the region
         *   to be chopped between query and reference.
         */
        void chopBack(size_t len, int numMismatch);
        /**
         * Remove fuzzy end from the front of the alignment
         * @return true if trimming happened.
         */
        bool trimFront();
        /**
         * Remove fuzzy end from the back of the alignment
         * @return true if trimming happened.
         */
        bool trimBack();
        /**
         * @return the trimming status for front and back.
         *   The first if for front, second for back.
         */
        pair<bool,bool> trim();
        /**
         * Move deletion from oldloc to newloc
         * the sequence in between must be repeat sequence.
         */
        void moveDeletion(int oldloc, int newloc);
        void moveInsertion(int oldloc, int newloc);
        /**
         * Only change the bases near the end < 7 nt
         * if the Query base is different from the reference base
         * on both ends. MD tag will also needs to be updated.
         * This might be a bad idea and cause trouble.
         */
        void patchEnd();
        /**
         * Recalculate the value for tag NM
         * This is needed after taking a subsequence or 
         * NM tag is missing from the object.
         * @param refseq the reference sequence whole string
         *  not just the subsequence for this object.
         */
         void updateNMTag(const string& refseq);
         void reduceNMTag(int diff);
         /**
         * friendly wrapper for better code.
         * @return -1 if no tag
         */
         int getASValue() const;
         /**
          * Not sure limiting to uisnged int16 is good or not.
          * This is OK for short reads, but maybe too small.
          * Performance may not be better if using smaller types.
          * @return the NM tag value or -1 if not found NM tag
          */
         uint16_t getNMValue() const;
         /**
          * @return the length of the reference sequence for this alignment.
          */
         int getReferenceLength() const {
            return rsname[getReferenceId()].second;
         }
         bool nearReferenceEnd(int d) const {
            return abs(getReferenceLength() - getEndPosition()) < d;
         }
         bool nearReferenceBegin(int d) const {
            return getPosition() < d;
         }
         /**
          * @return a positive number of how long does the mate cover
          *   on the reference sequence based on the MC (mate cigar)
          *   tag value. 
          * Note: this method is not very efficient since it has to 
          * get the tag and parse the tag.
          */
         int getMateRefwidth() const;
         /**
            * Try to calucate mate reference width. If could not then return -1
            */
         int getMateReferenceWidth() {
            return getMateRefwidth();
         }
         /**
          * must call setRefvector() before using this
          * function.
          */
         const string& getReferenceName() const {
            if (rsname.empty()) {
               cerr << __FILE__ << ":" << __LINE__ << ": rsname not loaded need to call setRefvector(BamReader::getReferenceMetaData())\n";
               throw runtime_error(string(__func__) + ": empty rsname, may need to call setRefvector()");
            }
            return rsname.at(getReferenceId()).first;
         }
         /**
          * Make this alignment into unaligned status.
          * flag: ummapped, mate_unmapped, improper_pair
          * major field: refid, position, mate_refid, mate_position
          *              CigarData, AlignedBases
          *      either set to -1 or clear for string values.
          * tag: NM, MD, MC removed
          * If reverse strand, then need to revcomp both base and quality
          */
         void makeUnmapped();
         /**
          * Will do all the operation as makeUnmapped() except for
          * refid, and position. After this operation the refid,pos
          * are left unchanged. This is useful if you want to 
          * keep the two mates next to each other in position sorted
          * bam files.
          */
         void markUnmapped();
         void makeMateUnmapped();
         int32_t     Length() const {
           return SupportData.QuerySequenceLength;
         }             
         /**
          * @return the query sequence length
          */
         int32_t length() const {
            return SupportData.QuerySequenceLength;
         }             
         /**
          * @return the length of the query sequence.
          */
         int32_t getLength() const;
         /**
          * length setter function.
          * This should only be use in conjunction of query base 
          * insertion and deletion.
          */
         void setLength(int32_t len) {
            SupportData.QuerySequenceLength=len;
            if ((int)QueryBases.size() != len) {
               cerr << __FILE__ << ":" << __LINE__ << ":WARN QueryBase length being changed\n";
               QueryBases.resize(len);
            }
         }

         ////// debug function ///

         void showTagData(ostream& ous) const;

         ////// static methods ////

         static void setPolishMax(int len) {
            TRIMLEN_MAX = len;
         }
         static void setPolishGap(int gap) {
            GAP_CUT = gap;
         }
         /**
          * Helper function to be used by other functions to convert
          * string version of cigar to vector version.
          * @param cigarstr cigar string.
          * @return vector of cigar operations [code,length]
          */
         static vector<pair<char,int>> parseCigar(const string& cigarstr);
         /**
          * convert cigar data to string version.
          */
         static string cigarToString(const vector<pair<char,int>>& cg);
         /**
         * The header part of BamReader
         * for indexing readname,readlength.
         * Bam file refid, matid directly look up into this table
         * and you can get REFNAME (first) and seqlen (second)
         * You can switch during BamAlignment processing to a different
         * file header.
         * This method will load refname2id and rsname.
         * @param refvec is the refernece vector to be obtained from 
         *    BamReader::getReferenceMetaData()
         */
         static void setRefvector(vector<pair<string,int>>&& refvec);
         /**
          * Before using this function you must load the reference look up table
          * with: setRefvector(BamReader::getReferenceMetaData())
          * You need a specific object to use getReferenceMetaData().
          * @name is the referece string name such as 1 or chr1 depends on 
          * the version of reference genome used.
          * @return the reference id given name for looking into 
          *   the id used by Bamfile.
          */
         static int referenceIdFromName(const string& name) {
            assert(!refname2id.empty());
            return refname2id[name];
         }
         /**
          * Given a reference id, what's the string name.
          * This is used to extract the genomic sequence
          * from external sources.
          */
         static const pair<string,int>& getRefnameFromId(int refid) {
            return rsname[refid];
         }

    private:
         //// static members to be shared by derived class ////
         /**
          * [refname, reflength] indexed on refid
          */
         static vector<pair<string,int>> rsname;
         /**
          * a quick look up table from refname to refid
          */
         static map<string,int> refname2id;
         /** 
          * Helper function apply to generic situation where
          * both Ref or Query may be advanced.
          *
          * move i to position b if i is not at b
          * j to the start of the query index
          * @param i index in reference seq 0-based index
          * @param j index in query sequence 0-based index
          * @param x desired index to advance i to on refseq. If x happen
          *     to reside in a Deletion, then x will be automatically
          *     advanced to the first base of the next M. Not sure
          * @param ci cigar segment index [0, NumberOfCigar-1]
          */
         //void advanceIndex(int &i, int &j, int &x, unsigned int &cigarIdx, unsigned int &ci) const;
         // above new implementation did not pass test yet
         void advanceIndex(int &i, int &j, int &x, unsigned int &cigarIdx, unsigned int &ci, char &cigarState) const;

    // public data fields, TODO: these fileds should all become private in the future
    public:
        /**
         * The standard bases and their complements. A->T, C->G, G->C, T->A, ....,
         * N->N
         */
        static map<char,char> complementBase;
        /** 
         * read or query name 
         * Use getName() to read this one
         * */
        std::string Name;    
        /** 'original' sequence (contained in BAM file)
         */
        std::string QueryBases;         
        /** 
         * not sure this is useful 
         * 'aligned' sequence (QueryBases plus deletion, padding, clipping
         * chars) This field will be completely empty after reading from
         * BamReader/BamMultiReader when QueryBases is empty.
         * After any modifcation of the object on the query sequences
         * this field should be cleared.
         */
        std::string AlignedBases;       
        /** 
         * FASTQ qualities (ASCII characters, not numeric values)
         * String representation from char 33 to 93 in visible
         * range. Score = int value of char - 33
         * Qscore of 30 is char 63 '?'
         */
        std::string Qualities;          
        /** 
         * [TAG][T]{ data } for atomic type where data length is sizeof(T)
         * for array
         * [TAG][B][T][L]{ data } data length is sizeof(T)*L
         * TAG is two bytes long, BAM_TAG_TAGSIZE 2
         * T is a single char (byte) one of BAM_TAGTYPE_xxx
         *    length is defined by BAM_TAG_TYPESIZE 1
         * L is sizeof(int32_t) 4
         * BAM_TAG_ARRAYBASE_SIZE=8 is the array metadata length.
         * tag data (use provided methods to query/modify)
         * String encoding of everything. Should consider using
         * a proper data structure. 
         * TODO: rewrite this part.
         */
        std::string TagData;            
        /** 
         * ID number for reference sequence
         *  -1 for unmapped reads
         */
        int32_t     RefID;              
        /**
         * (0-based) where alignment starts on reference.
         * If unmapped then this value can be set to -1
         */
        int32_t     Position;
        /** 
         * BAM (standard) index bin number for this alignment
         */
        uint16_t    Bin;                
        uint16_t    MapQuality;         // mapping quality score
        /**
         * This is sam/bam file field #2 containing 12 bit information
         * alignment bit-flag. 
         * use the provided methods to query/modify.
         */
        uint32_t    AlignmentFlag;
       /**  
        * CIGAR operations for this alignment.
        * CigarOp has Type,Length public field
        * Defined in BamAux.h
        */
        std::vector<CigarOp> CigarData; 
        /** 
         * ID number for reference sequence where alignment's mate was aligned.
         * If unpaired -1. For sorting purpose, if mate is unmapped but the current
         * alignment is mapped, the mateRefId will be same as this one.
         */
        int32_t MateRefID;          
        /** position (0-based) where alignment's mate starts.
         *  SEELF position [a, b], MATE[c, d]
         *  If mate on the same reference, [a, b] and [c, d]
         *  could be identical.
         *  If single read or unampped then -1 to indicate the null situation.
         *  */
        int32_t MatePosition;       
        /**
         * Field 9: TLEN in bam/sam format.  signed observed Template
         *          LENgth.
         *     Observed template length equals the number of bases
         *     from the leftmost mapped base to the rightmost mapped base.
         *     It is set as zero for singl-segment template or when the
         *     information is unavailable.
         *
         * Note: even when I set this member to > 0, when the 
         * read is marked as unpaired, the BamWriter will still
         * output the filed as ZERO. Have not figured why.
         *
         * This value may be misleading in computing real fragment
         * length of the read. Need another function.
         * For both reads in the negative strand the insert size is wrong
         * <pre>
         * BWA problem:
         * <--R1-- <--R2--
         * |--I----| insert size is the number of bases in the range
         *       |--I----| should the the correct value
         * </pre>
         * Need some correction for the -/- configuration
         * Then you can get the range properly
         */
        int32_t InsertSize;        

    //! \internal
    // internal utility methods
    private:
        /**
         * If not type of char, int then will throw exception.
         * @return the elemental type length
         */
        short getAtomicTagLength(const char tagt) const;
        /**
         * Only work with Atomic, String, or Hex types.
         * [TAG][T]{ data } data length is the sizeof(T)
         *      |
         *      p
         * @param p poiner at tag type position, 1 before tag data string.
         *     p at 3rd char.
         * @return the length of the data length in number of chars.
         *    The type char is not included.
         */
        size_t getBasicTagLength(const char* p) const;
        /**
         * [TAG:2][B:1][T:1][L:4]{ data }
         *             |<--   result -->|
         *             p
         * @param p at the array-element type char position (4th char)
         * @return the number of char from the array-element type char to the end of array.
         *   including the element type char.
         */
        size_t getArrayTagLength(const char* p) const;
        /**
         * @param p is at the first byte for the Tag data which is 
         *    [TAG][TYEP]{  data } 
         *    [TAG][B][T][LENGTH]{ array data }  for array
         *    |<---                        -->|
         *    p points to the first byte
         * @return the full width of a tag pointed to by the pointer
         */
        size_t getTagWidth(const char* p) const;
        /**
         * This is much more efficient implementation and design
         * compared to FindTag().
         * @param tag is the two letter tag name such as MD
         * @return pointer to the first char of the desired tag
         *   if success otherwise return nullptr.
         */
        const char* findTag(const std::string& tag) const;
        char* findTag(const std::string& tag);
        /** 
         *  Searches for requested tag in BAM tag data.
         *  @param  tag            requested 2-character tag name
         *  @param  pTagData       pointer to current position in BamAlignment::TagData
         *     After the call it will be pointing to the first char for tag
         *     in TagDATA.
         *  @param  tagDataLength  length of BamAlignment::TagData
         *  @param  numBytesParsed number of bytes parsed so far before TAG
         *  @return true if found
         *  If tag is found, pTagData will point to the byte where the tag data begins.
         *        numBytesParsed will correspond to the position in the full TagData string.
        */
        bool FindTag(const std::string& tag, char*& pTagData, 
              const unsigned int& tagDataLength, unsigned int& numBytesParsed) const;
        bool IsValidSize(const std::string& tag, const std::string& type) const;
        /** 
         *  Moves to next available tag in tag data string such that
         *  pTagData will be at the first char of TAG.

         *  @param[in]     storageType    BAM tag type-code that determines how far to move cursor
         *  @param[in,out] pTagData       pointer to current position (cursor) in tag string (data).
         *                    |<-pTagData is here
         *      |TAG-2|Type-1|Data-determined by Type|
         *      in case of array tag, at the element type char; regardless at 4th char.
         *      After the operation, pTagData will be at first char of next TAG.
         *  @param[in,out] numBytesParsed report of how many bytes were parsed (cumulatively)
         *  @return \c if storageType was a recognized BAM tag type

         *  \post \a pTagData       will point to the byte where the next tag data begins.
         *        \a numBytesParsed will correspond to the cursor's position in the full TagData string.
         */
        //bool SkipToNextTag(const char storageType, char*& pTagData, unsigned int& numBytesParsed) const;
        bool SkipToNextTag(char*& pTagData, unsigned int& numBytesParsed) const;
         /**
          * Max length from either end to to polishing
          * default 6
          */
         static int TRIMLEN_MAX; // = 6;
         /**
          * Max length between mismatch to do polishing
          * default 3
          */
         static int GAP_CUT; // = 3;

    // internal data member
    private:
        // nested class TODO: simplify in future versions
        struct BamAlignmentSupportData {
            /** 
             * cigarop, tags are stored here
             * Method to populate this field:
             *    BamReaderPrivate::LoadNextAlignment(BamAlignment& alignment) 
             */
            std::string AllCharData;
            uint32_t    BlockLength;  // not sure what this is
            /**
             * TODO: make sure This is the same as Cigar.size()?
             */
            uint32_t    NumCigarOperations; // should be calculated on the fly
            uint32_t    QueryNameLength;  // duplicate data discard in the future
            /**
             * TODO: this is rudundant with QueryBases.size()
             * consider remove this field in the future.
             */
            uint32_t    QuerySequenceLength; // is this duplicate of QueryLength?
            bool        HasCoreOnly;
            
            // constructor
            BamAlignmentSupportData()
                : BlockLength(0) , NumCigarOperations(0)
                , QueryNameLength(0)
                , QuerySequenceLength(0)
                , HasCoreOnly(false)
            { }
            BamAlignmentSupportData(const BamAlignmentSupportData& o) 
               : AllCharData(o.AllCharData), BlockLength(o.BlockLength),
                 NumCigarOperations(o.NumCigarOperations),
                 QueryNameLength(o.QueryNameLength),
                 QuerySequenceLength(o.QuerySequenceLength),
                 HasCoreOnly(o.HasCoreOnly) { }
            BamAlignmentSupportData(BamAlignmentSupportData&& o) 
               : AllCharData(std::move(o.AllCharData)), BlockLength(o.BlockLength),
                 NumCigarOperations(o.NumCigarOperations),
                 QueryNameLength(o.QueryNameLength),
                 QuerySequenceLength(o.QuerySequenceLength),
                 HasCoreOnly(o.HasCoreOnly) { }
            BamAlignmentSupportData& operator=(const BamAlignmentSupportData& o) {
               if (this != &o) {
                 AllCharData=o.AllCharData;
                 BlockLength=o.BlockLength;
                 NumCigarOperations=o.NumCigarOperations;
                 QueryNameLength=o.QueryNameLength;
                 QuerySequenceLength=o.QuerySequenceLength;
                 HasCoreOnly=o.HasCoreOnly; 
               }
               return *this;
            }
            BamAlignmentSupportData& operator=(BamAlignmentSupportData&& o) {
               if (this != &o) {
                 AllCharData=std::move(o.AllCharData);
                 BlockLength=o.BlockLength;
                 NumCigarOperations=o.NumCigarOperations;
                 QueryNameLength=o.QueryNameLength;
                 QuerySequenceLength=o.QuerySequenceLength;
                 HasCoreOnly=o.HasCoreOnly; 
               }
               return *this;
            }
        };
        /**
         * Terrible data duplication! REMOVE.
         */
        BamAlignmentSupportData SupportData;
        friend class Internal::BamReaderPrivate;
        friend class Internal::BamWriterPrivate;
    //! \endinternal
};

// ---------------------------------------------------------
// BamAlignment tag access methods implementation

// should use a data structure other than playing with fire here
// The optimization is useless.  TODO: replace with a proper 
// data structure
template<typename T>
inline bool BamAlignment::AddTag(const std::string& tag, const std::string& type, const T& value) 
{
    // if char data not populated, do that first
    if (SupportData.HasCoreOnly)
        BuildCharData();
    // check tag/type size
    if (!IsValidSize(tag, type)) {
        // TODO: set error string?
        return false;
    }
    // check that storage type code is OK for T
    //if ( !TagTypeHelper<T>::CanConvertTo(type.at(0)) ) {
    if ( !TagTypeHelper<T>::CanConvertTo(type.front()) ) {
        // TODO: set error string?
        return false;
    }
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    // if tag already exists, return false
    // use EditTag explicitly instead
    if (FindTag(tag, pTagData, tagDataLength, numBytesParsed)) {
        // TODO: set error string?
        return false;
    }
    // otherwise, convert value to string
    union { T value; char valueBuffer[sizeof(T)]; } un;
    un.value = value;
    // copy original tag data to temp buffer
    const std::string newTag = tag + type;
    const size_t newTagDataLength = tagDataLength + newTag.size() + sizeof(T); // leave room for new T
    RaiiBuffer originalTagData(newTagDataLength); // will self destruct
    memcpy(originalTagData.Buffer, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
    // append newTag
    strcat(originalTagData.Buffer + tagDataLength, newTag.data());
    memcpy(originalTagData.Buffer + tagDataLength + newTag.size(), un.valueBuffer, sizeof(T));
    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData.Buffer;
    TagData.assign(newTagData, newTagDataLength); // will automatically add NULL char at end
    return true;
}

// type is string version
template<typename T>
inline void BamAlignment::addTag(const std::string& tag, const std::string& type, const T& value) 
{
    // if char data not populated, do that first
    if (SupportData.HasCoreOnly)
        BuildCharData();
    // check tag/type size
    if (!isValidTagName(tag)) {
       throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR Invalid tag name to add: " + tag);
    }
    // check that storage type code is OK for T
    //if ( !TagTypeHelper<T>::CanConvertTo(type.at(0)) ) {
    //if (!TagTypeHelper<T>::CanConvertTo(type.front())) {
    //if (!TagTypeHelper<T>::canStore(type.front(), value)) {
    if (!canStore<T>(type.front(), value)) {
       cerr << *this << endl;
       throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR type "
             + type + " cannot store value " + to_string(value) + " in tag: " + tag);
    }
    char *p = findTag(tag);
    if (p != nullptr) {
       throw logic_error("Tag " + tag + " already exists cannot add again");
    }
    size_t prevTDL = TagData.size(); // previous TagData length
    TagData.resize(prevTDL + Constants::BAM_TAG_TAGSIZE + Constants::BAM_TAG_TYPESIZE +sizeof(T));
    // p move to TAG NAME filed to write
    p = TagData.data() + prevTDL;
    memcpy(p, tag.c_str(), Constants::BAM_TAG_TAGSIZE);
    p += Constants::BAM_TAG_TAGSIZE;
    memcpy(p, type.c_str(), Constants::BAM_TAG_TYPESIZE);
    //short typeLen = getAtomicTagLength(*p);
    short typeLen = getAtomicTagLength(*p);
    if (typeLen != sizeof(T)) {
       //throw logic_error("write more code to deal with addTag type and T not same length");
       cerr << __FILE__ << ":" << __LINE__ << ":WARN T type " << typeid(T).name() << " " << value << " length differ from " << type 
          << " but small enough can store in stated type" << endl;
    }
    p += Constants::BAM_TAG_TYPESIZE;
    if (!TagTypeHelper<T>::CanConvertTo(type.front())) {
       //memcpy(p, &value, sizeof(T)); // canStore() already rulled out value too large 
       storeToAs(value, p, type.front());
    }
    else { // T can be converted to type
       storeToTag(p, value, type.front());
       /*
       if (type.front() == Constants::BAM_TAG_TYPE_INT8) {
          int8_t tmpv = value;
          memcpy(p, &tmpv, 1);
       }
       else if (type.front() == Constants::BAM_TAG_TYPE_UINT8) {
          uint8_t tmpv = value;
          memcpy(p, &tmpv, sizeof(uint8_t));
       }
       else if (type.front() == Constants::BAM_TAG_TYPE_INT16) {
          int16_t tmpv = value;
          memcpy(p, &tmpv, sizeof(int16_t));
       }
       else if (type.front() == Constants::BAM_TAG_TYPE_UINT16) {
          uint16_t tmpv = value;
          memcpy(p, &tmpv, sizeof(uint16_t));
       }
       else if (type.front() == Constants::BAM_TAG_TYPE_INT32) {
          int32_t tmpv = value;
          memcpy(p, &tmpv, sizeof(int32_t));
       }
       else if (type.front() == Constants::BAM_TAG_TYPE_UINT32) {
          uint32_t tmpv = value;
          memcpy(p, &tmpv, sizeof(uint32_t));
       }
       else if (type.front() == Constants::BAM_TAG_TYPE_FLOAT) {
          float tmpv = value;
          memcpy(p, &tmpv, sizeof(float));
       }
       else if (type.front() == Constants::BAM_TAG_TYPE_STRING) {
          throw logic_error("strin as template specializaiton should not get here");
       }
       else {
          throw logic_error("write more code for new tag type: " + type);
       }
       */
    }
}

template<typename T>
inline void BamAlignment::addTag(const std::string& tag, const char type, const T& value) 
{
    // if char data not populated, do that first
    if (SupportData.HasCoreOnly)
        BuildCharData();
    // check tag/type size
    if (!isValidTagName(tag)) {
       throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR Invalid tag name to add: " + tag);
    }
    if (!canStore<T>(type, value)) {
       cerr << *this << endl;
       throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR type "
             + string(1, type) + " cannot store value " + to_string(value) + " in tag: " + tag);
    }
    char *p = findTag(tag);
    if (p != nullptr) {
       throw logic_error("Tag " + tag + " already exists cannot add again");
    }
    size_t prevTDL = TagData.size(); // previous TagData length
    TagData.resize(prevTDL + Constants::BAM_TAG_TAGSIZE + Constants::BAM_TAG_TYPESIZE +sizeof(T));
    // p move to TAG NAME filed to write
    p = TagData.data() + prevTDL;
    memcpy(p, tag.c_str(), Constants::BAM_TAG_TAGSIZE);
    p += Constants::BAM_TAG_TAGSIZE;
    //memcpy(p, type.c_str(), Constants::BAM_TAG_TYPESIZE);
    short typeLen = getAtomicTagLength(type);
    *p = type; // ++p; 
    p += Constants::BAM_TAG_TYPESIZE;
    if (typeLen != sizeof(T)) {
       //throw logic_error("write more code to deal with addTag type and T not same length");
       cerr << __FILE__ << ":" << __LINE__ <<  ":WARN addTag() type " << type << " and T " << typeid(T).name() << " not same length. Will try the best" << endl;
    }
    if (!TagTypeHelper<T>::CanConvertTo(type)) {
       //cerr << "tryont to convert " << typeid(T).name() << " to " << type << endl;
       //memcpy(p, &value, sizeof(T));
       storeToAs(value, p, type);
    }
    else {
       storeToTag(p, value, type);
       /*
       if (type == Constants::BAM_TAG_TYPE_INT8) {
          int8_t tmpv = value;
          memcpy(p, &tmpv, 1);
       }
       else if (type == Constants::BAM_TAG_TYPE_UINT8) {
          uint8_t tmpv = value;
          memcpy(p, &tmpv, 1);
       }
       else if (type == Constants::BAM_TAG_TYPE_INT16) {
          int16_t tmpv = value;
          memcpy(p, &tmpv, 2);
       }
       else if (type == Constants::BAM_TAG_TYPE_UINT16) {
          uint16_t tmpv = value;
          memcpy(p, &tmpv, 2);
       }
       else if (type == Constants::BAM_TAG_TYPE_INT32) {
          int32_t tmpv = value;
          memcpy(p, &tmpv, 3);
       }
       else if (type == Constants::BAM_TAG_TYPE_UINT32) {
          uint32_t tmpv = value;
          memcpy(p, &tmpv, 3);
       }
       else if (type == Constants::BAM_TAG_TYPE_FLOAT) {
          float tmpv = value;
          memcpy(p, &tmpv, sizeof(float));
       }
       else if (type == Constants::BAM_TAG_TYPE_STRING) {
          throw logic_error("strin as template specializaiton should not get here");
       }
       else {
          throw logic_error("write more code for new tag type: " + string(1, type));
       }
       */
    }
}

// string specialization
template<> inline bool BamAlignment::AddTag<std::string>(const std::string& tag,
                              const std::string& type, const std::string& value)
{
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();
    // check tag/type size
    if ( !IsValidSize(tag, type) ) {
        // TODO: set error string?
        return false;
    }
    // check that storage type code is OK for string
    if ( !TagTypeHelper<std::string>::CanConvertTo(type.at(0)) ) {
        // TODO: set error string?
        return false;
    }
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    // if tag already exists, return false
    // use EditTag explicitly instead
    if (FindTag(tag, pTagData, tagDataLength, numBytesParsed)) {
        // TODO: set error string?
        return false;
    }
    // otherwise, copy tag data to temp buffer
    const std::string newTag = tag + type + value;
    const size_t newTagDataLength = tagDataLength + newTag.size() + 1; // leave room for null-term
    RaiiBuffer originalTagData(newTagDataLength);
    memcpy(originalTagData.Buffer, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
    // append newTag (removes original null-term, then appends newTag + null-term)
    strcat(originalTagData.Buffer + tagDataLength, newTag.data());
    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData.Buffer; // cast to const char*
    TagData.assign(newTagData, newTagDataLength); // TagData \0\0. last null char is at TagData.size()
    return true;
}

// string specialization
// can use this function to add Hex string, need to use type="H"
template<> inline void BamAlignment::addTag<std::string>(
      const std::string& tag, const std::string& type, const std::string& value)
{
    // if char data not populated, do that first
    if (SupportData.HasCoreOnly)
        BuildCharData();
    // check tag/type size
    if (!isValidTagName(tag)) {
       throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR Invalid tag name to add: " + tag);
    }
    if (!TagTypeHelper<string>::CanConvertTo(type.front())) {
       throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR type "
             + type + " is not Hex or string type when tag: " + tag);
    }
    char *p = findTag(tag);
    if (p != nullptr) {
       throw logic_error("Tag " + tag + " already exists cannot add again");
    }
    size_t prevTDL = TagData.size(); // previous TagData length
    TagData.resize(prevTDL + Constants::BAM_TAG_TAGSIZE + Constants::BAM_TAG_TYPESIZE + value.size() + 1);
    // p move to TAG NAME filed to write
    p = TagData.data() + prevTDL;
    memcpy(p, tag.c_str(), Constants::BAM_TAG_TAGSIZE);
    p += Constants::BAM_TAG_TAGSIZE;
    memcpy(p, type.c_str(), Constants::BAM_TAG_TYPESIZE);
    p += Constants::BAM_TAG_TYPESIZE;
    memcpy(p, value.c_str(), value.size()); // copy the trailing null char terminator of value
    *(p+value.size()) = '\0';
}

template<> inline void BamAlignment::addTag<std::string>(
      const std::string& tag, const char type, const std::string& value)
{
    // if char data not populated, do that first
    if (SupportData.HasCoreOnly)
        BuildCharData();
    // check tag/type size
    if (!isValidTagName(tag)) {
       throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR Invalid tag name to add: " + tag);
    }
    if (!TagTypeHelper<string>::CanConvertTo(type)) {
       throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR type "
             + string(1,type) + " is not Hex or string type when tag: " + tag);
    }
    char *p = findTag(tag);
    if (p != nullptr) {
       throw logic_error("Tag " + tag + " already exists cannot add again");
    }
    size_t prevTDL = TagData.size(); // previous TagData length
    TagData.resize(prevTDL + Constants::BAM_TAG_TAGSIZE + Constants::BAM_TAG_TYPESIZE + value.size() + 1);
    // p move to TAG NAME filed to write
    p = TagData.data() + prevTDL;
    memcpy(p, tag.c_str(), Constants::BAM_TAG_TAGSIZE);
    p += Constants::BAM_TAG_TAGSIZE;
    *p = type; 
    p += Constants::BAM_TAG_TYPESIZE;
    // the C++ string may not have the trailing \0 we need to add it 
    // explicitly
    memcpy(p, value.c_str(), value.size()); 
    *(p+value.size()) = '\0';
}

// implementation had a bug, fixed, will deprecate this version
template<typename T>
inline bool BamAlignment::AddTag(const std::string& tag, const std::vector<T>& values) {
    // if char data not populated, do that first
    if (SupportData.HasCoreOnly) BuildCharData();
    // check for valid tag name length
    if (tag.size() != Constants::BAM_TAG_TAGSIZE) {
        cerr << __FILE__ << ":" << __LINE__ << ":WARN tag " + tag + " TAGSIZE wrong\n";
        return false;
    }
    // localize the tag data, work as char pointer type
    char* pTagData = (char*)TagData.data();
    const unsigned int prevTDLen = TagData.size(); // previous TagData length
    unsigned int numBytesParsed = 0;

    // if tag already exists, return false
    // use EditTag explicitly instead
    if (FindTag(tag, pTagData, prevTDLen, numBytesParsed)) {
        cerr << __FILE__ << ":" << __LINE__ << ":WARN tag " + tag + " already exists. Cannot add another.\n";
        // TODO: set error string?
        return false;
    }
    // build new tag's base information, BAM_TAG_ARRAYBASE_SIZE=8
    char newTagBase[Constants::BAM_TAG_ARRAYBASE_SIZE];
    // c_str() is '\0' terminated
    std::memcpy(newTagBase, tag.c_str(), Constants::BAM_TAG_TAGSIZE); // BAM_TAG_TAGSIZE=2
    newTagBase[Constants::BAM_TAG_TAGSIZE] = Constants::BAM_TAG_TYPE_ARRAY; // 'B'
    newTagBase[Constants::BAM_TAG_TAGSIZE+1] = TagTypeHelper<T>::TypeCode(); // uint32_t => 'I'

    // add number of array elements to newTagBase
    const int32_t numElements  = values.size();
    std::memcpy(newTagBase + 4, &numElements, sizeof(int32_t)); // 32 bits is 4 bytes
    // copy current TagData string to temp buffer, leaving room for new tag's contents
    const size_t newTagDataLength = 
       prevTDLen + Constants::BAM_TAG_ARRAYBASE_SIZE + numElements*sizeof(T); // C++ string does not
       // need the last extra '\0', but this library manipulate data() as C null 
       // terminate string, so we have to add this extra. TODO: future versions
       // will eliminate the need to work in C space.
    RaiiBuffer originalTagData(newTagDataLength+1);
    //std::memcpy(originalTagData.Buffer, TagData.c_str(), prevTDLen+1); // '+1' for TagData's null-term
    std::memcpy(originalTagData.Buffer, TagData.c_str(), prevTDLen); // '+1' for TagData's null-term
    // write newTagBase (removes old null term), strcat need '\0' from originalTagData.Buffer
    //std::strcat(originalTagData.Buffer + prevTDLen, (const char*)newTagBase);
    std::memcpy(originalTagData.Buffer + prevTDLen, newTagBase, Constants::BAM_TAG_ARRAYBASE_SIZE);
    // add vector elements to tag
    int elementsBeginOffset = prevTDLen + Constants::BAM_TAG_ARRAYBASE_SIZE;
    for (int i = 0 ; i < numElements; ++i) {
        const T& value = values.at(i); // no need for at which slows down
        //cerr << "copying value: " << value << endl;
        memcpy(originalTagData.Buffer + elementsBeginOffset + i*sizeof(T), &value, sizeof(T));
    }
    // store temp buffer back in TagData
    //const char* newTagData = (const char*)originalTagData.Buffer;
    originalTagData.Buffer[newTagDataLength]='\0';
    //TagData.assign(newTagData, newTagDataLength); // not sure about the '\0' null terminator
    TagData.assign(originalTagData.Buffer, newTagDataLength); // not sure about the '\0' null terminator
    return true;
}

template<typename T>
inline void BamAlignment::addTag(const std::string& tag, const std::vector<T>& values) {
    // if char data not populated, do that first
    if (SupportData.HasCoreOnly) BuildCharData();
    // check for valid tag name length
    if (!isValidTagName(tag)) {
        cerr << __FILE__ << ":" << __LINE__ << ":WARN tag " + tag + " TAGSIZE wrong\n";
        throw logic_error("Invalid BamTag name: " + tag);
    }
    if (values.empty()) {
       throw logic_error("empty values when adding array tag");
    }
    // localize the tag data, work as char pointer type
    char* p = findTag(tag);
    if (p != nullptr) {
       throw logic_error("BamTag " + tag + " already exists");
    }
    const int32_t numElements  = values.size();
    size_t newtaglen = Constants::BAM_TAG_ARRAYBASE_SIZE + sizeof(T)*numElements;
    size_t prevTDL = TagData.size();
    TagData.resize(prevTDL + newtaglen);
    p = TagData.data() + prevTDL;
    memcpy(p, tag.c_str(), tag.size());
    p += Constants::BAM_TAG_TAGSIZE; // same as tag.size()
    *p = Constants::BAM_TAG_TYPE_ARRAY; // 'B'
    ++p;
    *p = TagTypeHelper<T>::TypeCode(); // int32_t => 'i', uint32_t => 'I'
    ++p;
    std::memcpy(p, &numElements, sizeof(int32_t)); // 32 bits is 4 bytes
    p += sizeof(int32_t); // now at first element data 
    memcpy(p, values.data(), sizeof(T)*numElements);
}

template<typename T>
inline bool BamAlignment::EditTag(const std::string& tag, const std::string& type, const T& value) 
{
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();
    // remove existing tag if present, then append tag with new value
    //if ( HasTag(tag) )
    if (hasTag(tag)) {
        //RemoveTag(tag);
        removeTag(tag);
    }
    return AddTag(tag, type, value);
}

template<typename T>
inline void BamAlignment::editTag(const std::string& tag, const std::string& type, const T& value) 
{
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();
    // remove existing tag if present, then append tag with new value
    if (hasTag(tag)) {
        removeTag(tag);
    }
    //return addTag(tag, type, value);
    addTag(tag, type, value);
}

// char type version
template<typename T>
inline void BamAlignment::editTag(const std::string& tag, const char type, const T& value) 
{
    //if (getName() == "S262055") {
    //    cerr << __LINE__ << ": before edit, TagData " << TagData << endl << " new value: "
    //        << value << endl;
    //}
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();
    // remove existing tag if present, then append tag with new value
    if (hasTag(tag)) {
        removeTag(tag);
    }
    addTag(tag, type, value); 
    //if (getName() == "S262055") {
    //    cerr << __LINE__ << ": after edit, TagData " << TagData << endl;
    //    auto x = getTag<string>("XA");
    //    cerr << x.first << endl;
    //}
}

// this versin has a bug need to fix
template<typename T>
inline bool BamAlignment::EditTag(const std::string& tag, const std::vector<T>& values) {
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();
    // remove existing tag if present, then append tag with new values
    if ( HasTag(tag) )
        //RemoveTag(tag);
        removeTag(tag);
    return AddTag(tag, values);
}

template<typename T>
inline void BamAlignment::editTag(const std::string& tag, const std::vector<T>& values) {
    // if char data not populated, do that first
    if ( SupportData.HasCoreOnly )
        BuildCharData();
    // remove existing tag if present, then append tag with new values
    if (hasTag(tag))
        removeTag(tag);
    //return AddTag(tag, values);
    addTag(tag, values);
}

/**
 * Obtain the tag value tag and put into destination
 * @param tag name of the tag such as NM
 * @param destination value of which will be updated from this 
 *    object.
 */
template<typename T>
inline bool BamAlignment::GetTag(const std::string& tag, T& destination) const {
    // skip if alignment is core-only
    if (SupportData.HasCoreOnly) {
        // TODO: set error string?
        return false;
    }
    // skip if no tags present
    if ( TagData.empty() ) {
        // TODO: set error string?
        return false;
    }
    // localize the tag data
    char* pTagData = (char*)TagData.data(); // string's low level representation
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    // return failure if tag not found
    if (!FindTag(tag, pTagData, tagDataLength, numBytesParsed)) {
        // TODO: set error string?
        return false;
    }
    const char type = *(pTagData - 1); // make sure type compatible
    if (!TagTypeHelper<T>::CanConvertFrom(type)) {
       //cerr << __FILE__ << ":" << __LINE__ << ":ERROR " << tag << " stored type " 
       //   << type << " cannot be converted to " << typeid(T).name() << endl;
       throw BamTypeException("Cannot convert from " + string(1, type) + " to " + string(typeid(T).name()));
       //return false;
    }
    size_t destinationLength = 0;
    try {
       destinationLength = getAtomicTagLength(*(pTagData-1));
    }
    catch (logic_error& err) {
       cerr << err.what() << endl;
       return false;
    }
    // store data in destination
    destination = 0;
    memcpy(&destination, pTagData, destinationLength);
    //cerr << "Just after memcpy destination value: " << destination << endl;
    return true;
}

// should not use this for vector version
// Will automatically convert from stored type to T if possible.
template<typename T> inline pair<T,bool> BamAlignment::getTag(const std::string& tag) const {
    // skip if alignment is core-only
    if (SupportData.HasCoreOnly || TagData.empty()) {
       throw BamNotagException("Only has core data");
       // return false;
    }
    const char* p = findTag(tag);
    if (p == nullptr) {
       //throw BamNotagException(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR Bam tag: " + tag + " not found");
       return make_pair(0, false);
    }
    // p point to first char of TAG
    //cerr << __LINE__ << ":DEBUG p is at " << p << endl;
    p += Constants::BAM_TAG_TAGSIZE;
    char typechar = *p;
    //cerr << "p is at " << p << endl;
    /*
    if (!TagTypeHelper<T>::CanConvertFrom(*p)) {
       //cerr << *this << endl;
       //cerr << __FILE__ << ":" << __LINE__ << ":ERROR " << tag << " stored type " 
       //   << *p << " cannot be converted to " << string(typeid(T).name()) << endl;
       throw BamTypeException(string(__FILE__) + ":" + to_string(__LINE__) 
             + ":ERROR Cannot convert from stored data type " + string(1, typechar) 
             + " in tag " + tag + " to requested type: " + string(typeid(T).name())
             + " consider using the type of the same signage");
    }
    */
    //size_t tagdatalen = getAtomicTagLength(*p);
    size_t tagdatalen = getAtomicTagLength(typechar); // for numeric type
    T res=0;
    if (sizeof(T) >= tagdatalen) { // return type is longer or same than stored type 
       //char typechar = *p;
       p += Constants::BAM_TAG_TYPESIZE;
       //cerr << "output data size is larger than stored value\n";
       if (typechar == Constants::BAM_TAG_TYPE_ASCII ||
             typechar == Constants::BAM_TAG_TYPE_INT8) 
       {
            int8_t x;
            memcpy(&x, p, sizeof(int8_t));
            res = x;
       }
       else if (typechar == Constants::BAM_TAG_TYPE_UINT8) {
            uint8_t x;
            memcpy(&x, p, sizeof(uint8_t));
            res = x;
       }
       else if (typechar == Constants::BAM_TAG_TYPE_INT16) {
            int16_t x;
            memcpy(&x, p, sizeof(int16_t));
            res = x;
       }
       else if (typechar == Constants::BAM_TAG_TYPE_UINT16) {
            uint16_t x;
            memcpy(&x, p, sizeof(uint16_t));
            res = x;
       }
       else if (typechar == Constants::BAM_TAG_TYPE_INT32) {
            int32_t x;
            memcpy(&x, p, sizeof(int32_t));
            res = x;
       }
       else if (typechar == Constants::BAM_TAG_TYPE_UINT32) {
            uint32_t x;
            memcpy(&x, p, sizeof(uint32_t));
            res = x;
       }
       else { //if (typechar == Constants::BAM_TAG_TYPE_FLOAT) {
          assert(typechar == Constants::BAM_TAG_TYPE_FLOAT);
            float x;
            memcpy(&x, p, sizeof(float));
            res = x;
       }
    }
    else {
       //cerr << "WARN: Tag type " << typechar << " longer than return type " << typeid(T).name() << endl;
       //if (!TagTypeHelper<T>::CanConvertFrom(typechar)) {
       //   //cerr << *this << endl;
       //   //cerr << __FILE__ << ":" << __LINE__ << ":ERROR " << tag << " stored type " 
       //   //   << *p << " cannot be converted to " << string(typeid(T).name()) << endl;
       //   throw BamTypeException(string(__FILE__) + ":" + to_string(__LINE__) 
       //         + ":ERROR Cannot convert from stored data type " + string(1, typechar) 
       //         + " in tag " + tag + " to requested type: " + string(typeid(T).name())
       //         + " consider using the type of the same signage");
       //}
      //cerr << "tagatalen=" << tagdatalen << endl;
      //p += Constants::BAM_TAG_TYPESIZE;
      //cerr << "p is at " << p << endl;
      //memcpy(&res, p, tagdatalen);

       p += Constants::BAM_TAG_TYPESIZE;
       //cerr << "output data size is larger than stored value\n";
       if (typechar == Constants::BAM_TAG_TYPE_ASCII || typechar == Constants::BAM_TAG_TYPE_INT8) 
       {
            int8_t x;
            memcpy(&x, p, sizeof(int8_t));
            if (x > static_cast<int8_t>(numeric_limits<T>::max())) {
               throw runtime_error("tag value " + to_string(x) + " cannot be store in user request type: " + string(typeid(T).name()));
            }
            res = x;
       }
       else if (typechar == Constants::BAM_TAG_TYPE_UINT8) {
            uint8_t x;
            memcpy(&x, p, sizeof(uint8_t));
            if (x > static_cast<uint8_t>(numeric_limits<T>::max())) {
               throw runtime_error("tag value " + to_string(x) + " cannot be store in user request type: " + string(typeid(T).name()));
            }
            res = x;
       }
       else if (typechar == Constants::BAM_TAG_TYPE_INT16) {
            int16_t x;
            memcpy(&x, p, sizeof(int16_t));
            if (x > static_cast<int16_t>(numeric_limits<T>::max())) {
               throw runtime_error("tag value " + to_string(x) + " cannot be store in user request type: " + string(typeid(T).name()));
            }
            res = x;
       }
       else if (typechar == Constants::BAM_TAG_TYPE_UINT16) {
            uint16_t x;
            memcpy(&x, p, sizeof(uint16_t));
            if (x > static_cast<uint16_t>(numeric_limits<T>::max())) {
               throw runtime_error("tag value " + to_string(x) + " cannot be store in user request type: " + string(typeid(T).name()));
            }
            res = x;
       }
       else if (typechar == Constants::BAM_TAG_TYPE_INT32) {
            int32_t x;
            memcpy(&x, p, sizeof(int32_t));
            if (x > static_cast<int32_t>(numeric_limits<T>::max())) {
               throw runtime_error("tag value " + to_string(x) + " cannot be store in user request type: " + string(typeid(T).name()));
            }
            res = x;
       }
       else if (typechar == Constants::BAM_TAG_TYPE_UINT32) {
            uint32_t x;
            memcpy(&x, p, sizeof(uint32_t));
            if (x > static_cast<uint32_t>(numeric_limits<T>::max())) {
               throw runtime_error("tag value " + to_string(x) + " cannot be store in user request type: " + string(typeid(T).name()));
            }
            res = x;
       }
       else { //if (typechar == Constants::BAM_TAG_TYPE_FLOAT) {
          assert(typechar == Constants::BAM_TAG_TYPE_FLOAT);
            float x;
            memcpy(&x, p, sizeof(float));
            if (x > static_cast<float>(numeric_limits<T>::max())) {
               throw runtime_error("tag value " + to_string(x) + " cannot be store in user request type: " + string(typeid(T).name()));
            }
            res = x;
       }
    }
    //cerr << "res=" << res << endl;
    return make_pair(res, true);
}

template<> inline pair<string,bool> BamAlignment::getTag<string>(const std::string& tag) const {
    //cerr << "getting tag " << tag << endl;
    // skip if alignment is core-only
    if (SupportData.HasCoreOnly || TagData.empty()) {
       throw BamNotagException("Only has core data");
       // return false;
    }
    const char* p = findTag(tag);
    //cerr << "now at tag: " << p << endl;
    if (p == nullptr) {
       //throw BamNotagException(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR Bam tag: " + tag + " not found");
       return make_pair(string(), false);
    }
    // p point to first char of TAG
    p += Constants::BAM_TAG_TAGSIZE;
    if (!TagTypeHelper<string>::CanConvertFrom(*p)) {
       //cerr << __FILE__ << ":" << __LINE__ << ":ERROR " << tag << " stored type " 
       //   << type << " cannot be converted to " << typeid(T).name() << endl;
       throw BamTypeException(string(__FILE__) + ":" + to_string(__LINE__) 
             + ":ERROR Cannot convert from stored data type " + string(1, *p) 
             + " in tag " + tag + " to string type");
    }
    // above code will make sure p point to H or Z
    //if (*p != Constants::BAM_TAG_TYPE_HEX && *p != Constants::BAM_TAG_TYPE_STRING) {
    //    throw BamTypeException("Tag: " + tag + " does not hold string data");
    //}
    p += Constants::BAM_TAG_TYPESIZE; // move to data
    return make_pair(string(p), true);
}

/**
 * Specialization for string
 */
template<>
inline bool BamAlignment::GetTag<std::string>(const std::string& tag, std::string& destination) const
{
    // skip if alignment is core-only
    if ( SupportData.HasCoreOnly ) {
        // TODO: set error string?
        return false;
    }
    // skip if no tags present
    if ( TagData.empty() ) {
        // TODO: set error string?
        return false;
    }
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    // return failure if tag not found
    if ( !FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        // TODO: set error string?
        return false;
    }

    // otherwise copy data into destination
    const unsigned int dataLength = strlen(pTagData); // len without null char
    //destination.clear();
    //destination.resize(dataLength);
    //memcpy( (char*)destination.data(), pTagData, dataLength );
    //assign is more efficient than above 3 statements
    destination.assign(pTagData, dataLength);
    return true;
}

/*! \fn template<typename T> bool GetTag(const std::string& tag, std::vector<T>& destination) const
    \brief Retrieves the numeric array associated with a BAM tag.

    \param tag[in]          2-character tag name
    \param destination[out] retrieved values
    \return \c true if found
*/
template<typename T>
inline bool BamAlignment::GetTag(const std::string& tag, std::vector<T>& destination) const {
    // skip if alignment is core-only
    if (SupportData.HasCoreOnly ) {
        cerr << __FILE__ << ":" << __LINE__ << ": bam has only core data failed to get tag: "
            << tag << "\n";
        // TODO: set error string?
        return false;
    }
    // skip if no tags present
    if (TagData.empty()) {
        cerr << __FILE__ << ":" << __LINE__ << ": empty tag data, failed to get tag value\n";
        return false;
    }
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    // return false if tag not found
    if (!FindTag(tag, pTagData, tagDataLength, numBytesParsed)) {
        //cerr << __FILE__ << ":" << __LINE__ << ": cannot find tag: " << tag 
        //   << ", failed to get tag value\n";
        return false;
    }
    // check that tag is array type
    const char tagType = *(pTagData - 1);
    if (tagType != Constants::BAM_TAG_TYPE_ARRAY) {
        cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ 
           << ":ERROR cannot store a non-array tag in array destination\n";
        throw logic_error("Tag: " + tag + " does not hold array data");
        //return false;
    }
    // fetch element type
    const char elementType = *pTagData;
    if (!TagTypeHelper<T>::CanConvertFrom(elementType) ) {
        cerr << __FILE__ << ":" << __LINE__ << ": cannot convert from "
            << elementType << endl;
        throw BamTypeException(string(typeid(T).name()) + " cannot hold Bam type: " + string(1, elementType));
        //return false;
    }
    if (!Constants::isAtomicBamTagType(elementType)) {
       throw logic_error(string(1, elementType) + " is not atomic bam type");
    }
    // now T and elementType are of the same data length is asserted
    /*
    ++pTagData;

    // calculate length of each element in tag's array
    //int elementLength = 0;
    switch (elementType ) {
        case (Constants::BAM_TAG_TYPE_ASCII) :
        case (Constants::BAM_TAG_TYPE_INT8)  :
        case (Constants::BAM_TAG_TYPE_UINT8) :
            //elementLength = sizeof(uint8_t);
            break;
        case (Constants::BAM_TAG_TYPE_INT16)  :
        case (Constants::BAM_TAG_TYPE_UINT16) :
            //elementLength = sizeof(uint16_t);
            break;
        case (Constants::BAM_TAG_TYPE_INT32)  :
        case (Constants::BAM_TAG_TYPE_UINT32) :
        case (Constants::BAM_TAG_TYPE_FLOAT)  :
            //elementLength = sizeof(uint32_t);
            break;
        // var-length types not supported for numeric destination
        case (Constants::BAM_TAG_TYPE_STRING) :
        case (Constants::BAM_TAG_TYPE_HEX)    :
        case (Constants::BAM_TAG_TYPE_ARRAY)  :
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ 
               << ":ERROR invalid array data, variable-length elements are not allowed\n";
            return false;
        // unknown tag type
        default:
            cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ 
               << ":ERROR invalid array element type: " << elementType << endl;
            return false;
    }
    */
    // not using elementLength
    // TODO: remove above code block
    //std::cerr << "BamTagData arrary element data width (Byte): " << elementLength << std::endl;
    // get number of elements
    ++pTagData;
    int32_t numElements;
    memcpy(&numElements, pTagData, sizeof(int32_t));
    pTagData += 4;
    destination.clear();
    destination.reserve(numElements);
    // read in elements
    T value;
    for (int i = 0 ; i < numElements; ++i ) {
        memcpy(&value, pTagData, sizeof(T));
        pTagData += sizeof(T);
        destination.push_back(value);
    }
    return true;
}

template<typename T>
inline bool BamAlignment::getTag(const std::string& tag, std::vector<T>& destination) const {
    // skip if alignment is core-only
    if (SupportData.HasCoreOnly) {
        cerr << __FILE__ << ":" << __LINE__ << ": bam has only core data failed to get tag: "
            << tag << "\n";
        // TODO: set error string?
        return false;
    }
    // skip if no tags present
    if (TagData.empty()) {
        cerr << __FILE__ << ":" << __LINE__ << ": empty tag data, failed to get tag value\n";
        return false;
    }
    const char* p = findTag(tag);
    if (!isValidArrayTag(p)) {
       throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR Invalid array tag: " + tag);
    }
    p += (Constants::BAM_TAG_TAGSIZE + Constants::BAM_TAG_TYPESIZE);
    const char elementType = *p;
    if (!TagTypeHelper<T>::CanConvertFrom(elementType) ) {
        cerr << __FILE__ << ":" << __LINE__ << ": cannot convert from "
            << elementType << endl;
        throw logic_error(string(typeid(T).name()) + " cannot hold Bam type: " + string(1, elementType));
    }
    int32_t numE;
    ++p;
    memcpy(&numE, p, sizeof(int32_t));
    p += sizeof(int32_t); // should be 4 bytes
    destination.resize(numE);
    //size_t len = numE * sizeof(T);
    memcpy(destination.data(), p, numE*sizeof(T));
    return true;
}

template<typename T>
inline vector<T> BamAlignment::getArrayTag(const std::string& tag) const {
    // skip if alignment is core-only
    if (SupportData.HasCoreOnly || TagData.empty()) {
        cerr << __FILE__ << ":" << __LINE__ << ": bam has only core data failed to get tag: "
            << tag << "\n";
        //return false;
        throw BamNotagException("Has only core data");
    }
    const char* p = findTag(tag);
    if (p == nullptr) {
       //throw BamNotagException("There is no BamTag: " + tag);
       return vector<T>(); // return empty vector object
    }
    if (!isValidArrayTag(p)) { // p at start of TAG
       throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR Invalid array tag: " + tag);
    }
    p += (Constants::BAM_TAG_TAGSIZE + Constants::BAM_TAG_TYPESIZE);
    const char elementType = *p;
    if (!TagTypeHelper<T>::CanConvertFrom(elementType) ) {
        cerr << __FILE__ << ":" << __LINE__ << ": cannot convert from "
            << elementType << endl;
        throw logic_error(string(typeid(T).name()) + " cannot hold Bam type: " + string(1, elementType));
    }
    int32_t numE;
    ++p;
    memcpy(&numE, p, sizeof(int32_t));
    p += sizeof(int32_t); // should be 4 bytes
    vector<T> res(numE);
    //destination.resize(numE);
    //size_t len = numE * sizeof(T);
    memcpy(res.data(), p, numE*sizeof(T));
    return res;
}

typedef std::vector<BamAlignment> BamAlignmentVector;

} // namespace BamTools

#endif // BAMALIGNMENT_H
