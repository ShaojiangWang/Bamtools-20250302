// ***************************************************************************
// BamAux.h (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 25 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides data structures & utility methods that are used throughout the API.
// ***************************************************************************
// This is a header-only file, no implementation

#ifndef BAMAUX_H
#define BAMAUX_H

#include "api/api_global.h"
#include <cstring>
#include <fstream> 
#include <iostream>
#include <string>
#include <vector>
#include <array>

using namespace std;

/*! \file BamAux.h

    Provides data structures & utility methods that are used throughout the API.
*/

/*! \namespace BamTools
    \brief Contains all BamTools classes & methods.

    The BamTools API contained in this namespace contains classes and methods
    for reading, writing, and manipulating BAM alignment files.
*/
namespace BamTools {

// ----------------------------------------------------------------
// CigarOp

/*! \struct BamTools::CigarOp
    \brief Represents a CIGAR alignment operation.

    \sa \samSpecURL for more details on using CIGAR operations.
*/
struct API_EXPORT CigarOp {
    char     Type;   //!< CIGAR operation type (MIDNSHPX=)
    uint32_t Length; //!< CIGAR operation length (number of bases)
    
    /** 
     * default constructor 
     */
    CigarOp() : Type('\0'), Length(0) { }
    /**
     * Constructor from known value
     */
    CigarOp(const char type, const uint32_t& length)
        : Type(type) , Length(length) { }
    bool operator==(const CigarOp& co) const {
       return Type == co.Type && Length == co.Length;
    }
    friend ostream& operator<<(ostream& ous, const CigarOp& co) {
       //ous << co.Type << co.Length; return ous; 
       ous << co.Length << co.Type; return ous; 
    }
    /**
     * Convert this to a more universal type
     * If this class does not provide any special operation,
     * then there is no use to use a special type. Just 
     * more noise.
     */
    pair<char, int> topair() const { return pair<char,int>(Type, Length); }
    void fromPair(const pair<char,int> &p) { Type=p.first; Length = p.second; }
    /**
     * @return the length as 32-bits integer
     */
    uint32_t getLength() const { return Length; }
    /**
     * @return the cigar operation one of MIDNSHPX=
     */
    char getType() const { return Type; }
    bool isDeletion() const {
       return Type == 'D';
    }
    bool isInsertion() const {
       return Type == 'I';
    }
    bool isMatch() const {
       return Type == 'M';
    }
    bool isSoft() const {
       return Type == 'S';
    }
    /**
     * change the length for the cigar segment to l.
     * @param l new length.
     */
    void setLength(uint32_t l) { Length = l; }
    /**
     * Reduce the length by l
     */
    void shrink(uint32_t l) { 
       if (Length <= l) 
          throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR Length "
                + to_string(Length) + " shrink by " + to_string(l) + " will be zeor or negative");
       Length -= l; 
    }
    /**
     * Expand length by l
     */
    void expand(uint32_t l) { Length += l; }
    void setType(char t) { Type = t; }
};

// ----------------------------------------------------------------
// RefData

/**
 *  Represents a brief reference sequence entry: [name, length].
*/
struct API_EXPORT RefData {
    std::string RefName;    //!< name of reference sequence
    int32_t     RefLength;  //!< length of reference sequence
    
    /** 
     * default constructor
     */
    RefData()
        : RefName(), RefLength(0)
    { }
    RefData(const std::string& name,
            const int32_t& length)
        : RefName(name), RefLength(length)
    { }
    pair<string, int32_t> asPair() const {
       return make_pair(RefName, RefLength);
    }
    /**
     * @return the reference name
     */
    const string& getRefname() const {
       return RefName;
    }
    const string& getName() const {
       return RefName;
    }
    int32_t getReflength() const {
       return RefLength;
    }
    int32_t getLength() const {
       return RefLength;
    }
    bool operator==(const string& refn) const {
       return RefName == refn;
    }
    friend bool operator==(const string& refn, const RefData& rd) {
       return refn == rd.RefName;
    }
    friend ostream& operator<<(ostream& ous, const RefData& rd) {
       ous << rd.RefName << '\t' << rd.RefLength;
       return ous;
    }
};

/**
 * Convenience typedef for vector of RefData entries.
 * Alias for vector of RefData.
 */
typedef std::vector<RefData> RefVector;

// ----------------------------------------------------------------
// BamRegion

/** 
 *  Represents a sequential genomic region
 *
 *  Allowed to span multiple (sequential) references.
 *
 *  Warning: BamRegion now represents a zero-based, HALF-OPEN interval.
 *  In previous versions of BamTools (0.x & 1.x) all intervals were treated
 *  as zero-based, CLOSED.
 *
 *  Region can span multiple references if the references are
 *  sorted in some order.
*/
struct API_EXPORT BamRegion {
    int LeftRefID;      //!< reference ID for region's left boundary
    /**
     * position for region's left boundary. 0-based index.
     */
    int LeftPosition;   
    int RightRefID;     //!< reference ID for region's right boundary
    /**
     * [leftPos, rightPos)
     * Right position is the desired end + 1
     * -1 means the end of the entire chromosome.
     */
    int RightPosition; 
    
    /**
     * constructor from full information.
     * @param leftId left reference id, zero-based index
     * @param leftPos left position zero-indexed position. -1 means to the end.
     * @param rightId -1 means all the way to the end of the file
     *     if rightId == leftId, then means the entire leftId.
     */
    BamRegion(const int& leftID   = -1, const int& leftPos  = -1,
              const int& rightID  = -1, const int& rightPos = -1)
        : LeftRefID(leftID), LeftPosition(leftPos)
        , RightRefID(rightID), RightPosition(rightPos)
    { }
    /**
     * Constructor for only one reference
     * @param refid reference id, chr1 is 1, chrm is 0, chrX is 23
     *    Should use the refid from the header of the Bam file.
     * @param leftp left position
     * @param rightp right position
     */
    //explicit BamRegion(int refid, int leftp, int rightp) 
    explicit BamRegion(int refid, const pair<int,int>& reg) 
       : LeftRefID(refid), LeftPosition(reg.first), RightRefID(refid), RightPosition(reg.second)
    { }
    BamRegion(const array<int,3>& refid_b_e) 
       : LeftRefID(refid_b_e[0]), LeftPosition(refid_b_e[1]), RightRefID(refid_b_e[0]), RightPosition(refid_b_e[2])
    { }
    /**
     * Region on a single chromosome.
     */
    BamRegion(int refid, int b, int e) 
       : LeftRefID(refid), LeftPosition(b), RightRefID(refid), RightPosition(e)
    { }
    
    /** 
      * copy constructor
    */
    BamRegion(const BamRegion& other)
        : LeftRefID(other.LeftRefID), LeftPosition(other.LeftPosition)
        , RightRefID(other.RightRefID), RightPosition(other.RightPosition)
    { }

    BamRegion& operator=(const BamRegion& o) {
      if (this != &o) {
         LeftRefID=o.LeftRefID;
         LeftPosition=o.LeftPosition;
         RightRefID=o.RightRefID;
         RightPosition=o.RightPosition;
      }
      return *this;
    }
    /**
     * Updating interval without changing the reference
     */
    void setInterval(const pair<int,int>& itv) {
       LeftPosition = itv.first;
       RightPosition = itv.second;
    }
    void set(const array<int,4>& rawreg) {
       LeftRefID = rawreg[0];
       LeftPosition = rawreg[1];
       RightRefID = rawreg[2];
       RightPosition = rawreg[3];
    }
    //! Clears region boundaries
    void clear(void) {
        LeftRefID  = -1; LeftPosition  = -1;
        RightRefID = -1; RightPosition = -1;
    }

    //! Returns true if region has a left boundary
    bool isLeftBoundSpecified(void) const {
        return ( LeftRefID >= 0 && LeftPosition >= 0 );
    }

    //! Returns true if region boundaries are not defined
    bool isNull(void) const {
        return ( !isLeftBoundSpecified() && !isRightBoundSpecified() );
    }

    //! Returns true if region has a right boundary
    bool isRightBoundSpecified(void) const {
        return ( RightRefID >= 0 && RightPosition >= 1 );
    }
    bool isSingleReference() const {
       return LeftRefID == RightRefID && LeftPosition != -1;
    }
    friend ostream& operator<<(ostream& ous, const BamRegion& reg) {
       ous << reg.LeftRefID << ":" << reg.LeftPosition << "-"   
         << reg.RightRefID << ":" << reg.RightPosition; 
       return ous;
    }
};

struct CustomHeaderTag {
  std::string TagName;
  std::string TagValue;
};


// ----------------------------------------------------------------
// General utility methods

/*! \fn bool FileExists(const std::string& filename)
    \brief returns true if the file exists
*/
API_EXPORT inline bool FileExists(const std::string& filename) {
    std::ifstream f(filename.c_str(), std::ifstream::in);
    return !f.fail();
}

/*! \fn void SwapEndian_16(int16_t& x)
    \brief swaps endianness of signed 16-bit integer, in place
*/
API_EXPORT inline void SwapEndian_16(int16_t& x) {
    x = ((x >> 8) | (x << 8));
}

/*! \fn void SwapEndian_16(uint16_t& x)
    \brief swaps endianness of unsigned 16-bit integer, in place
*/
API_EXPORT inline void SwapEndian_16(uint16_t& x) {
    x = ((x >> 8) | (x << 8));
}

/*! \fn void SwapEndian_32(int32_t& x)
    \brief swaps endianness of signed 32-bit integer, in place
*/
API_EXPORT inline void SwapEndian_32(int32_t& x) {
    x = ( (x >> 24) | 
         ((x << 8) & 0x00FF0000) | 
         ((x >> 8) & 0x0000FF00) | 
          (x << 24)
        );
}

/*! \fn void SwapEndian_32(uint32_t& x)
    \brief swaps endianness of unsigned 32-bit integer, in place
*/
API_EXPORT inline void SwapEndian_32(uint32_t& x) {
    x = ( (x >> 24) | 
         ((x << 8) & 0x00FF0000) | 
         ((x >> 8) & 0x0000FF00) | 
          (x << 24)
        );
}

/*! \fn void SwapEndian_64(int64_t& x)
    \brief swaps endianness of signed 64-bit integer, in place
*/
API_EXPORT inline void SwapEndian_64(int64_t& x) {
    x = ( (x >> 56) | 
         ((x << 40) & 0x00FF000000000000ll) |
         ((x << 24) & 0x0000FF0000000000ll) |
         ((x << 8)  & 0x000000FF00000000ll) |
         ((x >> 8)  & 0x00000000FF000000ll) |
         ((x >> 24) & 0x0000000000FF0000ll) |
         ((x >> 40) & 0x000000000000FF00ll) |
          (x << 56)
        );
}

/*! \fn void SwapEndian_64(uint64_t& x)
    \brief swaps endianness of unsigned 64-bit integer, in place
*/
API_EXPORT inline void SwapEndian_64(uint64_t& x) {
    x = ( (x >> 56) | 
         ((x << 40) & 0x00FF000000000000ll) |
         ((x << 24) & 0x0000FF0000000000ll) |
         ((x << 8)  & 0x000000FF00000000ll) |
         ((x >> 8)  & 0x00000000FF000000ll) |
         ((x >> 24) & 0x0000000000FF0000ll) |
         ((x >> 40) & 0x000000000000FF00ll) |
          (x << 56)
        );
}

/*! \fn void SwapEndian_16p(char* data)
    \brief swaps endianness of the next 2 bytes in a buffer, in place
*/
API_EXPORT inline void SwapEndian_16p(char* data) {
    uint16_t& value = (uint16_t&)*data; 
    SwapEndian_16(value);
}

/*! \fn void SwapEndian_32p(char* data)
    \brief swaps endianness of the next 4 bytes in a buffer, in place
*/
API_EXPORT inline void SwapEndian_32p(char* data) {
    uint32_t& value = (uint32_t&)*data; 
    SwapEndian_32(value);
}

/*! \fn void SwapEndian_64p(char* data)
    \brief swaps endianness of the next 8 bytes in a buffer, in place
*/
API_EXPORT inline void SwapEndian_64p(char* data) {
    uint64_t& value = (uint64_t&)*data; 
    SwapEndian_64(value);
}

/*! \fn bool SystemIsBigEndian(void)
    \brief checks host architecture's byte order
    \return \c true if system uses big-endian ordering
*/
API_EXPORT inline bool SystemIsBigEndian(void) {
   const uint16_t one = 0x0001;
   return ((*(char*) &one) == 0 );
}

/*! \fn void PackUnsignedInt(char* buffer, unsigned int value)
    \brief stores unsigned integer value in a byte buffer

    \param[out] buffer destination buffer
    \param[in]  value  value to 'pack' in buffer
*/
API_EXPORT inline void PackUnsignedInt(char* buffer, unsigned int value) {
    buffer[0] = (char)value;
    buffer[1] = (char)(value >> 8);
    buffer[2] = (char)(value >> 16);
    buffer[3] = (char)(value >> 24);
}

/*! \fn void PackUnsignedShort(char* buffer, unsigned short value)
    \brief stores unsigned short integer value in a byte buffer

    \param[out] buffer destination buffer
    \param[in]  value  value to 'pack' in buffer
*/
API_EXPORT inline void PackUnsignedShort(char* buffer, unsigned short value) {
    buffer[0] = (char)value;
    buffer[1] = (char)(value >> 8);
}

/*! \fn double UnpackDouble(const char* buffer)
    \brief reads a double value from byte buffer

    \param[in] buffer source byte buffer
    \return the (double) value read from the buffer
*/
API_EXPORT inline double UnpackDouble(const char* buffer) {
    union { double value; unsigned char valueBuffer[sizeof(double)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    un.valueBuffer[4] = buffer[4];
    un.valueBuffer[5] = buffer[5];
    un.valueBuffer[6] = buffer[6];
    un.valueBuffer[7] = buffer[7];
    return un.value;
}

/*! \fn double UnpackDouble(char* buffer)
    \brief reads a double value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (double) value read from the buffer
*/
API_EXPORT inline double UnpackDouble(char* buffer) {
    return UnpackDouble( (const char*)buffer );
}

/*! \fn double UnpackFloat(const char* buffer)
    \brief reads a float value from byte buffer

    \param[in] buffer source byte buffer
    \return the (float) value read from the buffer
*/
API_EXPORT inline float UnpackFloat(const char* buffer) {
    union { float value; unsigned char valueBuffer[sizeof(float)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

/*! \fn double UnpackFloat(char* buffer)
    \brief reads a float value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (float) value read from the buffer
*/
API_EXPORT inline float UnpackFloat(char* buffer) {
    return UnpackFloat( (const char*)buffer );
}

/*! \fn signed int UnpackSignedInt(const char* buffer)
    \brief reads a signed integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (signed int) value read from the buffer
*/
API_EXPORT inline signed int UnpackSignedInt(const char* buffer) {
    union { signed int value; unsigned char valueBuffer[sizeof(signed int)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

/*! \fn signed int UnpackSignedInt(char* buffer)
    \brief reads a signed integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (signed int) value read from the buffer
*/
API_EXPORT inline signed int UnpackSignedInt(char* buffer) {
    return UnpackSignedInt( (const char*) buffer );
}

/*! \fn signed short UnpackSignedShort(const char* buffer)
    \brief reads a signed short integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (signed short) value read from the buffer
*/
API_EXPORT inline signed short UnpackSignedShort(const char* buffer) {
    union { signed short value; unsigned char valueBuffer[sizeof(signed short)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    return un.value;
}

/*! \fn signed short UnpackSignedShort(char* buffer)
    \brief reads a signed short integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (signed short) value read from the buffer
*/
API_EXPORT inline signed short UnpackSignedShort(char* buffer) {
    return UnpackSignedShort( (const char*)buffer );
}

/*! \fn unsigned int UnpackUnsignedInt(const char* buffer)
    \brief reads an unsigned integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (unsigned int) value read from the buffer
*/
API_EXPORT inline unsigned int UnpackUnsignedInt(const char* buffer) {
    union { unsigned int value; unsigned char valueBuffer[sizeof(unsigned int)]; } un;
    un.value = 0;
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
    un.valueBuffer[2] = buffer[2];
    un.valueBuffer[3] = buffer[3];
    return un.value;
}

/*! \fn unsigned int UnpackUnsignedInt(char* buffer)
    \brief reads an unsigned integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (unsigned int) value read from the buffer
*/
API_EXPORT inline unsigned int UnpackUnsignedInt(char* buffer) {
    return UnpackUnsignedInt( (const char*)buffer );
}

/*! \fn unsigned short UnpackUnsignedShort(const char* buffer)
    \brief reads an unsigned short integer value from byte buffer

    \param[in] buffer source byte buffer
    \return the (unsigned short) value read from the buffer
*/
API_EXPORT inline unsigned short UnpackUnsignedShort(const char* buffer) {
    union { unsigned short value; unsigned char valueBuffer[sizeof(unsigned short)]; } un;
    un.value = 0;
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    un.valueBuffer[0] = buffer[0];
    un.valueBuffer[1] = buffer[1];
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    un.valueBuffer[0] = buffer[1];
    un.valueBuffer[1] = buffer[0];
#else
    #error "Unsupported hardware"
#endif
    return un.value;
}

/*! \fn unsigned short UnpackUnsignedShort(char* buffer)
    \brief reads an unsigned short integer value from byte buffer

    This is an overloaded function.

    \param[in] buffer source byte buffer
    \return the (unsigned short) value read from the buffer
*/
API_EXPORT inline unsigned short UnpackUnsignedShort(char* buffer) {
    return UnpackUnsignedShort( (const char*)buffer );
}

// ----------------------------------------------------------------
// 'internal' helper structs

/**
 * struct RaiiBuffer internal.
 *  Work with n allocated, not n+1.
 *  Insuffient method provided worse than managing your own memory.
 */
struct RaiiBuffer {
    // data members
   /**
    * Internal buffer
    */
    char* Buffer;
    /**
     * Number of bytes allocated
     */
    const size_t NumBytes;

    /** 
     * constructor of size n buffer
     */
    RaiiBuffer(const size_t n)
        : Buffer( new char[n]() ), NumBytes(n)
    { }
    ~RaiiBuffer(void) {
        delete[] Buffer;
    }
    /** 
     * Set all location to 0 char
     */
    void Clear(void) {
        memset(Buffer, 0, NumBytes);
    }
};

/**
 * ==M0==~~X0~~==M1===~~~X1~~===M2===
 * Match segment is always one more than mismatch segments.
 * If no match at start or end then the value is zero.
 * The match segment stores the number of exact matched 
 * length of bases between reference and query. the mismatched
 * segment stores the reference sequences either as mismatch
 * taht is represented by base character or insertion started
 * with the ^ character.
 */
class Matchdiff {
   public:
      Matchdiff() { }
      /**
       * BamAlignment::getMDArray() can parse the MD tag
       * and convert it to the two vectors used here.
       */
      Matchdiff(vector<int>&& m, vector<string>&& x)
        : mseg(m), xseg(x)
      { }
      /**
       * Can be used to update the MD tag value
       */
      string toString() const;
      /**
       * @param idx is 0 based index from 0 to the length of the
       * aligment minux one. [0, L-1] where L is the length
       * of the alignment. idx must > 0.
       * @return the number of mismatches in removed region.
       *   mismatch are query deletion and mismatch (excluding query insertion)
       */
      int removeBefore(int idx);
      int removeAfter(int idx);
      int getMseglen(int i) const {
         return mseg[i];
      }
      int getXseglen(int i) const {
         if (xseg[i][0] == '^') return xseg[i].size()-1;
         else return xseg[i].size();
      }
      /**
       * @param x is a valid index ins xseg
       * @return true If query is deleted for the xth xsegment.
       */
      bool isDeletion(int x) const {
         return xseg[x][0] == '^';
      }
      /**
       * Total legth of the mseg and xseg
       */
      int length() const;

   private:
      vector<int> mseg;
      vector<string> xseg; // mismatch of refseq
};

} // namespace BamTools

#endif // BAMAUX_H
