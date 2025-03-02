// ***************************************************************************
// bamtools_utilities.h (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 October 2011
// ---------------------------------------------------------------------------
// Provides general utilities used by BamTools sub-tools.
// ***************************************************************************

#ifndef BAMTOOLS_UTILITIES_H
#define BAMTOOLS_UTILITIES_H

#include <api/BamAux.h>
#include <utils/utils_global.h>
#include <string>
#include <vector>
#include <regex>

#define BAMTOOLS_ASSERT_UNREACHABLE BT_ASSERT_UNREACHABLE
#define BAMTOOLS_ASSERT_MESSAGE( condition, message ) BT_ASSERT_X( condition, message )

namespace BamTools {

class BamReader;
class BamMultiReader;

class UTILS_EXPORT Utilities {
  
    public: 
        // returns true if 'source' contains 'pattern' or 'c'
        static bool Contains(const std::string& source, const std::string& pattern);
        static bool Contains(const std::string& source, const char c);

        // returns true if 'source' ends with 'pattern' or 'c'
        static bool EndsWith(const std::string& source, const std::string& pattern);
        static bool EndsWith(const std::string& source, const char c);

        // check if a file exists
        static bool FileExists(const std::string& fname);
        
        /** 
         * Parses a region string, uses reader to do validation (valid ID's, positions), stores in Region struct
         * Returns success (true/false)
         * @param region is the result to be outputed by this function.
         */
        static bool ParseRegionString(const std::string& regionString,
                                      const BamReader& reader,
                                      BamRegion& region);
        // Same as above, but accepts a BamMultiReader
        static bool ParseRegionString(const std::string& regionString,
                                      const BamMultiReader& reader,
                                      BamRegion& region);

        // sequence utilities
        static void Reverse(std::string& sequence);
        static void ReverseComplement(std::string& sequence);

        // split string on delimiter character (or string of allowed delimiters)
        static std::vector<std::string> Split(const std::string& source, const char delim);
        static std::vector<std::string> Split(const std::string& source, const std::string& delims);

        // returns true if 'source' starts with 'pattern' or 'c'
        static bool StartsWith(const std::string& source, const std::string& pattern);
        static bool StartsWith(const std::string &source, const char c);
        /**
         * @param reader is a helper object to check for the validity of the refname and pos.
         * @return true if both refname is valid and pos is not outside the end of the genomic DNA
         */
        static bool validRefpos(const string& refname, int pos, const BamReader& reader);
        //bool validRefpos(int rid, int pos, const BamReader& reader);
        template<class T>
        static bool validRefpos(int rid, int pos, const T& reader) {
             if ( rid == -1 ) {
                cerr << __LINE__ << ": refid must be positive" << endl;
                return false;
             }
             const RefVector& refv = reader.GetReferenceData();
             // startPos cannot be greater than or equal to reference length
             const RefData& refinfo = refv.at(rid);
             if (pos >= refinfo.getLength()) return false;
             return true;
         }
        template<class T>
        static void setToEnd(int rid, int& pos, const T& reader) {
             const RefVector& refv = reader.GetReferenceData();
             const RefData& refinfo = refv.at(rid);
             pos = refinfo.getLength();
         }
        static tuple<string, int, string, int> extractRegion(const string& regstr);

        /**
         * @param br is only used for reading only.
         */
        template<class T>
        static array<int, 4> parseRegion(const string& regstr, const T& br) {
           cerr << __LINE__ << ": parsing region string: " << regstr << endl;
            regex singleGenomic("([_A-Za-z0-9]+):([[:digit:]]+)(?:\\.\\.|-)([[:digit:]]+)");
            regex doubleGenomic("([_A-Za-z0-9]+):([[:digit:]]+)(?:\\.\\.|-)([_A-Za-z0-9]+):([[:digit:]]+)");
            smatch matchResult;
            array<int,4> res{-1, 0, -1, -1};
            if (regex_match(regstr, matchResult, singleGenomic)) {
               cerr << __LINE__ << ": single genomic specified\n";
              res[0] = br.getReferenceId(matchResult[1]);
              if (res[0] == -1) {
                 throw runtime_error("invalide first genomic id: " + matchResult[1].str());
              }
              res[2] = res[0];
              res[1] = stoi(matchResult[2]);
              if (!validRefpos<T>(res[0], res[1], br)) {
                 throw runtime_error("invalide start position " + matchResult[2].str() + " on " + matchResult[1].str());
              }
              res[3] = stoi(matchResult[3]);
              if (res[3] == -1) { 
                 setToEnd<T>(res[0], res[3], br);
              }
              if (!validRefpos<T>(res[0], res[3], br)) {
                 cerr << __LINE__ << ": invalid end position " << res[3] << " on " << matchResult[1] << " we will assume the end of it\n";
                 //res[3] = -1;
                 // -1 means to the end of the bam file. So must be set to the end of the chromosome
                 setToEnd<T>(res[0], res[3], br);
              }
            }
            else if (regex_match(regstr, matchResult, doubleGenomic)) {
               cerr << __LINE__ << ": double genomic specified\n";
              res[0] = br.getReferenceId(matchResult[1]);
              if (res[0] == -1) {
                 throw runtime_error("invalide first genomic id: " + matchResult[1].str());
              }
              res[1] = stoi(matchResult[2]);
              if (!validRefpos<T>(res[0], res[1], br)) {
                 throw runtime_error("invalide start position " + matchResult[2].str() + " on " + matchResult[1].str());
              }
              res[2] = br.getReferenceId(matchResult[3]);
              if (res[2] == -1) {
                 throw runtime_error("invalide sedond genomic id: " + matchResult[3].str());
              }
              res[3] = stoi(matchResult[4]);
              if (!validRefpos<T>(res[2], res[3], br)) {
                 cerr << __LINE__ << ": invalid end position " << res[3] << " on " << matchResult[3] << " we will assume the end of it\n";
                 //res[3] = -1;
                 setToEnd<T>(res[0], res[3], br);
              }
            }
            else {
               auto i = regstr.find('-');
               if (i != string::npos) {
                  throw runtime_error("invalid region: " + regstr);
               }
               i = regstr.find("..");
               if (i != string::npos) {
                  throw runtime_error("invalid region: " + regstr);
               }
               i = regstr.find(':');
               if (i != string::npos) { // chrX:123
                    res[0] = br.getReferenceId(regstr.substr(0,i));
                    if (res[0] == -1) {
                       throw runtime_error("invalide first genomic name: " + regstr.substr(0,i));
                    }
                    res[1] = stoi(regstr.substr(i+1));
                    if (!validRefpos<T>(res[0], res[1], br)) {
                       throw runtime_error("invalid reference position " + regstr.substr(i+1));
                    }
                    res[2] = res[0];
                    //res[3] = -1;
                    setToEnd(res[0], res[3], br);
               }
               else { // chrom3 no colon
                    res[0] = br.getReferenceId(regstr);
                    if (res[0] == -1) {
                       throw runtime_error("invalide first genomic name: " + regstr);
                    }
                    res[2] = res[0];
                    res[1]=0;
                    setToEnd<T>(res[0], res[3], br);
               }
            }
            cerr << __LINE__ << ": parseRegion result: " << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << endl;
            return res;
         }
};

} // namespace BamTools
  
#endif // BAMTOOLS_UTILITIES_H
