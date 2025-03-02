// ***************************************************************************
// bamtools_utilities.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 October 2011
// ---------------------------------------------------------------------------
// Provides general utilities used by BamTools sub-tools.
// ***************************************************************************

#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <utils/bamtools_utilities.h>
using namespace BamTools;

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
//#include <regex>

using namespace std;

namespace BamTools {
  
const char REVCOMP_LOOKUP[] = {'T',  0,  'G', 'H',
                                0,   0,  'C', 'D',
                                0,   0,   0,   0,
                               'K', 'N',  0,   0,
                                0,  'Y', 'W', 'A',
                               'A', 'B', 'S', 'X',
                               'R',  0 };
  
} // namespace BamTools 
  
// returns true if 'source' contains 'pattern'
bool Utilities::Contains(const string& source, const string& pattern) {
    return ( source.find(pattern) != string::npos );
}

// returns true if 'source' contains 'c'
bool Utilities::Contains(const std::string &source, const char c) {
    return ( source.find(c) != string::npos );
}

// returns true if 'source' ends with 'pattern'
bool Utilities::EndsWith(const string& source, const string& pattern) {
    return ( source.find(pattern) == (source.length() - pattern.length()) );
}

// returns true if 'source' ends with 'c'
bool Utilities::EndsWith(const std::string& source, const char c) {
    return ( source.find(c) == (source.length() - 1) );
}

// check if a file exists
bool Utilities::FileExists(const string& filename) {
    ifstream f(filename.c_str(), ifstream::in);
    return !f.fail();
}

bool Utilities::validRefpos(const string& refname, int pos, const BamReader& reader) {
    int rid = reader.getReferenceID(refname);
    if ( rid == -1 ) {
       cerr << __LINE__ << ": genomic name=" << refname << " not in bam file\n";
       return false;
    }
    if (pos == -1) return true;
    const RefVector refv = reader.getReferenceData();
    // startPos cannot be greater than or equal to reference length
    const RefData& refinfo = refv.at(rid);
    if (pos >= refinfo.getLength()) return false;
    return true;
}

/*
bool Utilities::validRefpos(int rid, int pos, const BamReader& reader) {
    if ( rid == -1 ) {
       cerr << __LINE__ << ": refid must be positive" << endl;
       return false;
    }
    const RefVector& refv = reader.getReferenceData();
    // startPos cannot be greater than or equal to reference length
    const RefData& refinfo = refv.at(rid);
    if (pos >= refinfo.getLength()) return false;
    return true;
}
*/

tuple<string, int, string, int> Utilities::extractRegion(const string& regstr) {
   regex singleGenomic("([_A-Za-z0-9]+):(\\d+)(?:\\.\\.|-)(\\d+)");
   regex doubleGenomic("([_A-Za-z0-9]+):(\\d+)(?:\\.\\.|-)([_A-Za-z0-9]+):(\\d+)");
   smatch matchResult;
   tuple<string,int,string,int> res{"", -1, "", -1};
   if (regex_match(regstr, matchResult, singleGenomic)) {
     get<0>(res) = matchResult[1];
     get<1>(res) = stoi(matchResult[2]);
     get<2>(res) = matchResult[1];
     get<3>(res) = stoi(matchResult[3]);
   }
   else if (regex_match(regstr, matchResult, doubleGenomic)) {
     get<0>(res) = matchResult[1];
     get<1>(res) = stoi(matchResult[2]);
     get<2>(res) = matchResult[3];
     get<3>(res) = stoi(matchResult[4]);
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
      if (i != string::npos) {
           get<0>(res) = regstr.substr(0, i);
           get<1>(res) = stoi(regstr.substr(i+1));
           get<2>(res) = get<0>(res);
      }
      else {
           get<0>(res) = regstr;
      }
   }
   return res;
}

/* converted to template; now in the header file
array<int, 4> Utilities::parseRegion(const string& regstr, const BamReader& br) {
   regex singleGenomic("([_A-Za-z0-9]+):(\d+)(?:\.\.|-)(\d+)");
   regex doubleGenomic("([_A-Za-z0-9]+):(\d+)(?:\.\.|-)([_A-Za-z0-9]+):(\d+)");
   smatch matchResult;
   array<int,4> res{-1, -1, -1, -1};
   if (regex_match(regstr, matchResult, singleGenomic)) {
     res[0] = br.getReferenceId(matchResult[1]);
     if (res[0] == -1) {
        throw runtime_error("invalide first genomic id: " + matchResult[1]);
     }
     res[2] == res[0];
     res[1] = stoi(matchResult[2]);
     if (!validRefpos(res[0], res[1])) {
        throw runtime_error("invalide start position " + matchResult[2] + " on " + matchResult[1]);
     }
     res[3] = stoi(matchResult[3]);
     if (!validRefpos(res[0], res[3])) {
        cerr << __LINE__ << ": invalid end position " << res[3] << " on " << matchResult[1] << " we will assume the end of it\n";
        res[3] = -1;
     }
   }
   else if (regex_match(regstr, matchResult, doubleGenomic)) {
     res[0] = br.getReferenceId(matchResult[1]);
     if (res[0] == -1) {
        throw runtime_error("invalide first genomic id: " + matchResult[1]);
     }
     res[1] = stoi(matchResult[2]);
     if (!validRefpos(res[0], res[1])) {
        throw runtime_error("invalide start position " + matchResult[2] + " on " + matchResult[1]);
     }
     res[2] == br.getReferenceId(matchResult[3]);
     if (res[2] == -1) {
        throw runtime_error("invalide sedond genomic id: " + matchResult[3]);
     }
     res[3] = stoi(matchResult[4]);
     if (!validRefpos(res[2], res[3])) {
        cerr << __LINE__ << ": invalid end position " << res[3] << " on " << matchResult[3] << " we will assume the end of it\n";
        res[3] = -1;
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
      if (i != string::npos) {
           res[0] = br.getReferenceId(regstr.substr(0,i));
           if (res[0] == -1) {
              throw runtime_error("invalide first genomic name: " + regstr.substr(0,i));
           }
           res[1] = stoi(regstr.substr(i+1));
           if (!validRefpos(res[0], res[1])) {
              throw runtime_error("invalid reference position " + regstr.substr(i+1));
           }
           res[2] = res[0];
      }
      else {
           res[0] = br.getReferenceId(regstr);
           if (res[0] == -1) {
              throw runtime_error("invalide first genomic name: " + regstr;
           }
      }
   }
   return res;
}
*/

// Parses a region string, does validation (valid ID's, positions), stores in Region struct
// Returns success (true/false)
bool Utilities::ParseRegionString(
      const string& regionString, const BamReader& reader, BamRegion& region)
{
    // -------------------------------
    // parse region string
    // check first for empty string
    if ( regionString.empty() ) {
        cerr << __LINE__ << ":ERROR restionString is empty\n";
        return false;   
    }
    cerr << __LINE__ << ": parsing region: " << regionString << " ...\n";
    try {
       array<int,4> rawreg = parseRegion<BamReader>(regionString, reader);
       region.set(rawreg);
       return true;
    }
    catch (const runtime_error& err) {
       cerr << __LINE__ << ": failed to parse region string\n";
       return false;
    }
    /*

    // non-empty string, look for a colom
    size_t foundFirstColon = regionString.find(':');
    // store chrom strings, and numeric positions
    string startChrom;
    string stopChrom;
    int startPos;
    int stopPos;
    
    // no colon found
    // going to use entire contents of requested chromosome 
    // just store entire region string as startChrom name
    // use BamReader methods to check if its valid for current BAM file
    if ( foundFirstColon == string::npos ) {
        startChrom = regionString;
        startPos   = 0;
        stopChrom  = regionString;
        stopPos    = 0;
    }
    // colon found, so we at least have some sort of startPos requested
    else {
        // store start chrom from beginning to first colon
        startChrom = regionString.substr(0,foundFirstColon);
        // look for ".." after the colon chr:BEG..END
        size_t foundRangeDots = regionString.find("..", foundFirstColon+1);
        // no dots found
        // so we have a startPos but no range
        // store contents before colon as startChrom, after as startPos
        if ( foundRangeDots == string::npos ) { // chr8:9891
            startPos   = atoi( regionString.substr(foundFirstColon+1).c_str() ); 
            stopChrom  = startChrom;
            stopPos    = -1;
        } 
        // ".." found, so we have some sort of range selected
        else {
            // store startPos between first colon and range dots ".."
            startPos = atoi( regionString.substr(foundFirstColon+1, foundRangeDots-foundFirstColon-1).c_str() );
            // look for second colon
            size_t foundSecondColon = regionString.find(':', foundRangeDots+1);
            // no second colon found
            // so we have a "standard" chrom:start..stop input format (on single chrom)
            if ( foundSecondColon == string::npos ) {
                stopChrom  = startChrom;
                stopPos    = atoi( regionString.substr(foundRangeDots+2).c_str() );
            }
            // second colon found
            // so we have a range requested across 2 chrom's
            else {
                stopChrom  = regionString.substr(foundRangeDots+2, foundSecondColon-(foundRangeDots+2));
                stopPos    = atoi( regionString.substr(foundSecondColon+1).c_str() );
            }
        }
    }
    // -------------------------------
    // validate reference IDs & genomic positions
    const RefVector references = reader.GetReferenceData();
    
    // if startRefID not found, return false
    int startRefID = reader.GetReferenceID(startChrom);
    if ( startRefID == -1 ) {
       cerr << __LINE__ << ": genomic name=" << startChrom << " not in bam file\n";
       return false;
    }
    
    // startPos cannot be greater than or equal to reference length
    const RefData& startReference = references.at(startRefID);
    if ( startPos >= startReference.RefLength ) {
       cerr << __LINE__ << ": genomic name=" << startChrom << " start position after end " << startPos << endl;
       return false;
    }
    
    // if stopRefID not found, return false
    int stopRefID = reader.GetReferenceID(stopChrom);
    if ( stopRefID == -1 ) {
       cerr << __LINE__ << ": second genomic name=" << stopChrom << " not in bam file" << endl;
       return false;
    }
    // stopPosition cannot be larger than reference length
    const RefData& stopReference = references.at(stopRefID);
    if ( stopPos > stopReference.RefLength ) {
       cerr << __LINE__ << ": genomic name=" << startChrom << " stop position " << stopPos << " after " << stopReference.RefLength << endl;
       return false;
    }
    // if no stopPosition specified, set to reference end
    if ( stopPos == -1 ) stopPos = stopReference.RefLength;  
    // -------------------------------
    // set up Region struct & return
    
    region.LeftRefID     = startRefID;
    region.LeftPosition  = startPos;
    region.RightRefID    = stopRefID;;
    region.RightPosition = stopPos;
    return true;
    */
}

// Same as ParseRegionString() above, but accepts a BamMultiReader
bool Utilities::ParseRegionString(const string& regionString, 
      const BamMultiReader& reader, BamRegion& region)
{
    // -------------------------------
    // parse region string
  
    // check first for empty string
    if (regionString.empty()) {
        cerr << __LINE__ << ":ERROR regionString is empty\n";
        return false;   
    }
    /*
    // non-empty string, look for a colom
    size_t foundFirstColon = regionString.find(':');
    
    // store chrom strings, and numeric positions
    string startChrom;
    string stopChrom;
    int startPos;
    int stopPos;

    // no colon found
    // going to use entire contents of requested chromosome 
    // just store entire region string as startChrom name
    // use BamReader methods to check if its valid for current BAM file
    if ( foundFirstColon == string::npos ) {
        startChrom = regionString;
        startPos   = 0;
        stopChrom  = regionString;
        stopPos    = -1;
    }
    // colon found, so we at least have some sort of startPos requested
    else {
        regex singleGenomic("([_A-Za-z0-9]+):(\d+)(?:\.\.|-)(\d+)");
        regex doubleGenomic("([_A-Za-z0-9]+):(\d+)(?:\.\.|-)([_A-Za-z0-9]+):(\d+)");
        smatch matchResult;
        if (regex_match(regionString, matchResult, singleGenomic)) {
           startChrom = matchResult[1];
           startPos = matchResult[2];
           stopChrom = startChrom;
           stopPos = matchResult[3];
        }
        else if (regex_match(regionString, matchResult, doubleGenomic)) {
           startChrom = matchResult[1];
           startPos = matchResult[2];
           stopChrom = matchResult[3];
           stopPos = matchResult[4];
        }
        else {
           // store start chrom from beginning to first colon
           startChrom = regionString.substr(0,foundFirstColon);
           // look for ".." after the colon
           size_t foundRangeDots = regionString.find("..", foundFirstColon+1);
           // no dots found
           // so we have a startPos but no range
           // store contents before colon as startChrom, after as startPos
           if ( foundRangeDots == string::npos ) {
               startPos   = atoi( regionString.substr(foundFirstColon+1).c_str() ); 
               stopChrom  = startChrom;
               stopPos    = -1;
           } 
           // ".." found, so we have some sort of range selected
           else {
               // store startPos between first colon and range dots ".."
               startPos = atoi( regionString.substr(foundFirstColon+1, foundRangeDots-foundFirstColon-1).c_str() );
               // look for second colon
               size_t foundSecondColon = regionString.find(':', foundRangeDots+1);
               // no second colon found
               // so we have a "standard" chrom:start..stop input format (on single chrom)
               if ( foundSecondColon == string::npos ) {
                   stopChrom  = startChrom;
                   stopPos    = atoi( regionString.substr(foundRangeDots+2).c_str() );
               }
               // second colon found
               // so we have a range requested across 2 chrom's
               else {
                   stopChrom  = regionString.substr(foundRangeDots+2, foundSecondColon-(foundRangeDots+2));
                   stopPos    = atoi( regionString.substr(foundSecondColon+1).c_str() );
               }
           }
       }
    }

    // -------------------------------
    // validate reference IDs & genomic positions

    const RefVector references = reader.GetReferenceData();

    // if startRefID not found, return false
    int startRefID = reader.GetReferenceID(startChrom);
    if ( startRefID == -1 ) {
       cerr << __LINE__ << ": genomic name=" << startChrom << " not in bam file\n";
       return false;
    }
    // startPos cannot be greater than or equal to reference length
    const RefData& startReference = references.at(startRefID);
    if ( startPos >= startReference.RefLength ) return false;

    // if stopRefID not found, return false
    int stopRefID = reader.GetReferenceID(stopChrom);
    if ( stopRefID == -1 ) {
       cerr << __LINE__ << ": second genomic name=" << stopChrom << " not in bam file" << endl;
       return false;
    }

    // stopPosition cannot be larger than reference length
    const RefData& stopReference = references.at(stopRefID);
    if ( stopPos > stopReference.RefLength ) {
       cerr << __LINE__ << ": stop position " << stopPos << " after chromosome " << stopReference.RefLength << endl;
       return false;
    }
    // if no stopPosition specified, set to reference end
    if ( stopPos == -1 ) stopPos = stopReference.RefLength;

    // -------------------------------
    // set up Region struct & return

    region.LeftRefID     = startRefID;
    region.LeftPosition  = startPos;
    region.RightRefID    = stopRefID;;
    region.RightPosition = stopPos;
    return true;
    */
    try {
       array<int,4> rawreg = parseRegion<BamMultiReader>(regionString, reader);
       region.set(rawreg);
       return true;
    }
    catch (const runtime_error& err) {
       cerr << __LINE__ << ": failed to parse region string: " << regionString << endl;
       return false;
    }
}

void Utilities::Reverse(string& sequence) {
    reverse(sequence.begin(), sequence.end());
}

void Utilities::ReverseComplement(string& sequence) {
    
    // do complement, in-place
    size_t seqLength = sequence.length();
    for ( size_t i = 0; i < seqLength; ++i )
        sequence.replace(i, 1, 1, REVCOMP_LOOKUP[(int)sequence.at(i) - 65]);
    
    // reverse it
    Reverse(sequence);
}

vector<string> Utilities::Split(const string& source, const char delim) {

    stringstream ss(source);
    string field;
    vector<string> fields;

    while ( getline(ss, field, delim) )
        fields.push_back(field);
    return fields;
}

vector<string> Utilities::Split(const string& source, const string& delims) {

    vector<string> fields;

    char* tok;
    char* cchars = new char[source.size()+1];
    char* cstr = &cchars[0];
    strcpy(cstr, source.c_str());
    tok = strtok(cstr, delims.c_str());
    while (tok != NULL) {
        fields.push_back(tok);
        tok = strtok(NULL, delims.c_str());
    }

    delete[] cchars;

    return fields;
}

// returns true if 'source' starts with 'pattern'
bool Utilities::StartsWith(const string& source, const string& pattern) {
    return ( source.find(pattern) == 0 );
}

// returns true if 'source' starts with 'c'
bool Utilities::StartsWith(const std::string &source, const char c) {
    return ( source.find(c) == 0 );
}
