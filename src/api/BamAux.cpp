#include "api/BamAux.h"
#include <cassert>

/// Matchdiff Class ////

using namespace BamTools;

string Matchdiff::toString() const {
   if (mseg.empty()) {
      throw logic_error("Matchdiff::mseg is empty");
   }
   string res = to_string(mseg[0]);
   for (size_t i=0; i<xseg.size(); ++i) {
      res += (xseg[i] + to_string(mseg[i+1]));
   }
   return res;
}

int Matchdiff::length() const {
   int res=0;
   for (auto l : mseg) res += l;
   for (size_t i=0; i < xseg.size(); ++i)
      res += getXseglen(i);
   return res;
}

int Matchdiff::removeBefore(int idx) {
   //cerr << " idx at start of removeBefore: " << idx << endl;
   assert(idx > 0);
   unsigned int d = 0; // index into mseg
   int miscnt = 0; // mismatch count in trimmed region
   int i = 0;
   while (d < mseg.size() && i < idx) {
      if (i+getMseglen(d) > idx) { // idx between i, i+mdvec.first[d]
         // i idx   i+matchlen
         // |  |    |
         // ---=====
         // stop here
         if (d > 0) {
            mseg.erase(mseg.begin(), mseg.begin()+d);
            xseg.erase(xseg.begin(), xseg.begin()+d);
         }
         mseg.front() -= (idx-i);
         break; // job done
      }
      else { // check mismatch segment
         i += getMseglen(d);
         int E = i + getXseglen(d);
         //cerr << " i=" << i << " E=" << E << endl;
         if (isDeletion(d)) { // ref insertion query deletion
            miscnt += getXseglen(d);
            if (E >= idx) { // ^ is meta character need to remove
               //            idx
               //     i  idx E
               //     |  |   |
               // ====~~~~~~~======
               mseg.erase(mseg.begin(), mseg.begin()+d+1);
               xseg.erase(xseg.begin(), xseg.begin()+d+1);
               //idx = E;
               break;
            }
            /*
            else if (E == idx) { // ^ is meta character need to remove
               //            idx
               //     i      E 
               //     |      |
               // ====~~~~~~~======
               mseg.erase(mseg.begin(), mseg.begin()+d+1);
               xseg.erase(xseg.begin(), xseg.begin()+d+1);
               break;
            }
            */
            else {
               i += getXseglen(d);
            }
         }
         else { // in mismatch
            if (E > idx) {
               // i  idx  E
               // |  |    |
               // ---=====
               if (d > 0) {
                  mseg.erase(mseg.begin(), mseg.begin()+d);
                  xseg.erase(xseg.begin(), xseg.begin()+d);
               }
               mseg.front() = 0; // for padding last must be a number
               miscnt += (idx-i);
               xseg.front() = xseg.front().substr(idx-i);
               break; // job done
            }
            else if (E == idx) {
               //         idx
               // i       E
               // |       |
               // ========
               miscnt += getXseglen(d);
               mseg.erase(mseg.begin(), mseg.begin()+d+1);
               xseg.erase(xseg.begin(), xseg.begin()+d+1);
               break;
            }
            else {
               miscnt += getXseglen(d);
               i += getXseglen(d);
            }
         }
      }
      ++d;
   }
   return miscnt;
}

int Matchdiff::removeAfter(int idx) {
   //assert(idx > 0 && idx < length()-1);
   if (idx < 0 || idx > length()-1) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR idx out of range\n";
      throw out_of_range(to_string(idx) + " is outof range in " + string(__func__));
   }
   int i = length()-1;
   int d = mseg.size() - 1; 
   int miscnt = 0; // mismatch count in trimmed region
   while (d > -1 && i > idx) {
      //cerr << "d=" << d << " i=" << i << " idx=" << idx << endl;
      //if (d-1 < 0) { // can be single M segment then ther eis no mismatch
      //   throw logic_error("d out of bound in " + string(__func__));
      //}
      if (i-getMseglen(d) < idx) { // idx between (i-mseglen, i]
         // i-ml   idx  i
         // |      |    |
         //  =======-----
         if (d < static_cast<int>(mseg.size()-1)) {
            mseg.resize(d+1);
            xseg.resize(d);
         }
         mseg.back() -= (i - idx);
         break;
      }
      else { // check mismatch segment
         i -= getMseglen(d); // i now in the mismatch segment
         int B = i - static_cast<int>(getXseglen(d-1)); 
         if (isDeletion(d-1)) { // query deletion
            miscnt += getXseglen(d-1);
            if (B <= idx) {
               mseg.resize(d);
               xseg.resize(d-1);
               break; // done
            }
            /*
            else if (B == idx) {
               mseg.resize(d);
               xseg.resize(d-1);
               break; // done
            }
            */
            else {
               i -= getXseglen(d-1);
            }
         }
         else { // in mismatch
            if (B < idx) {
               // B      idx  i
               // |      |    |
               //  =======-----
               if (d < static_cast<int>(mseg.size()-1)) {
                  mseg.resize(d+1);
                  xseg.resize(d);
               }
               mseg.back() = 0; // for padding last must be a number
               miscnt += (i - idx);
               int new_xlen = xseg.back().size() - i + idx;
               xseg.back().resize(new_xlen);
               break; // job done
            }
            else if (B == idx) {
               // idx
               //  B        i
               //  |        |
               //   =========
               miscnt += getXseglen(d-1);
               mseg.resize(d);
               xseg.resize(d-1);
               break;
            }
            else {
               miscnt += xseg[d-1].size();
               i -= getXseglen(d-1);
            }
         }
      }
      --d;
   }
   return miscnt;
}

