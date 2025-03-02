// ***************************************************************************
// BamReader.cpp (c) 2009 Derek Barnett, Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 29 July 2013 (DB)
// ---------------------------------------------------------------------------
// Provides read access to BAM files.
// ***************************************************************************

#include "api/BamReader.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace std;
using namespace BamTools;

BamReader::BamReader() 
// swang : d(new BamTools::Internal::BamReaderPrivate(this)) { }
: d(new BamTools::Internal::BamReaderPrivate()) { }
// swang: why do we need "this" passed to Internal::BamReaderPrivate for its m_parent which is not used in the class ?

BamReader::~BamReader(void) {
   delete d;
   d = 0;
}

bool BamReader::Close(void) {
   return d->Close();
}

const string& BamReader::GetFilename(void) const {
   return d->Filename();
}

bool BamReader::IsOpen(void) const {
   return d->IsOpen();
}

bool BamReader::Jump(int refID, int position) {
   return d->SetRegion( BamRegion(refID, position) );
}
bool BamReader::Open(const std::string& filename) {
   return d->Open(filename);
}
bool BamReader::Rewind(void) {
   return d->Rewind();
}
bool BamReader::SetRegion(const BamRegion& region) { 
   return d->SetRegion(region);
}
bool BamReader::SetRegion(const int& leftRefID, const int& leftPos,
               const int& rightRefID, const int& rightPos) 
{
   return d->SetRegion(BamRegion(leftRefID, leftPos, rightRefID, rightPos));
}

bool BamReader::GetNextAlignment(BamAlignment& alignment) {
   // uses Internal::BamReaderPrivate to do the work
   return d->GetNextAlignment(alignment);
}

BamAlignment* BamReader::next() {
   BamAlignment* ptr = new BamAlignment();
   if (nextAlignment(*ptr)) {
      return ptr;
   }
   delete ptr;
   return nullptr;
}

BamAlignment* BamReader::nextCore() {
   BamAlignment* ptr = new BamAlignment();
   if (nextAlignmentCore(*ptr)) {
      return ptr;
   }
   delete ptr;
   return nullptr;
}

bool BamReader::GetNextAlignmentCore(BamAlignment& alignment) {
   return d->GetNextAlignmentCore(alignment);
}
std::string BamReader::GetHeaderText(void) const {
   return d->GetHeaderText();
}
int BamReader::GetReferenceCount(void) const {
   return d->GetReferenceCount();
}

int BamReader::GetReferenceID(const std::string& refName) const {
   return d->GetReferenceID(refName);
}

bool BamReader::CreateIndex(const BamIndex::IndexType& type) {
   return d->CreateIndex(type);
}
bool BamReader::HasIndex(void) const {
   return d->HasIndex();
}
bool BamReader::LocateIndex(const BamIndex::IndexType& preferredType)
{
   return d->LocateIndex(preferredType);
}
bool BamReader::OpenIndex(const std::string& indexFilename) {
   return d->OpenIndex(indexFilename);
}
void BamReader::SetIndex(BamIndex* index) {
   d->SetIndex(index);
}

//std::string BamReader::GetErrorString(void) const {
//   return d->GetErrorString();
//}
        
vector<pair<string,int>> BamReader::getReferenceMetaData() const {
   const RefVector& tmp = GetReferenceData();
   //vector<pair<string, int> > res(tmp.size());
   vector<pair<string, int> > res;
   res.reserve(tmp.size());
   for (size_t i=0; i<tmp.size(); ++i) {
      //res[i]=make_pair(tmp[i].getRefname(), tmp[i].getReflength());
      res.push_back(make_pair(tmp[i].getRefname(), tmp[i].getReflength()));
   }
   return res;
}

const string& BamReader::getReferenceName(int refid) const {
   const RefVector& rd = GetReferenceData();
   if (rd.empty()) {
      throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR BamReader::GetReferenceData() returned empty object");
   }
   if (refid < 0 || refid >= static_cast<int>(rd.size())) {
      cerr << __FILE__ << ":" << __LINE__ << ":ERROR refid:" << refid << " is invalid\n";
      throw logic_error(string(__FILE__) + ":" + to_string(__LINE__) + ":ERROR reference id " + to_string(refid) + " is invalid");
   }
   return rd[refid].getRefname();
}
