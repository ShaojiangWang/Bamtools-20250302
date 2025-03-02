#include <gtest/gtest.h>
#include "api/BamAlignment.h"

using namespace std;
using namespace BamTools;

class BamtoolTest : public testing::Test {
   protected:
      BamtoolTest() 
         : Test(),
           ba("S123745075", 53, 129951, 99, 53, 130061, 
               "GCTCATGTATGCTTGAACGACAAATAAAAGTTCGGGGGGGAGAAGAGAGGAGAGAGAGAGAGCGAAGGGGAGAGAGGGGGGAGAGGGGGGGGGGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA", 
               "LZjjddZjjZjjZjd]jjjjdETEEjdZZjZEjZZ8ZjdZETjjLdEd]ZZjjjjSjj]Ljj@0L-AT%LELdTE\\IL8\\TIT!\\8\\LD\\\\\\,L8I\\\\\\T\\8888HHH$$8$HHH8HH8H8$H$HHH8HH88H", 
               "83M1D50M")
      {
         ba.addTag<string>("MD", BamTools::Constants::BAM_TAG_TYPE_STRING, "62A2C27^A40");
         //ba.addTag<uint8_t>("NM", BamTools::Constants::BAM_TAG_TYPE_UINT8, 3);
         ba.addTag("NM", BamTools::Constants::BAM_TAG_TYPE_UINT8, 3);
      }
      ~BamtoolTest() override { }
      BamAlignment ba;
};

TEST_F(BamtoolTest, cutBefore) {
    EXPECT_FALSE(ba.valid());
    cerr << "before cut\n" << ba << endl;
    ba.chopBefore(130041);
    cerr << "after cut\n" << ba << endl;
    ASSERT_EQ("44M", ba.getCigarString());
}


