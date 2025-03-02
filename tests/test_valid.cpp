#include <gtest/gtest.h>
#include "api/BamAlignment.h"

using namespace std;
using namespace BamTools;

class BamtoolTestValid : public testing::Test {
   protected:
      BamtoolTestValid() 
         : Test(),
           ba("S618", 53, 159479, 16, -1, -1, 
                 "GGCGGCGGTGGTGGGGGTGGGGGGGGTCCTCCCCCGCCCCCCCCCCCCACGCCTCCTCCCCTCCTCCCGCCCACGCCCCGCTCCCCGCCCCCGGAGCCCCGCGGACGCGACGCCGCGACGAGTAGG",
                 "9II-IIII9II-IIIII-9-99IIII9II--I99I9-I9IIIII9IIIII9II9-I9-IIII9I--9I-I9II9I9III-I9II99-9IIII999I9II---III-I--9II99I9II9IIIIIII",
                 "18M2I106M")
      {
         ba.addTag<string>("MD", BamTools::Constants::BAM_TAG_TYPE_STRING, "8C2C7T8T4C0G71T17");
         //ba.addTag<uint8_t>("NM", BamTools::Constants::BAM_TAG_TYPE_UINT8, 3);
         //ba.addTag("NM", BamTools::Constants::BAM_TAG_TYPE_UINT8, 9);
         ba.addTag<uint8_t>("NM", BamTools::Constants::BAM_TAG_TYPE_UINT8, 9);
      }
      ~BamtoolTestValid() override { }
      BamAlignment ba;
};

TEST_F(BamtoolTestValid, valid) {
    cerr << "bam object\n" << ba << endl;
    EXPECT_TRUE(ba.valid());
}


