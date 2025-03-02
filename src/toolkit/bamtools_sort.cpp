// ***************************************************************************
// bamtools_sort.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 27 March 2012 (DB)
// ---------------------------------------------------------------------------
// Sorts an input BAM file
// ***************************************************************************

#include "bamtools_sort.h"

#include <api/SamConstants.h>
#include <api/BamMultiReader.h>
#include <api/BamWriter.h>
#include <api/algorithms/Sort.h>
#include <utils/bamtools_options.h>
using namespace BamTools;
using namespace BamTools::Algorithms;

#include <cstdio>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <thread> // swang
#include <atomic> // swang
#include <mutex> // swang
#include <cmath> // swang
#include <unistd.h> // swang
using namespace std;

namespace BamTools {
  
// defaults
//
// ** These defaults should be tweaked & 'optimized' per testing ** //
//
//    I say 'optimized' because each system will naturally perform
//    differently.  We will attempt to determine a sensible
//    compromise that should perform well on average.
const unsigned int SORT_DEFAULT_MAX_BUFFER_COUNT  = 10000000;  // max numberOfAlignments for buffer
const unsigned int SORT_DEFAULT_MAX_BUFFER_MEMORY = 409600;    // Mb
    
const unsigned int SORT_DEFAULT_MAX_BUFFER_COUNT_PER_THREAD  = 500000;  // should be bigger max numberOfAlignments for buffer, swang 34298/2 = 17149

} // namespace BamTools

// ---------------------------------------------
// SortSettings implementation

struct SortTool::SortSettings {

    // flags
    bool HasInputBamFilename;
    bool HasMaxBufferCount;
    bool HasMaxBufferMemory;
    bool HasOutputBamFilename;
    bool IsSortingByName;

    // filenames
    string InputBamFilename;
    string OutputBamFilename;

    // parameters
    unsigned int MaxBufferCount;
    unsigned int MaxBufferCountPerThrerad; // swang
    unsigned int MaxBufferMemory;

    // constructor
    SortSettings(void)
        : HasInputBamFilename(false)
        , HasMaxBufferCount(false)
        , HasMaxBufferMemory(false)
        , HasOutputBamFilename(false)
        , IsSortingByName(false)
        , InputBamFilename(Options::StandardIn())
        , OutputBamFilename(Options::StandardOut())
        , MaxBufferCount(SORT_DEFAULT_MAX_BUFFER_COUNT)
        , MaxBufferCountPerThrerad(SORT_DEFAULT_MAX_BUFFER_COUNT_PER_THREAD)  // swang
        , MaxBufferMemory(SORT_DEFAULT_MAX_BUFFER_MEMORY)
    { }
};

// ---------------------------------------------
// SortToolPrivate implementation

class SortTool::SortToolPrivate {
    // ctor & dtor
    public:
        SortToolPrivate(SortTool::SortSettings* settings);
        ~SortToolPrivate(void) { }
        
    // 'public' interface
    public:
        bool Run(void);
        
    // internal methods
    private:
        bool CreateSortedTempFile(vector<BamAlignment>& buffer);
        bool GenerateSortedRuns(void);
        bool MergeSortedRuns(void);
        bool WriteTempFile(const vector<BamAlignment>& buffer, const string& tempFilename);
        void SortBuffer(vector<BamAlignment>& buffer);

        // swang begin
        static const bool usingMultipThread = false; // or set it false and let BamWriter do multiThread; Here we can also do it to speed up
        bool GenerateSortedRunsNew(void);
        // swang end

    // data members
    private:
        SortTool::SortSettings* m_settings;
        string m_tempFilenameStub;
        //int m_numberOfRuns;
        atomic<int> m_numberOfRuns {0}; // swang
        //std::mutex m_mutex_create; // swang
        std::mutex m_mutex_write; // swang
        string m_headerText;
        RefVector m_references;
        /**
         * Hold the names of the tmpfiles
         */
        vector<string> m_tempFilenames;
};

// constructor
SortTool::SortToolPrivate::SortToolPrivate(SortTool::SortSettings* settings) 
    : m_settings(settings)
    //, m_numberOfRuns(0)
{
    // set filename stub depending on inputfile path
    // that way multiple sort runs don't trip on each other's temp files
    if ( m_settings) {
        //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ <<  " InputBamFilename " << m_settings->InputBamFilename << endl; // Col0_C1.100k.bam
        size_t extensionFound = m_settings->InputBamFilename.find(".bam");
        if ( extensionFound != string::npos )
            m_tempFilenameStub = m_settings->InputBamFilename.substr(0,extensionFound);
        m_tempFilenameStub.append(".sort.temp.");
    }
}

// generates mutiple sorted temp BAM files from single unsorted BAM file
bool SortTool::SortToolPrivate::GenerateSortedRuns(void) {
    // open input BAM file
    BamReader reader;
    if ( !reader.Open(m_settings->InputBamFilename) ) {
        cerr << "bamtools sort ERROR: could not open " << m_settings->InputBamFilename
             << " for reading... Aborting." << endl;
        return false;
    }

    // swang begin
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ <<  " InputBamFilename " << m_settings->InputBamFilename << endl; // Col0_C1.100k.bam
    // swang end

    // get basic data that will be shared by all temp/output files 
    SamHeader header = reader.GetHeader();
    if ( !header.HasVersion() )
        header.Version = Constants::SAM_CURRENT_VERSION;
    header.SortOrder = ( m_settings->IsSortingByName
                       ? Constants::SAM_HD_SORTORDER_QUERYNAME
                       : Constants::SAM_HD_SORTORDER_COORDINATE );
    m_headerText = header.ToString();
    m_references = reader.GetReferenceData();
    

    // set up alignments buffer
    BamAlignment al;
    vector<BamAlignment> buffer;

    buffer.reserve( (size_t)(m_settings->MaxBufferCount*1.1) );
    bool bufferFull = false;

    // if sorting by name, we need to generate full char data
    // so can't use GetNextAlignmentCore()
    if ( m_settings->IsSortingByName ) {
        // iterate through file
        while ( reader.GetNextAlignment(al)) {
            // check buffer's usage
            bufferFull = ( buffer.size() >= m_settings->MaxBufferCount );
            // store alignments until buffer is "full"
            if ( !bufferFull )
                buffer.push_back(al);
            // if buffer is "full"
            else {
                // so create a sorted temp file with current buffer contents
                // then push "al" into fresh buffer
                cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " call CreateSortedTempFile" << endl; // swang
                CreateSortedTempFile(buffer);
                buffer.push_back(al);
            }
        }
        cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " swang IsSortingByName GenerateSortedRuns buffer size" << buffer.size() << endl; // swang
    }
    // sorting by position, can take advantage of GNACore() speedup
    else {

        // iterate through file
        while ( reader.GetNextAlignmentCore(al) ) {
            // check buffer's usage
            bufferFull = ( buffer.size() >= m_settings->MaxBufferCount );

            // store alignments until buffer is "full"
            if ( !bufferFull )
            {
                buffer.push_back(al);
            }
            else {
                // create a sorted temp file with current buffer contents
                // then push "al" into fresh buffer
                CreateSortedTempFile(buffer);
                buffer.push_back(al);
            }
        }
    }
    // handle any leftover buffer contents
    if ( !buffer.empty() )
    {
        CreateSortedTempFile(buffer); // hit here
    }
    
    // close reader & return success
    reader.Close();
    return true;
}


bool SortTool::SortToolPrivate::CreateSortedTempFile(vector<BamAlignment>& buffer) {
    // do sorting
    //sleep(1);
    m_mutex_write.lock();
    SortBuffer(buffer);
    // write sorted contents to temp file, store success/fail
    stringstream tempStr;
    tempStr << m_tempFilenameStub << m_numberOfRuns++; // swang atomic<int>m_numberOfRuns
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << ":mtempStr.str() = " << tempStr.str() << endl; // Col0_C1.100k.sort.temp.0
    bool success = WriteTempFile( buffer, tempStr.str() );
    // save temp filename for merging later
    m_tempFilenames.push_back(tempStr.str());
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << ": m_tempFilenames.size() = " << m_tempFilenames.size() << endl;
    // clear buffer contents & update run counter
    buffer.clear();
    //++m_numberOfRuns; // swang
    // return success/fail of writing to temp file
    // TODO: a failure returned here is not actually caught and handled anywhere
    m_mutex_write.unlock();
    return success;
}

// merges sorted temp BAM files into single sorted output BAM file
bool SortTool::SortToolPrivate::MergeSortedRuns(void) {
    // open up multi reader for all of our temp files
    // this might get broken up if we do a multi-pass system later ??
    BamMultiReader multiReader;
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << ":m_settings->OutputBamFilename = " << m_settings->OutputBamFilename << endl; // myOutput3B.bam swang

    //if ( !multiReader.Open(m_tempFilenames) ) {
    if ( !multiReader.CreateThreads(m_tempFilenames) ) { // swang we can try this

        cerr << __FILE__ << ":" << __LINE__ << ":" << __func__ 
           << ":ERROR could not open BamMultiReader for merging temp files... Aborting.\n";
        return false;
    }

    // open writer for our completely sorted output BAM file
    BamWriter mergedWriter;
    if ( !mergedWriter.Open(m_settings->OutputBamFilename, m_headerText, m_references) ) {
        cerr << "bamtools sort ERROR: could not open " << m_settings->OutputBamFilename
             << " for writing... Aborting." << endl;
        multiReader.Close();
        return false;
    }
    // while data available in temp files
    BamAlignment al;
    while ( multiReader.GetNextAlignmentCore(al) )
        mergedWriter.SaveAlignment(al);
  
    // close files
    multiReader.Close();
    mergedWriter.Close();
    
    // delete all temp files
    for (auto& tf : m_tempFilenames) {
       remove(tf.c_str());
    }
  
    // return success
    return true;
}

bool SortTool::SortToolPrivate::Run(void) {
    // this does a single pass, chunking up the input file into smaller sorted temp files, 
    // then write out using BamMultiReader to handle merging
    // swang begin
    if(usingMultipThread)
    {
        if ( GenerateSortedRunsNew() )
            return MergeSortedRuns();
        else
            return false;
    }
    else{
        if ( GenerateSortedRuns() )
            return MergeSortedRuns();
        else
            return false;
    }

    // swang end
} 
    
void SortTool::SortToolPrivate::SortBuffer(vector<BamAlignment>& buffer) {
    // ** add further custom sort options later ?? **
    // sort buffer by desired method
    if ( m_settings->IsSortingByName )
        std::stable_sort( buffer.begin(), buffer.end(), Sort::ByName() );
    else
        std::stable_sort( buffer.begin(), buffer.end(), Sort::ByPosition() );
}

// swang begin:
// generates mutiple sorted temp BAM files from single unsorted BAM file
bool SortTool::SortToolPrivate::GenerateSortedRunsNew(void) {
    // open input BAM file
    BamReader reader;
    if ( !reader.Open(m_settings->InputBamFilename) ) {
        cerr << "bamtools sort ERROR: could not open " << m_settings->InputBamFilename
             << " for reading... Aborting." << endl;
        return false;
    }

    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ <<  " InputBamFilename " << m_settings->InputBamFilename << endl; // Col0_C1.100k.bam

    // get basic data that will be shared by all temp/output files
    SamHeader header = reader.GetHeader();
    if ( !header.HasVersion() )
        header.Version = Constants::SAM_CURRENT_VERSION;
    header.SortOrder = ( m_settings->IsSortingByName
                       ? Constants::SAM_HD_SORTORDER_QUERYNAME
                       : Constants::SAM_HD_SORTORDER_COORDINATE );
    m_headerText = header.ToString();
    m_references = reader.GetReferenceData();

    // set up alignments buffer
    BamAlignment al;
    //vector<BamAlignment> buffer;
    vector<BamAlignment> sub_buffer;
    vector<vector<BamAlignment>> sub_bufferVec;

    int NumOfThread = 0;
    //buffer.reserve( (size_t)(m_settings->MaxBufferCount*1.1) );
    sub_buffer.reserve( (size_t)(m_settings->MaxBufferCountPerThrerad*1.1) );
    bool bufferFull = false;

    // if sorting by name, we need to generate full char data
    // so can't use GetNextAlignmentCore()
    if ( m_settings->IsSortingByName ) {
        // iterate through file
        //u_int32_t count_al = 0;
        while ( reader.GetNextAlignment(al)) {
            //count_al++;
            // check buffer's usage
            // bufferFull = ( buffer.size() >= m_settings->MaxBufferCount );

            bufferFull = ( sub_buffer.size() >= m_settings->MaxBufferCountPerThrerad );
            // store alignments until buffer is "full"
            if ( !bufferFull )
            {
                //buffer.push_back(al);
                sub_buffer.push_back(al);
            }
            // if buffer is "full"
            else {
                // so create a sorted temp file with current buffer contents
                // then push "al" into fresh buffer
                //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " call CreateSortedTempFile" << endl;
                //CreateSortedTempFile(buffer);
                //buffer.push_back(al);

                // ??? CreateSortedTempFile(sub_buffer);
                sub_buffer.push_back(al);
                sub_bufferVec.push_back(sub_buffer);
                //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " sub_bufferVec.size() = " << sub_bufferVec.size() << endl;
                NumOfThread++;
                sub_buffer.clear();
                bufferFull = false;
            }
        }

        if (!sub_buffer.empty())
        {
            sub_bufferVec.push_back(sub_buffer);
            NumOfThread++;
        }
        //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " swang IsSortingByName GenerateSortedRuns buffer size" << buffer.size()
        //     << " count_al " << count_al << endl; // swang
    }
    // sorting by position, can take advantage of GNACore() speedup
    else {

        // iterate through file
        while ( reader.GetNextAlignmentCore(al) ) {
            // check buffer's usage
            // bufferFull = ( buffer.size() >= m_settings->MaxBufferCount );
            bufferFull = ( sub_buffer.size() >= m_settings->MaxBufferCountPerThrerad );

            // store alignments until buffer is "full"
            if ( !bufferFull )
            {
                //buffer.push_back(al);
                sub_buffer.push_back(al);
            }
            else {
                // create a sorted temp file with current buffer contents
                // then push "al" into fresh buffer
                //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " call CreateSortedTempFile " << endl;//for local_al_name " << local_al_name << endl;
                // ??? CreateSortedTempFile(sub_buffer);
                sub_buffer.push_back(al);
                sub_bufferVec.push_back(sub_buffer);
                //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " sub_bufferVec.size() = " << sub_bufferVec.size() << endl; // swang
                NumOfThread++;
                sub_buffer.clear();
                bufferFull = false;
            }
        }
        if (!sub_buffer.empty())
        {
            sub_bufferVec.push_back(sub_buffer);
            NumOfThread++;
        }
        // buffer.size() = 34298  count_al = 34298
        //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " sub_buffer size " << sub_buffer.size() << endl; // 34298
    }
    // handle any leftover buffer contents
    //if ( !buffer.empty() )
    if ( !sub_bufferVec.empty() )
    {
        //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " sub_bufferVec.size() = " << sub_bufferVec.size() << endl;
        // Launch a thread for each reader
        std::vector<std::thread> threads;

        for (int ix = 0; ix < NumOfThread; ix++)
        {
            vector<BamAlignment>& buffer_local = sub_bufferVec[ix];
            threads.emplace_back([&]() {
                        CreateSortedTempFile(buffer_local);
                });

        }
        // Join all threads
        for (auto& th : threads) {
            if (th.joinable()) {
                //cout  << __FILE__ << ":" << __LINE__ << ":" << __func__ << " join thread. " << endl;
                th.join();
                //cout  << __FILE__ << ":" << __LINE__ << ":" << __func__ << " done join thread. " << endl;
            }
        }
    }

    // close reader & return success
    reader.Close();
    return true;
}
// swang end

bool SortTool::SortToolPrivate::WriteTempFile(const vector<BamAlignment>& buffer,
                                              const string& tempFilename)
{
    // open temp file for writing
    BamWriter tempWriter;
    if ( !tempWriter.Open(tempFilename, m_headerText, m_references) ) {
        cerr << "bamtools sort ERROR: could not open " << tempFilename
             << " for writing." << endl;
        return false;
    }
  
    //cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " tempFilename :" << tempFilename << ", buffer.size() = " << buffer.size() << endl; // 34298
    // write data
    vector<BamAlignment>::const_iterator buffIter = buffer.begin();
    vector<BamAlignment>::const_iterator buffEnd  = buffer.end();
    for ( ; buffIter != buffEnd; ++buffIter )  {
        const BamAlignment& al = (*buffIter);
        tempWriter.SaveAlignment(al);
    }
  
    // close temp file & return success
    tempWriter.Close();
    return true;
}

// ---------------------------------------------
// SortTool implementation

SortTool::SortTool(void)
    : AbstractTool()
    , m_settings(new SortSettings)
    , m_impl(0)
{
    // set program details
    Options::SetProgramInfo("bamtools sort", "sorts a BAM file", "[-in <filename>] [-out <filename>] [sortOptions]");
    // set up options
    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");
    Options::AddValueOption("-in",  "BAM filename", "the input BAM file",  "",
                            m_settings->HasInputBamFilename,  m_settings->InputBamFilename,
                            IO_Opts, Options::StandardIn());
    Options::AddValueOption("-out", "BAM filename", "the output BAM file", "",
                            m_settings->HasOutputBamFilename, m_settings->OutputBamFilename,
                            IO_Opts, Options::StandardOut());

    OptionGroup* SortOpts = Options::CreateOptionGroup("Sorting Methods");
    Options::AddOption("-byname", "sort by alignment name", m_settings->IsSortingByName, SortOpts);

    OptionGroup* MemOpts = Options::CreateOptionGroup("Memory Settings");
    Options::AddValueOption("-n",   "count", "max number of alignments per tempfile", "",
                            m_settings->HasMaxBufferCount,  m_settings->MaxBufferCount,
                            MemOpts, SORT_DEFAULT_MAX_BUFFER_COUNT);
    Options::AddValueOption("-mem", "Mb", "max memory to use", "",
                            m_settings->HasMaxBufferMemory, m_settings->MaxBufferMemory,
                            MemOpts, SORT_DEFAULT_MAX_BUFFER_MEMORY);
}

SortTool::~SortTool(void) {

    delete m_settings;
    m_settings = 0;

    delete m_impl;
    m_impl = 0;
}

int SortTool::Help(void) {
    Options::DisplayHelp();
    return 0;
}

int SortTool::Run(int argc, char* argv[]) {

    // parse command line arguments
    Options::Parse(argc, argv, 1);

    // initialize SortTool with settings
    m_impl = new SortToolPrivate(m_settings);

    // run SortTool, return success/fail
    if ( m_impl->Run() )
        return 0;
    else
        return 1;
}
