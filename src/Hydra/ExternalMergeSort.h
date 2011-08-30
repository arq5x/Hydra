/****************************************************************************
   ExternalMergeSort.h (c) 2009 Aaron Quinlan

   Hall Lab
   Department of Biochemistry and Molecular Genetics
   University of Virginia

   All rights reserved.

   Object for storing the BamAlignments (N=2) for
   a readpair mapping combination
****************************************************************************/
#ifndef EXTMERGESORT_H
#define EXTMERGESORT_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <queue>
#include <cstdio>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <libgen.h> //for basename()
using namespace std;

// "local" includes
#include "gzstream.h"


bool isRegularFile(const string& filename);
bool isGzipFile(const string& filename);
//STLized version of basename()
// (because POSIX basename() modifies the input string pointer)
// Additionally: removes any extension the basename might have.
std::string stl_basename(const std::string& path);


template <class T>
class MERGE_DATA {

public:
    // data
    T data;
    istream *stream;
    bool (*compFunc)(const T &a, const T &b);

    // constructor
    MERGE_DATA (const T &data,
                istream *stream,
                bool (*compFunc)(const T &a, const T &b))
    :
        data(data),
        stream(stream),
        compFunc(compFunc)
    {}

    // comparison operator for maps keyed on this structure
    bool operator < (const MERGE_DATA &a) const
    {
        // recall that priority queues try to sort from
        // highest to lowest. thus, we need to negate.
        return !(compFunc(data, a.data));
    }
};


//************************************************
// DECLARATION
// Class methods and elements
//************************************************
template <class T>
class ExtMergeSort {

public:
    // constructor
    ExtMergeSort(const string &inFile,
                 ostream *out,
                 bool (*compareFunction)(const T &a, const T &b),
                 int  maxBufferSize  = 10000,
                 bool compressOutput = false,
                 string tempPath     = "");
    // destructor
    ~ExtMergeSort(void);

    // drives the creation of sorted sub-files stored on disk.
    void DivideAndSort();

    // drives the merging of the sorted temp files.
    // final, sorted and merged output is written to "out".
    void Merge();

private:
    string _inFile;
    bool (*_compareFunction)(const T &a, const T &b);
    string _tempPath;
    vector<string>    _vTempFileNames;
    vector<ifstream*>  _vTempFiles;
    unsigned int _maxBufferSize;
    unsigned int _runCounter;
    bool _compressOutput;
    bool _tempFileUsed;
    ostream *_out;

    void WriteToTempFile(const vector<T> &lines);
    void OpenTempFiles();
    void CloseTempFiles();
};



//************************************************
// IMPLEMENTATION
// Class methods and elements
//************************************************
// constructor
template <class T>
ExtMergeSort<T>::ExtMergeSort (const string &inFile,
                               ostream *out,
                               bool (*compareFunction)(const T &a, const T &b),
                               int maxBufferSize,
                               bool compressOutput,
                               string tempPath)
    : _inFile(inFile)
    , _out(out)
    , _compareFunction(compareFunction)
    , _tempPath(tempPath)
    , _maxBufferSize(maxBufferSize)
    , _runCounter(0)
    , _compressOutput(compressOutput)
{}


// destructor
template <class T>
ExtMergeSort<T>::~ExtMergeSort(void)
{}


// Sort
template <class T>
void ExtMergeSort<T>::DivideAndSort() {

    istream *input = new ifstream(_inFile.c_str(), ios::in);
    // gzipped
    if ((isGzipFile(_inFile) == true) && (isRegularFile(_inFile) == true)) {
        delete input;
        input = new igzstream(_inFile.c_str(), ios::in);
    }

    // bail unless the file is legit
    if ( input->good() == false ) {
        cerr << "Error: The requested input file (" << _inFile << ") could not be opened. Exiting!" << endl;
        exit (1);
    }
    vector<T> lineBuffer;
    lineBuffer.reserve(100000);
    unsigned int totalBytes = 0;  // track the number of bytes consumed so far.

    // track whether or not we actually had to use a temp
    // file based on the memory that was allocated
    _tempFileUsed = false;

    // keep reading until there is no more input data
    T line;
    while (*input >> line) {

        // add the current line to the buffer
        // and track the memory used.
        lineBuffer.push_back(line);
        totalBytes += sizeof(line);  // buggy?

        // sort the buffer and write to a temp file
        // if we have filled up our quota
        if (totalBytes > _maxBufferSize) {
            sort(lineBuffer.begin(), lineBuffer.end(), *_compareFunction);
            WriteToTempFile(lineBuffer);
            // clear the buffer for the next run
            lineBuffer.clear();
            _tempFileUsed = true;
            totalBytes = 0;
        }
    }

    // handle the run (if any) from the last chunk of the input file.
    if (lineBuffer.empty() == false) {
        // write the last "chunk" to the tempfile if
        // a temp file had to be used (i.e., we exceeded the memory)
        if (_tempFileUsed == true) {
            sort(lineBuffer.begin(), lineBuffer.end(), *_compareFunction);
            WriteToTempFile(lineBuffer);
        }
        // otherwise, the entire file fit in the memory given,
        // so we can just dumpto stdout.
        else {
            sort(lineBuffer.begin(), lineBuffer.end(), *_compareFunction);
            for (size_t i = 0; i < lineBuffer.size(); ++i)
                *_out << lineBuffer[i];
        }
    }
}


template <class T>
void ExtMergeSort<T>::WriteToTempFile(const vector<T> &lineBuffer) {
    // name the current tempfile
    stringstream tempFileSS;
    if (_tempPath.size() == 0)
        tempFileSS << _inFile << "." << _runCounter;
    else
        tempFileSS << _tempPath << "/" << stl_basename(_inFile) << "." << _runCounter;
    string tempFileName = tempFileSS.str();

    // do we want a regular or a gzipped tempfile?
    ofstream *output;
    //if (_compressOutput == true)
        //output = new ogzstream(tempFileName.c_str(), ios::out);
    //else
    output = new ofstream(tempFileName.c_str(), ios::out);

    // write the contents of the current buffer to the temp file
    for (size_t i = 0; i < lineBuffer.size(); ++i)
        *output << lineBuffer[i];

    // update the tempFile number and add the tempFile to the list of tempFiles
    ++_runCounter;
    output->close();
    delete output;
    _vTempFileNames.push_back(tempFileName);
}


//---------------------------------------------------------
// MergeDriver()
//
// Merge the sorted temp files.
// uses a priority queue, with the values being a pair of
// the record from the file, and the stream from which the record came.
// SEE: http://stackoverflow.com/questions/2290518/c-n-way-merge-for-external-sort, post from Eric Lippert.
//----------------------------------------------------------
template <class T>
void ExtMergeSort<T>::Merge() {

    // we can skip this step if there are no temp files to
    // merge.  That is, the entire inout file fit in memory
    // and thus we just dumped to stdout.
    if (_tempFileUsed == false)
        return;


    // open the sorted temp files up for merging.
    // loads ifstream pointers into _vTempFiles
    OpenTempFiles();

    // priority queue for the buffer.
    priority_queue< MERGE_DATA<T> > outQueue;

    // extract the first line from each temp file
    T line;
    for (size_t i = 0; i < _vTempFiles.size(); ++i) {
        *_vTempFiles[i] >> line;  // overloaded >> for PAIR struct
        outQueue.push( MERGE_DATA<T>(line, _vTempFiles[i], _compareFunction) );
    }

    // keep working until the queue is empty
    while (outQueue.empty() == false) {
        // grab the lowest element, print it, then ditch it.
        MERGE_DATA<T> lowest = outQueue.top();

        // write the entry from the top of the queue
        // and then delete that entry forever.
        *_out << lowest.data;
        outQueue.pop();

        // add the next line from the lowest stream (above) to the queue
        // as long as it's not EOF.
        *(lowest.stream) >> line;  // overloaded >> for PAIR struct
        if (*(lowest.stream))
            outQueue.push( MERGE_DATA<T>(line, lowest.stream, _compareFunction) );
    }

    // clean up the gazillion temp files.
    CloseTempFiles();
}


template <class T>
void ExtMergeSort<T>::OpenTempFiles() {
    for (size_t i=0; i < _vTempFileNames.size(); ++i) {

        ifstream *file;

        // not gzipped
        if ((isGzipFile(_vTempFileNames[i]) == false) && (isRegularFile(_vTempFileNames[i]) == true)) {
            file = new ifstream(_vTempFileNames[i].c_str(), ios::in);
        }
        // gzipped
        //else if ((isGzipFile(_vTempFileNames[i]) == true) && (isRegularFile(_vTempFileNames[i]) == true)) {
        //    file = new igzstream(_vTempFileNames[i].c_str(), ios::in);
        //}

        if (file->good() == true) {
            // add a pointer to the opened temp file to the list
            _vTempFiles.push_back(file);
        }
        else {
            cerr << "Unable to open temp file (" << _vTempFileNames[i]
                 << ").  I suspect a limit on number of open file handles.  Exiting."
                 << endl;
             exit(1);
        }
    }
}


template <class T>
void ExtMergeSort<T>::CloseTempFiles() {
    // delete the pointers to the temp files.
    for (size_t i=0; i < _vTempFiles.size(); ++i) {
        _vTempFiles[i]->close();
        delete _vTempFiles[i];
    }
    // delete the temp files from the file system.
    for (size_t i=0; i < _vTempFileNames.size(); ++i) {
        remove(_vTempFileNames[i].c_str());  // remove = UNIX "rm"
    }
}


/*
   returns TRUE if the file is a regular file:
     not a pipe/device.

   This implies that the file can be opened/closed/seek'd multiple times without losing information
 */
bool isRegularFile(const string& filename) {
       struct stat buf ;
       int i;

       i = stat(filename.c_str(), &buf);
       if (i!=0) {
               cerr << "Error: can't determine file type of '" << filename << "': " << strerror(errno) << endl;
               exit(1);
       }
       if (S_ISREG(buf.st_mode))
               return true;

       return false;
}


/*
   returns TRUE if the file has a GZIP header.
   Should only be run on regular files.
 */
bool isGzipFile(const string& filename) {
       //see http://www.gzip.org/zlib/rfc-gzip.html#file-format
       struct  {
               unsigned char id1;
               unsigned char id2;
               unsigned char cm;
       } gzip_header;
       ifstream f(filename.c_str(), ios::in|ios::binary);
       if (!f)
               return false;

       if (!f.read((char*)&gzip_header, sizeof(gzip_header)))
               return false;

       if ( gzip_header.id1 == 0x1f
                       &&
                       gzip_header.id2 == 0x8b
                       &&
                       gzip_header.cm == 8 )
               return true;

       return false;
}


string stl_basename(const string &path) {
    string result;

    char* path_dup = strdup(path.c_str());
    char* basename_part = basename(path_dup);
    result = basename_part;
    free(path_dup);

    size_t pos = result.find_last_of('.');
    if (pos != string::npos )
        result = result.substr(0,pos);

    return result;
}


#endif /* EXTMERGESORT_H */

