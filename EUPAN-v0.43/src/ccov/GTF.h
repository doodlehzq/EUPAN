#ifndef GTF_H
#define GTF_H

#include <vector>
#include <string>
using namespace std;

typedef enum GeneStructureType{
  START_CODON,
  CDS,
  STOP_CODON,
  EXON,
  THREEUTR,
  FIVEUTR,
  GENE,
  UNDEF,
} GeneStructureType;

typedef enum SequenceDirection{
  PLUS,
  MINUS,
} SequenceDirection;

class GeneItem{
 public:
  GeneStructureType type;
  string chr;
  string source;
  string phase;
  string score;
  SequenceDirection symbol;
  unsigned int begin;
  unsigned int end;
  string gene_id;
  string transcript_id;
 GeneItem():type(UNDEF),chr("undef"),source("undef"),phase("."),score(""),begin(0),end(0),gene_id(""),transcript_id(""){}; 
  //GeneStructure(const GeneStructureType t, const unsigned int b, const unsigned int e): type(t), begin(b), end(e){}; 
};

class GeneStructure{
 public:
  string name;
  string chr;
  unsigned int begin;
  unsigned int end;
  int complete;
  int single_exon;
  SequenceDirection direction;
  vector<GeneItem> itemseq;
  double cdsCov;
  double cdsDep;
  double cdsEven;
  double geneCov;
  double geneDep;
  double geneEven;
  void update();
  void clear();
};

class GTF{
 public:
  int size;
  string file;
  vector<GeneStructure> GeneStructureSequence;
  void ReadGTF(const string &file);
  void ReadGFF(const string &file);
  void WriteGTF();
  void clear();
};




class GTFFileList{
 public:
  vector<string> list;
  void CheckFileState();
};

class GTFList{
 public:
  vector<GTF> gtflist;
  void ReadGTFList(const GTFFileList &filelist);
};

#endif
