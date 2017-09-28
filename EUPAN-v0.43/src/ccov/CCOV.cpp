#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include <zconf.h>
#include <time.h>
#include <stdio.h>
#include <map>
#include "GTF.h"
#include <api/BamReader.h>

using namespace std;
#define BUFFER 2048

void PrintUsage(){
  cerr << "ccov <annotation.gtf/.gff> <mapping.bam>" << endl;
  cerr << "anntation file is suffux-sensitive and only .gtf and .gff can be used." << endl;
}
void initMapper(map<string, vector<unsigned int> > &mapper,string genomeFile);
void checkCDSCoverage(map<string,vector<unsigned int> > &mapper, GTF &gtf);
void checkGeneCoverage(map<string,vector<unsigned int> > &mapper, GTF &gtf);
void checkGeneEvenness(map<string,vector<unsigned int> > &mapper, GTF &gtf);
void checkCDSEvenness(map<string,vector<unsigned int> > &mapper, GTF &gtf);
void PrintChrSize(map<string,vector<unsigned int> > mapper);

int main(int argc, char* argv[]){
  string of = "output";
  if(argc<3){
    PrintUsage();
    exit(-1);
  }
  string annotationFile=argv[1];
  string bamFile=argv[2];

  GTF gtf;
  cerr <<"Read annotation file: "<<annotationFile<<" ......"<<endl;
  string suffix=annotationFile.substr(annotationFile.size()-3,3);
  //  cout<< suffix << "\t" << suffix.size()<<endl;
  if(suffix=="gtf"){
    gtf.ReadGTF(annotationFile);
  }
  else if(suffix=="gff"){
    gtf.ReadGFF(annotationFile);
  }
  else{
    gtf.ReadGTF(annotationFile);
  }

  cerr << "Transcript number: " << gtf.size <<endl;


  cerr<<"Init mapper ..." << endl;

  BamTools::BamReader reader;
  if(!reader.Open(bamFile)){
    cerr<<"Couldn't open bam file:" <<bamFile<<endl;
    exit(-1);
  }
  BamTools::RefVector refvector=reader.GetReferenceData();
    map<string, vector<unsigned int> > mapper;
    vector<BamTools::RefData>::iterator iter=refvector.begin();
    while(iter!=refvector.end()){
    vector<unsigned int> tmp(iter->RefLength);
    mapper[iter->RefName]=tmp;
    ++iter;
    }

  //reading SAM
  cerr <<"Read mapping file ..."<<endl;

  BamTools::BamAlignment al;
  while(reader.GetNextAlignment(al)){
    vector<BamTools::CigarOp> cigar=al.CigarData;
    vector<BamTools::CigarOp>::iterator citer=cigar.begin();
    int start=al.Position; 
    while(citer!=cigar.end()){
      switch(citer->Type){
      case 'M':
	{
	  for(int i=0;i<citer->Length;i++){
	    //	    cout<<start<<endl;
	    if(start>=refvector[al.RefID].RefLength){
	      cerr<<al.RefID<<"\t"<<refvector[al.RefID].RefName<<"\t"<<refvector[al.RefID].RefLength<<"\t"<<al.Name<<"\t"<<al.Position<<"\t"<<i<<mapper[refvector[al.RefID].RefName].size()<<"\t"<<start<<"\t"<<citer->Length<<citer->Type<<endl;
	      exit(-1);
	    }
	    mapper[refvector[al.RefID].RefName][start]++;
	    start++;
	  }
	  break;
	}
      case 'N':
	{
	  start+=citer->Length;
	  break;
	}

      case 'I':
	{
	  break;
	}
      case 'D':
	{
	  for(int i=0;i<citer->Length;i++){
	    if(start>=refvector[al.RefID].RefLength){
	      cerr<<al.RefID<<"\t"<<refvector[al.RefID].RefName<<"\t"<<refvector[al.RefID].RefLength<<"\t"<<al.Name<<"\t"<<al.Position<<"\t"<<i<<mapper[refvector[al.RefID].RefName].size()<<"\t"<<start<<"\t"<<citer->Length<<citer->Type<<endl;
	      exit(-1);
	    }
	    mapper[refvector[al.RefID].RefName][start]++;
	    start++;
	  }
	  break;
	}
      case 'S':
	{
	  break;
	}
      case 'H':
	{
	  break;
	}
      default:
	{
	  cerr<<"unknown cigar:"<<citer->Type<<endl;
	  exit(-1);
	}
      }
      ++citer;
    }
  }
  reader.Close();

  //Check coverage
  cerr <<"Calculate coverage and evenness" <<" ......"<<endl;
  checkCDSCoverage(mapper,gtf);
  checkGeneCoverage(mapper,gtf);
  checkCDSEvenness(mapper,gtf);
  checkGeneEvenness(mapper,gtf);

  //Output coverage
 
  cout<<"#name\tchr\tbegin\tend\tgene_coverage\tgene_depth\tgene_evenness\ttranscript_coverage\ttranscript_depth\ttranscript_evenness"<<endl;

  vector<GeneStructure>::iterator oiter=gtf.GeneStructureSequence.begin();
  while(oiter!=gtf.GeneStructureSequence.end()){
    cout<<oiter->name <<"\t"<<oiter->chr<<"\t"<<oiter->begin<<"\t"<<oiter->end;
    cout<< setprecision(4)
	<<"\t"<< oiter->geneCov <<"\t"<< oiter->geneDep <<"\t"<< oiter->geneEven 
	<<"\t"<< oiter->cdsCov  <<"\t"<< oiter->cdsDep  <<"\t"<< oiter->cdsEven 
	<<endl;
    ++oiter;
  }

  return 0;
}

void checkCDSCoverage(map<string,vector<unsigned int> > &mapper, GTF &gtf){
  vector<GeneStructure>::iterator giter=gtf.GeneStructureSequence.begin();
  while(giter!=gtf.GeneStructureSequence.end()){
    int tlength=0;
    int slength=0;
    int alength=0;
    vector<GeneItem>::iterator iiter=giter->itemseq.begin();
    while(iiter!=giter->itemseq.end()){
      tlength+=(iiter->end-iiter->begin+1);
      for(int i=iiter->begin;i<=iiter->end;++i){
	alength+=mapper[iiter->chr][i-1];
	if(mapper[iiter->chr][i-1]>0){
	  slength++;
	}
      }
      ++iiter;
    }
    giter->cdsDep=(double)alength/(double)tlength;
    giter->cdsCov=(double)slength/(double)tlength;
    ++giter;
  }
}

void checkCDSEvenness(map<string,vector<unsigned int> > &mapper, GTF &gtf){
  vector<GeneStructure>::iterator giter=gtf.GeneStructureSequence.begin();
  while(giter!=gtf.GeneStructureSequence.end()){
    int tlength=0;
    vector<GeneItem>::iterator iiter=giter->itemseq.begin();
    int even=0;
    while(iiter!=giter->itemseq.end()){
      tlength+=(iiter->end-iiter->begin);
      for(int i=iiter->begin+1;i<=iiter->end;i++){
	if( (mapper[iiter->chr][i-1]>0 && mapper[iiter->chr][i-2]>0)
	    ||  (mapper[iiter->chr][i-1]==0 && mapper[iiter->chr][i-2]==0)){
	  even++;
	}
      }
      ++iiter;
    }
    giter->cdsEven=(double)even/(double)tlength;
    ++giter;
  }
}


void checkGeneEvenness(map<string,vector<unsigned int> > &mapper, GTF &gtf){
  vector<GeneStructure>::iterator giter=gtf.GeneStructureSequence.begin();
  while(giter!=gtf.GeneStructureSequence.end()){
    int tlength=giter->end - giter->begin;
    int even=0;
    for(int i=giter->begin+1;i<=giter->end;i++){
      if( (mapper[giter->chr][i-1]>0 && mapper[giter->chr][i-2]>0)
	  ||  (mapper[giter->chr][i-1]==0 && mapper[giter->chr][i-2]==0)){
	even++;
      }
    }
    giter->geneEven=(double)even/(double)tlength;
    ++giter;
  }
}

void checkGeneCoverage(map<string,vector<unsigned int> > &mapper, GTF &gtf){
  vector<GeneStructure>::iterator giter=gtf.GeneStructureSequence.begin();
  while(giter!=gtf.GeneStructureSequence.end()){
    int tlength=giter->end - giter->begin + 1;
    int slength=0;
    int alength=0;
    for(int i=giter->begin;i<=giter->end;i++){
      //      if(i-1>=mapper[giter->chr].size()){
      //	cerr << i <<"\t**"<<mapper[giter->chr].size()<<endl;
      //	exit(-1);
      //      }  
    alength+=mapper[giter->chr][i-1];
      if(mapper[giter->chr][i-1]>0){
	slength++;
	alength+=mapper[giter->chr][i-1];
      }
    }
    giter->geneDep=(double)alength/(double)tlength;
    giter->geneCov=(double)slength/(double)tlength;
    ++giter;
  }
}


void initMapper(map<string,vector<unsigned int> > &mapper,string genomeFile){
  ifstream genomeStream;
  if(genomeFile.substr(genomeFile.size()-3)==".gz"){
    genomeFile = "zcat " + genomeFile +"|";
  }
  genomeStream.open(genomeFile.c_str());
  if(!genomeStream){
    cerr << "Unable to open reference genome file: "
         << genomeFile << endl;
    exit(-1);
  }
  string seq;
  string line;
  string chr="";
  while(getline(genomeStream, line)){
    if(line[0]=='>'){
      if(chr!=""){
	vector<unsigned int> tmp(seq.size(),0);
	mapper[chr]=tmp;
      }
      chr=line.substr(1,line.size()-1);
      seq="";
    }
    else{
      seq+=line;
    }
  }
  genomeStream.close();
  vector<unsigned int> tmp(seq.size(),0);
  mapper[chr]=tmp;
}

void PrintChrSize(map<string,vector<unsigned int> > mapper){
  map<string,vector<unsigned int> >::iterator iter=mapper.begin();
  while(iter!=mapper.end()){
    cerr<<iter->first << "\t" << iter->second.size()<<endl;
    iter++;
  }
}
