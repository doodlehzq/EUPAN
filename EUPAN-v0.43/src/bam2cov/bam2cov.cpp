#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <api/BamReader.h>

using namespace std;

void PrintUsage(){
  cerr << "bam2cov <bam> " << endl;
}
void OutputSegments(vector<vector<unsigned int> > mapper,BamTools::RefVector rvector);

int main(int argc, char* argv[]){
  string of = "output";
  if(argc<2){
    PrintUsage();
    exit(-1);
  }
  string alignfile=argv[1];
  BamTools::BamReader reader;
  if(!reader.Open(alignfile)){
    cerr<<"Couldn't open bam file:" <<alignfile<<endl;
    exit(-1);
  }
  BamTools::RefVector refvector=reader.GetReferenceData();
  vector<vector<unsigned int> > mapper;
  vector<BamTools::RefData>::iterator iter=refvector.begin();
  while(iter!=refvector.end()){
    vector<unsigned int> tmp(iter->RefLength);
    mapper.push_back(tmp);
    //    cerr<<iter->RefName<<"\t"<<mapper.back().size()<<endl;
    ++iter;
  }



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
	      cerr<<al.RefID<<"\t"<<refvector[al.RefID].RefName<<"\t"<<refvector[al.RefID].RefLength<<"\t"<<al.Name<<"\t"<<al.Position<<"\t"<<i<<mapper[0].size()<<"\t"<<start<<"\t"<<citer->Length<<citer->Type<<endl;
	      exit(-1);
	    }
	    mapper[al.RefID][start]++;
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
	      cerr<<al.RefID<<"\t"<<refvector[al.RefID].RefName<<"\t"<<refvector[al.RefID].RefLength<<"\t"<<al.Name<<"\t"<<al.Position<<"\t"<<i<<mapper[0].size()<<"\t"<<start<<"\t"<<citer->Length<<citer->Type<<endl;
	      exit(-1);
	    }
	    mapper[al.RefID][start]++;
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

  OutputSegments(mapper,refvector);
  return 0;
}


void OutputSegments(vector<vector<unsigned int> > mapper,BamTools::RefVector rvector){

  for(int j=0;j<mapper.size();j++){
    int s=0;
    int flag=0;
    for(int i=0;i<mapper[j].size();i++){
      if(mapper[j][i]==0 && !flag){
	continue;
      }
      else if(mapper[j][i]>0 && !flag){
	s=i;
	flag=1;
	continue;
      }
      else if(mapper[j][i]>0 && flag){
	continue;
      }
      else{
	cout<<rvector[j].RefName<<"\t"<<s+1<<"\t"<<i<<endl;
	flag=0;
      }
    }
    if(flag){
      cout<<rvector[j].RefName<<"\t"<<s+1<<"\t"<<rvector[j].RefLength<<endl;
    }
  }
}



