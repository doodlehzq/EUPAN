#ifndef GTF_CPP
#define GTF_CPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <unistd.h>

#include "GTF.h"

void GTFList::ReadGTFList(const GTFFileList &filelist){
  unsigned int i;
  GTF gtf;
  for(i=0;i<filelist.list.size();i++){
    gtf.clear();
    gtf.file=filelist.list[i];
    gtf.ReadGTF(filelist.list[i]);
    this->gtflist.push_back(gtf);
  }
}

void GTFFileList::CheckFileState(){
  unsigned int i;
  for(i=0;i<this->list.size();i++){
    if(-1==access(this->list[i].c_str(),R_OK)){
      cerr << "Error: GTF-formated file"
	   << this->list[i]
	   << "doesn't exist or cannot be read!" 
	   <<endl;
      exit(-1);
    }
  }
} 

void GTF::clear(){
  this->file="";
  this->GeneStructureSequence.clear();
}

void GTF::ReadGTF(const string &file){
  ifstream infile(file.c_str());
  if(!infile){
    cerr << "Error: Cannot open input GTF-formated file: "
	 << file << endl; 
    exit(-1);
  }
  string gene_id_pre="";
  string line;
  string tmp;
  GeneItem item;
  int initflag=1;
  GeneStructure gene;
  if(infile.peek()==EOF){
    size=0;
  }
  else{
    while(getline(infile,line)){
      istringstream iline(line);
      iline >> item.chr;      iline >> item.source;
      iline >> tmp;
      if(tmp == "start_codon"){
	item.type = START_CODON;
      }
      else if(tmp == "stop_codon"){
	item.type = STOP_CODON;
      }
      else if(tmp == "CDS"){
	item.type = CDS;
      }
      else if(tmp == "exon"){
	item.type = EXON;
      }
      else if(tmp == "gene"){
	item.type = GENE;
      }
      else if(tmp == "transcript"){
	continue;
      }
      else if(tmp == "five_prime_UTR" || tmp == "5UTR"){
	item.type = FIVEUTR;
      }
      else if(tmp == "three_prime_UTR" || tmp == "3UTR"){
	item.type = THREEUTR;
      }
      else{
	cerr << "Unknown type in " << file
	     << ": " << tmp <<endl;
      }
      iline >> item.begin;
      iline >> item.end;
      iline >> item.score;
      iline >> tmp;
      if(tmp == "+"){
	item.symbol = PLUS;
      }
      else if(tmp == "-"){
	item.symbol = MINUS;
      }
      else{
	cerr << "Unknown symbol(helix direction) in " << file
	     << ": " << tmp <<endl;
      }
      iline >> item.phase;
      item.gene_id = "";
      iline >> tmp;
      iline >> tmp;
      iline >> tmp;
      iline >> item.gene_id;
      item.gene_id = item.gene_id.substr(1,item.gene_id.size()-3);
      item.transcript_id = "";
      if(item.gene_id == gene_id_pre){
	gene.itemseq.push_back(item);
      }
      else{
	gene_id_pre=item.gene_id;
	if(initflag==1){
	  gene.itemseq.push_back(item);
	  initflag=0;
	}
	else{
	  gene.update();
	  GeneStructureSequence.push_back(gene);
	  gene.clear();
	  gene.itemseq.push_back(item);
	}
      }
    }
    gene.update();
    GeneStructureSequence.push_back(gene);
    gene.clear();
  }
  infile.close();
  this->size=GeneStructureSequence.size();
}


void GTF::ReadGFF(const string &file){
  ifstream infile(file.c_str());
  if(!infile){
    cerr << "Error: Cannot open input GFF-formated file: "
	 << file << endl; 
    exit(-1);
  }
  string gene_id_pre="";
  string line;
  string tmp;
  GeneItem item;
  int initflag=1;
  GeneStructure gene;
  if(infile.peek()==EOF){
    size=0;
  }
  else{
    while(getline(infile,line)){
      istringstream iline(line);
      iline >> item.chr;      iline >> item.source;
      iline >> tmp;
      if(tmp == "start_codon"){
	item.type = START_CODON;
      }
      else if(tmp == "stop_codon"){
	item.type = STOP_CODON;
      }
      else if(tmp == "CDS"){
	item.type = CDS;
      }
      else if(tmp == "exon"){
	item.type = EXON;
      }
      else if(tmp == "gene"){
	item.type = GENE;
	continue;
      }
      else if(tmp == "transcript"){
	continue;
      }

      else if(tmp == "five_prime_UTR" || tmp == "5UTR"){
	item.type = FIVEUTR;
      }
      else if(tmp == "three_prime_UTR" || tmp == "3UTR"){
	item.type = THREEUTR;
      }
      else{
	cerr << "Unknown type in " << file
	     << ": " << tmp <<endl;
      }
      iline >> item.begin;
      iline >> item.end;
      iline >> item.score;
      iline >> tmp;
      if(tmp == "+"){
	item.symbol = PLUS;
      }
      else if(tmp == "-"){
	item.symbol = MINUS;
      }
      else{
	cerr << "Unknown symbol(helix direction) in " << file
	     << ": " << tmp <<endl;
      }
      iline >> item.phase;
      item.gene_id = "";
      while(iline >> tmp){
	if(tmp.substr(0,7)=="Parent="){
	  item.gene_id = tmp.substr(7,item.gene_id.size()-7);
	  break;
	}
      }
      item.transcript_id = "";
      if(item.gene_id == gene_id_pre){
	gene.itemseq.push_back(item);
      }
      else{
	gene_id_pre=item.gene_id;
	if(initflag==1){
	  gene.itemseq.push_back(item);
	  initflag=0;
	}
	else{
	  gene.update();
	  GeneStructureSequence.push_back(gene);
	  gene.clear();
	  gene.itemseq.push_back(item);
	}
      }
    }
    gene.update();
    GeneStructureSequence.push_back(gene);
    gene.clear();
  }
  infile.close();
  this->size=GeneStructureSequence.size();
}

void GeneStructure::clear(){
  itemseq.clear();
}

void GeneStructure::update(){
  this->chr=this->itemseq[0].chr;
  this->direction=this->itemseq[0].symbol;
  this->name=this->itemseq[0].gene_id;
  int size=this->itemseq.size();
  this->begin=this->itemseq[0].begin;
  this->end=this->itemseq[size-1].end;
  if(this->direction==PLUS && this->itemseq[0].type==START_CODON && this->itemseq[size-1].type==STOP_CODON){
    this->complete=1;
  }
  else if(this->direction==MINUS && this->itemseq[0].type==STOP_CODON && this->itemseq[size-1].type==START_CODON){
    this->complete=1;
  }
  else{
    this->complete=0;
  }
  if(size==3 && this->complete){
    single_exon=1;
  }
  else{
    single_exon=0;
  }
}

void GTF::WriteGTF(){
  ofstream outfile(file.c_str());
  unsigned int i=0;
  while(i<GeneStructureSequence.size()){
    unsigned int j=0;
    while(j<GeneStructureSequence[i].itemseq.size()){
      outfile << GeneStructureSequence[i].itemseq[j].chr           << "\t"
	      << GeneStructureSequence[i].itemseq[j].source        << "\t";
      if(GeneStructureSequence[i].itemseq[j].type == START_CODON){
	outfile << "start_codon\t";
      }
      else if(GeneStructureSequence[i].itemseq[j].type == STOP_CODON){
	outfile << "stop_codon\t";
      }
      else{
	outfile << "CDS\t";
      }
      outfile << GeneStructureSequence[i].itemseq[j].begin         << "\t"
	      << GeneStructureSequence[i].itemseq[j].end           << "\t"
	      << GeneStructureSequence[i].itemseq[j].score         << "\t";
      if(GeneStructureSequence[i].itemseq[j].symbol == PLUS){
	outfile << "+\t";
      }
      else{
	outfile << "-\t";
      }
      outfile << GeneStructureSequence[i].itemseq[j].phase         << "\t"
	      << GeneStructureSequence[i].itemseq[j].gene_id
	      << GeneStructureSequence[i].itemseq[j].transcript_id << endl;
      j++;
    }
    i++;
  }
  outfile.close();
}

#endif
