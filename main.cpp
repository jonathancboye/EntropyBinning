//File:main.cpp
//by: Jonathan Carpenter
//carpenter.102@wright.edu

#include<iostream>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<string>
#include<utility>
#include<cstdlib>
#include<algorithm>
#include<iomanip>
#include<list>
#include<map>
#include<typeinfo>
#include "gene.h"
using namespace std;

//Load the genes from a file into a vector
void loadGenes(vector<Gene>& genes, string inputFile);

//Return ture if the lhs gene has a lower IS value than rhs
bool greaterThan(const Gene lhs, const Gene rhs){return lhs > rhs;}

//########## The following functions create the required files for task 1 ##########
//Create file GenesEntropyOrder.txt
//with the following format: gene number, split value,
//distribution (number of elements in lower set, number of elements in higher set),
//size-weighted entropy
void createFile_GenesEntropyOrder(const vector<Gene>& genes);

//Create file DataEntropyMethod.txt
//In decending order of information gain we have the genes and their orginial values
void createFile_DataEntropyMethod(const vector<Gene>& genes, int k);

//Create file DiscretizedDataEntropyMethod.txt
//In descending order of information gain this file lists the discretized values for each of
//the k highest ranked genes. Where a and b represent lower and higher intervalues respectivly
void createFile_DiscretizedDataEntropyMethod(const vector<Gene>& genes, int k);

//########## These functions create the required files for task 2 ##########
//create file TopKgeneSetSFS.txt
//Print top k genes
void createFile_TopKgeneSetSFS(const vector<Gene>& genes);

//create file DataTopKgeneSetSFS.txt
//Print top k genes in SFS order with their original values
void createFile_DataTopKgeneSetSFS(const vector<Gene>& genes);

//create file DiscretizedDataTopKgeneSetSFS
//Print top k genes in SFS order with their discretized values
void createFile_DiscretizedDataTopKgeneSetSFS(const vector<Gene>& genes);

//Preform Sequential Forward Search to select genes 
//and return the a vector of selected genes 
vector<Gene> SFS(const vector<Gene>& genes, int k);

int main(int argc, char* argv[]){
  string inputFile = "";
  vector<Gene> genes; //Vector of genes
  vector<Gene> selectedGenes; //Vector of genes selected by sequential forward search
  int k = 10; //number of genes to print out

  if(argc != 3){
    cout << "Program needs two command line arguments:" << endl;
    cout << "<int> <file name>" << endl;
    return EXIT_FAILURE;
  }

  if(!(k = atoi(argv[1])) ){
    cout << "Please enter in a integer for the first command line argument.\n"
	 << "Must be greater than zero."
	 << endl;
    return EXIT_FAILURE;
  }

  inputFile = argv[2];
    

  //load genes from file
  loadGenes(genes, inputFile);
  
  //preform entropy based binning on all genes
  for(int i=0;i<genes.size(); ++i){
    genes[i].entropyBinning();
    //genes[i].printGene();//#################### remove comment for debugging
    genes[i].sortByBiopsy();
  }

  sort(genes.begin(), genes.end(), greaterThan); //sort genes according to IS value
  selectedGenes = SFS(genes, k); //select k genes by using sequential forward search

  createFile_GenesEntropyOrder(genes); //create file GenesEntropyOrder.txt

  createFile_DataEntropyMethod(genes, k); //create file DataEntropyMethod.txt
  
  createFile_DiscretizedDataEntropyMethod(genes,k); //create file DiscretizedDataEntropyMethod.txt

  createFile_TopKgeneSetSFS(selectedGenes); //create file TopKgeneSetSFS.txt

  createFile_DataTopKgeneSetSFS(selectedGenes); //create file DataTopKgeneSetSFS.txt
  
  createFile_DiscretizedDataTopKgeneSetSFS(selectedGenes); //create file DiscretizedDataTopKgeneSetSFS

  
  return EXIT_SUCCESS;
}

void loadGenes(vector<Gene>& genes, string inputFile){
  ifstream ifs;
  string tuple; //A tuple read in from the file
  string strBuilder = ""; //Builds up a string to be converted to a number
  char c; //biopsie result
  geneData data;
  string bioResult = ""; //negative or positve depending on the result of the biopsie
  int count;
  int biopsy = 0; //Number of current biopsy

  ifs.open(inputFile.c_str());
  if(!ifs.fail()){    
    try{
      while(getline(ifs,tuple)){
	count = 0;
	strBuilder = "";
	bioResult = tuple.substr(tuple.size() - 8);
	c = (bioResult == "positive") ? 'P' : 'N';
      
	for(int i=0;i<tuple.size();++i){
	  if(tuple[i] != ',' && tuple[i+1] != ' ')
	    strBuilder += tuple[i];
	  else{
	    data.value = atof(strBuilder.c_str());
	    data.classifier = c;
	    data.biopsy = biopsy;
 
	    if(count  >= genes.size()){
	      Gene g(count+1);
	      genes.push_back(g);
	    }
	    genes[count].insert(data);
	    i++;
	    count++;
	    strBuilder = "";
	  }
	  biopsy++;
	}
      }
    }
    catch(...){
      cout << "Bad input file: needs to be of the exact format of p1data.txt,\n"
	   << "but can have more genes and biopsies."
	   << endl;
      ifs.close();
      exit(-1);
    }
    ifs.close();
  }
  else{
    cout << "Could not open file: " << inputFile << endl;
    return;
  }
}

void createFile_GenesEntropyOrder(const vector<Gene>& genes){
  ofstream ofs;

  ofs.open("GenesEntropyOrder.txt");
  if(!ofs.fail()){
    int pad = 10;
    map<char,pair<int,int> > dist;

    for(int i=0;i<genes.size(); ++i){
      string str;
      dist = genes[i].getDistribution();
      
      ofs <<"g"
	  << setw(5) << setfill('0') << genes[i].getGeneNumber()
	  << ", "
	  << setw(5) << setfill(' ') << genes[i].getSplitValue() << ", "
	  << setw(5) << '('
	  << setw(2) << setfill('0') << dist['a'].first 
	  << setfill(' ') << ',' 
	  << setw(2)<< setfill('0') << dist['a'].second 
	  << setfill(' ') << ')'
	  << ','
	  << ' '<< '('
	  << setw(2) << setfill('0') << dist['b'].first 
	  << setfill(' ') << ',' 
	  << setw(2)<< setfill('0') << dist['b'].second 
	  << setfill(' ') << ')'
	  << ','
	  << setw(5) << ' '
	  << setprecision(9) << genes[i].getISValue()
    	  << endl;
    }
  }
  else
    cout << "Failed to open file: GenesEntropyOrder.txt" << endl;
  ofs.close();
}

//Helper function for formatting output
void writeGeneValues(ofstream& ofs, const vector<Gene>& genes, int k){
  int lastValue = genes[0].size();
  int lastGene = (k < genes.size()) ? k : genes.size();
  for(int i=0;i < lastGene;++i){
      
      ofs << setw(8) << setfill (' ') << "g"
	  << setw(5) << setfill('0')
	  << genes[i].getGeneNumber()
	  << ','
	  << setw(2) << setfill(' ') << ' ';
    }
    ofs << endl;
    for(int i=0;i < lastValue;i++){
      for(int j=0;j < lastGene;++j){
	ofs << setw(13) << setfill(' ')
	    << scientific 
	    << genes[j].getBiopsyValue(i)
	    << ','
	    << setw(2) << setfill(' ')<< ' ';
      }
      ofs << genes[lastGene - 1].getGeneData(i).classifier
          << endl;
    }
}

//Helper function for formatting output
void writeGeneDiscreteValues(ofstream& ofs, const vector<Gene>& genes, int k){
  int lastValue = genes[0].size();
  int lastGene = (k < genes.size()) ? k : genes.size();

  for(int i=0;i < lastGene;++i){
    ofs << setw(1) << setfill (' ') << "g"
	<< setw(5) << setfill('0')
	<< genes[i].getGeneNumber()
	<< ','
	<< setw(3) << setfill(' ') << ' ';
  }
  ofs << endl;
  for(int i=0;i< lastValue ;i++){
    for(int j=0;j < lastGene;++j){
      ofs << setw(6) << setfill(' ')
	  << genes[j].getDiscreteValue(i)
	  << ','
	  << setw(3) << setfill(' ')<< ' ';
    }
    ofs << genes[lastGene - 1].getGeneData(i).classifier
	<< endl;
  }
}

void createFile_DataEntropyMethod(const vector<Gene>& genes, int k){
  ofstream ofs;
  ofs.open("DataEntropyMethod.txt");
  if(!ofs.fail()){
    writeGeneValues(ofs, genes, k);
  }
  else
    cout << "Error: can't open file DataEntropyMethod.txt" << endl;
  ofs.close();	
}

void createFile_DiscretizedDataEntropyMethod(const vector<Gene>& genes, int k){
  ofstream ofs;
  ofs.open("DiscretizedDataEntropyMethod.txt");
  if(!ofs.fail()){
    writeGeneDiscreteValues(ofs, genes, k);
  }
  else
    cout << "Error: can't open file DiscretizedDataEntropyMethod.txt" << endl;
  ofs.close();	
} 

void createFile_TopKgeneSetSFS(const vector<Gene>& genes){
  ofstream ofs;
  ofs.open("TopKgeneSetSFS.txt");
  
  if(!ofs.fail()){
    for(int i=0;i < genes.size();++i){
      ofs << "g"
	  << setw(5) << setfill('0')
	  << genes[i].getGeneNumber()
	  << endl;
    }
  }
  else
    cout << "Error: can't open file TopKgeneSetSFS.txt" << endl;
  ofs.close();
}

void createFile_DataTopKgeneSetSFS(const vector<Gene>& genes){
  ofstream ofs;
  ofs.open("DataTopKgeneSetSFS.txt");
  if(!ofs.fail()){
    writeGeneValues(ofs, genes, genes.size());
  }
  else
    cout << "Error: can't open file DataTopKgeneSetSFS.txt" << endl;
}

void createFile_DiscretizedDataTopKgeneSetSFS(const vector<Gene>& genes){
    ofstream ofs;
  ofs.open("DiscretizedDataTopKgeneSetSFS.txt");
  if(!ofs.fail()){
    writeGeneDiscreteValues(ofs, genes, genes.size());
  }
  else
    cout << "Error: can't open file DiscretizedDataTopKgeneSetSFS.txt" << endl;
  ofs.close();	
}

//struct for holding information about how
//many times a pattern occurs with a negative
//or positive biopsy result
struct sfsData{
  string pattern; //sequence of a's and b's
  map<char,int> NPCount; //keeps track of how many times a pattern
                         //occured with the a N or P
  sfsData(){
    pattern = "";
    NPCount['N'] = 0;
    NPCount['P'] = 0;
  }
};

list<Gene>::iterator
selectGene(const vector<string>& patterns, list<Gene>& genes){
  
  list<Gene>::iterator  selectedGene; //pointer to the Gene with the highest consistency Rate
  list<Gene>::iterator it;
  double maxConsistencyRate = 0; //holds the value of the highest consistency Rate

  for(it = genes.begin(); it!=genes.end(); ++it){
    string pattern = ""; //holds the current pattern
    geneData gdata; 
    map<string,sfsData> m; //map that keeps track the number of unique patterns and how many N's and P's occur
    int predicted = 0; //The number of correct predictions
    double consistencyRate; // The consistencyRate for a gene
    int total = it -> size(); //Total number of values in a gene
    
    //generate map 
    for(int j=0;j < total;++j){
      gdata = it -> getGeneData(j);
      pattern = patterns[j] + gdata.discreteValue;

      if(m.count(pattern))
	m[pattern].NPCount[gdata.classifier] += 1;
      else{
	sfsData sdata;
	sdata.pattern = pattern;
	sdata.NPCount[gdata.classifier] += 1;
	m[pattern] = sdata;
      }
    }
    
    map<string,sfsData>::iterator mit;
    
    for(mit = m.begin();mit != m.end(); ++mit){
      sfsData d = mit -> second;
      predicted += max(d.NPCount['N'],
		       d.NPCount['P']);
    }

    //calculate consistency rate
    consistencyRate = (2*predicted - total)/(double)total;
    
    //keep track of the current max consistency rate in a list
    if(consistencyRate > maxConsistencyRate){
      maxConsistencyRate = consistencyRate;
      selectedGene = it;
    }
  }
  return selectedGene;
}
      
vector<Gene> SFS(const vector<Gene>& genes, int k){
  list<Gene> g(genes.begin(), genes.end());
  vector<string> patterns(g.front().size(), "");
  vector<Gene> selectedGenes;
  
  for(int j =0;j < k; ++j){
    list<Gene>::iterator it = selectGene(patterns,g);
    selectedGenes.push_back(*it); //add selected gene to the set of selected genes
    for(int i=0;i < patterns.size(); ++i){
      patterns[i] += it -> getDiscreteValue(i);
    }
    g.erase(it); //remove the selected gene from the linked list of genes
  }
  return selectedGenes;
}
 
