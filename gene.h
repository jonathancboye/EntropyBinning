#ifndef GENE_H
#define GENE_H

#include <vector>
#include <utility>
#include <map>

struct geneData{
  double value; //value of gene
  int biopsy; //biopsy gene appears in
  char classifier; //either N or P depending on biopsy result

  //either a or b depending on if gene is in
  //upper or lower split of biopsy respectively
  char discreteValue;
                    
};

class Gene{
 public:
  
 Gene(int num): geneNumber(num){};
  
  Gene(const Gene& g);

  //Add gene to values
  void insert(const geneData data);
  
  //Print out  the elements of sorted_values
  void printGene() const;

  //Preform entropy based binning on the gene
  void entropyBinning();

  //Return the gene's number
  int getGeneNumber() const {return geneNumber;}

  //Return the split value
  int getSplitValue() const {return splitValue;}

  //Return the max IS value
  double getISValue() const {return ISValue;}

  //Return the Distribution of the split
  //Keys are 'a' and 'b'
  //Values are pair(number of N's, number of P's)
  std::map<char, std::pair<int, int> >getDistribution() const {return distribution;} 

  //Return gene's value of a particluar tuple
  //input is a integer that represents the tuple number
  double getBiopsyValue(int i) const;

  //Return a if gene's value is below the split point else b
  char getDiscreteValue(int i) const;

  //Return a gene's data in a particular biopsy
  geneData getGeneData(int i) const{return values[i];}

 

  //Return gene values
  std::vector<geneData> getValues(){return values;}

  //Return the size of values
  int size()const {return values.size();}
  
  //Return true if lhs has a higher info gain than rhs
  bool operator>(const Gene& rhs) const;

  //sort gene by value in ascending order
  void sortByValues();

  //sort gene by biopsy in ascending order
  void sortByBiopsy();
  
 private:

  int geneNumber; //Number of the gene

  std::vector<geneData> values; //vector of gene data
  
  //Keys are 'a' and 'b'
  //Values are pair(number of N's, number of P's)
  std::map<char,std::pair<int,int> > distribution; // Shows how may items are in the lower and higher orderset respectivley
 
  int splitValue; //split value determined by entropy based binning

  double ISValue; //size-weighted entropy, smaller equals more information gain aka higher rank  

  //set the value of discreteValue for each gene after entropyBinning is preformed
  void discretizeValues();

  void setDistribution();

};

#endif //GENE_H



