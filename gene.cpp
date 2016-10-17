#include "gene.h"
#include <iterator>
#include <algorithm>
#include <functional>
#include <iostream>
#include <cmath>
#include <map>
#include <utility>
#include <iomanip>
using namespace std;

bool lessValue(geneData g1, geneData g2){
  return g1.value < g2.value;
}

bool lessBiopsy(geneData g1, geneData g2){
  return g1.biopsy < g2.biopsy;
}

Gene::Gene(const Gene& g){
  geneNumber = g.geneNumber;
  values = g.values;
  distribution = g.distribution;
  splitValue = g.splitValue;
  ISValue = g.ISValue;
}

void Gene::sortByValues(){
  sort(values.begin(), values.end(), lessValue);
}

void Gene::sortByBiopsy(){
  sort(values.begin(), values.end(), lessBiopsy);
}

void Gene::insert(const geneData data){
  geneData v = data;  
  values.push_back(v);  
}

//Helper function for entropyBinning
//Caluculates the entropy of a set of gene values
double entropy(const vector<geneData>& v, int start, int end){
  int S = end - start;
  if(S == 0){return 0;}
  int p = 0, n = 0;
  double P,N;

  for(int i=start;i<end;++i){
    if(v[i].classifier == 'P'){p++;}
    else{n++;}
  }
  P = (double)p/S;
  N = (double)n/S;

  double lvalue = (P == 0 ) ? 0 : -(P * (log(P)/log(2)));
  double rvalue = (N == 0) ? 0 : -(N * (log(N)/log(2)));

  return lvalue + rvalue;
}

void Gene::entropyBinning(){
  sortByValues();

  int S = values.size();
  double entropyOf_S = entropy(values, 0, S);
  double maxGain = 0;

  for(int i=1;i<size();++i){
    int S1 = i - 0;
    int S2 = S - i;
    double entropyOf_S1 = entropy(values,0,i);
    double entropyOf_S2 = entropy(values,i,S);
    double IS = ((double)S1/S) * entropyOf_S1 + ((double)S2/S) * entropyOf_S2;
    double Gain = entropyOf_S - IS;
    
    if(Gain > maxGain){
      maxGain = Gain;
      splitValue = i;
      ISValue = IS;
    }
  }  

  discretizeValues();
  setDistribution();
}

//Prints out a gene's values,
//number, size
void Gene::printGene() const{
  map<char, pair<int,int> > dist;
  dist = getDistribution();
  cout << 'g' << geneNumber << endl;
  cout << "[ ";
  for(int i=0;i<values.size();++i){
    if(i%5 ==0){cout << endl;}
    cout << "(" << values[i].value << "," << values[i].classifier
	 << "," << values[i].discreteValue
	 << ")" << " ";
  }
  cout << "]" <<  endl;
  cout << "split point: " << splitValue << endl;
  cout  << '('
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
	<< endl;
}

double Gene::getBiopsyValue(int i) const{
  if(i >= 0 && i<values.size()){return values[i].value;}
  cout << "Value out of range";
  return -1;
}

char Gene::getDiscreteValue(int i) const{
  if(i >= 0 && i < values.size()){return values[i].discreteValue;}
  cout << "Value out of range";
  return 'e';
}

bool Gene::operator>(const Gene& rhs) const {
  return (ISValue < rhs.getISValue()) ? true : false;
}

void Gene::discretizeValues(){
  for(int i=0;i<values.size();++i){
    values[i].discreteValue = (i < splitValue) ? 'a':'b';
  }
}

void Gene::setDistribution(){
  map<char,pair<int,int> > m;
  m['a'] = make_pair(0,0); //pair(number of N's, number of P's)
  m['b'] = make_pair(0,0); //pair(number of N's, number of P's)
  for(int i=0;i < values.size();++i){
    if(values[i].classifier == 'N')
      m[values[i].discreteValue].first += 1;
    else
      m[values[i].discreteValue].second += 1;
  }
  distribution = m;
}
