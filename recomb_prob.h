//Copyright (C) 2015 Thomas Brown, Xavier Didelot, Daniel J. Wilson, Nicola De Maio
//
// This file is part of SimBac.
//
// SimBac is free software: you can redistribute it and/or modify
// it under the terms of the

// General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SimBac is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SimBac.  If not, see <http://www.gnu.org/licenses/>.

#ifndef RECOMB_PROB
#define RECOMB_PROB
#include <math.h>

void calc_nonClonalRecomb(const int G, const double delta, vector<double> &prob, double &recombRate, const list<int> &starts, const list<int> &ends, const double noStop, const double siteRecomb, int &totMaterial){
  //Total amount of ancestral material in the node
  list<int>::const_iterator itStart = starts.begin(), itEnd = ends.begin();
  totMaterial = 0;
  int b = starts.size();
  //Reset recombination rate and probability of recombination along the genome
  prob.clear();
  recombRate = 0.0;
  if (starts.size() == 0) return;
  while (itStart != starts.end()){
      totMaterial += (*itEnd - *itStart) + 1;
      ++itStart;
      ++itEnd;
  }
  if (totMaterial <= 1) return;
  if ((starts.front() == 0) && (ends.back() == G-1)){
    //Interval wraps around the end of the genome
    if (b > 1){
      prob.push_back(0.0);
      itStart = starts.begin();
      itEnd = ends.begin();
      ++itStart;
      while (itStart != starts.end()){
        prob.push_back(siteRecomb*delta*(1-pow(noStop,*itStart-*itEnd))*(1-pow(noStop,G-*itStart+*itEnd)));
        recombRate += prob.back();
        ++itStart;
        ++itEnd;
      }
      recombRate += siteRecomb*(totMaterial-(b-1))*(1-pow(noStop,G-1));
      prob[0] = siteRecomb*(1-pow(noStop,G-1));
    }else{
      //Interval is the entire genome
      recombRate = siteRecomb*totMaterial*(1-pow(noStop,G-1));
      prob.push_back(siteRecomb*(1-pow(noStop,G-1)));
    }
  }else{
    //Intervals are contained within the genome
    prob.push_back(siteRecomb*delta*(1-pow(noStop,G+starts.front()-ends.back()))*(1-pow(noStop,ends.back()-starts.front())));
    recombRate += prob[0];
    itStart = starts.begin();
    itEnd = ends.begin();
    ++itStart;
    while (itStart != starts.end()){
      prob.push_back(siteRecomb*delta*(1-pow(noStop,*itStart-*itEnd))*(1-pow(noStop,G-*itStart+*itEnd)));
      recombRate += prob.back();
      ++itStart;
      ++itEnd;
    }
    recombRate += siteRecomb*(totMaterial-b)*(1-pow(noStop,G-1));
  }
}

void calc_clonalRecomb(const int G, const double delta, vector<double> &prob, double &recombRate, const list<int> &starts, const list<int> &ends, const double noStop, const double siteRecomb, int &totMaterial){
  //Total amount of ancestral material in the node
  list<int>::const_iterator itStart = starts.begin(), itEnd = ends.begin();
  totMaterial = 0;
  int b = starts.size();
  //Reset recombination rate and probability of recombination along the genome
  prob.clear();
  recombRate = 0;
  if (starts.size() == 0) return;
  while (itStart != starts.end()){
      totMaterial += (*itEnd - *itStart) + 1;
      ++itStart;
      ++itEnd;
  }
  if (totMaterial <= 1) return;
  if ((starts.front() == 0) && (ends.back() == G-1)){
    //Interval wraps around the end of the genome
    if (b > 1){
      prob.push_back(0.0);
      itStart = starts.begin();
      itEnd = ends.begin();
      ++itStart;
      while (itStart != starts.end()){
        prob.push_back(siteRecomb*delta*(1-pow(noStop,*itStart-*itEnd)));
        recombRate += prob.back();
        ++itStart;
        ++itEnd;
      }
      recombRate += siteRecomb*(totMaterial-(b-1));
      prob[0] = siteRecomb;
    }else{
      //Interval is the entire genome
      recombRate = siteRecomb*totMaterial;
      prob.push_back(siteRecomb);
    }
  }else{
    //Intervals are contained within the genome
    prob.push_back(siteRecomb*delta*(1-pow(noStop,G+starts.front()-ends.back())));
    recombRate += prob[0];
    itStart = starts.begin();
    itEnd = ends.begin();
    ++itStart;
    while (itStart != starts.end()){
      prob.push_back(siteRecomb*delta*(1-pow(noStop,*itStart-*itEnd)));
      recombRate += prob.back();
      ++itStart;
      ++itEnd;
    }
    recombRate += siteRecomb*(totMaterial-b);
  }
}
#endif
