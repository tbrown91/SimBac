//Copyright (C) 2015 Thomas Brown, Xavier Didelot, Daniel J. Wilson, Nicola De Maio
//
// This file is part of SimBac.
//
// SimBac is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
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

#ifndef RECOMB_EVENT
#define RECOMB_EVENT
#include <math.h>

void split_ancestries(list<int> &starts_1, list<int> &ends_1, list<int> &starts_2, list<int> &ends_2, const int beg, const int end){
  //Split the ancestry of the recombinant node into two, extracting the recombinant interval
  list<int>::iterator itStart1=starts_1.begin(), itEnd1=ends_1.begin();
  if (beg < end){
    //Recombinant break does not wrap around the end of the genome
      while (itStart1!=starts_1.end()){
        //Find the interval containing the start and end points of the recombinant interval
        if (end < *itStart1) return; //No more intervals affected by recombinant interval
        else if (beg > *itEnd1){
          //Ancestral interval is not in recombinant break
          ++itStart1;
          ++itEnd1;
        }else if ((beg <= *itStart1) && (end > *itEnd1)){
          //Recombinant break takes entire ancestral interval
          starts_2.push_back(*itStart1);
          ends_2.push_back(*itEnd1);
          itStart1 = starts_1.erase(itStart1);
          itEnd1 = ends_1.erase(itEnd1);
        }else if ((beg <= *itStart1) && (end > *itStart1)){
          //Break falls at start of ancestral interval
          starts_2.push_back(*itStart1);
          ends_2.push_back(end-1);
          *itStart1 = end;
          return; //No more ancestral inervals to check
        }else if ((beg <= *itEnd1) && (end > *itEnd1)){
          //Break falls at end of ancestral interval
          starts_2.push_back(beg);
          ends_2.push_back(*itEnd1);
          *itEnd1 = beg - 1;
          ++itStart1;
          ++itEnd1;
        }else{
          //Break falls inside ancestral interval
          starts_2.push_back(beg);
          ends_2.push_back(end-1);
          ends_1.insert(itEnd1,beg-1);
          ++itStart1;
          starts_1.insert(itStart1,end);
          return;//No more ancestral inervals to check
        }
      }
  }else{
    //Recombinant break wraps around the end of the genome
    if (beg == end){
      //Recombinant interval takes entire genome
      starts_2 = starts_1;
      ends_2 = ends_1;
      starts_1.clear();
      ends_1.clear();
    }else{
      while (itStart1 != starts_1.end()){
        if ((end > *itEnd1) || (beg <= *itStart1)){
          //Break covers the ancestral interval
          starts_2.push_back(*itStart1);
          ends_2.push_back(*itEnd1);
          itStart1 = starts_1.erase(itStart1);
          itEnd1 = ends_1.erase(itEnd1);
        }else if ((end > *itStart1) && (beg <= *itEnd1)){
          //Break covers start and end of interval
          starts_2.push_back(*itStart1);
          ends_2.push_back(end-1);
          starts_2.push_back(beg);
          ends_2.push_back(*itEnd1);
          *itStart1 = end;
          *itEnd1 = beg-1;
          ++itStart1;
          ++itEnd1;
        }else if ((end > *itStart1) && (end < *itEnd1)){
          //Break falls across an interval at the start
          starts_2.push_back(*itStart1);
          ends_2.push_back(end-1);
          *itStart1 = end;
          ++itStart1;
          ++itEnd1;
        }else if ((beg <= *itEnd1) && (beg > *itStart1)){
          //Break falls across an interval at the end
          starts_2.push_back(beg);
          ends_2.push_back(*itEnd1);
          *itEnd1 = beg - 1;
          ++itStart1;
          ++itEnd1;
        }else{
          ++itEnd1;
          ++itStart1;
        }
      }
    }
  }
}

void choose_nonClonalRecomb(const vector<double> &prob, const int G, const list<int> &starts, const list<int> &ends, int &beg, int &end, const double noStop, const int totMaterial, const double recombRate){
  //Choose a recombination interval for the chosen lineage which is non-clonal
  int b=starts.size();
  double r_1 = (gsl_rng_uniform(rng)*recombRate);
  int index = 0;
  list<int>::const_iterator itStart = starts.begin(), itEnd = ends.begin();
  while (r_1 > prob[index]){
    r_1 -= prob[index];
    ++index;
    if (index == b) break;
  }
  if (index == b){
    //Choose a start point within the ancestral intervals
    beg = (int)floor(gsl_rng_uniform(rng)*(totMaterial-b));
    beg += *itStart+1;
    ++itStart;
    for (itEnd=ends.begin();itEnd!=ends.end();++itEnd){
      if (beg > *itEnd){
        beg = beg - *itEnd + *itStart;
      }else break;
      ++itStart;
    }
    double r_2 = gsl_rng_uniform(rng);
    int len = (int)ceil(log(1-r_2*(1-pow(noStop,G-1)))/log(noStop));
    end = (beg + len) % G;
  }else{
    //Start of recombination occurs at beginning of an ancestral interval
    if (index == 0){
      if ((starts.front() == 0) && (ends.back() == G-1)){
        //Start site is at beginning of genome inside an ancestral interval
        beg = 0;
        double r_2 = gsl_rng_uniform(rng);
        int len = (int)ceil(log(1-r_2*(1-pow(noStop,G-1)))/log(noStop));
        end = (beg + len) % G;
      }else{
        //Interval starts at first interval
        beg = starts.front();
        double r_2 = gsl_rng_uniform(rng);
        int len = (int)ceil(log(1-r_2*(1-pow(noStop,ends.back()-starts.front())))/log(noStop));
        end = (beg + len) % G;
      }
    }else{
      itStart = starts.begin();
      itEnd = ends.begin();
      advance(itEnd,index-1);
      advance(itStart,index);
      beg = *itStart;
      //Simulate recombinant break length via a truncated geometric distribution
      double r_2 = gsl_rng_uniform(rng);
      int len = (int)ceil(log(1-r_2*(1-pow(noStop,G+(*itEnd)-(*itStart))))/log(noStop));
      end = (beg + len) % G;
    }
  }
  itStart = starts.begin();
  itEnd = ends.begin();
  if ((starts.front() == 0) && (ends.back() == G-1)){
    //Ancestral material wraps around end of genome
    if (b>1){
      //Check if end of recombinant interval falls between ancestral intervals
      ++itStart;
      while (itStart!=starts.end()){
        if ((end > *itEnd) && (end <= *itStart)){
          end = *itEnd+1;
          break;
        }else if (end <= *itEnd){break;}
        ++itEnd;
        ++itStart;
      }
    }
  }else{
    //Check if end of recombinant interval falls between ancestral intervals
    if ((end <= starts.front()) || (end > ends.back())) end = ends.back()+1;
    else{
      ++itStart;
      while(itStart!=starts.end()){
        if ((end > *itEnd) && (end <= *itStart)){
          end = *itEnd+1;
          break;
        }else if (end <= *itEnd){break;}
        ++itEnd;
        ++itStart;
      }
    }
  }
}

void choose_clonalRecomb(const vector<double> &prob, const int G, const list<int> &starts, const list<int> &ends, int &beg, int &end, const double delta, const int totMaterial, const double recombRate){
  //Choose a recombination interval for the chosen lineage which is non-clonal
  int b=starts.size();
  double r_1 = (gsl_rng_uniform(rng)*recombRate);
  int index = 0;
  list<int>::const_iterator itStart = starts.begin(), itEnd = ends.begin();
  while (r_1 >= prob[index]){
    r_1 -= prob[index];
    ++index;
    if (index == b) break;
  }
  if (index == b){
    //Choose a start point within the ancestral intervals
    beg = (int)floor(gsl_rng_uniform(rng)*(totMaterial-b));
    beg += starts.front()+1;
    itStart = starts.begin();
    ++itStart;
    for (itEnd=ends.begin();itEnd!=ends.end();++itEnd){
      if (beg > *itEnd){
        beg = beg - *itEnd + *itStart;
      }else break;
      ++itStart;
    }
  }else{
    advance(itStart,index);
    beg = *itStart;
  }
  int len = gsl_ran_geometric(rng,1.0/delta);
  if (len > G){
    len = G;
  }
  end = (beg + len) % G;
  itStart = starts.begin();
  itEnd = ends.begin();
  if ((starts.front() == 0) && (ends.back() == G-1)){
    //Ancestral material wraps around end of genome
    if (b>1){
      //Check if end of recombinant interval falls between ancestral intervals
      ++itStart;
      while (itStart != starts.end()){
        if ((end > *itEnd) && (end <= *itStart)){
          end = *itEnd+1;
          break;
        }else if (end <= *itEnd){break;}
        ++itEnd;
        ++itStart;
      }
    }
  }else{
    //Check if end of recombinant interval falls between ancestral intervals
    if ((end <= starts.front()) || (end > ends.back())) end = ends.back()+1;
    else{
      ++itStart;
      while (itStart != starts.end()){
        if ((end > *itEnd) && (end <= *itStart)){
          end = *itEnd+1;
          break;
        }else if (end <= *itEnd){break;}
        ++itEnd;
        ++itStart;
      }
    }
  }
}



#endif
