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

#ifndef MRCA_CHECKS
#define MRCA_CHECKS

void removeAncMat(const int start, const int end, list<int> &starts, list<int> &ends){
  //Reove MRCA material from chosen nodes
  list<int>::iterator itStart=starts.begin(), itEnd = ends.begin();
  while (itEnd!=ends.end()){
    if (end < *itStart) return; //No more ancestral material to check
    else if (start <= *itStart){
      if (end >= *itEnd){
        //Entire interval to be removed
        itStart = starts.erase(itStart);
        itEnd = ends.erase(itEnd);
      }else{
        //Start of interval to be removed
        *itStart = end + 1;
        return; //No more ancestral material to check
      }
    }else if (start <= *itEnd){
      if (end >= *itEnd){
        //End of interval to be removed
        *itEnd = start - 1;
        ++itStart;
        ++itEnd;
      }else{
        //Middle of interval to be removed
        ends.insert(itEnd,start - 1);
        ++itStart;
        starts.insert(itStart,end + 1);
        return;//No more material to check;
      }
    }else{
        ++itStart;
        ++itEnd;
    }
  }
}

void modifyMRCA(Arg::MRCA &M, const int start, const int end, int &MRCA_check){
//update the list of MRCA intervals by subtracting an interval overlapping in the two coalesceing lineages
  M.itStart=M.starts.begin();
  M.itEnd=M.ends.begin();
  M.itValue=M.values.begin();
	if (M.itStart==M.starts.end()){
		cout << "Error : this part should be in the MRCA vector, but it isn't (1) " << *M.itStart << " " << *M.starts.end() << endl;
	 	return;
	}
	while (*(M.itEnd) < start) {
		++(M.itStart);
		++(M.itEnd);
		++(M.itValue);
		if (M.itStart==M.starts.end()){
			cout << "Error : this part should be in the MRCA vector, but it isn't (2)" << *M.itStart << " " << *M.starts.end() << endl;
      return;
		}
	}
	if (start<*(M.itStart)){
		cout << "Error : this part should be in the MRCA vector, but it isn't (3)" << start << " " << *M.itStart << endl;
    return;
	}
	//first element in MRCA list that overlaps, but starting values do not coincide
	if ((M.itStart!=M.starts.end()) && (*M.itStart<start)){
    ++(M.itStart);
    M.starts.insert(M.itStart,start);
    --(M.itStart);
		M.ends.insert(M.itEnd,start-1);
		M.values.insert(M.itValue,*M.itValue);
	}
	//iteratively look at all intervals in the MRCA structure overlapping completely the given interval
	while ((M.itStart!=M.starts.end()) && (*M.itEnd<=end)) {
    --*(M.itValue);
    if (*(M.itValue) == 1) MRCA_check = 1;//Need to remove some ancestral material from new node
		++(M.itStart);
		++(M.itEnd);
		++(M.itValue);
	}
  // look at the final interval in the MRCA, if present, that overlaps the given interval but has end after the end of the given interval
  if ((M.itStart!=M.starts.end()) && (*M.itStart<=end)){
    ++(M.itStart);
    M.starts.insert(M.itStart,end+1);
    M.ends.insert(M.itEnd,end);
    M.values.insert(M.itValue,*M.itValue);
    --M.itValue;
    --*(M.itValue);
    if (*(M.itValue) == 1) MRCA_check = 1;//Need to remove some ancestral material from new node
  }
}

void merge_MRCA(Arg::MRCA &M){
  //Merge any consecutive intervals that have the same MRCA value to reduce the size of the struct for searching purposes
  M.itStart=M.starts.begin();
  M.itEnd=M.ends.begin();
  M.itValue=M.values.begin();
  while (M.itValue != (M.values).end()){
    int tempVal = *(M.itValue);
    ++M.itValue;
    ++M.itStart;
    if (*(M.itValue) == tempVal){
      M.itStart = (M.starts).erase(M.itStart);
      M.itEnd = (M.ends).erase(M.itEnd);
      ++(M.itEnd);
      M.itValue = (M.values).erase(M.itValue);
    }else ++M.itEnd;
  }
}

void update_MRCA(Arg::MRCA &M, const list<int> &starts_1, const list<int> &ends_1, const list<int> &starts_2, const list<int> &ends_2, int &MRCA_check){
  //Find overlapping regions of ancestral material in the two children and update the MRCA lists
  if ((starts_1.size() == 0) || (starts_2.size() == 0)) return; //No ancestral material in one or both nodes, therefore MRCA struct does not need to be updated
  list<int>::const_iterator itStart1 = starts_1.begin(), itEnd1 = ends_1.begin();
  list<int>::const_iterator itStart2 = starts_2.begin(), itEnd2 = ends_2.begin();
  int currentStart1=0, currentStart2=0;
  currentStart1=*itStart1, currentStart2=*itStart2;
  int interval_check = 0;
  while ((itStart1 != starts_1.end()) && (itStart2 != starts_2.end())){
    if ((currentStart1 < currentStart2) && (*itEnd1>=currentStart2)){ //overlap
      currentStart1=currentStart2;
    }else if ((currentStart2 < currentStart1) && (*itEnd2>=currentStart1)){ //overlap
      currentStart2=currentStart1;
    }else if ((currentStart1 < currentStart2) && (*itEnd1<currentStart2)){ //no overlap
      ++itStart1;
      ++itEnd1;
      if (itStart1!=starts_1.end()) currentStart1=*itStart1;
    }else if ((currentStart2 < currentStart1) && (*itEnd2<currentStart1)){ //no overlap
      ++itStart2;
      ++itEnd2;
      if (itStart2 != starts_2.end()) currentStart2=*itStart2;
    }else if ((currentStart1 == currentStart2) && (*itEnd1<*itEnd2)){ //overlap, same start
      modifyMRCA(M, currentStart1, *itEnd1, MRCA_check);
      interval_check = 1;//Potential to merge MRCA struct
      currentStart2=*itEnd1+1;
      ++itStart1;
      ++itEnd1;
      if (itStart1 != starts_1.end()) currentStart1=*itStart1;
    }else if ((currentStart1 == currentStart2) && (*itEnd2<*itEnd1)){ //overlap, same start
      modifyMRCA(M, currentStart1, *itEnd2, MRCA_check);
      interval_check = 1;//Potential to merge MRCA struct
      currentStart1=*itEnd2+1;
      ++itStart2;
      ++itEnd2;
      if (itStart2 != starts_2.end()) currentStart2=*itStart2;
    }else if ((currentStart1 == currentStart2) && (*itEnd2==*itEnd1)){ //overlap, same start, same end
      modifyMRCA(M, currentStart1, *itEnd1, MRCA_check);
      interval_check = 1;//Potential to merge MRCA struct
      ++itStart2;
      ++itEnd2;
      if (itStart2 != starts_2.end()) currentStart2=*itStart2;
      ++itStart1;
      ++itEnd1;
      if (itStart1 != starts_1.end()) currentStart1=*itStart1;
    }else{
      cout << "Error : this case should not happen in updating MRCA vector" << endl;
      return;
    }
  }
  if (interval_check == 1){
    //If the intervals of the MRCA struct have been altered, check if any can be merged together
    merge_MRCA(M);
  }
}
#endif
