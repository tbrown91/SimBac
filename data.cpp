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

#include "data.h"

Data::Data(int n,vector<int> blocks) {
	this->n=n;
	this->L=blocks.back();
	this->blocks=blocks;
	data=vector<string>(n,string(L,'N'));
	poly=vector<bool>(L,true);
	begEnd=vector<bool>(L+1,false);
	for (unsigned int i=0;i<blocks.size();i++) begEnd[blocks[i]]=true;
}

Data::Data(string filename) {
    string line;
    ifstream file;
    file.open(filename.data());//Open file
    unsigned int which=0;
    while (1) {
        getline(file,line);//Read next line from file
        if (file.eof()) break;//Stop if end of file
        if (line.size()==0 || line[0]=='#' || line[0]=='=') continue;//Ignore empty lines, comments, and end of block lines
        if (line[0]=='>') {//Header line
            line.erase(0,1);
            istringstream iss(line);
            iss>>which;
            while (which>=data.size()) data.push_back("");
            if (which==0) {
                blocks.push_back(data[0].size());
            };
            continue;
        }
        //Sequence data line
        data[which].append(line);
    }
    file.close();//Close file
    n=data.size();
    L=data[0].size();
    blocks.push_back(L);
    for (unsigned int i=1;i<n;i++) if (data[0].size()!=data[i].size()) {
            cerr<<"Data is inconsistent: "<<data[0].size()<<"!="<<data[i].size()<<endl;
            break;
        }
    cout<<"Read input file with "<<n<<" isolates and "<<getB()<<" blocks for a total of "<<L<<" sites."<<endl;
    for (unsigned int i=0;i<n;i++) for (unsigned int j=0;j<L;j++) data[i][j]=convert(data[i][j]);
    poly=vector<bool>(L,true);
    for (unsigned int i=0;i<L;i++) makePoly(i);
	begEnd=vector<bool>(L+1,false);
	for (unsigned int i=0;i<blocks.size();i++) begEnd[blocks[i]]=true;
}

char Data::convert(char in) {
    switch (in) {
    case 'a':
    case 'A':
        return 0;
    case 't':
    case 'T':
        return 1;
    case 'c':
    case 'C':
        return 2;
    case 'g':
    case 'G':
        return 3;
    default:
        return 'N';
    }
}

char Data::convertBack(char in) {
    switch (in) {
    case 0:
        return 'A';
    case 1:
        return 'T';
    case 2:
        return 'C';
    case 3:
        return 'G';
    default:
        return 'N';
    }
}

void Data::output(ostream * out)
{
	for (int i=0;i<getB();++i)
	{
		for (unsigned int j=0;j<n;++j)
		{
			*out<<">"<<j<<endl;
			for (int k=blocks[i];k<blocks[i+1];++k) *out<<convertBack(data[j][k]);
			*out<<endl;
		}
		if (blocks.size() != 2){
			*out<<"="<<endl;
		}
	}
}

double Data::watterson()
{
	int p=0;
	for (unsigned int i=0;i<poly.size();++i) if (poly[i]) p++;//bad with missing data
	double s=0;
	for (unsigned int i=1;i<=n;++i) s+=1.0/i;
	return 1.0*p/s;
}

Data::~Data() {}
