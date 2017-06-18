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

#include "arg.h"
#include "rng.h"

static const char * help=
    "\
    Usage: SimBac [OPTIONS]\n\
    \n\
    Options:\n\
    -N NUM   Sets the number of isolates (default is 100)\n\
    -T NUM   Sets the value of theta, between 0 and 1 (default is 0.01))\n\
    -m NUM   Sets the minimum probability of mutation in an interval of external recombination between 0 & 1 (default is 0)\n\
    -M NUM   Sets the maximum probability of mutation in an interval of external recombination between 0 & 1 (default is 0)\n\
    -R NUM   Sets the value of R, the site-specific internal recombination rate (default is 0.01)\n\
    -r NUM   Sets the rate of R external, the site-specific rate of external recombination (default is 0)\n\
    -D NUM   Sets the value of delta (default is 500)\n\
    -e NUM   Sets the average length of external recombinant interval (default is 500)\n\
    -B NUM,...,NUM Sets the number and length of the fragments\n\
             (default is 10000)\n\
    -G NUM   Sets the gap between each fragment(default is 0)\n\
    -s NUM   Use given seed to initiate random number generator\n\
    -o FILE  Export data to given file\n\
    -c FILE  Export clonal genealogy to given file\n\
    -l FILE  Export local trees to given file\n\
    -b FILE  Write log file of internal recombinant break interval locations\n\
    -f FILE  Write log file of external recombinant break interval locations\n\
    -g FILE  Write log file of recombinant break interval locations and relevant taxa (Use only recommended for small ARGs)\n\
    -d FILE  Export DOT graph to given file\n\
    -a       Include ancestral material in the DOT graph\n\
    ";

int main(int argc, char *argv[]) {

    //Initialization
    int n=100;
    double theta_site = 0.01;//Site-specific rate of mutation
    double theta_extMin=0.5;
    double theta_extMax=1.0;
    double rho_site = 0.01;//Site-specific rate of internal recombination
    double rho_extSite = 0.0;//Site-specific rate of external recombination
    double delta=500.0;
    double delta_ext=500.0;
    int seed=-1;
    bool am=false;
    vector<int> blocks;
    blocks.push_back(0);
    blocks.push_back(10000);
    vector<int> gaps;
    int c;
    string cfile="";//File to which the clonal genealogy is exported
    string lfile="";//File to which the local trees are exported
    string dfile="";//File to which the DOT graph is exported
    string ofile="";//File to which the data is exported
    string bfile="";//Log file of internal recombinant intervals
    string efile="";//Log file of external recombinant intervals
    string ffile="";//Log file of origins of recombinant intervals
    double p1,p2;
    char * pch;
    if (argc == 1){
      cout<<help<<endl;
      return 1;
    }
    while ((c = getopt (argc, argv, "ahN:T:m:M:R:r:D:e:s:B:G:c:l:d:o:b:f:g")) != -1)
    switch (c)
    {
        case('N'):n=atoi(optarg);break;
        case('m'):theta_extMin=atof(optarg);break;
        case('M'):theta_extMax=atof(optarg);break;
        case('T'):theta_site=atof(optarg);break;
        case('R'):rho_site=atof(optarg);break;
        case('r'):rho_extSite=atof(optarg);break;
        case('D'):delta=atof(optarg);break;
        case('e'):delta_ext=atof(optarg);break;
        case('s'):seed=atoi(optarg);break;
        case('B'):blocks=Arg::makeBlocks(optarg);break;
        case('G'):gaps=Arg::makeGaps(optarg);break;
        case('h'):cout<<help<<endl;return 1;break;
        case('a'):am=true;break;
        case('c'):cfile=optarg;break;
        case('l'):lfile=optarg;break;
        case('d'):dfile=optarg;break;
        case('o'):ofile=optarg;break;
        case('b'):bfile=optarg;break;
        case('f'):efile=optarg;break;
        case('g'):ffile=optarg;break;
        case '?':cout<<"Wrong arguments: did not recognise "<<c<<" "<<optarg<<endl<<help<<endl;return 1;
        default:abort();
    }


    if (seed==-1) makerng(); else {rng=gsl_rng_alloc(gsl_rng_default);gsl_rng_set(rng,seed);};
    if (gaps.size() == 0){
      for (size_t i=0;i<blocks.size()-1;++i) gaps.push_back(0);
    }
    if (rho_extSite == 0) delta_ext = 1;//Avoid dividing by zero

    if (gaps.size() != blocks.size()-1) {cout << "Wrong number of gaps given (must be same as number of blocks)" << endl;return 1;}

    double theta = theta_site * blocks.back();

    if (argc-optind>0) {cout<<"Wrong arguments."<<endl<<help<<endl;return 1;}
    //Build the ARG
    Arg * arg=new Arg(n,rho_site,rho_extSite,delta,delta_ext,blocks,gaps);
    //Build the data and export it
    if (ofile.length()>0) {
    Data * data=arg->drawData(theta,theta_extMin,theta_extMax);
    ofstream dat;
    dat.open(ofile.data());
    data->output(&dat);
    dat.close();
    delete(data);}
    //Extract the clonal genealogy and export it
    if (cfile.length()>0) {
    string truth=arg->extractCG();
    ofstream tru;
    tru.open(cfile.data());
    tru<<truth<<endl;
    tru.close();}
    //Extract the local trees and export them
    if (lfile.length()>0) {
    ofstream lf;
    lf.open(lfile.data());
    arg->outputLOCAL(&lf);
    lf.close();}
    //Export to DOT format
    if (dfile.length()>0) {
    ofstream dot;
    dot.open(dfile.data());
    arg->outputDOT(&dot,am);
    dot.close();}
    //Write internal recombinant breaks to file
    if (bfile.length()>0) {
    ofstream ibreaks;
    ibreaks.open(bfile.data());
    arg->outputIBREAKS(&ibreaks);
    ibreaks.close();}
    //Write extenal recombinant breaks to file
    if (efile.length()>0) {
    ofstream ebreaks;
    ebreaks.open(efile.data());
    arg->outputEBREAKS(&ebreaks);
    ebreaks.close();}
    //Write recombinant breaks and origins to file
    if (ffile.length()>0) {
    ofstream breaks;
    breaks.open(ffile.data());
    arg->outputBREAKS(&breaks);
    breaks.close();}
    delete(arg);
    gsl_rng_free(rng);
}
