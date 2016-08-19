#ifndef TRACK_RECOMBMATERIAL
#define TRACK_RECOMBMATERIAL

void find_taxa(int r,int n,vector<int> &poss_taxa,vector<vector<int> > s,vector<vector<bool> > mat,vector<bool> clonal){
  // Check if current node contains the relevant interval
  if (r<n){
    poss_taxa[r]=1;
  }else{
    if (s[r][1]==-1){
      find_taxa(s[r][0],n,poss_taxa,s,mat,clonal);
    }else{
      if ((clonal[s[r][0]])&&(clonal[s[r][1]])){
        find_taxa(s[r][0],n,poss_taxa,s,mat,clonal);
        find_taxa(s[r][1],n,poss_taxa,s,mat,clonal);
      }else if (clonal[s[r][1]]){
        find_taxa(s[r][1],n,poss_taxa,s,mat,clonal);
      }else if (clonal[s[r][0]]){
        find_taxa(s[r][0],n,poss_taxa,s,mat,clonal);
      }else{
        find_taxa(s[r][0],n,poss_taxa,s,mat,clonal);
        find_taxa(s[r][1],n,poss_taxa,s,mat,clonal);
      }
    }
  }
}

void find_originTaxa(int r,int n,vector<int> &poss_taxa,vector<vector<int> > s,vector<vector<bool> > mat,vector<bool> clonal){
  if (s[r][3]==-1){//Recombinant edge coalesces
    if (clonal[s[r][2]]){
      find_taxa(s[r][2],n,poss_taxa,s,mat,clonal);
    }else{
      find_originTaxa(s[r][2],n,poss_taxa,s,mat,clonal);
    }
  }else{//Recombinant edge undergoes another recombination event
    find_originTaxa(s[r][2],n,poss_taxa,s,mat,clonal);//Follow donor edge
  }
}

#endif
