//Input files: emaize_5M.add.gz merbin_infw.txt snpfunc_5M_emaize 200kpos.txt
//Usage: Put input files and summary.cpp at the same directory, then:
//             g++ -O3 summary.cpp -o summary
//             gzip -dc emaize_5M.add.gz | ./summary merbin_infw.txt snpfunc_5M_emaize 200kpos.txt 1 1 1 1 1 1 > out
//Then "out" is a space separated snp-index file and we can use various ML methods to analyze it.

#include<cstdio>
#include<cstdlib>
#include<cstring>

double weights[6]={1,1,1,1,1,1};
typedef struct snp_id{int chr; double pos;}sid;

int compare(sid a, sid b){
  if(a.chr<b.chr)return -1;
  else if(a.chr>b.chr)return 1;
  else if(a.pos<b.pos)return -1;
  else if(a.pos>b.pos)return 1;
  else return 0;
}

int s2n(char* s){
  int r=0;
  while(*s>='0'&&*s<='9'){r=r*10+*s-'0';s++;}
  return r;
}

sid split_name(char* s){
  s+=3;sid r;
  r.chr=s[0]-'0';
  if(s[1]=='.')s+=4;
  else{r.chr*=10;s+=5;}
  r.pos=(double)s2n(s);
  return r;
}

sid int_end(char* s){
  sid r;
  strtok(s,"\t");
  r.chr=s2n(strtok(NULL,"\t"));
  strtok(NULL,"\t");
  sscanf(strtok(NULL,"\t"),"%lf",&(r.pos));
  return r;
}

int func(char* s,sid* r){
  *r=split_name(strtok(s,","));
  strtok(NULL,",");
  return s2n(strtok(NULL,","));
}

sid exclude(char* s){
  sid r;
  char* a=strtok(s,"\t");
  char* b=strtok(NULL,"\t");
  sscanf(a,"%d",&(r.chr));
  sscanf(b,"%lf",&(r.pos));
  return r;
}

char buf1[102400];
char buf2[102400];
char buf3[102400];
char buf4[102400];
double cursum[6210];
int main(int argc, char* argv[]){
  char* intervals_file_name=argv[1];
  char* comments_file_name=argv[2];
  char* exclude_file_name=argv[3];
  for(int i=0;i<6;i++)
    sscanf(argv[4+i],"%lf",weights+i);
  
  fgets(buf1,102400,stdin);
  FILE* interval_f=fopen(intervals_file_name,"r");
  fgets(buf2,102400,interval_f);
  fgets(buf2,102400,interval_f);
  sid curint=int_end(buf2);
  FILE* func_f=fopen(comments_file_name,"r");
  fgets(buf3,102400,func_f);
  fgets(buf3,102400,func_f);
  sid cur_func;
  int cat=func(buf3, &cur_func);
  for(int i=0;i<6210;i++)cursum[i]=0;
  FILE* exclude_f=fopen(exclude_file_name,"r");
  fgets(buf4,102400,exclude_f);
  sid curex=exclude(buf4);
  while(!feof(stdin)){
    fgets(buf1,102400,stdin);
    //fputs(buf1,stdout);
    int snp_l[6210];
    char* snp_c=strtok(buf1,",");
    strtok(NULL,",");
    strtok(NULL,",");
    for(int i=0;i<6210;i++){
      snp_l[i]=s2n(strtok(NULL,","));
    }
    //fputs(snp_c,stdout);
    sid snp_name=split_name(snp_c);
    //printf("\n%d,%g\n",cur_func.chr,cur_func.pos);
      if(compare(snp_name,curex)==0)continue;
      while(compare(snp_name,curex)>0 && !feof(exclude_f)){
        fgets(buf4,102400,exclude_f);
        curex=exclude(buf4);
    }
    double snp_weight=0;
    while(compare(cur_func,snp_name)<0 && !feof(func_f)){
      fgets(buf3,102400,func_f);
      cat=func(buf3, &cur_func);
      //printf("%d,%g\n",cat,weights[cat-1]);
    }
    if(compare(cur_func,snp_name)==0)
      snp_weight=weights[cat-1];
    else
      continue;
    if(compare(snp_name, curint)<0){
      //printf("%g\n",snp_weight);
      for(int i=0;i<6210;i++)
	cursum[i]+=snp_weight*snp_l[i];
    }
    else{
      for(int i=0;i<6210-1;i++)printf("%g ",cursum[i]);
      printf("%g\n",cursum[6209]);
      while(compare(curint,snp_name)<=0){
	fgets(buf2,102400,interval_f);
	curint=int_end(buf2);
      }
      if(compare(snp_name,curint)<0)
	for(int i=0;i<6210;i++)
	  cursum[i]=snp_weight*snp_l[i];
      else
	break;
    }
  }  
  return 0;
}
