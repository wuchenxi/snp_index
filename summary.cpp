//Input files: emaize_5M.add.gz merbin_infw.txt snpfunc_5M_emaize 200kpos.txt wlst
//Usage: Put input files and summary.cpp at the same directory, then:
//             g++ -O3 summary.cpp -o summary
//             gzip -dc emaize_5M.add.gz | ./summary merbin_infw.txt snpfunc_5M_emaize 200kpos.txt wlst > out
//Then "out" is a space separated snp-index file and we can use various ML methods to analyze it.

#include<cstdio>
#include<cstdlib>
#include<cstring>
//#include <execinfo.h>
//#include <signal.h>
//#include <unistd.h>
//#include <cmath>


//void handler(int sig) {
//  void *array[10];
//  size_t size;

  //   get void*'s for all entries on the stack
//  size = backtrace(array, 10);

  //print out all the frames to stderr
//  fprintf(stderr, "Error: signal %d:\n", sig);
//  backtrace_symbols_fd(array, size, STDERR_FILENO);
//  exit(1);
//}

double weights[6]={1,1,1,1,1,1};
typedef struct snp_id{int chr; double pos;}sid;

int compare(sid a, sid b){
  if(a.chr<b.chr)return -1;
  else if(a.chr>b.chr)return 1;
  else if(a.pos<b.pos)return -1;
  else if(a.pos>b.pos)return 1;
  else return 0;
}

sid split_name(char* s){
  if(s==NULL||strlen(s)<5)return {1,0};
  s+=3;sid r;
  r.chr=s[0]-'0';
  if(s[1]=='.')s+=4;
  else{r.chr*=10;s+=5;}
  sscanf(s,"%lf",&(r.pos));
  return r;
}

sid int_end(char* s){
  sid r;
  strtok(s,"\t");
  char* a=strtok(NULL,"\t");
  strtok(NULL,"\t");
  char* b=strtok(NULL,"\t");
  if(a==NULL||b==NULL)return {0,0};
  sscanf(a,"%d",&(r.chr));
  sscanf(b,"%lf",&(r.pos));
  return r;
}

int func(char* s,sid* r){
  *r=split_name(strtok(s,","));
  strtok(NULL,",");
  long c;
  char* a=strtok(NULL,",");
  if(a==NULL)return 0;
  sscanf(a,"%ld",&c);
  return c;
}

sid exclude(char* s){
  //fputs(s,stdout);
  //printf("hi");
  sid r;
  char* a=strtok(s,"\t");
  char* b=strtok(NULL,"\t");
  if(a==NULL||b==NULL)return {11,0};
  sscanf(a,"%d",&(r.chr));
  sscanf(b,"%lf",&(r.pos));
  return r;
}

char buf1[102400];
char buf2[102400];
char buf3[102400];
char buf4[102400];
double cursum[6210];
int snp_l[6210];

int main(int argc, char* argv[]){
  //  signal(SIGSEGV, handler);   // install our handler
 
  char* intervals_file_name=argv[1];
  char* comments_file_name=argv[2];
  char* exclude_file_name=argv[3];
  FILE* weight_file=fopen(argv[4],"r");
  for(int i=0;i<6;i++)
    fscanf(weight_file,"%lf",weights+i);
  
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
  //int snpcount=0;
  while(!feof(stdin)){
    fgets(buf1,102400,stdin);
    if(strlen(buf1)<2)break;
    //if(rand()%1024!=23)continue;
    //fputs(buf1,stdout);
    char* snp_c=strtok(buf1,",");
    strtok(NULL,",");
    strtok(NULL,",");
    for(int i=0;i<6210;i++)
      snp_l[i]=(int)(strtok(NULL,",")[0]-'0');
    //fputs(snp_c,stdout);
    sid snp_name=split_name(snp_c);
    //printf("\n%d,%g\n",cur_func.chr,cur_func.pos);
      if(compare(snp_name,curex)==0)continue;
      while(compare(snp_name,curex)>0 && !feof(exclude_f)){
        fgets(buf4,102400,exclude_f);
        if(strlen(buf4)>2)curex=exclude(buf4);
	else break;}
    double snp_weight=0;
    while(compare(cur_func,snp_name)<0 && !feof(func_f)){
      fgets(buf3,102400,func_f);
      if(strlen(buf3)>2)cat=func(buf3, &cur_func);
      else break;
      //printf("%d,%g\n",cat,weights[cat-1]);
    }
    if(compare(cur_func,snp_name)==0 && cat<7 && cat>0)
      snp_weight=weights[cat-1];
    else
      continue;
    if(compare(snp_name, curint)<0){
      //snpcount++;
      for(int i=0;i<6210;i++)
	cursum[i]+=snp_weight*snp_l[i];
    }
    else{//if(snpcount){
      for(int i=0;i<6210-1;i++)printf("%g ",cursum[i]);
      printf("%g\n",cursum[6209]);//}
      while(compare(curint,snp_name)<=0 && !feof(interval_f)){
	fgets(buf2,102400,interval_f);
	if(strlen(buf2)>2)curint=int_end(buf2);
	else return 0;
	if(curint.chr==0)return 0;
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
