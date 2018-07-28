//Input files: emaize_5M.add.gz merbin_infw.txt snpfunc_5M_emaize chr_name
//Usage: Put input files and summary.cpp at the same directory, then:
//             g++ -O3 summary_new.cpp -o summary
//             gzip -dc chr1_Ncii4M.num.gz | ./summary merbin_infw.txt chr1_anno.txt 1 > out
//Each interval will result in 19 rows in the file "out", corresponding to the summary of
//each type of SNP in said interval. The first entry of each row is the number of such SNPs,
//followed by the sums of such SNPs for each of the 6210 samples.
//Do it for 10 chromosomes then cat the results together :)

#include<cstdio>
#include<cstdlib>
#include<cstring>

const char* annos[]={"anno", "downstream_gene_variant", "3_prime_UTR_variant", "intron_variant", "missense_variant", "5_prime_UTR_variant", "synonymous_variant", "upstream_gene_variant", "splice_region_variant", "intergenic_variant", "splice_donor_variant", "stop_gained", "stop_lost", "start_lost", "splice_acceptor_variant", "stop_retained_variant", "non_coding_transcript_exon_variant", "coding_sequence_variant", "incomplete_terminal_codon_variant"};

char buf1[1024000];
char buf2[102400];
char buf3[102400];
double cursum[19][6211];
int snp_l[6210];
char* chr_name;

double int_end(char* s){
  double r;
  strtok(s,"\t");
  char* a=strtok(NULL,"\t");
  strtok(NULL,"\t");
  char* b=strtok(NULL,"\t");
  if(a==NULL||b==NULL)return -2;
  if(strcmp(a,chr_name))return -1;
  sscanf(b,"%lf",&r);
  return r;
}

double split_name(char* s){
  if(s==NULL||strlen(s)<5)return -1;
  int i=0;
  for(;i<strlen(s);i++)if(s[i]=='.'){s[i]=0;break;}
  if(strcmp(s+3,chr_name))return -1;
  double r;
  sscanf(s+i+3,"%lf",&r);
  return r;
}

int func(char* s,double* r){
  *r=split_name(strtok(s,"\t"));
  char* f=strtok(NULL,"\t");
  if(f[strlen(f)-1]=='\n')f[strlen(f)-1]=0;
  for(int i=0;i<19;i++)
    if(strcmp(f,annos[i])==0)
      return i;
}


int main(int argc, char* argv[]){
  char* intervals_file_name=argv[1];
  char* comments_file_name=argv[2];
  chr_name=argv[3];
  
  fgets(buf1,1024000,stdin);
  FILE* interval_f=fopen(intervals_file_name,"r");
  fgets(buf2,102400,interval_f);
  fgets(buf2,102400,interval_f);
  double curint=int_end(buf2);
  FILE* func_f=fopen(comments_file_name,"r");
  fgets(buf3,102400,func_f);
  fgets(buf3,102400,func_f);
  double cur_func;
  int cat=func(buf3, &cur_func);
  for(int i0=0;i0<19;i0++){for(int i=0;i<6211;i++)cursum[i0][i]=0;}
  double sumweight=0;
  while(!feof(stdin)){
    fgets(buf1,102400,stdin);
    if(strlen(buf1)<2)break;
    char* snp_c=strtok(buf1,"\t");
    //printf("%s\n",snp_c);
    for(int i=0;i<6210;i++)
      snp_l[i]=(int)(strtok(NULL,"\t")[0]-'0');
    double snp_name=split_name(snp_c);
    if(snp_name<0){continue;}
    //printf("%g\n",snp_name);
    
    while(cur_func<snp_name && !feof(func_f)){
      fgets(buf3,102400,func_f);
      if(strlen(buf3)>2)
	cat=func(buf3, &cur_func);
      else break;
    }
  
    if(cur_func==snp_name && cat<19 && cat>=0)
      ;
    else
      continue;
    if(snp_name<curint){
      sumweight+=1;
      cursum[cat][0]+=1;
      for(int i=0;i<6210;i++)
	cursum[cat][i+1]+=snp_l[i];
    }
     else{
       if(sumweight){
	for(int i0=0;i0<19;i0++){
	  for(int i=0;i<6210;i++)
	    printf("%g ",cursum[i0][i]);
          printf("%g\n",cursum[i0][6210]);}
      }
      while(curint<=snp_name && !feof(interval_f)){
	fgets(buf2,102400,interval_f);
	if(strlen(buf2)>2)curint=int_end(buf2);
	else return 0;
	if(curint==-2)return 0;
      }
      if(snp_name<curint){
	cursum[cat][0]++;
	for(int i=0;i<6210;i++)
	  cursum[cat][i+1]=snp_l[i];
	sumweight++;
      }
      else
	break;
	}
    }  
  return 0;
}
