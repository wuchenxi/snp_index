#include <cstdio>
#include <cstdlib>
#include <cstring>

const int chr=7;

char buf[100000];
char buf2[100000];
double avg[6210];

int main(int argc, char* argv[]){
  FILE* idx=fopen("merbin_infw.txt","rb");
  FILE* in=fopen("bin7.geno","rb");
  FILE* out=fopen("sum7.geno","wb");
  fgets(buf,100000,in);
  for(int i=0;i<30350;i++){
    int j; int c; double l; double r;
    fscanf(idx, "%d", &j);
    fscanf(idx, "%d", &c);
    fscanf(idx, "%lf", &l);
    fscanf(idx, "%lf", &r);
    if(c==chr){
    printf("%d: (%lf, %lf)\n", c, l, r);
    for(int i=0;i<6210;i++)avg[i]=0;
    int len=0;
    int id;
    while(true){
      strcpy(buf2,buf);
	char* p=strtok(buf2,"\t");
	sscanf(p,"%d",&id);
	//printf("%d\n",id);
	if(id>=l&&id<=r){
	  len++;
	  for(int j=0;j<6210;j++){
	    p=strtok(NULL,"\t");avg[j]+=(int)(p[0]-'0');
	  }
	}
	if(id>=r)break;
      if(!fgets(buf,100000,in))break; 
    }
    if(len)
      {printf("%d, %d\n",len,id);
	for(int j=0;j<6209;j++)fprintf(out,"%g\t",avg[j]/len);
	fprintf(out,"%g\n",avg[6209]/len);}
    }
  }
  return 0;
}
