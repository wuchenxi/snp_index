#include <cstdio>
#include <cstdlib>
#include <cstring>

const int chr=3;

char buf[100000];
double avg[6210];

int main(int argc, char* argv[]){
  FILE* idx=fopen("merbin_infw.txt","rb");
  FILE* in=fopen("bin3.geno","rb");
  FILE* out=fopen("sum3.geno","wb");
  for(int i=0;i<30350;i++){
    int j; int c; double l; double r;
    fscanf(idx, "%d", &j);
    fscanf(idx, "%d", &c);
    fscanf(idx, "%lf", &l);
    fscanf(idx, "%lf", &r);
    if(c==chr){
    printf("%d: (%g, %g)\n", c, l, r);
    for(int i=0;i<6210;i++)avg[i]=0;
    int len=0;
    while(true){
        if(!fgets(buf,100000,in))break;
	int id; 
	char* p=strtok(buf,"\t");
	sscanf(p,"%d",&id);
	if(id>=l&&id<=r){
	  len++;
	  for(int j=0;j<6210;j++){
	    p=strtok(NULL,"\t");avg[j]+=(int)(p[0]-'0');
	  }
	}
    }
    for(int j=0;j<6209;j++)fprintf(out,"%g\t",avg[j]/len);
    fprintf(out,"%g\n",avg[6209]/len);
    }
  }
  return 0;
}
