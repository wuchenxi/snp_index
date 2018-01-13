#include<stdio.h>
#include<string.h>
#include<stdlib.h>
char buf[100000];
int buf2[6210];
int main(){
  int count=0;
  FILE* in=fopen("chr3_emaize.genoMat","r");
  FILE* out=fopen("bin3.geno","w");
  fgets(buf,100000,in);
  for(int i=0;1;i++){
    if(i%1000==0)printf("%d\n",i);
    if(!fgets(buf,100000,in))break;
    int j;
    char*p=strtok(buf,"\t");
    p=strtok(NULL,"\t");
    char maj=p[0];
    char min=p[2];
    p=strtok(NULL,"\t");
    p=strtok(NULL,"\t");
    int pos;
    sscanf(p,"%d", &pos);
    for(j=0;j<6210;j++){
      p=strtok(NULL,"\t");
      buf2[j]=0;
      if(p[0]==maj)buf2[j]--; else if(p[0]==min) buf2[j]++; else puts(p);
      if(p[1]==maj)buf2[j]--; else if(p[1]==min) buf2[j]++; else puts(p);
    }
    int sum=0;for(j=0;j<6210;j++)sum+=buf2[j];
    if(sum<=6210*2*0.9&&sum>=-6210*2*0.9){
      count++;
      fprintf(out,"%d\t",pos);
      for(j=0;j<6209;j++){
	if(buf2[j]==2)fputs("2\t",out);
        else if(buf2[j]==0)fputs("1\t",out);
        else fputs("0\t",out);}
      if(buf2[j]==2)fputs("2\n",out);
      else if(buf2[j]==0)fputs("1\n",out);
      else fputs("0\n",out);
    }
  }
  printf("%d\n",count);
  return 0;
}
