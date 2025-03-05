//CODE TO FIND HOW TEMPERATURE IS VARYING IN A SQUARE SECTION
//GIVEN: BOUNDRY CONDITIONS FOR A STEADY STATE PROBLEM
//INPUT REQUIRED:IN A PLATE OF X(0 TO 1) AND Y(0 TO 1) SECTION SIZE OR MESH 
//USING POISONS EQUATION AS SOURCE ARE TAKEN INTO ACCOUNT AND JACOBI ITERATION METHOD IS ADOPTED
//POISONS EQUATION IS SOLVED BY CENTRAL DIFFERENCE SCHEME
//GIVE THE SOURCE OR SINK VALUES AND THEIR COORDINATE LOCATIONS IN THE DOMAIN
//TEMPERATURE VALUES AT TOP,LEFT,RIGHT AND BOTTOM
//(IF NEUMANN BOUNDRY CONDITIONS ARE USED THEN SELECT ACCORDINGLY)
//EXAMPLE INPUT IS GIVEN IN BRACKETS
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
int main()
{
    double s,s2,top,left,right,bottom;
    int sm,s1,sms,t=-1,t1=-1,so;
    printf("Enter the size of mesh(0.1)  :");
    scanf("%lf",&s);
    if (s <= 0) {
    printf("Invalid mesh size. It must be greater than zero.\n");
    return 1;
    }
    printf("[ MESH SIZE SHOULD BE LESS THAN SOURCE COORDINATES]\nEnter number of sources(2) :");
    scanf("%d",&so);
    double source[so][2];  // 2D array for source positions
    double sourcev[so]; 
    printf("Enter the source coordinates and value(x y source) :\n(Example Q1 :.4 .5 1)\n(Example Q2 :.5 .5 -1) :\n");
    for(int i=0;i<so;i++){
        printf(" Q %d   :",i+1);
        scanf("%lf %lf %lf",&source[i][0],&source[i][1],&sourcev[i]);
    }
    printf("Enter Boundry conditions:\nOn top(0) :");
    scanf("%lf",&top);
    printf("On left (10):");
    scanf("%lf",&left);
    printf("On Right (10):");
    scanf("%lf",&right);
    printf("On Bottom(0) :");
    scanf("%lf",&bottom);
    s1=(int)(1.0/s)+1;//no. of nodes
    double *frame = (double *)malloc(s1 * s1 * sizeof(double));
    memset(frame, 0, s1 * s1 * sizeof(double));
    sm=s1-2;
    sms=(s1-2)*(s1-2);//no. of unknown nodes in coefficient matrix
    double *bmat=(double*)calloc(sms,sizeof(double));

    double *coeff = (double *)malloc(sms * sms * sizeof(double));
    memset(coeff, 0, sms * sms * sizeof(double));
    
//////////////////////////////////////////////////////////////////  
    
    printf("howmany sides have Neumann boundry condition (1):");
    int w;
    double p,k,l,r;
    scanf("%d",&w);
    int grad[w];
    printf("Top-1\nLeft-2\nRight-3\nBottom-4\nAll are Dirichlet-0\nEnter the sides from above which has gradient (Neumann condition)(3) :");
    for(int i=0;i<w;i++){
        scanf("%d",&grad[i]);
    }
/////////////////////////////////////////////////////////////////////  
    for(int i=0;i<s1;i++){//frame temp variation initially
        for(int j=0;j<s1;j++){
            if(i==0){// SETTING TOP SIDE OF FRAME TEMPERATURE
                frame[i*s1+j]=top;
            }
            if(i==s1-1){// SETTING BOTTOM SIDE OF FRAME TEMPERATURE
                frame[i*s1+j]=bottom;
            }
            if(j==0){// SETTING LEFT SIDE OF FRAME TEMPERATURE
                frame[i*s1+j]=left;
                if(i==0 && j==0){//LEFT TOP
                    if(top>left){
                      frame[i*s1+j]=top;  
                    }
                }
                if(i==(s1-1) && j==0){//LEFT BOTTOM
                    if(bottom>left){
                      frame[i*s1+j]=bottom;  
                    }
                }
            }
            if(j==(s1-1)){// SETTING RIGHT SIDE OF FRAME TEMPERATURE
                frame[i*s1+j]=right;
                if(i==0 &&j==s1-1){//RIGHT TOP
                    if(top>right){
                        frame[i*s1+j]=top;
                    }
                }
                if(i==s1-1 && j==s1-1){//RIGHT BOTTOM
                    if(bottom>right){
                        frame[i*s1+j]=bottom;
                    }
                }
            }
            
        }
    }//closed
//////////////////////////////////////////////////////////////////////////////////

//SETTING NEUMANN BOUNDRY CONDITION IF APPLICABLEBY OVERRIDING FRAME MATRIX
    for(int q=0;q<w;q++){
        if(grad[q]==1){//Top
            printf("Enter Top Gradient value(Neumann BC) :");
            scanf("%lf",&p);
            for(int q1=0;q1<s1;q1++){
                if(q1==0){
                    frame[0*s1+q1]=top;
                }
                else{
                    frame[0*s1+q1] = frame[0*s1+q1-1]+p;  
                }  
            }
        }
        if(grad[q]==2){//left
            printf("Enter left Gradient value(Neumann BC) :");
            scanf("%lf",&l);
            for(int q1=s1-1;q1>0;q1--){
                if(q1==s1-1){
                    frame[q1*s1+0]=left;
                }
                else{
                    frame[q1*s1+0] = frame[(q1+1)*s1+0]+l;  
                }
            }
        }
        if(grad[q]==3){//right
            printf("Enter right Gradient value(Neumann BC) (-1):");
            scanf("%lf",&r);
            for(int q1=s1-1;q1>0;q1--){
                if(q1==s1-1){
                frame[q1*s1+(s1-1)]=right;
                }
                else{
                frame[q1*s1+(s1-1)] = frame[(q1+1)*s1+(s1-1)]+r;  
                }
            }
        }
        if(grad[q]==4){//bottom
            printf("Enter bottom Gradient value(Neumann BC):");
            scanf("%lf",&k);
            for(int q1=0;q1<s1;q1++){
                if(q1==0){
                    frame[(s1-1)*s1+0]=bottom;
                }
                else{
                    frame[(s1-1)*s1+q1] = frame[(s1-1)*s1+(q1-1)]+k;  
                }
            }
                
        }
        if(grad[q]==0){
            break;
        }
    }
////////////////////////preparing coefficient matrix(coeff)coefficients of poisons equation//////////////////////////////////////
////////////////////////AND (bmat) VECTOR///////////////////////////////////////////////////////////////////////////////////////
////////////////////////COEFF*T=BMAT/////////////////////////////
//(T) IS VECTOR OF UNKNOWN TEMPERATURE AT EACH NODE IN THE FRAME EXCLUDING BOUNDRIES
    for(int i=0;i<sm;i++){
        for(int j=0;j<sm;j++){
            t+=1;
            t1+=1;
            if(j==0){///left side node
                if(i==0){///top left corner
                    bmat[t1]=-frame[i*s1+(j+1)]-frame[(i+1)*s1+j];
                    coeff[t*sms+t]=-4;
                    coeff[t*sms+t+1]=1;
                    coeff[t*sms+t+sm]=1;
                }
                else if(i>0 && i<sm-1){//leftside middle
                    coeff[t*sms+t]=-4;
                    coeff[t*sms+t+1]=1;
                    coeff[(t)*sms+t+sm]=1;
                    coeff[(t)*sms+t-sm]=1;
                    bmat[t1]=-frame[(i+1)*s1+(j)];
                }
                else if(i==sm-1){//leftside bottom
                    coeff[t*sms+t]=-4;  
                    coeff[t*sms+t+1]=1;
                    coeff[(t)*sms+t-sm]=1;
                    bmat[t1]=-frame[(i+2)*s1+(j+1)]-frame[(i+1)*s1+(j)];
                }
            }
            else if(j==sm-1){//right side
                if(i==0){///top right corner
                    coeff[t*sms+t]=-4;  
                    coeff[t*sms+t-1]=1;
                    if(t+sm<sms) coeff[t*sms+t+sm]=1;
                    bmat[t1]=-frame[(i)*s1+(j+1)]-frame[(i+1)*s1+(j+2)];
                }
                else if(i<sm-1 && i>0){
                    coeff[t*sms+t]=-4;
                    coeff[t*sms+t-1]=1;
                    if(t+sm<sms) coeff[t*sms+t+sm]=1;
                    coeff[(t)*sms+t-sm]=1;
                    bmat[t1]=-frame[(i+1)*s1+(j+2)];
                }
                else if(i==sm-1){//bottom right corner
                    coeff[t*sms+t]=-4;
                    coeff[t*sms+t-1]=1;
                    coeff[(t)*sms+t-sm]=1;  
                    bmat[t1]=-frame[(i+2)*s1+(j+1)]-frame[(i+1)*s1+(j+2)];
                }
            }
            else if(i==0 && j>0 && j<sm-1){//top nodes
                coeff[t*sms+t]=-4;
                coeff[t*sms+t+1]=1;
                if(t+sm<sms) coeff[t*sms+t+sm]=1;
                coeff[t*sms+t-1]=1;
                bmat[t1]=-frame[(i)*s1+(j+1)];
            }
            else if(i==sm-1 && j>0 && j<sm-1){//bottom nodes
                coeff[t*sms+t]=-4;
                coeff[t*sms+t+1]=1;
                coeff[(t)*sms+t-sm]=1;
                coeff[t*sms+t-1]=1;
                bmat[t1]=-frame[(i+2)*s1+(j+1)];
            }
            else{
                coeff[t*sms+t]=-4;
                coeff[t*sms+t+1]=1;
                coeff[(t)*sms+t-sm]=1;
                coeff[t*sms+t-1]=1; 
                if(t+sm<sms) coeff[t*sms+t+sm]=1;    
            }
            for(int y=0;y<so;y++){//ADDING SOURCE TERMS TO bmat VECTOR
                if(fabs(source[y][0] - i)<1e-6 && fabs(source[y][1] - j)<1e-6){
                    bmat[t1]=bmat[t1]+(s*s*sourcev[y]);
                }
            }
        }
    }//closed          **** coefficient matrix(coeff) is ready****B matrix(bmat) is ready*****
        
    ///setting initial temperature values
    double *temp=(double*)calloc(sms,sizeof(double));
    double *T=(double*)calloc(sms,sizeof(double));
    int count=0;
    double max=1;
//////////////////////////////Jacobi method///////////////////////////////////////
    while(max > 0.01){//IF diff bw previous and current valof temp >0.01 enter the loop we want that diff to be very small
        max=0;//INITIALISING
        for(int i=0;i<sms;i++){
            double sum=0;
            for(int j=0;j<sms;j++){
                if(i!=j){
                sum=sum+(coeff[i*sms+j]*temp[j]); 
                }
            }
            double val=((bmat[i]-sum)/coeff[i*sms+i]);//Jacobi iteration
            max = fmax(max,fabs(val-temp[i])); // Track the largest change(IF PREVIOUS VAL AND CURRENT VAL DIFF IS>0 STORE IN MAX)
            T[i] = val;
        }
        memcpy(temp, T, sms* sizeof(double));
        count+=1;
    }
    printf("%d Iterations\n\n",count);
    ////printing area
    printf("The temperature at each unknown node excluding boundaries   :\n\n");
    int c=0;
    for(int z=0;z<sms;z++){
        printf("%.8lf  ",T[z]);
        
        if(c==sm-1){
            printf("\n");
            c=-1;
        }
        c+=1;
    }
  /*for(int i=0;i<s1;i++){
        for(int j=0;j<s1;j++){
            printf("%lf\t",frame[i*s1+j]);
        }
        printf("\n");
    }*/
    free(frame);
    free(coeff);
    free(bmat);
    free(temp);
    free(T);
    return 0;
}
