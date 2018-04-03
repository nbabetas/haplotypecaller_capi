#include <jni.h>
#include <stdio.h>
#include "Main.h"
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>


JNIEXPORT void JNICALL
Java_Main_nativePrint(JNIEnv *env, jobject obj, jlong *memory,jlong *results, jlong size_x, jlong size_y) {
    MAIN(memory,results,size_x,size_y);
        
     
}
void MAIN (long *memory,long *results,long size_x,long size_y){
    
    char *XDATA=(char*) malloc(size_x * sizeof (char));
    char *YDATA=(char*) malloc(size_y * sizeof (char));
    int j;
    
    Read_DATA(memory, XDATA, size_x); //read first sequence
    Read_DATA((char*) memory + size_x, YDATA, size_y); // read second sequence (but first move pointer)
	   

	


   for (j=0;j<size_x;j++){
        printf("%c",XDATA[j]); 
         *results = (int)j;
         results = (int*) results + 1;
    } 
    
    printf("\n");
    
    for (j=0;j<size_y;j++){
        printf("%c",YDATA[j]); 
         *results = (int)j;
         results = (int*) results + 1;
    }	
    
    free(XDATA);
    free(YDATA);
  
}
 void Read_DATA(long* memory, char* DATA, long size) {
    int i;
    for (i = 0; i < size; i++) {
        DATA[i] = (*(char*) (memory));
        memory = (char*) memory + 1; // Load next char
    }
    
}




 



