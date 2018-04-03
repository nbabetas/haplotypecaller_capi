javac HaplotypeCaller.java
gcc -shared  -Wall -c  -g -std=c99 -m64 -lrt -lpthread -fopenmp  -lgomp -I/usr/lib/jvm/java-8-oracle/include -I/usr/lib/jvm/java-8-oracle/include/linux -Ipslse/common -Ipslse/libcxl   -fPIC  -MMD -MP -MF "/home/bampetas/Documents/JNI/pairhmm.o.d"  -o pairhmm.o  pairhmm.c  
gcc -shared  -Wall  -g -std=c99 -m64 -lrt -lpthread -fopenmp  -lgomp -I/usr/lib/jvm/java-8-oracle/include -I/usr/lib/jvm/java-8-oracle/include/linux  -I./pslse/common -I./pslse/libcxl  -fPIC  -MMD -MP -o pairhmm.so  pairhmm.o utils.c batch.c ./pslse/libcxl/libcxl.a -fPIC 
