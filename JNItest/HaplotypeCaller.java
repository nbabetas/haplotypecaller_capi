/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Scanner;
import sun.misc.Unsafe;
import static sun.misc.Unsafe.getUnsafe;
import static sun.nio.ch.IOStatus.normalize; //ayti

/**
 *
 * @author bampetas
 */
public class HaplotypeCaller{
    static ArrayList<String> hapl;
    static ArrayList<String> read;
    static ArrayList<Integer> prob;
	static ArrayList<Integer> Size;
	static int haplChars;
  	static int readChars;
	static int choise = 0;
    private final static long INT_SIZE_IN_BYTES = 1;

    static {
		 System.load("/home/nbampetas/bulk/Simulation/sim/JNI test/pairhmm.so");
    }		

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException  {
        Read_file();	 
        showBytes();
        showAllocateTooMuch();

    }



    private native void nativePrint(long memoryHapl,long memoryRead,long memoryProb,long memorySize,long memoryResults,long done, long afu,int prob_size, int size_size,int readChars,int haplChars,int choise);
	private native void freeAccel(long afu);


    public static void Read_file() throws FileNotFoundException {
		hapl =new ArrayList<>();
		read =new ArrayList<>(); 
		prob =new ArrayList<>();
		Size =new ArrayList<>();
		haplChars=0;
		readChars=0;

        Scanner s = new Scanner(new File("/shares/bulk/nbampetas/Simulation/sim/JNI test/m.txt")).useDelimiter("\n");        
        while (s.hasNext()){
			Size.add(Integer.parseInt(s.next()));
       		String content = s.next();
			String[] parts = content.split(" ");
			
			if(parts.length==1){				
				read.add(content);
				readChars=readChars+content.length();

			}else{				
				hapl.add(parts[0]);
				haplChars=haplChars+parts[0].length();

				for (int i=1; i<parts.length;i++){			
					prob.add(Integer.parseInt(parts[i]));
				}
					
			}
		
		}
    }


    public static void showBytes() {
        long memoryHapl = 0;
		long memoryRead = 0;
		long memorySize = 0;
		long memoryProb	= 0;
        long memoryResults = 0;
		long done = 0 ;
		long afu =0;
		
        try {
            Unsafe unsafe = getUnsafe();
            
            // Allocate given memory size
			int batches = (Size.size()/2)/16;

            memoryHapl = unsafe.allocateMemory(haplChars);
			memoryRead = unsafe.allocateMemory(readChars+10);
			memorySize = unsafe.allocateMemory(Size.size()*4);
			memoryProb = unsafe.allocateMemory(prob.size()*4);
			memoryResults=  unsafe.allocateMemory(((7*4)*batches*16)*2);

			done = unsafe.allocateMemory(32);
			afu  = unsafe.allocateMemory(8);
			
            
            // Write value to the allocated memory   
			int ch=0;
			for (int i= 0;i< hapl.size();i++){
				for(int j=0; j<hapl.get(i).length();j++){
					 unsafe.putAddress(memoryHapl+ (ch * INT_SIZE_IN_BYTES), hapl.get(i).charAt(j));
					 ch++;
				}
			} 
			
		    ch=0;   
			for (int i= 0;i< read.size();i++){
				for(int j=0; j<read.get(i).length();j++){
					 unsafe.putAddress(memoryRead+ (ch * INT_SIZE_IN_BYTES), read.get(i).charAt(j));
					 ch++;
				}
			}
			readChars=ch;
				for (int i= 0;i< prob.size();i++){
					 unsafe.putAddress(memoryProb+ (i * 4), prob.get(i));
			} 
			for (int i= 0;i< Size.size();i++){
					 unsafe.putAddress(memorySize+ (i * 4), Size.get(i));
					 
			}        
        	

            // Call C function
            new HaplotypeCaller().nativePrint(memoryRead,memoryHapl,memoryProb,memorySize,memoryResults,done,afu,prob.size(),Size.size(),haplChars,readChars,choise);
			choise ++;
			
            // Read value from memory
            	           
			//lint readValue = (int) unsafe.getAddress(memoryResults);
			//	System.out.println(readValue);  
           //  Take results from accelerator


			// Get the address of all WEDs and check when you have a ready result
			long wed0 =(long)unsafe.getAddress(done);
			long wed1 =(long)unsafe.getAddress(done+8); 
			long wed2 =(long)unsafe.getAddress(done+16);
			long wed3 =(long)unsafe.getAddress(done+24);
		   
			byte kernel0 = (byte)unsafe.getAddress(wed0); 
			byte kernel1 = (byte)unsafe.getAddress(wed1);
			byte kernel2 = (byte)unsafe.getAddress(wed2);	
			byte kernel3 = (byte)unsafe.getAddress(wed3);
			System.out.println(kernel0);
			System.out.println(kernel1);
			System.out.println(kernel2);
			System.out.println(kernel3);
			byte done1 =0;
			byte done2 =0;
			byte done3 =0;
			byte done4 =0;
			while (done1 !=1 || done2 !=1 || done3 !=1 || done4 !=1) {
				if(done1==0){
					 kernel0 = (byte)unsafe.getAddress(wed0); 
				}
				if(done2==0){
				 	 kernel1 = (byte)unsafe.getAddress(wed1);
				}
				if(done3==0){
					 kernel2 = (byte)unsafe.getAddress(wed2);	
				}
				if(done4==0){
				 	kernel3 = (byte)unsafe.getAddress(wed3);
				}
				
				if (kernel0 ==1){
					System.out.println("done_1");
					done1=1;
					kernel0 =0;
				}
				if (kernel1 ==1){
					System.out.println("done_2");
					done2=1;
					kernel1 =0;
				}
				if (kernel2 ==1){
					System.out.println("done_3");
					done3=1;
					kernel2 =0;
				}
				if (kernel3 ==1){
					System.out.println("done_4");
					done4=1;
					kernel3 =0;
				}					
			}
			System.out.println("done ALL");
			long afu_address =(long)unsafe.getAddress(afu);// get address of accelerator
			new HaplotypeCaller().freeAccel(afu_address); // free and unmap accelerator
			
            for (int i =0 ; i<7*batches*16*2;i++){
                 int readValue =(int) unsafe.getAddress(memoryResults+i*4);
                System.out.println( Integer.toHexString(readValue));
           }
            //Free memory
            unsafe.freeMemory(memoryHapl);
			unsafe.freeMemory(memoryRead);
			unsafe.freeMemory(memoryProb);
			unsafe.freeMemory(memorySize);
            unsafe.freeMemory(memoryResults);
			unsafe.freeMemory(done);
			unsafe.freeMemory(afu);
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }

    public static void showAllocateTooMuch() {
        try {
            Unsafe unsafe = getUnsafe();

            long bytes = Integer.MAX_VALUE;	// It's way too much memory!!		
            // Allocate given memory size
            long memoryAddress = unsafe.allocateMemory(bytes);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static Unsafe getUnsafe() throws Exception {
        // Get the Unsafe object instance
        Field field = sun.misc.Unsafe.class.getDeclaredField("theUnsafe");
        field.setAccessible(true);
        return (sun.misc.Unsafe) field.get(null);
    }

    public static long sizeOf(Object object) throws Exception {
        return getUnsafe().getAddress(
                normalize(getUnsafe().getInt(object, 4L)) + 12L);
    }

}

