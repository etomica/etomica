package etomica.zeolite;

import java.io.*;
import java.io.IOException;

public class DataCutter {
	public DataCutter(String filename,int cuts){
		this.cut = cuts;
		//Try reading the file
		try {
            fileReader = new FileReader(filename);
            System.out.println("Successivly opened inputFile");
        }
        catch(IOException e) {
            throw new RuntimeException("Cannot open "+filename+", caught IOException: " + e.getMessage());
        }
		//Read in the first line - numAtoms
        try {
            buffReader = new BufferedReader(fileReader);
            numAtoms = Integer.parseInt(buffReader.readLine());
        } catch(IOException e) {
            throw new RuntimeException("Problem reading "+filename+", caught IOException: " + e.getMessage());
        }
		//getNumLines
        try {
            numLines=0;
            while (buffReader.readLine()!=null){
              numLines++;
            }
        }catch(IOException e) {
            throw new RuntimeException("Problem reading "+filename+", caught IOException: " + e.getMessage());
        }
        //cut the number of lines into parts
        int blocks = numLines/(numAtoms+1);
        int blocksPerCut = blocks/cut+1;
        int linesPerCut = blocksPerCut*(numAtoms+1);
        System.out.println(blocksPerCut+" blocks per cut");
        //Print to new files
        try{
        	fileReader = new FileReader(filename);
        	buffReader = new BufferedReader(fileReader);
        	//Loop over all cuts
        	for(int i=0;i<cut;i++){
        		//System.out.println("Iteration "+i);
        		//Create fileWriter
        		String output = (filename+"_"+i);
        		fileWriter = new FileWriter(output,false);
        		//Print numAtoms
        		fileWriter.write(numAtoms+"\n");
        		if(i==0){
        			buffReader.readLine();
        		}
        		for(int j=0;j<linesPerCut;j++){
        			//System.out.println(j);
        			String out = buffReader.readLine();
        			if(out==null){
        				break;
        			}
        			fileWriter.write(out+"\n");
        		}
        		fileWriter.close();
        		fileWriter = null;
        	}
        }catch(IOException e) {
            System.err.println("Cannot close a file, caught IOException: " + e.getMessage());
        }   
	}
	
	public static void main(String[] args) {
		DataCutter cuts = new DataCutter("32_1000000_0.00611_2000",25);
		System.out.println("here goes");
	}
	
	private FileReader fileReader;
    private BufferedReader buffReader;
	private FileWriter fileWriter;
	private int numAtoms;
	private int cut;
	private int numLines;
}
