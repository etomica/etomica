package etomica.zeolite;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.space.IVector;
import etomica.space.Space;

public class MSDProcessor {

	private int methane;
	
    public MSDProcessor(Space space, String inputFile, String outputFile){
    	methane =1;
        msdOutput = outputFile;
        msdInput = inputFile;
        try {
            fileReader = new FileReader(inputFile);
            System.out.println("Successivly opened inputFile");
        }
        catch(IOException e) {
            throw new RuntimeException("Cannot open "+inputFile+", caught IOException: " + e.getMessage());
        }
		
        try {
            buffReader = new BufferedReader(fileReader);
            numAtoms = Integer.parseInt(buffReader.readLine())+1;
        } catch(IOException e) {
            throw new RuntimeException("Problem reading "+inputFile+", caught IOException: " + e.getMessage());
        }
        
        try {
            int numLines=0;
            while (buffReader.readLine()!=null){
              numLines++;
            }
            //Sanity check! numAtoms should be a multiple of numLines.            
            if (numLines%numAtoms!=0){
                throw new RuntimeException("Read "+numLines+" line(s) and there are "+numAtoms+" atom(s) in the simulation.");
            }
            
            numBlocks = numLines/numAtoms;
            deltaTmax = numBlocks/3;
            System.out.println("deltaTmax= "+deltaTmax);            
        } catch(IOException e) {
            throw new RuntimeException("Problem reading "+inputFile+", caught IOException: " + e.getMessage());
        }
        try {
        	buffReader.close();
        	fileReader.close();
        } catch(IOException e){
        	throw new RuntimeException("Couldn't shut down readers, caught IOException: " +e.getMessage());
        }
        
        coordBlock1 = new IVector[numAtoms];
        coordVector2 = space.makeVector();
        
        for (int j=0; j<numAtoms; j++){
            coordBlock1[j] = space.makeVector();
        }
       
	}
    
    /**
     * Allows for a new value of deltaTmax.  This variable defaults to 1/3 of numBlocks.
     * deltaTmax is the largest difference in block number for which the Mean Squared
     *   Displacement will be calculated.
     * @param newDeltaTmax - The new value of deltaTmax.
     */
    public void setDeltaTmax(int newDeltaTmax){
        deltaTmax = newDeltaTmax;
    }
    public void setMethane(int m){
    	methane=m;
    }
    public void fillArrays(){
    	//Total RMS displacement
        double[] totalRsquared = new double[deltaTmax];  
        //XYZ Components
        double[][] RsquaredXYZ = new double[deltaTmax][3];
        double[] RsquaredComp = new double[3];
        double[] Temperature = new double[numBlocks];
        //Pull out the temperatures first
        try{
        	fileReader = new FileReader(msdInput);
        	buffReader = new BufferedReader(fileReader);
        	//Get's rid of numAtoms
        	buffReader.readLine();
        	//Get temperatures
        	int counter = 0;
        	for(int i=0;i<numBlocks*numAtoms;i++){
        		String positionLine = buffReader.readLine();
                String[] coordString = positionLine.split("\t");
                if(coordString.length==1){
                	Temperature[counter]=Double.valueOf(coordString[0]).doubleValue();
                	counter++;
                }
        	}
        }catch(IOException e) {
        	
        }
        
        //Fills Block1 and 2, subtracts, and fills totalRsquared.  Repeat.
        for (int i=1; i<deltaTmax+1; i++){
            System.out.println("Solving for iteration "+i);
        	try{
            	fileReader = new FileReader(msdInput);
            	buffReader = new BufferedReader(fileReader);
            //Gets buffReader to start of block 1 in question
            for (int j=0; j<(i-1)*(numAtoms)+1; j++){
                buffReader.readLine();
            }
            
            //Get temperature
            //Block 1 Loop - Adds XYZ lines from block 1
            for (int k=0; k<numAtoms-1; k++){
                String positionLine = buffReader.readLine();
                String[] coordString = positionLine.split("\t");
                for (int l=0; l<coordString.length; l++) {
                    double coord = Double.valueOf(coordString[l]).doubleValue();
                    coordBlock1[k].setX(l,coord);
                }
            }
            //Block 2 Loop - Restricts number of block pairs subtracted
            for (int deltaT=1; deltaT<deltaTmax+1; deltaT++){
                //Get rid of temperature in this block
            	buffReader.readLine();
                //Block 2 Loop - Adds XYZ lines from block 2
                for (int iatom=0; iatom<numAtoms-1; iatom++){
                    String positionLine = buffReader.readLine();
                    String [] coordString = positionLine.split("\t");
               
                    for (int icoord=0; icoord<coordString.length; icoord++) {
                        double coord = Double.valueOf(coordString[icoord]).doubleValue();
                        coordVector2.setX(icoord,coord);
                    }
                    
                    coordVector2.ME(coordBlock1[iatom]);
                    totalRsquared[deltaT-1] += coordVector2.squared();
                    for(int j=0;j<RsquaredComp.length;j++){
                    	RsquaredComp[j] +=Math.pow(coordVector2.x(j),2.0);
                    	RsquaredXYZ[deltaT-1][j] += Math.pow(coordVector2.x(j),2.0);
                    }
                }
            }
            
            } catch(IOException e) {
                throw new RuntimeException("Problem creating array of positions, caught IOException: " + e.getMessage());
            }
        }
    
        /*
         * Each row of totalRsquared initially contains the summation of all
         * the atoms position differences for a specific deltaT (i.e. row 
         * three is the difference of block1 and block4, plus the difference 
         * of block2 and block5, etc.)
         * 
         * These sums are being divided by the number of atoms, and the respective
         * deltaT
         */
        for (int ideltaT=0; ideltaT<deltaTmax; ideltaT++){
            totalRsquared[ideltaT] /= ((numAtoms-1)*(numBlocks-ideltaT+1));
            for(int j=0;j<3;j++){
            	RsquaredXYZ[ideltaT][j] /= ((numAtoms-1)*(numBlocks-ideltaT+1));
            }
        }
        
        //Writes totalRsquared to file
        try{
        	System.out.println("Created new output file");
            fileWriter = new FileWriter(msdOutput, false);
            fileWriter.write((numAtoms-1)+"\n");
            fileWriter.write(numBlocks+"\n");
            
            int temp = 0;
            for(int i=0;i<Temperature.length;i++){
            	temp+=Temperature[i];
            }
            temp /=Temperature.length;
            
            fileWriter.write(temp+"\n");
            
            /*
            for(int i=0;i<Temperature.length;i++){
            	fileWriter.write(Temperature[i]+"\n");
            }*/
            
            for (int irow=0; irow<deltaTmax; irow++){
                fileWriter.write(irow+"\t"+totalRsquared[irow]+"\n");
                //fileWriter.write(irow+"\t"+RsquaredTotal[irow][0]+"\n");
                //fileWriter.write(irow+"\t"+RsquaredTotal[irow][1]+"\n");
                //fileWriter.write(irow+"\t"+RsquaredTotal[irow][2]+"\n");
                //fileWriter.write("\n"); 
            }
            
            fileWriter.write("Time dependent data for X,Y,Z\n");
            for(int j=0;j<3;j++){
            	for(int i=0;i<deltaTmax;i++){
            		fileWriter.write(RsquaredXYZ[i][j]+"\n");
            	}
            	fileWriter.write("\n");
            }
            fileWriter.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    
	private IVector [] coordBlock1;
	private IVector coordVector2;
	private int numAtoms;
    private int numBlocks;
    private int deltaTmax;
	private FileReader fileReader;
    private BufferedReader buffReader;
	private FileWriter fileWriter;
    private String msdOutput;
    private String msdInput;

}
