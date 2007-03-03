package etomica.meam;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import etomica.space.IVectorRandom;
import etomica.space.Space;

public class MSDProcessor {

    public MSDProcessor(Space space, String inputFile, String outputFile){
        
        msdOutput = outputFile;
        
        try {
            fileReader = new FileReader(inputFile);
        }
        catch(IOException e) {
            throw new RuntimeException("Cannot open "+inputFile+", caught IOException: " + e.getMessage());
        }
		
        try {
            buffReader = new BufferedReader(fileReader);
            numAtoms = Integer.parseInt(buffReader.readLine());
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
                        
        } catch(IOException e) {
            throw new RuntimeException("Problem reading "+inputFile+", caught IOException: " + e.getMessage());
        }
        
        coordBlock1 = new IVectorRandom[numAtoms];
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
    
    public void fillArrays(){
        double[] totalRsquared = new double[deltaTmax];        
        //Fills Block1 and 2, subtracts, and fills totalRsquared.  Repeat.
        for (int i=1; i<numBlocks; i++){
            try{
            
            buffReader.reset();
            
            //Gets buffReader to start of block 1 in question
            for (int j=0; j<(i-1)*numAtoms+1; j++){
                buffReader.readLine();
            }
            
            //Block 1 Loop - Adds XYZ lines from block 1
            for (int k=0; k<numAtoms; k++){
                String positionLine = buffReader.readLine();
                String [] coordString = positionLine.split(" +");
           
                for (int l=0; l<coordString.length; l++) {
                    double coord = Double.valueOf(coordString[l]).doubleValue();
                    coordBlock1[k].setX(l,coord);
                }
            }
            //Block 2 Loop - Restricts number of block pairs subtracted
            for (int deltaT=1; deltaT<deltaTmax+1; deltaT++){
                
                //Block 2 Loop - Adds XYZ lines from block 2
                for (int iatom=0; iatom<numAtoms; iatom++){
                    String positionLine = buffReader.readLine();
                    String [] coordString = positionLine.split(" +");
               
                    for (int icoord=0; icoord<coordString.length; icoord++) {
                        double coord = Double.valueOf(coordString[icoord]).doubleValue();
                        coordVector2.setX(icoord,coord);
                    }
                    
                    coordVector2.ME(coordBlock1[iatom]);
                    totalRsquared[deltaT-1] += coordVector2.squared();
                }
                
            }
            
            } catch(IOException e) {
                throw new RuntimeException("Problem creating array of positions, caught IOException: " + e.getMessage());
            }
        }
    
        /*
         * Each row of totalRsquared initially contains the summation of all
         * the atoms position differences for a specific deltaT (e.g. row 
         * three is the difference of block1 and block4, plus the difference 
         * of block2 and block5, etc.)
         * 
         * These sums are being divided by the number of atoms, and the respective
         * deltaT
         */
        for (int ideltaT=0; ideltaT<deltaTmax; ideltaT++){
            totalRsquared[ideltaT] /= (numAtoms*(numBlocks-ideltaT+1));
        }
        
        //Writes totalRsquared to file
        try{
            fileWriter = new FileWriter(msdOutput, false);
            fileWriter.write(numAtoms+"/n");
            fileWriter.write(numBlocks+"/n");
        
            for (int irow=0; irow<deltaTmax; irow++){
                fileWriter.write(irow+"     "+totalRsquared[irow]+"/n");
            }
            
            fileWriter.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
    
	private IVectorRandom [] coordBlock1;
	private IVectorRandom coordVector2;
	private int numAtoms;
    private int numBlocks;
    private int deltaTmax;
	private FileReader fileReader;
    private BufferedReader buffReader;
	private FileWriter fileWriter;
    private String msdOutput;

}
