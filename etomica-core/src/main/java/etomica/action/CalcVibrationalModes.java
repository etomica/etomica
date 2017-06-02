/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.box.Box;
import etomica.api.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.math.numerical.CalcGradientDifferentiable;

/**
 * Calculates the eigenvalues of an NxN matrix.  Assumes these are actually changes in
 * a group of molecules' force with respect to changes in other molecules' positions, and
 * returns eigenvalues.
 * 
 * @author msellers and ajschultz
 *
 */

public class CalcVibrationalModes implements IAction, Serializable {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;
    double [] modes;
    double [] frequencies;
    double [] positions;
    int[] d;
    double [][] dForces;
    double prodFreq;
    int[] modeSigns;
    Matrix fC;
    double mass;
    int writeCount;
    IMoleculeList ms;
    EigenvalueDecomposition eigenDecomp;
    CalcGradientDifferentiable cgd;
    
    public CalcVibrationalModes(){
        
        writeCount = 0;
       
        }

    public void setup(Box aBox, PotentialMaster aPotentialMaster, IMoleculeList movableSet, Space _space){
        ms = movableSet; 
        mass = ms.getMolecule(0).getType().getAtomType(0).getMass();
        cgd = new CalcGradientDifferentiable(aBox, aPotentialMaster, ms, _space);
        d = new int[movableSet.getMoleculeCount()*3];
        positions = new double[d.length];
        dForces = new double[positions.length][positions.length];
        modeSigns = new int[3];
    }
    
    public void actionPerformed() {
        // setup position array
        for(int i=0; i<ms.getMoleculeCount(); i++){
            for(int j=0; j<3; j++){
                positions[(3*i)+j] = ms.getMolecule(i).getChildList().getAtom(0).getPosition().getX(j);
            }
        }
        // fill dForces array
        for(int l=0; l<d.length; l++){
            d[l] = 1;
            System.arraycopy(cgd.df2(d, positions), 0, dForces[l], 0, d.length);
            //System.out.println("  -Calculating force constant row "+l+"...");
            d[l] = 0;
        }
        fC = new Matrix(dForces);
        eigenDecomp = new EigenvalueDecomposition(fC);
        modes = eigenDecomp.getRealEigenvalues();
        for(int i=0; i<modes.length; i++){   
            if(modes[i]>0.0){modeSigns[0]++;}
            else {modeSigns[1]++;}
        }
        modeSigns[2] = modes.length;
        frequencies = new double[modes.length];
        for(int i=0; i<frequencies.length; i++){
            //Negative mode catch
            if(modes[i]<0.0){
                frequencies[i] = 0.0;
                continue;
            }
            frequencies[i] = Math.sqrt(modes[i]) / (2*Math.PI);
        }
        prodFreq = 1;      
        for(int i=0; i<frequencies.length; i++){
            if(frequencies[i] == 0.0){
                //System.out.println("ZERO TEST");
                continue;
                
            }
            prodFreq = prodFreq * frequencies[i];
        }
        
        System.out.println("Normal mode vibrational data calculated.");
        //System.out.println(prodFreq);
        
    }
    
    /**
     * Performs an eigenvalue decomposition of the NxN matrix, fC.
     * 
     * @return
     */
    public double[] getLambdas(){
        return modes;
    }
    
    /**
     * Returns a one dimensional array of length 3 with the total number of
     * positive, negative, and imaginary modes.
     * 
     * @return modeSigns one-dimensional array
     */
    public int[] getModeSigns(){
        return modeSigns;
    }
    
    /**
     * Calculates the frequencies (omegas) of wave vectors described by the eigenvalues (lambdas) of
     * our system.
     * 
     * lambda = 4 * pi^2 * omega^2.
     * 
     * @return frequencies one-dimensional array of doubles
     */
    public double[] getFrequencies(){
        // where ;
        modes = getLambdas();
        

        return frequencies;
    }
    
    public double getProductOfFrequencies(){
        return prodFreq;
    }
    
    public void writeDataToFile(String file){
        FileWriter writer;
        //LAMBDAS
        try { 
            writer = new FileWriter(file);
            writer.write("Output "+writeCount+"--------"+"\n");
            writer.write(modeSigns[0]+" positive modes"+"\n");
            writer.write(modeSigns[1]+" negative modes"+"\n");
            writer.write(modeSigns[2]+" total modes"+"\n");
            writer.write("-modes"+"\n");
            for(int i=0; i<modes.length; i++){
                writer.write(modes[i]+"\n");
            }
            writer.write("-frequencies"+"\n");
            for(int i=0; i<frequencies.length; i++){
                writer.write(frequencies[i]+"\n");
            }
            writer.write("-frequency product"+"\n");
            writer.write(prodFreq+"\n");
            writer.close();
            writeCount++;
        }catch(IOException e) {
            System.err.println("Cannot open file, caught IOException: " + e.getMessage());
            return;
        }
    }

    
}
