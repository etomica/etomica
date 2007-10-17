package etomica.action;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Formatter;

import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.BoundaryPeriodic;
import etomica.space.IVector;
import etomica.space3d.BoundaryTruncatedOctahedron;
import etomica.util.Arrays;

/**
 * Dumps a box's configuration to a file.  The coordinates are written in a 
 * format that can be read in by ConfigurationFile.  The output file has a 
 * "pos_new" extension, which should be renamed to "pos" for use with
 * ConfigurationFile.
 */
public class WriteConfigurationDLPOLY implements Action {

    /**
     * Sets the configuration name.  The file written to is newConfName.pos_new
     */
    public void setConfName(String newConfName) {
        confName = newConfName;
    }
    
    /**
     * Returns the configuration name.  The file written to is confName.pos_new
     */
    public String getConfName() {
        return confName;
    }
    
    /**
     * Sets the box whose atom coordinates get written to the file.
     */
    public void setBox(Box newBox) {
        box = newBox;
        setDoApplyPBC(true);
    }
    
    /**
     * Returns the box whose atom coordinates get written to the file.
     */
    public Box getBox() {
        return box;
    }
    
    /**
     * Directs the writer to apply periodic boundary conditions or not (true 
     * by default).
     */
    public void setDoApplyPBC(boolean newDoApplyPBC) {
        doApplyPBC = newDoApplyPBC;
    }
    
    /**
     * Returns true if PBC are applied to coordinates written to the file.
     */
    public boolean getDoApplyPBC() {
        return doApplyPBC;
    }
    
    public void setWriteVelocity(boolean newWriteVelocity){
    	writeVelocity = newWriteVelocity;
    }
    
    public boolean getWriteVelocity(){
    	return writeVelocity;
    }
    
    /**
     * Writes the leaf Atom coordinates to the file confName.pos_new.  If the
     * file exists, it is overwritten.
     */
    public void actionPerformed() {
        FileWriter fileWriter = null;
        String fileName = confName + ".dlp";
        
//        try { 
//            fileWriter = null; //new FileWriter(fileName);
//        }catch(IOException e) {
//            System.err.println("Cannot open "+fileName+", caught IOException: " + e.getMessage());
//            return;
//        }
        try {
        	
        	Formatter formatter = new Formatter(new File(fileName));
        	
        	Boundary boundary = box.getBoundary();
        	int boundaryType = -1;
        	if (boundary instanceof BoundaryTruncatedOctahedron){
        		boundaryType = 4;
        	} else if(boundary instanceof BoundaryDeformablePeriodic){
        		boundaryType = 3;
        	} else if (boundary instanceof BoundaryPeriodic){
        		boolean[] periodicity = ((BoundaryPeriodic)boundary).getPeriodicity();
        		
        		/*
        		 * This handles: 
        		 *  0: non-periodic
        		 *  2: cubic/ orthorhombic periodic
        		 *  6: slit boundary with periodicity in x and y
        		 */
        		
        		if (periodicity[2]){
        			boundaryType = 2;
        		} else if (!periodicity[0]){
        			boundaryType = 0;
        		} else {
        			boundaryType = 6;
        		}
        	} else {
        		throw new RuntimeException("Input a correct boundary type");
        	}
        	
        	formatter.format("\n%10d%10d\n", new Object[]{new Integer(writeVelocity? 1:0), boundaryType});
        	
        	
        	IVector[] cell = boundary.getPeriodicVectors();
        	for (int i=0; i<cell.length; i++){
        		for (int j=0; j<3; j++){
        			formatter.format("%20f",new Object[]{cell[i].x(j)});
        		}
        		formatter.format("\n");
        	}
        	
            int[] elementNum = new int[0]; 
            String[] elementString = new String[0];
            int atomCount = 1;    
            
            for (int iMolec=0; iMolec<box.moleculeCount(); iMolec++) {
            	for (int num=0; num<elementNum.length; num++ ){
            		elementNum[num] = 0;
            		elementString[num] = null;
            	}
            	
            	IAtomGroup molecule = (IAtomGroup)box.molecule(iMolec);
                for (int iLeaf=0; iLeaf<molecule.getChildList().getAtomCount(); iLeaf++){
                	IAtomPositioned atom = (IAtomPositioned)molecule.getChildList().getAtom(iLeaf);
                	String atomName = ((AtomTypeLeaf)atom.getType()).getElement().getSymbol();
                	
                	boolean success = false;
                	for (int num=0; num<elementNum.length; num++){
                		if (elementString[num]==null){
                			success = true;
                			elementString[num] = atomName;
                			elementNum[num] = 1;
                			atomName = atomName+"1";
                			break;
                		}
                		if (elementString[num].equals(atomName)){
                			success = true;
                			elementNum[num]++;
                			atomName = atomName+elementNum[num];
                			break;
                		}
                	}
                	if (!success){
                		elementString = (String[])Arrays.addObject(elementString, atomName);
                		elementNum = Arrays.resizeArray(elementNum, elementNum.length+1);
                		elementNum[elementNum.length-1] = 1;
                		atomName = atomName+1;
                	}
                	
                	formatter.format("%8s%10d\n", new Object[]{atomName, atomCount});
                	atomCount++;
                	
                	IVector atomPos = atom.getPosition();
                	formatter.format("%20f%20f%20f\n", new Object[]{atomPos.x(0), atomPos.x(1), atomPos.x(2)});
                	
                	if (writeVelocity){
                		IVector atomVelocity = ((IAtomKinetic)atom).getVelocity();
                		formatter.format("%20f%20f%20f\n", 
                				new Object[]{atomVelocity.x(0), atomVelocity.x(1), atomVelocity.x(2)});
                	}
                }
                
            }
            formatter.close();
        } catch(IOException e) {
            System.err.println("Problem writing to "+fileName+", caught IOException: " + e.getMessage());
        }
    }

    private String confName;
    private Box box;
    private boolean doApplyPBC;
    private boolean writeVelocity;
}
