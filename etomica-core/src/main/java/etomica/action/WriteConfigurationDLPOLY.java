/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.chem.elements.*;
import etomica.molecule.IMolecule;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Vector;
import etomica.space3d.BoundaryTruncatedOctahedron;

import java.io.File;
import java.io.IOException;
import java.util.Formatter;
import java.util.HashMap;


/**
 * 
 * Dumps a box's configuration to a file.  The coordinates are written in a 
 * format that can be read in by ConfigurationFile.  The output file is normally 
 * named 'CONFIG', which is read by the DL_MULTI package as an input file
 * 
 */
public class WriteConfigurationDLPOLY implements IAction {
	
	public WriteConfigurationDLPOLY(){
		elementHash = new HashMap<IElement,String>();
		elementHash.put(Carbon.INSTANCE, "CA");
		elementHash.put(Hydrogen.INSTANCE, "HY");
		elementHash.put(Nitrogen.INSTANCE, "NI");
		elementHash.put(Oxygen.INSTANCE, "OX");
		
	}
	
	public HashMap<IElement,String> getElementHash(){
		return elementHash;
	}
	
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
    	if(!(box.getBoundary() instanceof Boundary)) {
    		throw new RuntimeException("The boundary within the box in WriteConfigurationDLPOLY MUST be derived from etomica.space.Boundary.");
        }
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
        String fileName = confName;
        
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
        	} else {
        		/*
        		 * This handles: 
        		 *  0: non-periodic
        		 *  2: cubic/ orthorhombic periodic
        		 *  6: slit boundary with periodicity in x and y
        		 */
        		
        		if (boundary.getPeriodicity(2)){
        			boundaryType = 2;
        		} else if (!boundary.getPeriodicity(0)){
        			boundaryType = 0;
        		} else {
        			boundaryType = 6;
        		}
        	}
        	
        	formatter.format("\n%10d%10d\n", new Object[]{new Integer(writeVelocity? 1:0), boundaryType});


            for (int i=0; i<3; i++){
                Vector cell = boundary.getEdgeVector(i);
        		for (int j=0; j<3; j++){
        			formatter.format("%20f",new Object[]{cell.getX(j)});
        		}
        		formatter.format("\n");
        	}
        	
            int atomCount =1;
            
            for (int flipIndex=0; flipIndex<2; flipIndex++){
            	
	            for (int iMolec=0; iMolec<box.getMoleculeList().getMoleculeCount(); iMolec++) {
	            	
	            	//for Orthorhombic Paracetamol
	            	if(boundaryType==2){
		            	if ((iMolec/4)%2 != flipIndex){
		            		continue;
		            	}
	            	}
	            	
	            	//for Monoclinic Paracetamol
	            	if(boundaryType==3){
		            	if ((iMolec/2)%2 != flipIndex){
		            		continue;
		            	}
	            	}
	            	IMolecule molecule = box.getMoleculeList().getMolecule(iMolec);
	                for (int iLeaf = 0; iLeaf<molecule.getChildList().size(); iLeaf++){
	                	IAtom atom = molecule.getChildList().get(iLeaf);
	                	String atomName = elementHash.get(atom.getType().getElement());
	       
	                	formatter.format("%8s%10d\n", new Object[]{atomName, atomCount});
	                	atomCount++;
	                	Vector atomPos = atom.getPosition();
	                	formatter.format("%20.12f%20.12f%20.12f\n", new Object[]{atomPos.getX(0), atomPos.getX(1), atomPos.getX(2)});
	                		                	
	                	
	                	if (writeVelocity){
	                		Vector atomVelocity = ((IAtomKinetic)atom).getVelocity();
	                		formatter.format("%20f%20f%20f\n", 
	                				new Object[]{atomVelocity.getX(0), atomVelocity.getX(1), atomVelocity.getX(2)});
	                	}
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
    private HashMap<IElement,String> elementHash;
}
