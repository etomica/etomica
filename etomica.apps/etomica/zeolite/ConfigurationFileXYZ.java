package etomica.zeolite;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.config.Configuration;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.space.Space;

public class ConfigurationFileXYZ extends Configuration{

		public ConfigurationFileXYZ(String aConfName){
			super();
			confName = aConfName;
		}
		
		public void initializeCoordinates(Phase phase) {
            min = phase.getSpace().makeVector();
            max = phase.getSpace().makeVector();
            dim = phase.getSpace().makeVector();
	        String fileName = confName+".xyz";
	        FileReader fileReader;
	        try {
	            fileReader = new FileReader(fileName);
	        }catch(IOException e) {
	            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
	        }
            AtomIteratorLeafAtoms atomIterator = new AtomIteratorLeafAtoms(phase);
	        try {
	            BufferedReader bufReader = new BufferedReader(fileReader);
	            atomIterator.reset();
	            //Skips the first line, which contains numAtoms
	            bufReader.readLine();
	            while (atomIterator.hasNext()) {
	                AtomLeaf atom = (AtomLeaf)atomIterator.nextAtom();
	                setPosition(atom,bufReader.readLine());
	            }
	            fileReader.close();
	        } catch(IOException e) {
	            throw new RuntimeException("Problem writing to "+fileName+", caught IOException: " + e.getMessage());
	        }
	        atomIterator.reset();
	        
	        while(atomIterator.hasNext()){
	        	AtomLeaf molecule = (AtomLeaf)atomIterator.nextAtom();
	        	translatePosition(molecule);
	        }
	        
		}
		
		private void setPosition(AtomLeaf atom, String string) {
	        String[] coordStr = string.split(" +");
            IVector pos = atom.getPosition();
            for (int i=0; i<pos.getD(); i++) {
	            double coord = Double.valueOf(coordStr[i]).doubleValue();
	            if(coord<min.x(i)) min.setX(i,coord);
	            if(coord>max.x(i)) max.setX(i,coord);
                pos.setX(i, coord);
	        }
            dim.E(max);
            dim.ME(min);
	    }

        private void translatePosition(AtomLeaf atom){
            atom.getPosition().ME(min);
            atom.getPosition().PEa1Tv1(-0.5,dim);
		}
        
		public int[] getNumAtoms(){
			String fileName = confName+".xyz";
	        FileReader fileReader;
	        try {
	            fileReader = new FileReader(fileName);
	        }catch(IOException e) {
	            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
	        }
	        BufferedReader bufReader = new BufferedReader(fileReader);
	        String numLine;
	        try{
	        	numLine = bufReader.readLine();
	        }
	        catch(IOException e) {
	            throw new RuntimeException("Problem writing to "+fileName+", caught IOException: " + e.getMessage());
	        }
	        String[] coordStr = numLine.split(" +");
	        double[] coord = new double[coordStr.length];
	        nAtomsList = new int[coordStr.length];
	        for (int i=0; i<coord.length; i++) {
	            nAtomsList[i] = Integer.parseInt(coordStr[i]);
	        }
	        return nAtomsList;
		}
		
		public IVector getUpdatedDimensions(){
			updatedDimensions = Space.makeVector(dim.getD());
			
			updatedDimensions.TE(1.01);
			return updatedDimensions;
		}
		
        private static final long serialVersionUID = 2L;
		private IVector updatedDimensions;
		private IVector min;
		private IVector max;
		private IVector dim;
		private int[] nAtomsList;
		private String confName;

}
