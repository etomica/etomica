package etomica.zeolite;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.config.Configuration;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.space.Vector;

public class ConfigurationFileXYZ extends Configuration{

		public ConfigurationFileXYZ(String aConfName){
			super();
			confName = aConfName;
	        min = new double[3];
	        max = new double[3];
	        dim = new double[3];
	        for(int i=0;i<3;i++){
	        	min[i]=0;
	        	max[i]=0;
	        	dim[i]=0;
	        }
	                        
		}
		
		public void initializeCoordinates(Phase phase) {
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
            Vector pos = atom.getCoord().getPosition();
	        for (int i=0; i<pos.D(); i++) {
	            double coord = Double.valueOf(coordStr[i]).doubleValue();
	            if(coord<min[i]) min[i] = coord;
	            if(coord>max[i]) max[i] = coord;
                pos.setX(i, coord);
	        }
	        for(int i=0;i<3;i++){
	        	dim[i] = max[i]-min[i];
	        }
	    }
		private void translatePosition(AtomLeaf atom){
			for(int i=0;i<min.length;i++){
				atom.getCoord().getPosition().PE(i,-1*min[i]);	
				atom.getCoord().getPosition().PE(i,-0.5*dim[i]);
			}
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
		
		public Vector getUpdatedDimensions(){
			updatedDimensions = Space.makeVector(dim);
			
			updatedDimensions.TE(0,1.01);
			updatedDimensions.TE(1,1.01);
			updatedDimensions.TE(2,1.01);
			return updatedDimensions;
		}
		
        private static final long serialVersionUID = 2L;
		private Vector updatedDimensions;
		private double[] min;
		private double[] max;
		private double[] dim;
		private int[] nAtomsList;
		private String confName;

}
