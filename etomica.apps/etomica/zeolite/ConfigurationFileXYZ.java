package etomica.zeolite;

import etomica.config.Configuration;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorArrayListCompound;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.space.Space;
import etomica.space.Vector;

public class ConfigurationFileXYZ extends Configuration{

		public ConfigurationFileXYZ(Space space, String aConfName){
			super(space);
			confName = aConfName;
	        newPos = space.makeVector();
	        atomIterator = new AtomIteratorArrayListCompound();
	        min = new double[3];
	        max = new double[3];
	        dim = new double[3];
	        for(int i=0;i<3;i++){
	        	min[i]=0;
	        	max[i]=0;
	        	dim[i]=0;
	        }
	                        
		}
		
		public void initializePositions(AtomArrayList[] atomLists) {
	        if (atomLists.length == 0) return;
	        String fileName = confName+".xyz";
	        FileReader fileReader;
	        try {
	            fileReader = new FileReader(fileName);
	        }catch(IOException e) {
	            throw new RuntimeException("Cannot open "+fileName+", caught IOException: " + e.getMessage());
	        }
	        try {
	            BufferedReader bufReader = new BufferedReader(fileReader);
	            atomIterator.setLists(atomLists);
	            AtomIteratorTree leafIterator = new AtomIteratorTree();
	            atomIterator.reset();
	            //Skips the first line, which contains numAtoms
	            bufReader.readLine();
	            while (atomIterator.hasNext()) {
	                Atom molecule = atomIterator.nextAtom();
	                if (molecule.getNode().isLeaf()) {
	                    setPosition((AtomLeaf)molecule,bufReader.readLine());
	                }
	                else {
	                    leafIterator.setRoot(molecule);
	                    leafIterator.reset();
	                    while (leafIterator.hasNext()) {
	                        setPosition((AtomLeaf)leafIterator.nextAtom(),bufReader.readLine());
	                    }
	                }
	            }
	            fileReader.close();
	        } catch(IOException e) {
	            throw new RuntimeException("Problem writing to "+fileName+", caught IOException: " + e.getMessage());
	        }
	        atomIterator.reset();
	        
	        while(atomIterator.hasNext()){
	        	Atom molecule = atomIterator.nextAtom();
	        	if(molecule.getNode().isLeaf()){
	        		translatePosition((AtomLeaf)molecule);
	        	}	
	        }
	        
		}
		
		private void setPosition(AtomLeaf atom, String string) {
	        String[] coordStr = string.split(" +");
	        double[] coord = new double[coordStr.length];
	        for (int i=0; i<coord.length; i++) {
	            coord[i] = Double.valueOf(coordStr[i]).doubleValue();
	            if(coord[i]<min[i]) min[i] = coord[i];
	            if(coord[i]>max[i]) max[i] = coord[i];
	        }
	        newPos.E(coord);
	        atom.coord.position().E(newPos);
	        for(int i=0;i<3;i++){
	        	dim[i] = max[i]-min[i];
	        }
	    }
		private void translatePosition(AtomLeaf atom){
			for(int i=0;i<min.length;i++){
				atom.coord.position().PE(i,-1*min[i]);	
				atom.coord.position().PE(i,-0.5*dim[i]);
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
			updatedDimensions = space.makeVector();
			updatedDimensions.E(dim);
			
			updatedDimensions.TE(0,1.01);
			updatedDimensions.TE(1,1.01);
			updatedDimensions.TE(2,1.01);
			return updatedDimensions;
		}
		
        private static final long serialVersionUID = 1L;
		private Vector updatedDimensions;
		private double[] min;
		private double[] max;
		private double[] dim;
		private int[] nAtomsList;
		private String confName;
	    private Vector newPos;
	    private final AtomIteratorArrayListCompound atomIterator;

}
