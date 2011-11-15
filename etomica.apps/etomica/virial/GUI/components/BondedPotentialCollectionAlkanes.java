package etomica.virial.GUI.components;

import java.util.HashMap;

import etomica.api.ISpecies;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.potential.P3BondAngle;
import etomica.potential.P4BondTorsion;
import etomica.potential.Potential;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;
import etomica.units.Kelvin;
import etomica.virial.MCMoveClusterTorsionMulti;
import etomica.virial.SpeciesAlkane;
import etomica.virial.simulations.SimulationVirialOverlap2;

public class BondedPotentialCollectionAlkanes extends BondedPotentialCollectionFactory{
	
	public int nSpheres;

	
	public BondedPotentialCollectionAlkanes(int NSpheres){
		nSpheres = NSpheres;
	}
	
	
	@SuppressWarnings("unchecked")
	public void AddBondedPotentialSets(PotentialCollectionFactory pObject, ISpace space, int speciesIndex){
		
		boolean IndexFlag = false;
		if(AtomicPotentialsCollection.class.isAssignableFrom(pObject.getClass())){
			IndexFlag = true;
		}
		
		
		if(nSpheres > 2){
			P3BondAngle P3 = new P3BondAngle(space);
			P3.setAngle(Math.PI*114.0/180.0);
			P3.setEpsilon(Kelvin.UNIT.toSim(62500));
			if(IndexFlag){
				pObject.getBondedPotentialSets(speciesIndex).put(new String[]{"P3-"+Integer.toString(speciesIndex),"P3-"+Integer.toString(speciesIndex)}, P3);}
			else{
				pObject.getBondedPotentialSets().put(new String[]{"P3-"+Integer.toString(speciesIndex),"P3-"+ Integer.toString(speciesIndex)}, P3);
			}
			
		
			int[][] triplets = new int[nSpheres-2][3];
            for (int i=0; i<nSpheres-2; i++) {
                triplets[i][0] = i;
                triplets[i][1] = i+1;
                triplets[i][2] = i+2;
            }
            if(IndexFlag){
            	pObject.getIteratorSets(speciesIndex).put(new String[]{"P3-"+Integer.toString(speciesIndex),
            														"P3-"+Integer.toString(speciesIndex)}, new Atomset3IteratorIndexList(triplets));
            }
            else{
            	pObject.getIteratorSets().put(new String[]{"P3-"+Integer.toString(speciesIndex),
												"P3-"+Integer.toString(speciesIndex)}, new Atomset3IteratorIndexList(triplets));
            }
			
		}
		if (nSpheres > 3) {
	        P4BondTorsion P4 = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
	        if(IndexFlag){
	        	pObject.getBondedPotentialSets(speciesIndex).put(new String[]{"P4-"+Integer.toString(speciesIndex),"P4-"+Integer.toString(speciesIndex)}, P4);}
	        else{
	        	pObject.getBondedPotentialSets().put(new String[]{"P4-"+Integer.toString(speciesIndex),"P4-"+Integer.toString(speciesIndex)}, P4);
	        }
	        int[][] quads = new int[nSpheres-3][4];
	        for (int i=0; i<nSpheres-3; i++) {
	            quads[i][0] = i;
	            quads[i][1] = i+1;
	            quads[i][2] = i+2;
	            quads[i][3] = i+3;
	        }
	        if(IndexFlag){
	        	pObject.getIteratorSets(speciesIndex).put(new String[]{"P4-"+Integer.toString(speciesIndex),"P4-"+Integer.toString(speciesIndex)},new Atomset4IteratorIndexList(quads));
	        }
	        else{
	        	pObject.getIteratorSets().put(new String[]{"P4-"+Integer.toString(speciesIndex),"P4-"+Integer.toString(speciesIndex)},new Atomset4IteratorIndexList(quads));
	        }
	    }

	    if (nSpheres > 4) {
	    	
	    	if(IndexFlag){
	    		pObject.getIteratorSets(speciesIndex).put(new String[]{"CH3","CH3"},new ApiIndexList(new int[][]{{0,nSpheres-1}}));
	    	}
	    	else{
	    		pObject.getIteratorSets().put(new String[]{"CH3","CH3"},new ApiIndexList(new int[][]{{0,nSpheres-1}}));
	    	}
	    }
	    if (nSpheres > 5) {
	        int[][] pairs = new int[2*(nSpheres-5)][2];
	        for (int i=0; i<nSpheres-5; i++) {
	            pairs[2*i][0] = 0;
	            pairs[2*i][1] = nSpheres-2-i;
	            pairs[2*i+1][0] = nSpheres-1;
	            pairs[2*i+1][1] = i+1;
	        }
	        if(IndexFlag){
	        	pObject.getIteratorSets(speciesIndex).put(new String[]{"CH2","CH3"},new ApiIndexList(pairs));
	        }
	        else{
	        	pObject.getIteratorSets().put(new String[]{"CH2","CH3"},new ApiIndexList(pairs));
	        }
	    }
	    
	    if (nSpheres > 6) {
	        int[][] pairs = new int[(nSpheres-6)*(nSpheres-5)/2][2];
	        int k = 0;
	        for (int i=1; i<nSpheres-5; i++) {
	            for (int j=i+4; j<nSpheres-1; j++) {
	                pairs[k][0] = i;
	                pairs[k][1] = j;
	                k++;
	            }
	        }
	        if(IndexFlag){
	        	pObject.getIteratorSets(speciesIndex).put(new String[]{"CH2","CH2"},new ApiIndexList(pairs));
	        }
	        else{
	        	pObject.getIteratorSets().put(new String[]{"CH2","CH2"},new ApiIndexList(pairs));
	        }
	    }
	
	}


}
