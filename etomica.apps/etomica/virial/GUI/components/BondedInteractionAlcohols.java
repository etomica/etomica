package etomica.virial.GUI.components;

import etomica.api.ISpecies;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.potential.P3BondAngle;
import etomica.potential.P4BondTorsion;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;
import etomica.units.Kelvin;
import etomica.virial.MCMoveClusterTorsionMulti;
import etomica.virial.simulations.SimulationVirialOverlap2;

public class BondedInteractionAlcohols extends BondedInteraction{
	
	
	@SuppressWarnings("unchecked")
	@Override
	public void AddBondedPotentialSets(PotentialObject pObject, ISpace space, int speciesIndex){
		// TODO Auto-generated method stub
		P3BondAngle P3 = new P3BondAngle(space);
        P3.setAngle(108.5*Math.PI/180);
        P3.setEpsilon(Kelvin.UNIT.toSim(55400));
        pObject.getBondedIIPotentialSetsSpeciesI(speciesIndex - 1).put(new String[]{"P3-"+Integer.toString(speciesIndex),"P3-"+Integer.toString(speciesIndex)}, P3);
        
        int[][] triplets = new int[][] {{0,1,2}};
        pObject.getIteratorSetsSpeciesI(speciesIndex - 1).put(new String[]{"P3-"+Integer.toString(speciesIndex),"P3-"+Integer.toString(speciesIndex)}, new Atomset3IteratorIndexList(triplets));
	}

}
