package etomica.threaded;

import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IPotential;

import etomica.potential.Potential;
import etomica.space.Space;

public class PotentialThreaded extends Potential {

	final protected IPotential[] potential;
	
	public PotentialThreaded(Space space, IPotential[] potential) {
		super(potential[0].nBody(), space);
		this.potential = potential;
	
	}

	public double energy(IAtomSet atoms) {
		//Only the energy from one thread (a partition of atoms)
		return potential[0].energy(atoms);
	}

	public double getRange() {
		
		return potential[0].getRange();
	}

	public void setBox(IBox box) {
		
		for(int i=0; i<potential.length; i++){
			potential[i].setBox(box);
		}

	}
	
	public IPotential[] getPotentials(){
		return potential;
	}

}
