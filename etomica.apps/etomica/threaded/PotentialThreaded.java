package etomica.threaded;

import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.potential.IPotential;
import etomica.potential.Potential;
import etomica.space.Space;

public class PotentialThreaded extends Potential {

	final protected IPotential[] potential;
	
	public PotentialThreaded(Space space, IPotential[] potential) {
		super(potential[0].nBody(), space);
		this.potential = potential;
	
	}

	public double energy(AtomSet atoms) {
		//Only the energy from one thread (a partition of atoms)
		return potential[0].energy(atoms);
	}

	public double getRange() {
		
		return potential[0].getRange();
	}

	public void setPhase(Phase phase) {
		
		for(int i=0; i<potential.length; i++){
			potential[i].setPhase(phase);
		}

	}
	
	public IPotential[] getPotentials(){
		return potential;
	}

}
