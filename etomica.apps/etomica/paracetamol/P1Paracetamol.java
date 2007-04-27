package etomica.paracetamol;
   
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2Dreiding;
import etomica.potential.Potential1;
import etomica.potential.PotentialGroup;
import etomica.space.Space;

public class P1Paracetamol extends Potential1 {


	public P1Paracetamol(Space space, PotentialMasterList potentialMaster) {
		super(space);
		
		PotentialGroup intramolecularPotential = potentialMaster.makePotentialGroup(1);
		
		/*
		 *  Bond Stretch Energy
		 */
		P2Dreiding potentialC4O1 = new P2Dreiding(space, 1.352, .17612581e6);
		intramolecularPotential.addPotential(potentialC4O1, 
				new ApiIndexList(new int[][]{{AtomParacetamol.indexC4, AtomParacetamol.indexO1}}));
		
		//bond length [amstrom]; parameter [K]
		
		potentialC1C2 = new P2Dreiding(space, 1.395, 0.264188175e6);
		potentialC2C3 = new P2Dreiding(space, 1.385, 0.264188175e6);
		potentialC1N1 = new P2Dreiding(space, 1.394, 0.17612581e6);
		potentialN1C7 = new P2Dreiding(space, 1.366, 0.17612581e6);
		potentialC7O2 = new P2Dreiding(space, 1.226, 0.35225162e6);
		potentialC7C8 = new P2Dreiding(space, 1.503, 0.17612581e6);
		
		//potentialMaster.addPotential(intramolecularPotential, new Species[] {speciesParacetamol});
		atomPair = new AtomPair();

	}

	@Override
	public double energy(AtomSet moleculeParacetamol) {
		
		AtomParacetamol atomNode = (AtomParacetamol)moleculeParacetamol;
		double sumEnergy = 0;
		
		atomPair.atom0 = atomNode.C4;
		atomPair.atom1 = atomNode.O1;
		sumEnergy += potentialC4O1.energy(atomPair);
		
		atomPair.atom0 = atomNode.C1;
		atomPair.atom1 = atomNode.C2;
		sumEnergy += potentialC1C2.energy(atomPair);
		
		atomPair.atom0 = atomNode.C3;
		atomPair.atom1 = atomNode.C4;
		sumEnergy += potentialC1C2.energy(atomPair);
		
		atomPair.atom0 = atomNode.C4;
		atomPair.atom1 = atomNode.C5;
		sumEnergy += potentialC1C2.energy(atomPair);
		
		atomPair.atom0 = atomNode.C6;
		atomPair.atom1 = atomNode.C1;
		sumEnergy += potentialC1C2.energy(atomPair);
		
		atomPair.atom0 = atomNode.C2;
		atomPair.atom1 = atomNode.C3;
		sumEnergy += potentialC2C3.energy(atomPair);
		
		atomPair.atom0 = atomNode.C5;
		atomPair.atom1 = atomNode.C6;
		sumEnergy += potentialC2C3.energy(atomPair);
		
		atomPair.atom0 = atomNode.C1;
		atomPair.atom1 = atomNode.N ;
		sumEnergy += potentialC1N1.energy(atomPair);
		
		atomPair.atom0 = atomNode.N;
		atomPair.atom1 = atomNode.C7;
		sumEnergy += potentialN1C7.energy(atomPair);
		
		atomPair.atom0 = atomNode.C7;
		atomPair.atom1 = atomNode.O2;
		sumEnergy += potentialC7O2.energy(atomPair);
		
		atomPair.atom0 = atomNode.C7;
		atomPair.atom1 = atomNode.C8;
		sumEnergy += potentialC7C8.energy(atomPair);
		
		return sumEnergy;
	}

	
	
	private AtomPair atomPair;
	
	private P2Dreiding potentialC4O1;
	private P2Dreiding potentialC1C2;
	private P2Dreiding potentialC2C3;
	private P2Dreiding potentialC1N1;
	private P2Dreiding potentialN1C7;
	private P2Dreiding potentialC7O2;
	private P2Dreiding potentialC7C8;	


	private static final long serialVersionUID = 1L;
}
