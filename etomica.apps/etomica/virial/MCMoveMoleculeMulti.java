package etomica.virial;

import etomica.*;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;

/**
 * @author kofke
 *
 * Extension of MCMoveAtom that does trial in which several atom positions are
 * perturbed.  However, position of first atom is never altered.  
 */
public class MCMoveMoleculeMulti extends MCMoveMolecule {

	/**
	 * Constructor for MCMoveAtomMulti.
	 * @param parentIntegrator
	 * @param nAtoms number of atoms to move in a trial.  Number of atoms in
	 * phase should be at least one greater than this value (greater
	 * because first atom is never moved)
	 */
	public MCMoveMoleculeMulti(IntegratorMC parentIntegrator, int nAtoms) {
		super(parentIntegrator);
		this.nAtoms = nAtoms;
		selectedList = new int[nAtoms];
		selectedAtoms = new Atom[nAtoms];
	}

	//note that total energy is calculated
	public boolean doTrial() {
		if(phase.moleculeCount()-1 < nAtoms) return false;
		selectAtoms();
		uOld = potential.calculate(phase, iteratorDirective.set(), energy.reset()).sum();
		for(int i=0; i<selectedAtoms.length; i++) selectedAtoms[i].coord.displaceWithin(stepSize);
		uNew = Double.NaN;
		return true;
	}
	
	// inefficient if selecting most of a large set of atoms
	public Atom[] selectAtoms() {
		resetSelectedList();
		for(int i=0; i<nAtoms; i++) {
			while(selectedAtoms[i] == null || !acceptableAtom(selectedAtoms[i])) selectedAtoms[i] = phase.randomMolecule();		
			selectedList[i] = selectedAtoms[i].node.index();
		}
		return selectedAtoms;
	}
	
	public void rejectNotify() {
		for(int i=0; i<selectedAtoms.length; i++) selectedAtoms[i].coord.replace();
	}

	public double lnProbabilityRatio() {
		uNew = potential.calculate(phase, iteratorDirective.set(), energy.reset()).sum();
//		if(phase.index == 3) System.out.println(phase.index+"  "+parentIntegrator.temperature()+"  "+(-(uNew-uOld)));
		return -(uNew - uOld)/parentIntegrator.temperature();
	}
	
	private boolean acceptableAtom(Atom a) {
		int index = a.node.index();
		int i=0;
		while(true) {
			if(index == selectedList[i]) return false;
			if(selectedList[i++] == 0) return true;//reached end of filled part of list without finding duplicate
		}
	}
	public void resetSelectedList() {
		for(int i=0; i<selectedList.length; i++) {
			selectedAtoms[i] = null;
			selectedList[i] = 0;
		}
	}
	
	
	
	private final int nAtoms;
	private final int[] selectedList;
	private final Atom[] selectedAtoms;
	
	public static void main(String arg[]) {
		
		Simulation sim = new Simulation();
		Species species = new SpeciesSpheres(10,3);
		species.setNMolecules(10);		
		Phase phase = new Phase();
		IntegratorMC integrator = new IntegratorMC();
		MCMoveMoleculeMulti mcmove = new MCMoveMoleculeMulti(integrator, 5);
		sim.elementCoordinator.go();
		
		for(int i=0; i<5; i++) {
			Atom[] atoms = mcmove.selectAtoms();
			for(int j=0; j<atoms.length; j++) System.out.print(atoms[j].node.index());
			System.out.println();
		}
	}//end of main
	/**
	 * @see etomica.integrator.MCMove#adjustStepSize()
	 */
	public void adjustStepSize() {
		super.adjustStepSize();
//		System.out.println("MCMoveAtomMulti nAtoms, stepsize, accept ratio: "+nAtoms+" "+stepSize+" "+acceptanceRatio());
	}

}
