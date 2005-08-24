package etomica.virial;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;

/**
 * @author kofke
 *
 * Extension of MCMoveAtom that does trial in which several atom positions are
 * perturbed.  However, position of first atom is never altered.  
 */
public class MCMoveClusterMoleculeMulti extends MCMoveMolecule implements MCMoveCluster {

    private final MeterClusterWeight weightMeter;

	/**
	 * Constructor for MCMoveAtomMulti.
	 * @param parentIntegrator
	 * @param nAtoms number of atoms to move in a trial.  Number of atoms in
	 * phase should be at least one greater than this value (greater
	 * because first atom is never moved)
	 */
	public MCMoveClusterMoleculeMulti(PotentialMaster potentialMaster, int nAtoms) {
		super(potentialMaster);
		this.nAtoms = nAtoms;
		selectedAtoms = new Atom[nAtoms];
        translationVectors = new Vector[nAtoms];
        for (int i=0; i<nAtoms; i++) {
            translationVectors[i] = potential.getSpace().makeVector();
        }
        weightMeter = new MeterClusterWeight(potential);
        setName("MCMoveClusterMolecule");
	}

    public void setPhase(Phase[] p) {
        super.setPhase(p);
        weightMeter.setPhase(p[0]);
    }
    
	//note that total energy is calculated
	public boolean doTrial() {
        if (selectedAtoms[0] == null) selectMolecules();
        uOld = weightMeter.getDataAsScalar();
        for(int i=0; i<selectedAtoms.length; i++) {
            translationVectors[i].setRandomCube();
            translationVectors[i].TE(stepSize);
            groupTranslationVector.E(translationVectors[i]);
            moveMoleculeAction.actionPerformed(selectedAtoms[i]);
        }
        ((PhaseCluster)phases[0]).trialNotify();
        uNew = Double.NaN;
        return true;
	}
	
	// inefficient if selecting most of a large set of atoms
	protected Atom[] selectMolecules() {
        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules(phases[0]);
        if (iterator.size()-1 != nAtoms) throw new IllegalStateException("move should work on number of molecules in phase-1");
        iterator.reset();
        int i=0;
        iterator.next();
        while (iterator.hasNext()) {
            selectedAtoms[i++] = iterator.nextAtom();
        }
		return selectedAtoms;
	}
	
	public void rejectNotify() {
        for(int i=0; i<selectedAtoms.length; i++) {
            groupTranslationVector.E(translationVectors[i]);
            groupTranslationVector.TE(-1);
            moveMoleculeAction.actionPerformed(selectedAtoms[i]);
        }
        ((PhaseCluster)phases[0]).rejectNotify();
	}

    public void acceptNotify() {
        ((PhaseCluster)phases[0]).acceptNotify();
    }
    public double trialRatio() {return 1.0;}
    
    public double lnProbabilityRatio() {
        return Math.log(probabilityRatio());
    }
    
    public double probabilityRatio() {
        uNew = weightMeter.getDataAsScalar();
//        if (Simulation.random.nextInt(200000) == 5) {
//            System.out.println("uOld "+uOld+" uNew "+uNew);
//        }
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }
    
	
	
	private final int nAtoms;
	private final Atom[] selectedAtoms;
    private final Vector[] translationVectors;

}
