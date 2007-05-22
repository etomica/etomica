package etomica.normalmode;

import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomMolecule;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMovePhaseStep;
import etomica.potential.PotentialGroup;
import etomica.potential.PotentialMaster;
import etomica.space.IVectorRandom;
import etomica.util.IRandom;

/**
 * Standard Monte Carlo molecule-displacement trial move.  Two molecules are moved at a
 * time in such a way that the geometric center of the system is not changed.
 *
 * @author Nancy Cribbin
 */
public class MCMoveMoleculeCoupled extends MCMovePhaseStep {

    private static final long serialVersionUID = 1L;
    protected final AtomGroupAction moveMoleculeAction;
    protected final IVectorRandom groupTransVect;
    protected IAtom molecule0, molecule1;
    protected final MeterPotentialEnergy energyMeter;
    protected AtomSource moleculeSource;
    protected double uOld, uNew;
    protected final IRandom random;
    protected final AtomIteratorArrayListSimple affectedMoleculeIterator;
    protected final AtomArrayList affectedMoleculeList;
    protected final AtomActionTranslateBy singleAction;
    protected final AtomPair pair;
    protected PotentialGroup potential;
    
    public MCMoveMoleculeCoupled(PotentialMaster potentialMaster, IRandom nRandom){
        super(potentialMaster);
        this.random = nRandom;
        moleculeSource = new AtomSourceRandomMolecule();
        moleculeSource.setPhase(potentialMaster.getSimulation().getPhases()[0]);
        ((AtomSourceRandomMolecule)moleculeSource).setRandom(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        energyMeter.setPhase(potentialMaster.getSimulation().getPhases()[0]);
        
        setName("McMoveMoleculeCoupled");
        affectedMoleculeList = new AtomArrayList(2);
        affectedMoleculeIterator = new AtomIteratorArrayListSimple(affectedMoleculeList);
        
        setStepSizeMax(2.0* potentialMaster.getSimulation().getDefaults().atomSize);
        setStepSizeMin(0.0);
        setStepSize(potentialMaster.getSimulation().getDefaults().atomSize);
        
        
        singleAction = new AtomActionTranslateBy(potentialMaster.getSpace());
        groupTransVect = (IVectorRandom)singleAction.getTranslationVector();
        
        moveMoleculeAction = new AtomGroupAction(singleAction);
        
        pair = new AtomPair();
        
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
    }
    
    public void setPotential(PotentialGroup newPotential){
        potential = newPotential;
    }
    
    public AtomIterator affectedAtoms() {
        affectedMoleculeList.clear();
        affectedMoleculeList.add(molecule0);
        affectedMoleculeList.add(molecule1);
        return affectedMoleculeIterator;
    }

    public double energyChange() {return uNew - uOld;}

    public void acceptNotify() {
        // I do believe nothing needs to happen here.
    }

    public boolean doTrial() {
//        System.out.println("doTrial MCMoveMoleculeCoupled called");
        
        molecule0 = moleculeSource.getAtom();
        molecule1 = moleculeSource.getAtom();
        if(molecule0==null || molecule1==null || molecule0==molecule1) return false;
        
        //make sure we don't double count the molecule0-molecule1 interaction
        pair.atom0 = molecule0;
        pair.atom1 = molecule1;
        energyMeter.setTarget(molecule0);
        uOld = energyMeter.getDataAsScalar();
        energyMeter.setTarget(molecule1);
        uOld += energyMeter.getDataAsScalar();
        uOld -= potential.energy(pair);
        
        
        if(uOld > 1e10){
            throw new RuntimeException(new ConfigurationOverlapException(phase));
        }
        
        groupTransVect.setRandomCube(random);
        groupTransVect.TE(stepSize);
        moveMoleculeAction.actionPerformed(molecule0);
        groupTransVect.TE(-1.0);
        moveMoleculeAction.actionPerformed(molecule1);
        
        energyMeter.setTarget(molecule0);
        uNew = energyMeter.getDataAsScalar();
        energyMeter.setTarget(molecule1);
        uNew += energyMeter.getDataAsScalar();
        /*
         * Because we have uNew is infinity, and we don't want to have to 
         * worry about the system subtracting infinity from infinity, and
         * setting uNew equal to zero, and accepting the move.
         */
        if(Double.isInfinite(uNew)) {return true;}  
        uNew -= potential.energy(pair);
        
        return true;
    }

    public double getA() {
        return 1.0;
    }

    public double getB() {
        return -(uNew - uOld);
    }

    public void rejectNotify() {
        moveMoleculeAction.actionPerformed(molecule0);
        groupTransVect.TE(-1.0);
        moveMoleculeAction.actionPerformed(molecule1);
    }
    
}