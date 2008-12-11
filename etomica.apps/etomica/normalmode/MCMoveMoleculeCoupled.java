package etomica.normalmode;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.atom.AtomArrayList;
import etomica.atom.MoleculePair;
import etomica.atom.MoleculeSource;
import etomica.atom.MoleculeSourceRandomMolecule;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

/**
 * Standard Monte Carlo molecule-displacement trial move.  Two molecules are moved at a
 * time in such a way that the geometric center of the system is not changed.
 *
 * @author Nancy Cribbin
 */
public class MCMoveMoleculeCoupled extends MCMoveBoxStep {

    private static final long serialVersionUID = 1L;
    protected final MoleculeChildAtomAction moveMoleculeAction;
    protected final IVectorRandom groupTransVect;
    protected IMolecule molecule0, molecule1;
    protected final MeterPotentialEnergy energyMeter;
    protected MoleculeSource moleculeSource;
    protected double uOld, uNew;
    protected final IRandom random;
    protected final AtomIteratorArrayListSimple affectedMoleculeIterator;
    protected final AtomArrayList affectedMoleculeList;
    protected final AtomActionTranslateBy singleAction;
    protected final MoleculePair pair;
    protected PotentialGroup potential;
    
    public MCMoveMoleculeCoupled(IPotentialMaster potentialMaster, IRandom nRandom,
    		                     ISpace _space){
        super(potentialMaster);
        this.random = nRandom;
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule)moleculeSource).setRandom(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        
        affectedMoleculeList = new AtomArrayList();
        affectedMoleculeIterator = new AtomIteratorArrayListSimple(affectedMoleculeList);
        
        singleAction = new AtomActionTranslateBy(_space);
        groupTransVect = (IVectorRandom)singleAction.getTranslationVector();
        
        moveMoleculeAction = new MoleculeChildAtomAction(singleAction);
        
        pair = new MoleculePair();
        
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
    }

    public void setBox(IBox newBox) {
        super.setBox(newBox);
        moleculeSource.setBox(newBox);
        energyMeter.setBox(newBox);
    }
    
    public void setPotential(PotentialGroup newPotential){
        potential = newPotential;
    }
    
    public AtomIterator affectedAtoms() {
        affectedMoleculeList.clear();
        affectedMoleculeList.addAll(molecule0.getChildList());
        affectedMoleculeList.addAll(molecule1.getChildList());
        return affectedMoleculeIterator;
    }

    public double energyChange() {return uNew - uOld;}

    public void acceptNotify() {
        // I do believe nothing needs to happen here.
    }

    public boolean doTrial() {
//        System.out.println("doTrial MCMoveMoleculeCoupled called");
        
        molecule0 = moleculeSource.getMolecule();
        molecule1 = moleculeSource.getMolecule();
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
            throw new RuntimeException(new ConfigurationOverlapException(box));
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