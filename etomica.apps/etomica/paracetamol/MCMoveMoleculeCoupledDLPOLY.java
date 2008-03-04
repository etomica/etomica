package etomica.paracetamol;

import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomMolecule;
import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.PotentialGroup;
import etomica.space.IVectorRandom;


/**
 * Standard Monte Carlo molecule-displacement trial move.  Two molecules are moved at a
 * time in such a way that the geometric center of the system is not changed.
 *
 * @author Nancy Cribbin
 */
public class MCMoveMoleculeCoupledDLPOLY extends MCMoveBoxStep {

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
    
    public MCMoveMoleculeCoupledDLPOLY(IPotentialMaster potentialMaster, IRandom nRandom){
        super(potentialMaster);
        this.random = nRandom;
        moleculeSource = new AtomSourceRandomMolecule();
        ((AtomSourceRandomMolecule)moleculeSource).setRandom(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        
        affectedMoleculeList = new AtomArrayList(2);
        affectedMoleculeIterator = new AtomIteratorArrayListSimple(affectedMoleculeList);
        
        singleAction = new AtomActionTranslateBy(potentialMaster.getSpace());
        groupTransVect = (IVectorRandom)singleAction.getTranslationVector();
        
        moveMoleculeAction = new AtomGroupAction(singleAction);
        
        pair = new AtomPair();
        
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
        uOld = energyMeter.getDataAsScalar();
        
        
        if(uOld > 1e10){
            throw new RuntimeException(new ConfigurationOverlapException(box));
        }
        
        groupTransVect.setRandomCube(random);
        groupTransVect.TE(stepSize);
        moveMoleculeAction.actionPerformed(molecule0);
        groupTransVect.TE(-1.0);
        moveMoleculeAction.actionPerformed(molecule1);
        
        uNew = energyMeter.getDataAsScalar();
        /*
         * Because we have uNew is infinity, and we don't want to have to 
         * worry about the system subtracting infinity from infinity, and
         * setting uNew equal to zero, and accepting the move.
         */
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