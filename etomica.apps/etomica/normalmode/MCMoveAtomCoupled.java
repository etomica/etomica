package etomica.normalmode;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.AtomSetSinglet;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

/**
 * Standard Monte Carlo atom-displacement trial move.  Two atoms are moved at a
 * time in such a way that the geometric center of the system is not changed.
 *
 * @author Andrew Schultz
 */
public class MCMoveAtomCoupled extends MCMoveBoxStep {
    
    private static final long serialVersionUID = 2L;
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final AtomArrayList affectedAtomList;
    protected final MeterPotentialEnergy energyMeter;
    protected final IVectorRandom translationVector;
    protected IAtom atom0, atom1;
    protected double uOld;
    protected double uNew;
    protected AtomSource atomSource;
    protected final IRandom random;
    protected IPotentialAtomic pairPotential;
    protected final AtomPair pair;
    protected final AtomSetSinglet atomSinglet;
    protected boolean doExcludeNonNeighbors, doIncludePair;
    protected IPotentialAtomic constraintPotential;

    public MCMoveAtomCoupled(IPotentialMaster potentialMaster, MeterPotentialEnergy energyMeter,
                             IRandom random, ISpace _space) {
        super(potentialMaster);
        this.random = random;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        this.energyMeter = energyMeter;
        translationVector = (IVectorRandom)_space.makeVector();
        setStepSizeMax(0.5);
        setStepSizeMin(0.0);
        setStepSize(0.1);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        affectedAtomList = new AtomArrayList(2);
        affectedAtomIterator = new AtomIteratorArrayListSimple(affectedAtomList);
        pair = new AtomPair();
        atomSinglet = new AtomSetSinglet();
    }
    
    public void setPotential(IPotentialAtomic newPotential) {
        if (newPotential.nBody() != 2) {
            throw new RuntimeException("must be a 2-body potential");
        }
        pairPotential = newPotential;
    }
    
    public void setConstraint(IPotentialAtomic newConstraintPotential) {
        if (newConstraintPotential.nBody() != 1) {
            throw new RuntimeException("must be a 1-body potential");
        }
        constraintPotential = newConstraintPotential;
    }

    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom0 = atomSource.getAtom();
        if (atom0 == null) return false;
        atom1 = atomSource.getAtom();
        if (atom1 == null || atom0 == atom1) return false;
        energyMeter.setTarget(atom0);
        uOld = energyMeter.getDataAsScalar();
        energyMeter.setTarget(atom1);
        uOld += energyMeter.getDataAsScalar();
        if(uOld > 1e10) {
            throw new ConfigurationOverlapException(box);
        }
        pair.atom0 = atom0;
        pair.atom1 = atom1;
        doIncludePair = true;
        if (doExcludeNonNeighbors && potential instanceof PotentialMasterList) {
            doIncludePair = false;
            IAtomList[] list0 = ((PotentialMasterList)potential).getNeighborManager(box).getDownList(atom0);
            for (int i=0; i<list0.length; i++) {
                if (((AtomArrayList)list0[i]).indexOf(atom1)>-1) {
                    doIncludePair = true;
                    break;
                }
            }
            if (!doIncludePair) {
                list0 = ((PotentialMasterList)potential).getNeighborManager(box).getUpList(atom0);
                for (int i=0; i<list0.length; i++) {
                    if (((AtomArrayList)list0[i]).indexOf(atom1)>-1) {
                        doIncludePair = true;
                        break;
                    }
                }
            }
        }
        if (doIncludePair) {
            uOld -= pairPotential.energy(pair);
        }

        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        atom0.getPosition().PE(translationVector);
        atom1.getPosition().ME(translationVector);

        return true;
    }

    public double getA() {return 1.0;}
    
    public double getB() {
        uNew = energyMeter.getDataAsScalar();
        energyMeter.setTarget(atom0);
        uNew += energyMeter.getDataAsScalar();
        if(!Double.isInfinite(uNew) && doIncludePair){
            uNew -= pairPotential.energy(pair);
        }
        if (constraintPotential != null) {
            // we could be double-counting the "energy" here (if the atoms are
            // neighbors), but the energy can only be 0 or +inf, so it's OK.
            atomSinglet.atom = atom0;
            constraintPotential.setBox(box);
            uNew += constraintPotential.energy(atomSinglet);
            atomSinglet.atom = atom1;
            uNew += constraintPotential.energy(atomSinglet);
        }

        return -(uNew - uOld);
    }
    
    public double energyChange() {return uNew - uOld;}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
    }
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.
     */
    public void rejectNotify() {
        atom0.getPosition().ME(translationVector);
        atom1.getPosition().PE(translationVector);
    }
        
    
    public AtomIterator affectedAtoms() {
        affectedAtomList.clear();
        affectedAtomList.add(atom0);
        affectedAtomList.add(atom1);
        return affectedAtomIterator;
    }
    
    public void setBox(IBox p) {
        super.setBox(p);
        energyMeter.setBox(p);
        atomSource.setBox(p);
        pairPotential.setBox(p);
    }
    
    /**
     * @return Returns the atomSource.
     */
    public AtomSource getAtomSource() {
        return atomSource;
    }
    /**
     * @param atomSource The atomSource to set.
     */
    public void setAtomSource(AtomSource source) {
        atomSource = source;
    }

    /**
     * Configures the move to not explicitly calculate the potential for atom
     * pairs that are not neighbors (as determined by the potentialMaster).
     * 
     * Setting this has an effect only if the potentialMaster is an instance of
     * PotentialMasterList
     */
    public void setDoExcludeNonNeighbors(boolean newDoExcludeNonNeighbors) {
        doExcludeNonNeighbors = newDoExcludeNonNeighbors;
    }

    /**
     * Returns true if the move does not explicitly calculate the potential for
     * atom pairs that are not neighbors (as determined by the potentialMaster).
     */
    public boolean getDoExcludeNonNeighbors() {
        return doExcludeNonNeighbors;
    }
}