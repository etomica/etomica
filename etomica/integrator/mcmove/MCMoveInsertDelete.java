package etomica.integrator.mcmove;

import etomica.action.MoleculeActionTranslateTo;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISpecies;
import etomica.atom.MoleculeArrayList;
import etomica.atom.MoleculeSetSinglet;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.space.ISpace;
import etomica.util.Debug;

/**
 * Elementary Monte Carlo move in which a molecule of a specified species is
 * inserted into or removed from a box.
 *
 * @author David Kofke
 */
public class MCMoveInsertDelete extends MCMoveBox {
    
    private static final long serialVersionUID = 2L;
    //chemical potential
    protected double mu;
    
    //directive must specify "BOTH" to get energy with all atom pairs
    protected final MeterPotentialEnergy energyMeter;
	protected ISpecies species;
	protected final AtomIteratorArrayListSimple affectedAtomIterator = new AtomIteratorArrayListSimple();
	protected IMolecule testMolecule;
	protected double uOld;
	protected double uNew = Double.NaN;
	protected boolean insert;
	protected final MoleculeArrayList reservoir;
    protected final MoleculeActionTranslateTo atomTranslator;
    protected IMoleculeList moleculeList;
    protected IRandom random;

    public MCMoveInsertDelete(IPotentialMaster potentialMaster, IRandom random,
    		                  ISpace _space) {
        super(potentialMaster);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setMu(0.0);
        energyMeter.setIncludeLrc(true);
        atomTranslator = new MoleculeActionTranslateTo(_space);
        reservoir = new MoleculeArrayList();
        this.random = random;
    }
    
//perhaps should have a way to ensure that two instances of this class aren't assigned the same species
    public void setSpecies(ISpecies s) {
        species = s;
        if(box != null) {
            moleculeList = box.getMoleculeList(species);
        }
    }
    public ISpecies getSpecies() {return species;}
    
    public void setBox(IBox p) {
        super.setBox(p);
        energyMeter.setBox(box);
        if(species != null) {
            moleculeList = box.getMoleculeList(species);
        }
    }
    
    /**
     * Chooses and performs with equal probability an elementary molecule insertion
     * or deletion.
     */
    public boolean doTrial() {
        insert = (random.nextInt(2) == 0);
        if(insert) {
            uOld = 0.0;
            
            if(!reservoir.isEmpty()) testMolecule = (IMolecule)reservoir.remove(reservoir.getMoleculeCount()-1);
            else testMolecule = species.makeMolecule();
            box.addMolecule(testMolecule);

            atomTranslator.setDestination(box.getBoundary().randomPosition());
            atomTranslator.actionPerformed(testMolecule);
        } else {//delete
            if(box.getNMolecules(species) == 0) {
                testMolecule = null;
                return false;
            }
            testMolecule = moleculeList.getMolecule(random.nextInt(moleculeList.getMoleculeCount()));
            //delete molecule only upon accepting trial
            energyMeter.setTarget(testMolecule);
            uOld = energyMeter.getDataAsScalar();
        } 
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double getA() {//note that moleculeCount() gives the number of molecules after the trial is attempted
        return insert ? box.getBoundary().volume()/box.getNMolecules(species) 
                      : (box.getNMolecules(species)+1)/box.getBoundary().volume();        
    }
    
    public double getB() {
        if(insert) {
            energyMeter.setTarget(testMolecule);
            uNew = energyMeter.getDataAsScalar();
            return (+mu - uNew);
        }
        uNew = 0.0;
        return (-mu + uOld);
    }
    
    public void acceptNotify() {
        if (insert && Debug.ON && Debug.DEBUG_NOW && Debug.allAtoms(new MoleculeSetSinglet(testMolecule))) System.out.println("accepted insertion of "+testMolecule);
        if(!insert) {
            // accepted deletion - remove from box and add to reservoir
            box.removeMolecule(testMolecule);
            reservoir.add(testMolecule);
        }
    }
    
    public void rejectNotify() {
        if(insert) {
            // rejected insertion - remove from box and return to reservoir
            box.removeMolecule(testMolecule);
            reservoir.add(testMolecule);
            // test molecule is no longer in the simulation and should not be 
            // returned by affectedAtoms
            testMolecule = null;
        }
    }
    
    public double energyChange() {return uNew - uOld;}

    /**
     * Returns an iterator giving molecule that is being added or deleted 
     * in the current or most recent trial.
     */
    public final AtomIterator affectedAtoms() {
        if(testMolecule == null) return AtomIteratorNull.INSTANCE;
        affectedAtomIterator.setList(testMolecule.getChildList());
        return affectedAtomIterator;
    }

    /**
     * Mutator method for the chemical potential of the insertion/deletion species.
     */
    public final void setMu(double mu) {this.mu = mu;}
    
    /**
     * Accessor method for the chemical potential of th insertion/deletion species.
     */
    public final double getMu() {return mu;}
    
    /**
     * Indicates that chemical potential has dimensions of energy.
     */
    public final etomica.units.Dimension getMuDimension() {return etomica.units.Energy.DIMENSION;}
  
}//end of MCMoveInsertDelete