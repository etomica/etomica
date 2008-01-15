package etomica.integrator.mcmove;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomFactory;
import etomica.atom.AtomSet;
import etomica.atom.AtomSetSinglet;
import etomica.atom.IMolecule;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.potential.PotentialMaster;
import etomica.species.Species;
import etomica.util.Debug;
import etomica.util.IRandom;

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
	protected Species species;
	protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
	protected IMolecule testMolecule;
	protected double uOld;
	protected double uNew = Double.NaN;
	protected boolean insert;
	protected final AtomArrayList reservoir;
    protected final AtomActionTranslateTo atomTranslator;
    protected AtomFactory moleculeFactory;
    protected AtomSet moleculeList;
    protected IRandom random;

    public MCMoveInsertDelete(PotentialMaster potentialMaster, IRandom random) {
        super(potentialMaster);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setMu(0.0);
        energyMeter.setIncludeLrc(true);
        atomTranslator = new AtomActionTranslateTo(potentialMaster.getSpace());
        reservoir = new AtomArrayList();
        this.random = random;
    }
    
//perhaps should have a way to ensure that two instances of this class aren't assigned the same species
    public void setSpecies(Species s) {
        species = s;
        if(box != null) {
            moleculeList = box.getMoleculeList(species);
        }
        moleculeFactory = species.getMoleculeFactory();
    }
    public Species getSpecies() {return species;}
    
    public void setBox(Box p) {
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
            
            if(!reservoir.isEmpty()) testMolecule = (IMolecule)reservoir.remove(reservoir.getAtomCount()-1);
            else testMolecule = (IMolecule)moleculeFactory.makeAtom();
            box.addMolecule(testMolecule);

            atomTranslator.setDestination(box.getBoundary().randomPosition());
            atomTranslator.actionPerformed(testMolecule);
        } else {//delete
            if(box.getNMolecules(species) == 0) {
                testMolecule = null;
                return false;
            }
            testMolecule = (IMolecule)moleculeList.getAtom(random.nextInt(moleculeList.getAtomCount()));
            //delete molecule only upon accepting trial
            energyMeter.setTarget(testMolecule);
            uOld = energyMeter.getDataAsScalar();
        } 
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double getA() {//note that moleculeCount() gives the number of molecules after the trial is attempted
        return insert ? box.volume()/box.getNMolecules(species) 
                      : (box.getNMolecules(species)+1)/box.volume();        
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
        if (insert && Debug.ON && Debug.DEBUG_NOW && Debug.allAtoms(new AtomSetSinglet(testMolecule))) System.out.println("accepted insertion of "+testMolecule);
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
        affectedAtomIterator.setAtom(testMolecule);
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