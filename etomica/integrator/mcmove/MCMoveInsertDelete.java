package etomica.integrator.mcmove;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.atom.AtomFactory;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.MCMove;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.species.Species;

/**
 * Elementary Monte Carlo move in which a molecule of a specified species is
 * inserted into or removed from a phase.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 07/09/02 (DAK) Added energyChange() method
  * 09/19/02 (DAK) Minor change in doTrial for case were deleting with N = 0
  */
public class MCMoveInsertDelete extends MCMove {
    
    //chemical potential
    protected double mu;
    
    //directive must specify "BOTH" to get energy with all atom pairs
    protected final MeterPotentialEnergy energyMeter;
	protected Species species;
	protected AtomTreeNodeGroup speciesAgentNode;
    protected SpeciesAgent speciesAgent;
	protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
	protected Atom testMolecule;
	protected double uOld;
	protected double uNew = Double.NaN;
	protected boolean insert;
	protected final AtomList reservoir;
    protected final AtomActionTranslateTo atomTranslator;
    protected AtomFactory moleculeFactory;

    public MCMoveInsertDelete(PotentialMaster potentialMaster) {
        super(potentialMaster, 1);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setMu(0.0);
        setTunable(false);
        energyMeter.setIncludeLrc(true);
        atomTranslator = new AtomActionTranslateTo(potentialMaster.getSpace());
        reservoir = new AtomList();
    }
    
//perhaps should have a way to ensure that two instances of this class aren't assigned the same species
    public void setSpecies(Species s) {
        species = s;
        if(phases[0] != null) {
            speciesAgent = species.getAgent(phases[0]);
            speciesAgentNode = (AtomTreeNodeGroup)speciesAgent.node; 
        }
        moleculeFactory = species.moleculeFactory();
    }
    public Species getSpecies() {return species;}
    
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        energyMeter.setPhase(p[0]);
        if(species != null) {
            speciesAgent = species.getAgent(phases[0]);
            speciesAgentNode = (AtomTreeNodeGroup)speciesAgent.node; 
        }
    }
    
    /**
     * Chooses and performs with equal probability an elementary molecule insertion
     * or deletion.
     */
    public boolean doTrial() {
        insert = Simulation.random.nextBoolean();
        if(insert) {
            uOld = 0.0;
            
            if(!reservoir.isEmpty()) testMolecule = reservoir.removeFirst();
            else testMolecule = moleculeFactory.makeAtom();
            phases[0].addMolecule(testMolecule, speciesAgent);

            atomTranslator.setDestination(phases[0].randomPosition());
            atomTranslator.actionPerformed(testMolecule);
        } else {//delete
            if(speciesAgent.moleculeCount() == 0) {
                testMolecule = null;
                return false;
            }
            testMolecule = speciesAgent.randomMolecule();
            //delete molecule only upon accepting trial
            energyMeter.setTarget(testMolecule);
            uOld = energyMeter.getDataAsScalar();
        } 
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {//note that moleculeCount() gives the number of molecules after the trial is attempted
        return insert ? Math.log(phases[0].volume()/speciesAgent.moleculeCount()) 
                      : Math.log((speciesAgent.moleculeCount()+1)/phases[0].volume());        
    }
    
    public double lnProbabilityRatio() {
        if(insert) {
            energyMeter.setTarget(testMolecule);
            uNew = energyMeter.getDataAsScalar();
            return (+mu - uNew)/temperature;
        }
        uNew = 0.0;
        return (-mu + uOld)/temperature;
    }
    
    public void acceptNotify() {
        //      accepted deletion - remove from phase and add to reservoir 
        if(!insert) {
            phases[0].removeMolecule(testMolecule);
            reservoir.add(testMolecule.seq);
        }
    }
    
    public void rejectNotify() {
        //      rejected insertion - remove from phase and return to reservoir
        if(insert) {
            phases[0].removeMolecule(testMolecule);
            reservoir.add(testMolecule.seq);
        }
    }
    
    public double energyChange(Phase phase) {return (this.phases[0] == phase) ? uNew - uOld : 0.0;}

    /**
     * Returns an iterator giving molecule that is being added or deleted 
     * in the current or most recent trial.
     */
    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.phases[0] != phase) return AtomIterator.NULL;
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
    public final etomica.units.Dimension getMuDimension() {return etomica.units.Dimension.ENERGY;}
  
}//end of MCMoveInsertDelete