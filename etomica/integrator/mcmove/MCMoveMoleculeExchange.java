package etomica.integrator.mcmove;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesAgent;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.MCMove;

/**
 * Performs a trial that results in the exchange of a molecule from one phase to another.
 * Primary use is as an elementary move in a Gibbs ensemble simulation
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 7/9/02 Added energyChange() method
  */
public final class MCMoveMoleculeExchange extends MCMove {
    
    private Phase firstPhase;
    private Phase secondPhase;
    private final double ROOT;
    private final MeterPotentialEnergy energyMeter;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();

    private transient Atom molecule;
    private transient Phase iPhase, dPhase;
    private transient SpeciesAgent iSpecies, dSpecies;
    private transient double uOld;
    private transient double uNew = Double.NaN;

    public MCMoveMoleculeExchange(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster, 2);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        ROOT = 1.0/space.D();
        setTunable(false);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(true);
    }
    
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        firstPhase = p[0];
        secondPhase = p[1];
    }
        
    public boolean doTrial() {
        if(Simulation.random.nextDouble() < 0.5) {
            iPhase = firstPhase;
            dPhase = secondPhase;
        }
        else {
            iPhase = secondPhase;
            dPhase = firstPhase;
        }
        if(dPhase.moleculeCount() == 0) { //no molecules to delete; trial is over
            uNew = uOld = 0.0;
            return false;
        }
        
        molecule = dPhase.randomMolecule();  //select random molecule to delete
        Species species = molecule.node.parentSpecies();
        
        iSpecies = species.getAgent(iPhase);  //insertion-phase speciesAgent
        dSpecies = species.getAgent(dPhase);  //deletion-phase species Agent
        
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar(dPhase);

        molecule.coord.displaceTo(iPhase.randomPosition());         //place at random in insertion phase
        molecule.node.setParent(iSpecies);
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {
        //note that dSpecies.nMolecules has been decremented
        //and iSpecies.nMolecules has been incremented
        return Math.log( (dSpecies.moleculeCount()+1)/dPhase.volume()
                         * iPhase.volume()/iSpecies.moleculeCount() ); 
    }
    
    public double lnProbabilityRatio() {
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar(iPhase);
        return -(uNew - uOld)/temperature;
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        molecule.coord.replace();
        molecule.node.setParent(dSpecies);
    }

    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.firstPhase != phase && this.secondPhase != phase) return AtomIterator.NULL;
        affectedAtomIterator.setAtom(molecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }
    
    public double energyChange(Phase phase) {
        if(phase == iPhase) return uNew;
        else if(phase == dPhase) return -uOld;
        else return 0.0;
    }

}//end of MCMoveMoleculeExchange