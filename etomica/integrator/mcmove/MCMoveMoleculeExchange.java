package etomica.integrator.mcmove;

import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.DataSourceCOM;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.MCMove;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.Species;

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
    private final MeterPotentialEnergy energyMeter;
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private final AtomActionTranslateTo moleculeTranslator;
    private final AtomActionTranslateBy moleculeReplacer;
    private final Vector translationVector;
    
    private transient Atom molecule;
    private transient Phase iPhase, dPhase;
    private transient SpeciesAgent iSpecies, dSpecies;
    private transient double uOld;
    private transient double uNew = Double.NaN;

    public MCMoveMoleculeExchange(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster, 2);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setTunable(false);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(true);
        moleculeReplacer = new AtomActionTranslateBy(space);
        moleculeTranslator = new AtomActionTranslateTo(space);
        translationVector = moleculeTranslator.getTranslationVector();
        setAtomPositionDefinition(new DataSourceCOM(space));
    }
    
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        firstPhase = p[0];
        secondPhase = p[1];
    }
        
    public boolean doTrial() {
        if(Simulation.random.nextBoolean()) {
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
        Species species = molecule.type.getSpecies();
        
        iSpecies = species.getAgent(iPhase);  //insertion-phase speciesAgent
        dSpecies = species.getAgent(dPhase);  //deletion-phase species Agent
        
        energyMeter.setTarget(molecule);
        energyMeter.setPhase(dPhase);
        uOld = energyMeter.getDataAsScalar();

        moleculeTranslator.setDestination(iPhase.randomPosition());         //place at random in insertion phase
        moleculeTranslator.actionPerformed(molecule);
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
        energyMeter.setPhase(iPhase);
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld)/temperature;
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        translationVector.TE(-1);
        moleculeReplacer.setTranslationVector(translationVector);
        moleculeReplacer.actionPerformed(molecule);
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

    
    /**
     * @return Returns the atomPositionDefinition.
     */
    public AtomPositionDefinition getAtomPositionDefinition() {
        return moleculeTranslator.getAtomPositionDefinition();
    }
    /**
     * @param atomPositionDefinition The atomPositionDefinition to set.
     */
    public void setAtomPositionDefinition(
            AtomPositionDefinition atomPositionDefinition) {
        moleculeTranslator.setAtomPositionDefinition(atomPositionDefinition);
    }
}//end of MCMoveMoleculeExchange