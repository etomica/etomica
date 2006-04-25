package etomica.integrator.mcmove;

import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.atom.AtomPositionCOM;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorPhase;
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
    private final IntegratorPhase integrator1, integrator2;
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
    

    public MCMoveMoleculeExchange(PotentialMaster potentialMaster, Space space,
            IntegratorPhase integrator1, IntegratorPhase integrator2) {
        super(potentialMaster);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        energyMeter.setIncludeLrc(true);
        moleculeReplacer = new AtomActionTranslateBy(space);
        moleculeTranslator = new AtomActionTranslateTo(space);
        translationVector = moleculeTranslator.getTranslationVector();
        setAtomPositionDefinition(new AtomPositionCOM(space));
        this.integrator1 = integrator1;
        this.integrator2 = integrator2;
        firstPhase = integrator1.getPhase();
        secondPhase = integrator2.getPhase();
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
        molecule.node.setParent((AtomTreeNodeGroup)iSpecies.node);
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double getA() {
        //note that dSpecies.nMolecules has been decremented
        //and iSpecies.nMolecules has been incremented
        return (dSpecies.getNMolecules()+1)/dPhase.volume()
               * iPhase.volume()/iSpecies.getNMolecules(); 
    }
    
    public double getB() {
        energyMeter.setPhase(iPhase);
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld);
    }
    
    public void acceptNotify() {
        try {
            integrator1.reset();
            integrator2.reset();
        } catch(ConfigurationOverlapException e) {
            throw new RuntimeException(e);
        }
    }
    
    public void rejectNotify() {
        translationVector.TE(-1);
        moleculeReplacer.setTranslationVector(translationVector);
        moleculeReplacer.actionPerformed(molecule);
        molecule.node.setParent((AtomTreeNodeGroup)dSpecies.node);
    }

    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.firstPhase != phase && this.secondPhase != phase) return AtomIteratorNull.INSTANCE;
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