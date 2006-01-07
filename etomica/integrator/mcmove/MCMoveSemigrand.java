package etomica.integrator.mcmove;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomList;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorCompound;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.DataSourceCOM;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.MCMove;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.species.Species;

/**
 * Basic Monte Carlo move for semigrand-ensemble simulations.  Move consists
 * of selecting a molecule at random and changing its species identity.  More precisely,
 * the molecule is removed and another molecule of a different species replaces it.
 * An arbitrary number of species may be designated as subject to these exchange moves.
 * Acceptance is regulated by a set of fugacity fractions that are specified at design time.
 *
 * @author Jhumpa Adhikari
 * @author David Kofke
 */
 
 /* History of changes
  * 7/9/02 Added energyChange() method
  */
  
public class MCMoveSemigrand extends MCMove {
    
    private Species[] speciesSet;
    private SpeciesAgent[] agentSet;
    private AtomArrayList[] reservoirs;
    private double[] fugacityFraction;
    private int nSpecies;
    private final AtomIteratorSinglet deleteAtomIterator;
    private final AtomIteratorSinglet insertAtomIterator;
    private final AtomIteratorCompound affectedAtomIterator; 
    private final MeterPotentialEnergy energyMeter;
    private final AtomActionTranslateTo moleculeTranslator;
    private AtomPositionDefinition atomPositionDefinition;
    
    private transient Atom deleteMolecule, insertMolecule;
    private transient double uOld;
    private transient double uNew = Double.NaN;
    private transient SpeciesAgent deleteAgent, insertAgent;
    private transient int iInsert, iDelete;

    
    public MCMoveSemigrand(PotentialMaster potentialMaster) {
        super(potentialMaster, 1);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        deleteAtomIterator = new AtomIteratorSinglet();
        insertAtomIterator = new AtomIteratorSinglet();
        affectedAtomIterator = new AtomIteratorCompound(new AtomIterator[] {deleteAtomIterator, insertAtomIterator});
        setTunable(false);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(true);
        moleculeTranslator = new AtomActionTranslateTo(potentialMaster.getSpace());
        setAtomPositionDefinition(new DataSourceCOM(potentialMaster.getSpace()));
    }
    
    /**
     * Extends the superclass method to initialize the exchange-set species agents for the phase.
     */
    public void setPhase(Phase p) {
        super.setPhase(p);
        energyMeter.setPhase(phase);
        if(speciesSet != null) {
            for(int i=0; i<nSpecies; i++) {
                agentSet[i] = speciesSet[i].getAgent(phase);
            }
        }
    }//end setPhase
    
    /**
     * Mutator method for the set of species that can participate in an exchange move.
     */
    public void setSpecies(Species[] species) {
        nSpecies = species.length;
        if(nSpecies < 2) throw new IllegalArgumentException("Wrong size of species array in MCMoveSemigrand");
        speciesSet = new Species[nSpecies];
        agentSet = new SpeciesAgent[nSpecies];
        fugacityFraction = new double[nSpecies];
        reservoirs = new AtomArrayList[nSpecies];
        for(int i=0; i<nSpecies; i++) {
            speciesSet[i] = species[i];
            if(phase != null) agentSet[i] = species[i].getAgent(phase);
            fugacityFraction[i] = 1.0/nSpecies;
            reservoirs[i] = new AtomArrayList();
        }
    }
    
    /**
     * Accessor method for the set of species that can participate in an exchange move.
     */
    public Species[] getSpecies() {return speciesSet;}
    
    /**
     * Specifies the fugacity fractions for the set of species that can participate in
     * an exchange move.  The given array must have the same dimension as the array of
     * species that was previously set in a call to setSpecies.  If the given set of "fractions"
     * does not sum to unity, the values will be normalized (e.g., sending the set {1.0, 1.0} 
     * leads to fugacity fractions of {0.5, 0.5}).
     */
    public void setFugacityFraction(double[] f) {
        if(f.length != nSpecies || speciesSet == null) 
            throw new IllegalArgumentException("Wrong size of fugacity-fraction array in MCMoveSemigrand");
            
        double sum = 0.0;
        for(int i=0; i<nSpecies; i++) {
            fugacityFraction[i] = f[i]; 
            sum += f[i];
            if(f[i] < 0.0) throw new IllegalArgumentException("Negative fugacity-fraction MCMoveSemigrand");
        }
        for(int i=0; i<nSpecies; i++) {fugacityFraction[i] /= sum;}//normalize to unity
    }
    /**
     * Sets fugacity fraction of the species corresponding to the given index.  Scales other
     * species fugacity fractions to normalize sum to unity.  If all other values were previously
     * zero (given species value was unity), they are all set to a uniform value that normalizes
     * the given new value for the species.
     */
    public void setFugacityFraction(int i, double f) {
        if(i < 0 || i >= nSpecies) 
            throw new IllegalArgumentException("Illegal fugacity-fraction index in MCMoveSemigrand");
            
        if(f > 1.0) f = 1.0;  //interpret any value greater than 1.0 as setting f[i] = 1.0
        else if(f < 0.0) f = 0.0; //interpret any value less than 0.0 as setting f[i] = 0.0
        
        if(fugacityFraction[i] == 1.0) { //old value is 1; set others uniformly
            double fNew = (1.0-f)/(nSpecies-1);
            for(int k=0; k<nSpecies; k++) fugacityFraction[k] = fNew;
        }
        else {
            double mult = (1.0 - f)/(1.0 - fugacityFraction[i]);
            for(int k=0; k<nSpecies; k++) fugacityFraction[k] *= mult;
        }
        fugacityFraction[i] = f;
    }
    public double getFugacityFraction(int i) {
        if(i < 0 || i >= nSpecies) 
            throw new IllegalArgumentException("Illegal fugacity-fraction index in MCMoveSemigrand");
        return fugacityFraction[i];
    }

    /**
     * Accessor method for the set of fugacity fractions.
     */
    public double[] getFugacityFraction() {return fugacityFraction;}
    
    public boolean doTrial() {
        //select species for deletion
        iDelete = Simulation.random.nextInt(nSpecies);//System.out.println("Random no. :"+randomNo);
        deleteAgent = agentSet[iDelete];
        if(deleteAgent.moleculeCount() == 0) {
            uNew = uOld = 0.0;
            return false;
        }

        //select species for insertion
        iInsert = iDelete;
        if(nSpecies == 2) iInsert = 1 - iDelete;
        else while(iInsert == iDelete) {iInsert = Simulation.random.nextInt(nSpecies);}
        insertAgent = agentSet[iInsert];
  
        AtomArrayList moleculeList = ((AtomTreeNodeGroup)deleteAgent.node).childList;
        deleteMolecule = moleculeList.get(Simulation.random.nextInt(moleculeList.size()));
        energyMeter.setTarget(deleteMolecule);
        uOld = energyMeter.getDataAsScalar();
        phase.removeMolecule(deleteMolecule);
        
        if(!reservoirs[iInsert].isEmpty()) insertMolecule = reservoirs[iInsert].remove(reservoir.size()-1);
        else insertMolecule = insertAgent.moleculeFactory().makeAtom();
        phase.addMolecule(insertMolecule, insertAgent);
        moleculeTranslator.setDestination(atomPositionDefinition.position(deleteMolecule));
        moleculeTranslator.actionPerformed(insertMolecule);
        //in general, should also randomize orintation and internal coordinates
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double getA() {
        return (double)(deleteAgent.moleculeCount()+1)/(double)insertAgent.moleculeCount()
                *(fugacityFraction[iInsert]/fugacityFraction[iDelete]);
    }
    
    public double getB() {
        energyMeter.setTarget(insertMolecule);
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld);
    }
    
    public void acceptNotify() {
        //put deleted molecule in reservoir
        reservoirs[iDelete].add(deleteMolecule);
    }

    public void rejectNotify() {
        //put deleted molecule back into phase
        phase.addMolecule(deleteMolecule, deleteAgent);
        //remove inserted molecule and put in reservoir
        phase.removeMolecule(insertMolecule);
        reservoirs[iInsert].add(insertMolecule);
    }
    
    

    public double energyChange(Phase p) {return (p == phase) ? uNew - uOld : 0.0;}
    
    public final AtomIterator affectedAtoms(Phase p) {
        if(p != phase) return AtomIterator.NULL;
        insertAtomIterator.setAtom(insertMolecule);
        deleteAtomIterator.setAtom(deleteMolecule);
        affectedAtomIterator.reset();
        return affectedAtomIterator;
    }

    /**
     * @return Returns the positionDefinition.
     */
    public AtomPositionDefinition geAtomPositionDefinition() {
        return atomPositionDefinition;
    }
    /**
     * @param positionDefinition The positionDefinition to set.
     */
    public void setAtomPositionDefinition(AtomPositionDefinition positionDefinition) {
        this.atomPositionDefinition = positionDefinition;
        moleculeTranslator.setAtomPositionDefinition(positionDefinition);
    }

}