package etomica.integrator.mcmove;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesAgent;
import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomList;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.iterator.AtomIteratorCompound;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.DataSourceCOM;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.MCMove;

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
    private AtomList[] reservoirs;
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
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        if(speciesSet != null) {
            for(int i=0; i<nSpecies; i++) {
                agentSet[i] = speciesSet[i].getAgent(phases[0]);
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
        reservoirs = new AtomList[nSpecies];
        for(int i=0; i<nSpecies; i++) {
            speciesSet[i] = species[i];
            if(phases[0] != null) agentSet[i] = species[i].getAgent(phases[0]);
            fugacityFraction[i] = 1.0/nSpecies;
            reservoirs[i] = new AtomList();
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
        iDelete = (int)(Simulation.random.nextDouble()*nSpecies);//System.out.println("Random no. :"+randomNo);
        deleteAgent = agentSet[iDelete];
        if(deleteAgent.moleculeCount() == 0) {
            uNew = uOld = 0.0;
            return false;
        }

        //select species for insertion
        iInsert = iDelete;
        if(nSpecies == 2) iInsert = 1 - iDelete;
        else while(iInsert == iDelete) {iInsert = (int)(Simulation.random.nextDouble()*nSpecies);}
        insertAgent = agentSet[iInsert];
  
        deleteMolecule = deleteAgent.randomMolecule();
        energyMeter.setTarget(deleteMolecule);
        uOld = energyMeter.getDataAsScalar(phases[0]);
        phases[0].removeMolecule(deleteMolecule);
        
        if(!reservoirs[iInsert].isEmpty()) insertMolecule = reservoirs[iInsert].removeFirst();
        else insertMolecule = insertAgent.moleculeFactory().makeAtom();
        phases[0].addMolecule(insertMolecule, insertAgent);
        moleculeTranslator.setDestination(atomPositionDefinition.position(deleteMolecule));
        moleculeTranslator.actionPerformed(insertMolecule);
        //in general, should also randomize orintation and internal coordinates
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double lnTrialRatio() {
        return Math.log((double)(deleteAgent.moleculeCount()+1)
                        /(double)insertAgent.moleculeCount());
    }
    
    public double lnProbabilityRatio() {
        energyMeter.setTarget(insertMolecule);
        uNew = energyMeter.getDataAsScalar(phases[0]);
        return -(uNew - uOld)/temperature +
                Math.log(fugacityFraction[iInsert]/fugacityFraction[iDelete]);
    }
    
    public void acceptNotify() {
        //put deleted molecule in reservoir
        reservoirs[iDelete].add(deleteMolecule.seq);
    }

    public void rejectNotify() {
        //put deleted molecule back into phase
        phases[0].addMolecule(deleteMolecule, deleteAgent);
        //remove inserted molecule and put in reservoir
        phases[0].removeMolecule(insertMolecule);
        reservoirs[iInsert].add(insertMolecule.seq);
    }

    public double energyChange(Phase phase) {return (this.phases[0] == phase) ? uNew - uOld : 0.0;}
    
    public final AtomIterator affectedAtoms(Phase phase) {
        if(this.phases[0] != phase) return AtomIterator.NULL;
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

//    public static void main(String[] args) {
//        
//        Simulation.instance = new etomica.graphics.SimulationGraphic(new Space2D());
//        Default.TRUNCATE_POTENTIALS = false;
//	    IntegratorMC integrator = new IntegratorMC();
//	    integrator.setDoSleep(true);
//		integrator.setSleepPeriod(2);
//	    MCMoveAtom mcMove = new MCMoveAtom(integrator);
//	    MCMoveMolecule mcMoveMolecule = new MCMoveMolecule(integrator);
//	    MCMoveRotateMolecule mcMoveRotate = new MCMoveRotateMolecule(integrator);
//	    final MCMoveSemigrand mcMoveSemi = new MCMoveSemigrand(integrator);
//	    //one species with 3 atoms per molecule
//	    SpeciesSpheres species0 = new SpeciesSpheres(1,3);
//	    //two species with 1 atom per molecule
//	    SpeciesSpheresMono species1 = new SpeciesSpheresMono(5);
//	    SpeciesSpheresMono species2 = new SpeciesSpheresMono(5);
//	    
//	    mcMoveSemi.setSpecies(new Species[] {species0, species1, species2});
//        etomica.graphics.ColorSchemeByType.setColor(species1, java.awt.Color.red);
//        etomica.graphics.ColorSchemeByType.setColor(species2, java.awt.Color.green);
//	    species0.setDiameter(3.0);
//	    species1.setDiameter(5.0);
//	    species2.setDiameter(1.0);
//	    Phase phase = new Phase();
//	    
//	    //intramolecular potential for 3-atom molecules
//	    P1TetheredHardSpheres potential0 = new P1TetheredHardSpheres();
//	    //hard-sphere intermolecular potential for 3-atom molecules
//	    PotentialGroup potential00 = new PotentialGroup(2);
//	        Potential2 potential00a = new P2HardSphere(potential00, 3.0);
//	    //0-1 species potential
//	    PotentialGroup potential01 = new PotentialGroup(2);
//	        Potential2 potential01a = new P2HardSphere(potential01, 4.0);
//	    //0-2 species potential
//	    PotentialGroup potential02 = new PotentialGroup(2);
//	        Potential2 potential02a =new P2HardSphere(potential02, 2.0);
//	    //1-1 species potential
//	    P2HardSphere potential11 = new P2HardSphere(5.0);
//	    //1-2 species potential
//	    P2HardSphere potential12 = new P2HardSphere(3.0);
//	    //2-2 species potential
//	    P2HardSphere potential22 = new P2HardSphere(1.0);
//	    potential0.setSpecies(new Species[] {species0});
//	    potential00.setSpecies(new Species[] {species0, species0});
//	    potential01.setSpecies(new Species[] {species0, species1});
//	    potential02.setSpecies(new Species[] {species0, species2});
//	    potential11.setSpecies(species1, species1);
//	    potential12.setSpecies(species1, species2);
//	    potential22.setSpecies(species2, species2);
//	    Controller controller = new Controller();
//	    etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase();
//
//		MeterMoleFraction density0 = new MeterMoleFraction();
//		density0.setPhase(new Phase[] {phase});
//		AccumulatorHistory density0Hist = new AccumulatorHistory(HistoryScrolling.FACTORY, 500);
//		AccumulatorManager density0HistManager = new AccumulatorManager(density0, new Accumulator[] {density0Hist});
//		etomica.graphics.DisplayPlot plot = new etomica.graphics.DisplayPlot();
//		plot.setLabel("Species0 density");
//		plot.setDataSources(density0Hist);
//		
//		MeterPotentialEnergy energyMeter = new MeterPotentialEnergy();
//		AccumulatorHistory energyHist = new AccumulatorHistory();
//		AccumulatorManager energyHistManager = new AccumulatorManager(energyMeter, new Accumulator[] {energyHist});
//		energyHistManager.setUpdateInterval(1);
//		etomica.graphics.DisplayPlot energyPlot = new etomica.graphics.DisplayPlot();
//		energyPlot.setDataSources(energyHist);
//		
//		etomica.graphics.DeviceTernarySelector selector = 
//		    new etomica.graphics.DeviceTernarySelector(Simulation.instance, 
//		                            new String[] {"Black", "Red", "Green"});
//		selector.addListener(new etomica.graphics.DeviceTernarySelector.Listener() {
//		    public void ternaryAction(double x0, double x1, double x2) {
//		        mcMoveSemi.setFugacityFraction(new double[] {x0, x1, x2});
//		    }
//		});
//		
//		Simulation.instance.elementCoordinator.go();
//		density0.setSpecies(species0);
//	    
//        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(Simulation.instance);
//        
//    }//end of main
// */

}//end of MCMoveSemigrand