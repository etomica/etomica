package etomica;

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
 
public class MCMoveSemigrand extends MCMove {
    
    private Species[] speciesSet;
    private SpeciesAgent[] agentSet;
    private double[] fugacityFraction;
    private int nSpecies;
    private transient Atom deleteMolecule, insertMolecule;
    private final AtomIteratorListSimple deleteAtomIterator;
    private final AtomIteratorListSimple insertAtomIterator;
    private final AtomIteratorCompound affectedAtomIterator; 

    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
    
    public MCMoveSemigrand(IntegratorMC parentIntegrator) {
        super(parentIntegrator);
        deleteAtomIterator = new AtomIteratorListSimple();
        insertAtomIterator = new AtomIteratorListSimple();
        affectedAtomIterator = new AtomIteratorCompound(new AtomIterator[] {deleteAtomIterator, insertAtomIterator});
        setTunable(false);
        setPerParticleFrequency(true);
    }
    
    /**
     * Extends the superclass method to initialize the exchange-set species agents for the phase.
     */
    public void setPhase(Phase p) {
        super.setPhase(p);
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
        for(int i=0; i<nSpecies; i++) {
            speciesSet[i] = species[i];
            if(phase != null) agentSet[i] = species[i].getAgent(phase);
            fugacityFraction[i] = 1.0/(double)nSpecies;
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
            double fNew = (1.0-f)/(double)(nSpecies-1);
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
    
        
    
    /**
     * Performs the trial and decides acceptance.
     */
    public final boolean thisTrial() {
        
        double uOld = 0.0;
        double uNew = 0.0;

        //select species for deletion
        int iDelete = (int)(Simulation.random.nextDouble()*nSpecies);//System.out.println("Random no. :"+randomNo);
        SpeciesAgent deleteAgent = agentSet[iDelete];
        if(deleteAgent.moleculeCount() == 0) return false;

        //select species for insertion
        int iInsert = iDelete;
        if(nSpecies == 2) iInsert = 1 - iDelete;
        else while(iInsert == iDelete) {iInsert = (int)(Simulation.random.nextDouble()*nSpecies);}
        SpeciesAgent insertAgent = agentSet[iInsert];
  
        deleteMolecule = deleteAgent.randomMolecule();
        uOld = potential.set(phase).calculate(iteratorDirective.set(deleteMolecule), energy.reset()).sum();
        deleteMolecule.sendToReservoir();
        
        insertMolecule = insertAgent.addNewAtom();
        insertMolecule.coord.translateTo(deleteMolecule.coord.position());
        //in general, should also randomize orintation and internal coordinates
        
        uNew = potential.calculate(iteratorDirective.set(insertMolecule), energy.reset()).sum();
//        if(iInsert==0) System.out.println(iDelete + "   "+ iInsert + "   "+uOld+"   "+uNew);

        if(uNew == Double.MAX_VALUE) {//reject
            deleteMolecule.node.setParent(deleteAgent);
            insertMolecule.sendToReservoir();
            return false;
        }
        double w = ((double)(deleteAgent.moleculeCount()+1)/(double)insertAgent.moleculeCount())
                * (fugacityFraction[iInsert]/fugacityFraction[iDelete])
                * Math.exp(-(uNew-uOld)/parentIntegrator.temperature);
        
        if(w > 1.0 || Simulation.random.nextDouble() < w) {//accept
            return true;
        }
        
        //reject
        deleteMolecule.node.setParent(deleteAgent);
        insertMolecule.sendToReservoir();
        return false;
        
    }//end of thisTrial
    
    public final AtomIterator affectedAtoms() {
        insertAtomIterator.setBasis(insertMolecule);
        deleteAtomIterator.setBasis(deleteMolecule);
        return affectedAtomIterator;
    }
/*
    public static void main(String[] args) {
        
        Simulation.instance = new etomica.graphics.SimulationGraphic(new Space2D());
        Default.TRUNCATE_POTENTIALS = false;
	    IntegratorMC integrator = new IntegratorMC();
	    MCMoveAtom mcMove = new MCMoveAtom(integrator);
	    MCMoveMolecule mcMoveMolecule = new MCMoveMolecule(integrator);
	    MCMoveRotateMolecule mcMoveRotate = new MCMoveRotateMolecule(integrator);
	    final MCMoveSemigrand mcMoveSemi = new MCMoveSemigrand(integrator);
	    //one species with 3 atoms per molecule
	    SpeciesSpheres species0 = new SpeciesSpheres(10,3);
	    //two species with 1 atom per molecule
	    SpeciesSpheresMono species1 = new SpeciesSpheresMono(2);
	    SpeciesSpheresMono species2 = new SpeciesSpheresMono(2);
	    
	    mcMoveSemi.setSpecies(new Species[] {species0, species1, species2});
        etomica.graphics.ColorSchemeByType.setColor(species1, java.awt.Color.red);
        etomica.graphics.ColorSchemeByType.setColor(species2, java.awt.Color.green);
	    species0.setDiameter(3.0);
	    species1.setDiameter(5.0);
	    species2.setDiameter(1.0);
	    Phase phase = new Phase();
	    
	    //intramolecular potential for 3-atom molecules
	    P1TetheredHardSpheres potential0 = new P1TetheredHardSpheres();
	    //hard-sphere intermolecular potential for 3-atom molecules
	    Potential2Group potential00 = new Potential2Group();
	        Potential2 potential00a = new P2HardSphere(potential00, 3.0);
	    //0-1 species potential
	    Potential2Group potential01 = new Potential2Group();
	        Potential2 potential01a = new P2HardSphere(potential01, 4.0);
	    //0-2 species potential
	    Potential2Group potential02 = new Potential2Group();
	        Potential2 potential02a =new P2HardSphere(potential02, 2.0);
	    //1-1 species potential
	    P2HardSphere potential11 = new P2HardSphere(5.0);
	    //1-2 species potential
	    P2HardSphere potential12 = new P2HardSphere(3.0);
	    //2-2 species potential
	    P2HardSphere potential22 = new P2HardSphere(1.0);
	    potential0.setSpecies(species0);
	    potential00.setSpecies(species0, species0);
	    potential01.setSpecies(species0, species1);
	    potential02.setSpecies(species0, species2);
	    potential11.setSpecies(species1, species1);
	    potential12.setSpecies(species1, species2);
	    potential22.setSpecies(species2, species2);
	    Controller controller = new Controller();
	    etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase();

		MeterMoleFraction density0 = new MeterMoleFraction();
		density0.setPhase(phase);
		density0.setHistorying(true);
		density0.setActive(true);		
		density0.getHistory().setNValues(500);		
		etomica.graphics.DisplayPlot plot = new etomica.graphics.DisplayPlot();
		plot.setLabel("Species0 density");
		plot.setDataSources(density0.getHistory());
		
		MeterPotentialEnergy energyMeter = new MeterPotentialEnergy();
		energyMeter.setUpdateInterval(1);
		energyMeter.setHistorying(true);
		etomica.graphics.DisplayPlot energyPlot = new etomica.graphics.DisplayPlot();
		energyPlot.setDataSources(energyMeter.getHistory());
		
		etomica.graphics.DeviceTernarySelector selector = 
		    new etomica.graphics.DeviceTernarySelector(Simulation.instance, 
		                            new String[] {"Black", "Red", "Green"});
		selector.addListener(new etomica.graphics.DeviceTernarySelector.Listener() {
		    public void ternaryAction(double x0, double x1, double x2) {
		        mcMoveSemi.setFugacityFraction(new double[] {x0, x1, x2});
		    }
		});
		
		integrator.setSleepPeriod(2);
		
		Simulation.instance.elementCoordinator.go();
		density0.setSpecies(species0);
	    
        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(Simulation.instance);
        
    }//end of main
    */
    
}//end of MCMoveSemigrand