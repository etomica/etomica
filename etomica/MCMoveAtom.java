package etomica;

import etomica.units.Dimension;

/**
 * Standard Monte Carlo atom-displacement trial move.
 *
 * @author David Kofke
 */
public class MCMoveAtom extends MCMove {
    
    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
    private final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private Atom atom;

    public MCMoveAtom(IntegratorMC parentIntegrator) {
        super(parentIntegrator);
        setStepSizeMax(Default.BOX_SIZE);
        setStepSizeMin(0.0);
        setStepSize(Default.ATOM_SIZE);
        setPerParticleFrequency(true);
    }
    
    public final Dimension getStepSizeDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMaxDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMinDimension() {return Dimension.LENGTH;}
    

    public boolean thisTrial() {
        double uOld, uNew;
        if(phase.atomCount()==0) {return false;}
        int i = (int)(Simulation.random.nextDouble()*phase.atomCount());
        atom = phase.firstAtom();
        // maybe try while(i-- >= 0) {}
        for(int j=i; --j>=0; ) {atom = atom.nextAtom();}  //get ith atom in list

        uOld = potential.set(phase).calculate(iteratorDirective.set(atom), energy.reset()).sum();
        atom.coord.displaceWithin(stepSize);
        uNew = potential.calculate(iteratorDirective.set(atom), energy.reset()).sum();//not thread safe for multiphase systems
        if(uNew < uOld) {   //accept
            return true;
        }
        if(uNew >= Double.MAX_VALUE ||  //reject
           Math.exp(-(uNew-uOld)/parentIntegrator.temperature) < Simulation.random.nextDouble()) {
             atom.coord.replace();
             return false;
        }
        //accept
        return true;
    }//end thisTrial
    
    public final AtomIterator affectedAtoms() {
        affectedAtomIterator.setBasis(atom);
        return affectedAtomIterator;
    }

/*    public static void main(String[] args) {
        
	    IntegratorMC integrator = new IntegratorMC();
	    MCMoveAtom mcMove = new MCMoveAtom(integrator);//comment this line to examine LRC by itself
	    SpeciesSpheres species = new SpeciesSpheres();
	    species.setNMolecules(72);
	    final Phase phase = new Phase();
	    P2LennardJones potential = new P2LennardJones();
	    Controller controller = new Controller();
	    DisplayPhase display = new DisplayPhase();

		MeterEnergy energy = new MeterEnergy();
		energy.setPhase(phase);
		energy.setHistorying(true);
		energy.setActive(true);		
		energy.getHistory().setNValues(500);		
		DisplayPlot plot = new DisplayPlot();
		plot.setLabel("Energy");
		plot.setDataSources(energy.getHistory());
		
		integrator.setSleepPeriod(2);
		
		DeviceToggleButton lrcToggle = new DeviceToggleButton(Simulation.instance,
		    new ModulatorBoolean() {
		        public void setBoolean(boolean b) {phase.setLrcEnabled(b);}
		        public boolean getBoolean() {return phase.isLrcEnabled();}
		    },"LRC enabled", "LRC disabled" );
		
		Simulation.instance.elementCoordinator.go();
	    
        potential.setIterator(new AtomPairIterator(phase));
        potential.set(species.getAgent(phase));
        
        Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main
   */ 
}