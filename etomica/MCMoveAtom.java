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
 //       System.out.println(phase.atomCount()+ "  "+i); System.exit(0);
        // maybe try while(i-- >= 0) {}
        for(int j=i; --j>=0; ) {
 //           System.out.println(atom.signature());  //get ith atom in list
            atom = atom.nextAtom();
        }
 //       if(atom.node.parentSpeciesAgent().node.index() != 2) {
 //           System.out.println();
 //       }   
 //       System.out.println();
        uOld = potential.set(phase).calculate(iteratorDirective.set(atom), energy.reset()).sum();
        atom.coord.displaceWithin(stepSize);
        uNew = potential.calculate(iteratorDirective.set(atom), energy.reset()).sum();//not thread safe for multiphase systems
 //       if(atom.node.parentSpeciesAgent().node.index() != 2) System.out.println(uOld+"  "+uNew);
        //System.out.println(uOld+"  "+uNew);
        if(uNew >= Double.MAX_VALUE) {
            atom.coord.replace();
            return false;
        } else if(uNew <= uOld) {   //accept
            return true;
        } else if(  //Metropolis test, reject
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
/*
    public static void main(String[] args) {
        
        Simulation.instance = new etomica.graphics.SimulationGraphic();
	    IntegratorMC integrator = new IntegratorMC();
	    MCMoveAtom mcMove = new MCMoveAtom(integrator);//comment this line to examine LRC by itself
	    SpeciesSpheresMono species = new SpeciesSpheresMono();
	    species.setNMolecules(72);
	    final Phase phase = new Phase();
	    P2LennardJones potential = new P2LennardJones();
	    Controller controller = new Controller();
	    etomica.graphics.DisplayPhase display = new etomica.graphics.DisplayPhase();

		MeterEnergy energy = new MeterEnergy();
		energy.setPhase(phase);
		energy.setHistorying(true);
		energy.setActive(true);		
		energy.getHistory().setNValues(500);		
		etomica.graphics.DisplayPlot plot = new etomica.graphics.DisplayPlot();
		plot.setLabel("Energy");
		plot.setDataSources(energy.getHistory());
		
		integrator.setSleepPeriod(2);
		
		etomica.graphics.DeviceToggleButton lrcToggle = new etomica.graphics.DeviceToggleButton(Simulation.instance,
		    new ModulatorBoolean() {
		        public void setBoolean(boolean b) {phase.setLrcEnabled(b);}
		        public boolean getBoolean() {return phase.isLrcEnabled();}
		    },"LRC enabled", "LRC disabled" );
		
		Simulation.instance.elementCoordinator.go();
	    
        potential.setIterator(new AtomPairIteratorGeneral(phase));
        potential.set(species.getAgent(phase));
        
        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(Simulation.instance);
    }//end of main
    */
}