package etomica;

import etomica.units.Dimension;

public class MCMoveAtom extends MCMove {
    
    private PotentialMaster.Agent phasePotential;
    private final IteratorDirective iteratorDirective = new IteratorDirective(IteratorDirective.BOTH);
    private final PotentialCalculation.EnergySum energy = new PotentialCalculation.EnergySum();

    public MCMoveAtom() {
        super();
        setStepSizeMax(Default.BOX_SIZE);
        setStepSizeMin(0.0);
        setStepSize(Default.ATOM_SIZE);
        setPerParticleFrequency(true);
    }
    
    public final Dimension getStepSizeDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMaxDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMinDimension() {return Dimension.LENGTH;}
    

    public void thisTrial() {
        double uOld, uNew;
        if(phase.atomCount()==0) {return;}
        int i = (int)(Simulation.random.nextDouble()*phase.atomCount());
        Atom a = phase.firstAtom();
        // maybe try while(i-- >= 0) {}
        for(int j=i; --j>=0; ) {a = a.nextAtom();}  //get ith atom in list
//        phasePotential.calculate(iteratorDirective.set(a), energy.reset());
//        uOld = energy.sum();
        uOld = phase.potential.calculate(iteratorDirective.set(a), energy.reset()).sum();
        a.coord.displaceWithin(stepSize);
//        phase.boundary().centralImage(a.coord.position());  //maybe a better way than this
        phase.iteratorFactory().moveNotify(a);
        uNew = phase.potential().calculate(iteratorDirective.set(a), energy.reset()).sum();
        if(uNew < uOld) {   //accept
            nAccept++;
            return;
        }
        if(uNew >= Double.MAX_VALUE ||  //reject
           Math.exp(-(uNew-uOld)/parentIntegrator.temperature) < Simulation.random.nextDouble()) {
             a.coord.replace();
             phase.iteratorFactory().moveNotify(a);
             return;
        }
        nAccept++;   //accept
    }//end thisTrial
    
    public static void main(String[] args) {
        
	    IntegratorMC integrator = new IntegratorMC();
	    MCMoveAtom mcMove = new MCMoveAtom();
	    SpeciesDisks species = new SpeciesDisks();
	    Phase phase = new Phase();
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
		plot.setDataSource(energy.getHistory());
		
		integrator.setSleepPeriod(2);
		
		Simulation.instance.elementCoordinator.go();
	    integrator.add(mcMove);
        Potential2.Agent potentialAgent = (Potential2.Agent)potential.getAgent(phase);
        potentialAgent.setIterator(new AtomPairIterator(phase));
        
        Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main
    
}