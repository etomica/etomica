package etomica;

import etomica.units.Dimension;

/**
 * Standard Monte Carlo atom-displacement trial move.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 07/09/02 Added energyChange() method 
  * 07/09/03 Changed fields from private to protected
  */
  
public class MCMoveAtom extends MCMove {
    
    protected final AtomIteratorSinglet affectedAtomIterator = new AtomIteratorSinglet();
    private final MeterPotentialEnergy energyMeter;
    protected Atom atom;
    protected double uOld;
    protected double uNew = Double.NaN;
/*debug* /    private PotentialCalculationEnergySumNearestPair energyDebug = 
             new PotentialCalculationEnergySumNearestPair();
    private int idx1 = 8;
    private int idx2 = 9;
    private int kmax = 3240;
    private int k = 0;
    private Atom atomIdx1, atomIdx2;
/* */
    public MCMoveAtom(PotentialMaster potentialMaster) {
        super(potentialMaster);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(Default.BOX_SIZE);
        setStepSizeMin(0.0);
        setStepSize(Default.ATOM_SIZE);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        setName("MCMoveAtom");
    }
    
    public final Dimension getStepSizeDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMaxDimension() {return Dimension.LENGTH;}
    public final Dimension getStepSizeMinDimension() {return Dimension.LENGTH;}
    
    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
      //  Atom atom0 = phase.speciesMaster.atomList.getFirst();
  //      System.out.println(((IteratorFactoryCell.NeighborSequencer)atom0.seq).nbrLink.previous.toString()+((IteratorFactoryCell.NeighborSequencer)atom0.seq).nbrLink.next.toString()+((AtomLinker.Tab[])((IteratorFactoryCell.NeighborSequencer)atom0.seq).cell.agents[0])[0].toString());
 /*debug* / atomIdx1 = phase.speciesMaster.atomList.get(idx1);
           atomIdx2 = phase.speciesMaster.atomList.get(idx2);  // */
        if(phase.atomCount()==0) return false;
 /*debug* /       k++;  // */
        atom = phase.speciesMaster.atomList.getRandom();
/*debug* /        if(k>3240 && (atom.node.index() == idx1 || atom.node.index() == idx2)) System.out.println(k + " " + atom.node.index()+atom.coord.position().toString());
/*debug* /        if(k == 1683) {
                     if(k>kmax && (atom.node.index() == 0)) System.out.println(((IteratorFactoryCell.NeighborSequencer)atom.seq).nbrLink.previous.toString()+((IteratorFactoryCell.NeighborSequencer)atom.seq).nbrLink.next.toString()+((AtomLinker.Tab[])((IteratorFactoryCell.NeighborSequencer)atom.seq).cell.agents[0])[0].toString());
    
            System.out.println("hi!");
        }// */
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar(phase);
        if(uOld > 1e10) {
            System.out.println("Uold: "+uOld);
   /*debug* /         uOld = potential.calculate(iteratorDirective.set(atom), energyDebug.reset()).sum();
            System.out.println(k + "  rMin = "+Math.sqrt(energyDebug.r2Min));
            System.out.println(energyDebug.atom1.node.index()+energyDebug.atom1.coord.position().toString());
            System.out.println(energyDebug.atom2.node.index()+energyDebug.atom2.coord.position().toString());
            System.out.println("Boundary: "+phase.boundary().dimensions().toString());
            //uOld = potential.set(phase).calculate(iteratorDirective, energy.reset()).sum();
            //System.out.println("Uold: "+uOld);
            //uOld = potential.calculate(iteratorDirective.set(), energyDebug.reset()).sum();
      //      AtomIteratorListSimple iterator = new AtomIteratorListSimple(phase.speciesMaster.atomList);
      //      while(iterator.hasNext()) {
      //          System.out.println(iterator.next().coord.position().toString());
      //      }
      //      System.out.println();
      //      */
            throw new RuntimeException("Overlap found in configuration");
        }    
        atom.coord.displaceWithin(stepSize);
        uNew = Double.NaN;
 //       phase.boundary().centralImage(atom.coord.position());
        return true;
    }//end of doTrial
    
    
    /**
     * Returns log of the ratio of the trial probabilities, ln(Tij/Tji) for the
     * states encountered before (i) and after (j) the most recent call to doTrial(). 
     * Tij is the probability that this move would generate state j from state i, and
     * Tji is the probability that a subsequent call to doTrial would return to state i
     * from state j.
     */
    public double lnTrialRatio() {return 0.0;}
    
    /**
     * Returns the log of the limiting-distribution probabilities of states, ln(Pj/Pi), 
     * for the states encountered before (i) and after (j) the most recent call to 
     * doTrial.
     */
    public double lnProbabilityRatio() {
//        energyMeter.setTarget(atom);
        uNew = energyMeter.getDataAsScalar(phase);
        return -(uNew - uOld)/temperature;
    }
    
    public double energyChange(Phase phase) {return (this.phase == phase) ? uNew - uOld : 0.0;}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
  /*debug* /      if(k>kmax && (atom.node.index() == idx1 || atom.node.index() == idx2)) System.out.println(k+"  acc2 " + atomIdx1.node.index()+atomIdx1.coord.position().toString() + atomIdx2.node.index() + atomIdx2.coord.position().toString() + Math.sqrt(parentIntegrator().parentSimulation().space.r2(atomIdx1.coord.position(),atomIdx2.coord.position(),phase.boundary())));
            //     if(k>kmax && (atom.node.index() == 0)) System.out.println(((IteratorFactoryCell.NeighborSequencer)atom.seq).nbrLink.previous.toString()+((IteratorFactoryCell.NeighborSequencer)atom.seq).nbrLink.next.toString()+((AtomLinker.Tab[])((IteratorFactoryCell.NeighborSequencer)atom.seq).cell.agents[0])[0].toString());
    // */  
    }
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.  This method should cause the system to be restored to the condition
     * before the most recent call to doTrial.
     */
    public void rejectNotify() {
  /*debug* /          if(k>kmax && (atom.node.index() == idx1 || atom.node.index() == idx2)) System.out.println(k+"  rej1 " + atomIdx1.node.index()+atomIdx1.coord.position().toString() + atomIdx2.node.index() + atomIdx2.coord.position().toString() + Math.sqrt(parentIntegrator().parentSimulation().space.r2(atomIdx1.coord.position(),atomIdx2.coord.position(),phase.boundary())));
            //         if(k>kmax && (atom.node.index() == 0)) System.out.println(((IteratorFactoryCell.NeighborSequencer)atom.seq).nbrLink.previous.toString()+((IteratorFactoryCell.NeighborSequencer)atom.seq).nbrLink.next.toString()+((AtomLinker.Tab[])((IteratorFactoryCell.NeighborSequencer)atom.seq).cell.agents[0])[0].toString());
  /* //  */          atom.coord.replace();
  /*debug* /          if(k>kmax && (atom.node.index() == idx1 || atom.node.index() == idx2)) System.out.println(k+"  rej2 " + atomIdx1.node.index()+atomIdx1.coord.position().toString() + atomIdx2.node.index() + atomIdx2.coord.position().toString() + Math.sqrt(parentIntegrator().parentSimulation().space.r2(atomIdx1.coord.position(),atomIdx2.coord.position(),phase.boundary())));
            //         if(k>kmax && (atom.node.index() == 0)) System.out.println(((IteratorFactoryCell.NeighborSequencer)atom.seq).nbrLink.previous.toString()+((IteratorFactoryCell.NeighborSequencer)atom.seq).nbrLink.next.toString()+((AtomLinker.Tab[])((IteratorFactoryCell.NeighborSequencer)atom.seq).cell.agents[0])[0].toString());
     // } // */
    }
        
    
    public AtomIterator affectedAtoms(Phase phase) {
        if(this.phase != phase) return AtomIterator.NULL;
        affectedAtomIterator.setAtom(atom);
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
    
    /**
     * Class used to examine configuration in cases when a large energy is found
     * before the trial is made.  A debugging tool.
     */
//      private class PotentialCalculationEnergySumNearestPair extends PotentialCalculationEnergySum {
//        public double r2Min;
//        Atom atom1, atom2;
//
//        public void calculate(AtomPairIterator iterator, Potential2 potential) {
//            r2Min = Double.MAX_VALUE;
//            atom1 = atom2 = null;
//            while(iterator.hasNext()) {
//                AtomPair pair = iterator.next();
//                sum += potential.energy(pair);
//                if(pair.r2() < r2Min) {
//                    r2Min = pair.r2();
//                    atom1 = pair.atom1();
//                    atom2 = pair.atom2();
//                }
//                if(sum >= Double.MAX_VALUE) return;
//            }//end while
//        }//end of calculate
//    }
 // */   
}