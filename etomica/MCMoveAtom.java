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
    private PotentialCalculationEnergySumNearestPair energyDebug = 
        new PotentialCalculationEnergySumNearestPair();

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
    
    int k = 0;
    public boolean thisTrial() {
        double uOld, uNew;
        if(phase.atomCount()==0) {return false;}
 //debug       k++;
        int i = (int)(Simulation.random.nextDouble()*phase.atomCount());
        atom = phase.speciesMaster.atomList.getRandom();
//debug        if(k>11000 && (atom.node.index() == 228 || atom.node.index() == 210)) System.out.println(k + " " + atom.node.index()+atom.coord.position().toString());
/*debug        if(k == 12597) {
            System.out.println("hi!");
        }*/
        uOld = potential.set(phase).calculate(iteratorDirective.set(atom), energy.reset()).sum();
        if(uOld > 1e10) {
            System.out.println("Uold: "+uOld);
            uOld = potential.calculate(iteratorDirective.set(atom), energyDebug.reset()).sum();
            System.out.println(k + "rMin = "+Math.sqrt(energyDebug.r2Min));
            System.out.println(energyDebug.atom1.node.index()+energyDebug.atom1.coord.position().toString());
            System.out.println(energyDebug.atom2.node.index()+energyDebug.atom2.coord.position().toString());
            System.out.println("Boundary: "+phase.boundary().dimensions().toString());
            //uOld = potential.set(phase).calculate(iteratorDirective, energy.reset()).sum();
            //System.out.println("Uold: "+uOld);
            //uOld = potential.calculate(iteratorDirective.set(), energyDebug.reset()).sum();
            //System.out.println("rMin = "+Math.sqrt(energyDebug.r2Min));
            //System.out.println(energyDebug.atom1.coord.position().toString());
            //System.out.println(energyDebug.atom2.coord.position().toString());
            AtomIteratorListSimple iterator = new AtomIteratorListSimple(phase.speciesMaster.atomList);
      //      while(iterator.hasNext()) {
      //          System.out.println(iterator.next().coord.position().toString());
      //      }
            System.out.println();
            throw new RuntimeException();
        }    
        atom.coord.displaceWithin(stepSize);
        phase.boundary().centralImage(atom.coord.position());
        uNew = potential.calculate(iteratorDirective.set(atom), energy.reset()).sum();//not thread safe for multiphase systems
        if(uNew >= Double.MAX_VALUE) {
  //debug          if(k>10000 && (atom.node.index() == 228 || atom.node.index() == 210)) System.out.println("  rej1 " + atom.node.index()+atom.coord.position().toString());
            atom.coord.replace();
  //debug          if(k>10000 && (atom.node.index() == 228 || atom.node.index() == 210)) System.out.println("  rej2 " + atom.node.index()+atom.coord.position().toString());
            return false;
        } else if(uNew <= uOld) {   //accept
  //debug          if(k>10000 && (atom.node.index() == 228 || atom.node.index() == 210)) System.out.println("  acc1 " + atom.node.index()+atom.coord.position().toString());
            return true;
        } else if(  //Metropolis test, reject
            Math.exp(-(uNew-uOld)/parentIntegrator.temperature) < Simulation.random.nextDouble()) {
  //debug          if(k>10000 && (atom.node.index() == 228 || atom.node.index() == 210)) System.out.println("  rej3 " + atom.node.index()+atom.coord.position().toString());
            atom.coord.replace();
   //debug         if(k>10000 && (atom.node.index() == 228 || atom.node.index() == 210)) System.out.println("  rej4 " + atom.node.index()+atom.coord.position().toString());
             return false;
        }
        //accept
  //debug      if(k>10000 && (atom.node.index() == 228 || atom.node.index() == 210)) System.out.println("  acc2 " + atom.node.index()+atom.coord.position().toString());
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
    
    /**
     * Class used to examine configuration in cases when a large energy is found
     * before the trial is made.  A debugging tool.
     */
    private class PotentialCalculationEnergySumNearestPair extends PotentialCalculationEnergySum {
        public double r2Min;
        Atom atom1, atom2;

        public void calculate(AtomPairIterator iterator, Potential2 potential) {
            r2Min = Double.MAX_VALUE;
            atom1 = atom2 = null;
            while(iterator.hasNext()) {
                AtomPair pair = iterator.next();
                sum += potential.energy(pair);
                if(pair.r2() < r2Min) {
                    r2Min = pair.r2();
                    atom1 = pair.atom1();
                    atom2 = pair.atom2();
                }
                if(sum >= Double.MAX_VALUE) return;
            }//end while
        }//end of calculate
    }
    
}