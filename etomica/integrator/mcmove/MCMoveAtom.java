package etomica.integrator.mcmove;

import etomica.Default;
import etomica.atom.Atom;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorSinglet;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.MCMove;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
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
    protected final MeterPotentialEnergy energyMeter;
    protected final Vector translationVector;
    protected Atom atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;

    public MCMoveAtom(PotentialMaster potentialMaster) {
        super(potentialMaster, 1);
        atomSource = new AtomSourceRandomLeaf();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        translationVector = potentialMaster.getSpace().makeVector();
        setStepSizeMax(Default.BOX_SIZE/2);
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
        atom = atomSource.getAtom();
        if (atom == null) return false;
        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        if(uOld > 1e10 && !Default.FIX_OVERLAP) {
            System.out.println("Uold: "+uOld);
            throw new RuntimeException("Overlap found in configuration");
        }
        translationVector.setRandomCube();
        translationVector.TE(stepSize);
        atom.coord.position().PE(translationVector);
        uNew = Double.NaN;
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
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld)/temperature;
    }
    
    public double energyChange(Phase phase) {return (this.phases[0] == phase) ? uNew - uOld : 0.0;}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
    }
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.  This method should cause the system to be restored to the condition
     * before the most recent call to doTrial.
     */
    public void rejectNotify() {
        translationVector.TE(-1);
        atom.coord.position().PE(translationVector);
    }
        
    
    public AtomIterator affectedAtoms(Phase phase) {
        if(this.phases[0] != phase) return AtomIterator.NULL;
        affectedAtomIterator.setAtom(atom);
        return affectedAtomIterator;
    }
    
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        energyMeter.setPhase(p[0]);
        atomSource.setPhase(p[0]);
    }
    
    /**
     * @return Returns the atomSource.
     */
    public AtomSource getAtomSource() {
        return atomSource;
    }
    /**
     * @param atomSource The atomSource to set.
     */
    public void setAtomSource(AtomSource source) {
        atomSource = source;
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
		    new ModifierBoolean() {
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