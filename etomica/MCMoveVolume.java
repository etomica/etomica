package etomica;

import etomica.units.Dimension;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *
 * @author David Kofke
 */

public class MCMoveVolume extends MCMove {
    
    protected double pressure;
    protected final PhaseAction.Inflate inflate = new PhaseAction.Inflate();
    private final IteratorDirective iteratorDirective = new IteratorDirective();
    private AtomIterator affectedAtomIterator;

    public MCMoveVolume(IntegratorMC parentIntegrator) {
        super(parentIntegrator);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(Default.PRESSURE);
    }
    
    public void setPhase(Phase p) {
        if(p == null) return;
        super.setPhase(p);
        inflate.setPhase(phase);
        affectedAtomIterator = phase.makeAtomIterator();
    }
    
    public boolean thisTrial() {
    //debug    System.out.print("volume");
        double hOld, hNew, vOld, vNew, uOld, uNew;
        vOld = phase.volume();
        uOld = potential.set(phase).calculate(iteratorDirective, energy.reset()).sum();
        hOld = uOld + pressure*vOld;
        double vScale = (2.*Simulation.random.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/(double)phase.parentSimulation().space().D());
        inflate.setScale(rScale);
        inflate.attempt();
        uNew = potential.set(phase).calculate(iteratorDirective, energy.reset()).sum();
        hNew = uNew + pressure*vNew;
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/parentIntegrator.temperature+(phase.moleculeCount()+1)*vScale)
                < Simulation.random.nextDouble()) 
            {  //reject
              inflate.undo();
     //debug         System.out.println("reject");
              return false;
            }
     //debug   System.out.println("accept");
        return true;   //accept
    }
    
    public AtomIterator affectedAtoms() {return affectedAtomIterator;}

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public final double pressure() {return pressure;}
    public Dimension getPressureDimension() {return Dimension.PRESSURE;}
    public final void setLogPressure(int lp) {setPressure(Math.pow(10.,(double)lp));}
    
    /**
     * main method to test and demonstrate this class
     */
/*    public static void main(String args[]) {
                    
        Simulation sim = new Simulation(new Space2D());
        Simulation.instance = sim;
        Species species = new SpeciesSpheresMono(sim);
//        species.setNMolecules(2);
        P2HardSphere potential = new P2HardSphere();
//        Potential potential = new P2LennardJones();
        IntegratorMC integrator = new IntegratorMC(sim);
        MCMove mcMoveAtom = new MCMoveAtom(integrator);
        MCMove mcMoveVolume = new MCMoveVolume(integrator);
        Controller controller = new Controller(sim);
        Phase phase = new Phase(sim);
        Display displayPhase = new DisplayPhase(sim);
        DeviceSlider slider = new DeviceSlider(mcMoveVolume, "pressure");
        slider.setMinimum(0);
        slider.setMaximum(100);

		Meter energy = new MeterPotentialEnergy();
		energy.setPhase(phase);
		energy.setHistorying(true);
		energy.setActive(true);		
		energy.getHistory().setNValues(500);		
		DisplayPlot plot = new DisplayPlot();
		plot.setLabel("Energy");
		plot.setDataSources(energy.getHistory());

	    Simulation.instance.elementCoordinator.go();
        
  //      phase.setDensity(0.1);
    		                                    
        Simulation.makeAndDisplayFrame(sim);
    }//end of main  
    */
}