package etomica;

import java.util.Random;
import etomica.units.Dimension;

public class MCMoveVolume extends MCMove {
    
    private final Random rand = new Random();
    protected double pressure;
    protected PhaseAction.Inflate inflate;

    public MCMoveVolume() {
        super();
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(Default.PRESSURE);
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        inflate = new PhaseAction.Inflate(phase);
    }
    
    public void thisTrial() {
        double hOld, hNew, vOld, vNew, uOld, uNew;
        vOld = phase.volume();
        uOld = phase.energy.potential();
        hOld = uOld + pressure*vOld;
        double vScale = (2.*Math.random()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/(double)phase.parentSimulation().space().D());
        inflate.actionPerformed(phase,rScale);
        uNew = phase.energy.potential();
        hNew = uNew + pressure*vNew;
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/parentIntegrator.temperature+(phase.moleculeCount+1)*vScale)
                < Math.random()) 
            {  //reject
              inflate.retractAction();
            }
        nAccept++;   //accept
    }
    
    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public final double pressure() {return pressure;}
    public Dimension getPressureDimension() {return Dimension.PRESSURE;}
    public final void setLogPressure(int lp) {setPressure(Math.pow(10.,(double)lp));}
    
    /**
     * main method to test and demonstrate this class
     */
    public static void main(String args[]) {
        
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
            
        Simulation sim = new Simulation(new Space2D());
        Simulation.instance = sim;
        Species species = new SpeciesDisks(sim);
        Potential2 potential = new P2SimpleWrapper(sim,new PotentialHardDisk(sim));
        IntegratorMC integrator = new IntegratorMC(sim);
        MCMove mcMoveAtom = new MCMoveAtom();
        MCMove mcMoveVolume = new MCMoveVolume();
        Controller controller = new Controller(sim);
        Phase phase = new Phase(sim);
        Display displayPhase = new DisplayPhase(sim);
        DeviceSlider slider = new DeviceSlider(mcMoveVolume, "pressure");
        slider.setMinimum(0);
        slider.setMaximum(100);
        
        mcMoveVolume.setPhase(phase);
        integrator.add(mcMoveAtom);
        integrator.add(mcMoveVolume);

	    Simulation.instance.elementCoordinator.go();
    		                                    
        f.add(Simulation.instance.panel()); //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }
     
}