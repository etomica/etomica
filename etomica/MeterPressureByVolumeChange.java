package etomica;
import etomica.units.*;

/**
 * Evaluates the pressure by examining the change in energy accompanying
 * small changes in volume.
 */
public class MeterPressureByVolumeChange extends MeterFunction implements EtomicaElement {
    
    PhaseAction.InflateAnisotropic inflater;
    Space.Vector[] scale;
    Space space;
    boolean[] inflateDimensions;
    private IteratorDirective iteratorDirective;
    private final PotentialCalculationEnergySum energy;
    private final PotentialMaster potential;
    private int nDimension;
    
    public MeterPressureByVolumeChange() {
        this(null);
    }
    
    /**
     * Constructor that indicates volume change should be performed anisotropically.
     * @param dimensions array indicating which dimensions should be inflated
     */
    public MeterPressureByVolumeChange(boolean[] dimensions) {
        this(Simulation.instance, dimensions);
    }
    
    public MeterPressureByVolumeChange(Simulation sim, boolean[] dimensions) {
        super(sim);
        space = sim.space;
        potential = sim.hamiltonian.potential;
        energy = sim.energySum(this);
        setX(-0.001, 0.001, 10);
        inflateDimensions = new boolean[space.D()];
        if(dimensions == null) {
            dimensions = new boolean[space.D()];
            for(int i=0; i<dimensions.length; i++) dimensions[i] = true;
        }
        setInflateDimensions(dimensions);
        iteratorDirective = new IteratorDirective();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pressure measured by trial volume perturbations");
        return info;
    }

    /**
     * For anisotropic volume change, indicates dimension in which volume is perturbed.
     */
    public final void setInflateDimensions(boolean[] directions) {
        if(directions.length != inflateDimensions.length) {
            throw new IllegalArgumentException();
        }
        nDimension = 0;
        for(int i=0; i<directions.length; i++) {
            inflateDimensions[i] = directions[i];
            if(inflateDimensions[i]) nDimension++;
            
        }
    }
    /**
     * Accessor method for setInflateDimension.
     */
    public boolean[] getInflateDimensions() {return inflateDimensions;}
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        inflater = new PhaseAction.InflateAnisotropic(p);
    }
    
    public void setX(double min, double max, int n) {
        super.setX(min, max, n);
        //x is scaling in volume if isotropic, but is linear scaling if not isotropic
        for(int i=0; i<nPoints; i++) { //disallow x = 0
            if(x[i] == 0.0) x[i] = 0.1*deltaX;
        }
        scale = new Space.Vector[nPoints];
        
        double mult = 1./(double)nDimension;
        for(int i=0; i<nPoints; i++) {
            scale[i] = space.makeVector();
            scale[i].E(Math.exp(mult*x[i]));
            for(int j=0; j<space.D(); j++) {
                if(!inflateDimensions[j]) scale[i].setComponent(j,1.0);
            }
        }
    }
    
    public double[] currentValue() {
        for(int i=0; i<nPoints; i++) {
            double uOld = potential.set(phase).calculate(iteratorDirective, energy.reset()).sum();
            inflater.setScale(scale[i]);
            inflater.attempt();
            double uNew = potential.calculate(iteratorDirective, energy.reset()).sum();
            y[i] = Math.exp(-(uNew-uOld)/phase.integrator().temperature()
                              + phase.moleculeCount()*x[i]);
            
            inflater.undo();
            //System.out.println( "  uNew " + uNew +" uOld " +uOld +" x " + x[i] +" scale" + scale[i]+ " y " +y[i] );
        }
        return y;
    }
    
    public double[] average() {
        average = super.average();
        for(int i=0; i<nPoints; i++) {
            average[i] = Math.log(average[i])/(x[i]*phase.moleculeCount());
        }
         //for(int i=0; i<nPoints; i++) {
           //     System.out.println("i " + i + "y "+ average[i]);
            //} 
        return average;
    }
    
    public Dimension getDimension() {return Dimension.NULL;}
    public Dimension getXDimension() {return Dimension.NULL;}

 /*     public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        
        Default.ATOM_SIZE =1.2;
        
        Simulation.instance = new Simulation(new Space2DCell());
        Phase phase1 = new Phase();
        MCMoveAtom mcmove= new MCMoveAtom();
        
        DisplayPlot plot1 = new DisplayPlot();

        MeterPressureByVolumeChange meterp = new MeterPressureByVolumeChange();
        meterp.setIsotropic(true);
        
        IntegratorMC integratorMC1 = new IntegratorMC();
	    integratorMC1.add(mcmove);

	    meterp.setPhase(phase1);
	    meterp.setActive(true);
        plot1.setDataSources(meterp);
	    plot1.setWhichValue(MeterAbstract.AVERAGE);
	    SpeciesSpheres speciesSphere1 = new SpeciesSpheres();
	    speciesSphere1.setNMolecules(200);
	    PotentialLJ potentialLJ = new PotentialLJ();
	    P2SimpleWrapper potential = new P2SimpleWrapper(potentialLJ);
	    
	    Controller controller1 = new Controller();
	  
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    MeterEnergy meterEnergy1 = new MeterEnergy();
	    
	    DisplayBox box1 = new DisplayBox();
	    
	    box1.setMeter(meterEnergy1);
	    box1.setWhichValue(MeterAbstract.AVERAGE);
		displayPhase1.setPhase(phase1);
		DeviceSlider temperatureSlider = new DeviceSlider(integratorMC1, "temperature");
	    temperatureSlider.setUnit(new Unit(Kelvin.UNIT));
	    temperatureSlider.setMinimum(50);
	    temperatureSlider.setMaximum(500);
	    
		Simulation.instance.elementCoordinator.go(); 
		
        integratorMC1.setTemperature(Kelvin.UNIT.toSim(10));
        
		Simulation.instance.setBackground(java.awt.Color.blue);		                                    
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
*/
}