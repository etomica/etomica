package etomica;
import etomica.units.*;

/**
 * Evaluates the pressure by examining the change in energy accompanying
 * small changes in volume.
 */
public class MeterPressureByVolumeChange extends MeterFunction implements EtomicaElement {
    
    PhaseAction.Inflate inflater;
    double[] scale;
    boolean isotropic = true;
    int inflateDimension = 0; //keys direction for inflation if not isotropic
    
    public MeterPressureByVolumeChange() {
        this(0);
    }
    
    /**
     * Constructor that indicates volume change should be performed anisotropically.
     * @param i index of the dimension in which volume should be expanded (e.g. i = 0 indicates x-dimension)
     */
    public MeterPressureByVolumeChange(int i) {
        this(Simulation.instance, i);
    }
    
    public MeterPressureByVolumeChange(Simulation sim, int i) {
        super(sim);
        setX(-0.001, 0.001, 10);
        setInflateDimension(i);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pressure measured by trial volume perturbations");
        return info;
    }

    /**
     * Method to indicate if volume change should or should not be performed isotropically.
     */
    public final void setIsotropic(boolean b) {
        isotropic = b;
        setX(xMin, xMax, nPoints);
    }
    /**
     * Accessor method accompanying setIsotropic.
     */
    public boolean getIsotropic() {return isotropic;}
    
    /**
     * For anisotropic volume change, indicates dimension in which volume is perturbed.
     */
    public final void setInflateDimension(int i) {
        inflateDimension = i;
        setIsotropic(false);
    }
    /**
     * Accessor method for setInflateDimension.
     */
    public int getInflateDimension() {return inflateDimension;}
    
    public boolean usesPhaseBoundary() {return false;}
    public boolean usesPhaseIteratorFactory() {return false;}
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        inflater = new PhaseAction.Inflate(p);
    }
    
    public void setX(double min, double max, int n) {
        super.setX(min, max, n);
        //x is scaling in volume if isotropic, but is linear scaling if not isotropic
        for(int i=0; i<nPoints; i++) { //disallow x = 0
            if(x[i] == 0.0) x[i] = 0.1*deltaX;
        }
        scale = new double[nPoints];
        
        double mult = isotropic ? 1./(double)parentSimulation().space().D() : 1.0;
        for(int i=0; i<nPoints; i++) {
            scale[i] = Math.exp(mult*x[i]);
        }
    }
    
    public double[] currentValue() {
        for(int i=0; i<nPoints; i++) {
            double uOld = phase.energy.potential();
            inflater.setScale(scale[i]);
            if(isotropic) inflater.attempt();
            else inflater.attempt(inflateDimension);
            double uNew = phase.energy.potential();
            y[i] = Math.exp(-(uNew-uOld)/phase.integrator().temperature()
                              + phase.moleculeCount()*x[i]);
            
            if (isotropic) inflater.undo();
            else inflater.undo(inflateDimension);
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
	    SpeciesDisks speciesDisk1 = new SpeciesDisks();
	    speciesDisk1.setNMolecules(200);
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