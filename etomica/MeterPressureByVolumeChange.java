package etomica;
import etomica.units.*;

public class MeterPressureByVolumeChange extends MeterFunction implements EtomicaElement {
    
    PhaseAction.Inflate inflater;
    double[] scale;
    boolean isotropic = true;
    int inflateDimension = 0; //keys direction for inflation if not isotropic
    
    public MeterPressureByVolumeChange() {
        this(0);
    }
    
    public MeterPressureByVolumeChange(int i) {
        this(Simulation.instance, i);
    }
    
    public MeterPressureByVolumeChange(Simulation sim, int i) {
        super(sim);
        setX(-0.001, 0.001, 10);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Pressure measured by trial volume perturbations");
        return info;
    }

    public void setIsotropic(boolean b) {
        isotropic = b;
        setX(xMin, xMax, nPoints);
        resetInflater();
    }
    public boolean getIsotropic() {return isotropic;}
    public void setInflateDimension(int i) {
        inflateDimension = i;
    }
    public int getInflateDimension() {return inflateDimension;}
    
    private void resetInflater() {
        //inflater = isotropic ? new PhaseAction.Inflate(phase) : new PhaseAction.InflateXY(phase,inflateDimension);
        //inflater = new PhaseAction.Inflate(phase);
    }    
    
    public boolean usesPhaseBoundary() {return false;}
    public boolean usesPhaseIteratorFactory() {return false;}
    
    public void setPhase(Phase p) {
        super.setPhase(p);
   
        inflater = new PhaseAction.Inflate(p);
        resetInflater();
    }
    
    public void setX(double min, double max, int n) {
        super.setX(min, max, n);
        //x is scaling in volume if isotropic, but is linear scaling if not isotropic
        for(int i=0; i<nPoints; i++) { //disallow x = 0
            if(x[i] == 0.0) x[i] = 0.1*deltaX;
        }
        scale = new double[nPoints];
        
        double mult = isotropic ? 1./(double)Simulation.instance.space.D() : 1.0;
        for(int i=0; i<nPoints; i++) {
            scale[i] = Math.exp(mult*x[i]);
        }
    }
    
    public double[] currentValue() {
           
        for(int i=0; i<nPoints; i++) {
            double uOld = phase.energy.potential();
            if(isotropic)
            inflater.actionPerformed(scale[i]);
            else inflater.actionPerformed(phase,scale[i],inflateDimension);
            double uNew = phase.energy.potential();
            y[i] = Math.exp(-(uNew-uOld)/phase.integrator().temperature()
                              + phase.moleculeCount()*x[i]);
            
            if ( isotropic) inflater.retractAction();
            else inflater.retractAction(inflateDimension);
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

      public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        Default.ATOM_SIZE =1.2;
        Phase phase1 = new Phase();
        MCMoveAtom mcmove= new MCMoveAtom();
        
        DisplayPlot plot1 = new DisplayPlot();

        MeterPressureByVolumeChange meterp = new MeterPressureByVolumeChange();
        IntegratorMC integratorMC1 = new IntegratorMC();
	    integratorMC1.add(mcmove);

	        meterp.setPhase(phase1);
	       meterp.setActive(true);
//	       plot1.setMeter(meterp);
           plot1.setMeterFunction(meterp);
	       plot1.setUseCurrentValue(false);
	     SpeciesDisks    speciesDisk1 = new SpeciesDisks();
	     speciesDisk1.setNMolecules(200);
	    
	    

        P2IdealGas p2IdealGas = new P2IdealGas();
	   
	    
	    Controller controller1 = new Controller();
	    
	  
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    MeterEnergy meterEnergy1 = new MeterEnergy();
	    
	    DisplayBox box1 = new DisplayBox();
	    
	    box1.setMeter(meterEnergy1);
	    box1.setUseCurrentValue(false);
		displayPhase1.setPhase(phase1);
		DeviceSlider temperatureSlider = new DeviceSlider(integratorMC1, "temperature");
	    temperatureSlider.setUnit(new Unit(Kelvin.UNIT));
	    temperatureSlider.setMinimum(50);
	    temperatureSlider.setMaximum(500);
	    
    	//Simulation.instance.add(displayPhase1);
    	//phase1.boundary().dimensions().setComponent(0,15);
       //Simulation.elementCoordinator = new Simulation.ElementCoordinator.Basic();
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
        //phase1.integrator().setTemperature();
        meterEnergy1.setPhase(phase1);
        //meterpxx.setPhase(phase1);
        phase1.setIntegrator(integratorMC1);
        integratorMC1.addIntervalListener(box1);
       integratorMC1.addIntervalListener(meterp);
       integratorMC1.setTemperature(Kelvin.UNIT.toSim(10));
        //integratorMC1.addIntervalListener(box2);
        integratorMC1.addIntervalListener(displayPhase1);
        controller1.add(integratorMC1);
        
		Simulation.instance.setBackground(java.awt.Color.blue);		                                    
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        //phase1.boundary().dimensions().setComponent(0,30);
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main

}