/*
 * Created on Apr 12, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package etomica.modules.dcvgcmd;
import etomica.atom.Atom;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.DisplayPhaseCanvas3DOpenGL;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorIntervalEvent;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.modifier.Modifier;
import etomica.nbr.PotentialMasterHybrid;
import etomica.potential.PotentialMaster;
import etomica.species.Species;
import etomica.units.Dimension;


/**
 * @author ecc4
 *
 * To change the template for this generated type comment go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
public class IntegratorDCVGCMD extends Integrator {
	
	IntegratorMC integratormc;
	IntegratorMD integratormd;
	double zFraction = 0.1;
	private MyMCMove mcMove1, mcMove2, mcMove3, mcMove4;
	private Species speciesA, speciesB;
	public DisplayPhaseCanvas3DOpenGL display;
	private double elapsedTime = 0.0;
    private final PotentialMasterHybrid potentialMasterHybrid;
	private int MDStepCount, MDStepRepetitions;
    
	public IntegratorDCVGCMD(PotentialMaster parent, double temperature, Species species1, Species species2) {
		super(parent, temperature);
		this.speciesA = species1;
		this.speciesB = species2;
		potentialMasterHybrid = (parent instanceof PotentialMasterHybrid)
                        ? (PotentialMasterHybrid)parent : null;
        setMDStepRepetitions(50);
    }
    
    public void setMDStepRepetitions(int interval) {
        MDStepRepetitions = interval;
        if (MDStepCount > interval || MDStepCount == 0) MDStepCount = interval;
    }
    
    protected void setup() throws ConfigurationOverlapException {
        integratormc.initialize();
        integratormd.initialize();
        super.setup();
    }
    
    public void setTemperature(double t) {
        super.setTemperature(t);
        if (integratormc != null) {
            integratormc.setTemperature(t);
        }
        if (integratormd != null) {
            integratormd.setTemperature(t);
        }
    }
	
	public void doStep() {
        if (potentialMasterHybrid != null) {
            potentialMasterHybrid.setUseNbrLists(MDStepCount > 0);
        }
		if(MDStepCount == 0){
		    MDStepCount = MDStepRepetitions;
			mcMove1.setupActiveAtoms();
			mcMove2.setupActiveAtoms();
			mcMove3.setupActiveAtoms();
			mcMove4.setupActiveAtoms();
			for(int i=0; i<50; i++) {
                integratormc.doStep();
                integratormc.fireIntervalEvent(intervalEventMC);
            }
            potentialMasterHybrid.setUseNbrLists(true);
            potentialMasterHybrid.getNeighborManager().setQuiet(true);
            potentialMasterHybrid.getNeighborManager().updateNbrsIfNeeded(integratormd);
            potentialMasterHybrid.getNeighborManager().setQuiet(false);
            try {
                integratormd.reset();
            } catch(ConfigurationOverlapException e) {}
	 	} else {
            MDStepCount--;
	 		integratormd.doStep();
            integratormd.fireIntervalEvent(intervalEventMD);
	 		elapsedTime += integratormd.getTimeStep();	
		} 
	}

//modulator for Mu's	
	public class Mu1Modulator implements Modifier {  
	   public void setValue(double x) {
		setMu(x, mcMove2.getMu());
	   }
	   public String getLabel() {return "mu1";}
	   public double getValue() {return mcMove1.getMu();}
	   public Dimension getDimension() {return Dimension.NULL;} } 
	
	public class Mu2Modulator implements Modifier {
	   public void setValue(double x) {
		setMu(mcMove1.getMu(), x);
	   }
	   public double getValue() {return mcMove2.getMu();}
	   public Dimension getDimension() {return Dimension.NULL;} 
	   public String getLabel() {return "mu2";}
	}
	
	
	public double elapsedTime() {
		return elapsedTime;
	}
	
	public void setIntegrators(IntegratorMC intmc, IntegratorMD intmd) {
		integratormc = intmc;
		integratormd = intmd;
        integratormc.setTemperature(temperature);
        integratormd.setTemperature(temperature);
		integratormd.addPhase(phase[0]);
		integratormc.addPhase(phase[0]);
		mcMove1 = new MyMCMove(potential, -zFraction);
		mcMove2 = new MyMCMove(potential, +zFraction);
		integratormc.addMCMove (mcMove1);
		integratormc.addMCMove (mcMove2);
		mcMove1.setSpecies(speciesA);
		mcMove2.setSpecies(speciesA);
		mcMove1.integrator = this;
		mcMove2.integrator = this;
		mcMove3 = new MyMCMove(potential, -zFraction);
		mcMove4 = new MyMCMove(potential, +zFraction);
		integratormc.addMCMove (mcMove3);
		integratormc.addMCMove (mcMove4);
		mcMove3.setSpecies(speciesB);
		mcMove4.setSpecies(speciesB);
		mcMove3.integrator = this;
		mcMove4.integrator = this;

        intervalEventMC = new IntegratorIntervalEvent(integratormc, 1);
        intervalEventMD = new IntegratorIntervalEvent(integratormd, 1);
	}
	
	public void setMu(double mu1, double mu2) {
		mcMove1.setMu(mu1);
		mcMove2.setMu(mu2);
//        mcMove3.setMu(mu2);
//        mcMove4.setMu(mu1);
        mcMove3.setMu(Double.NEGATIVE_INFINITY);
        mcMove4.setMu(Double.NEGATIVE_INFINITY);
    }
		
	/**
	 * @see etomica.integrator.Integrator#doReset()
	 */
	public void reset() throws ConfigurationOverlapException {
        if(!initialized) return;
        potentialMasterHybrid.setUseNbrLists(false);
		integratormc.reset();
        potentialMasterHybrid.setUseNbrLists(true);
		integratormd.reset();
		elapsedTime = 0.0;
	}

	/**
	 * @see etomica.integrator.Integrator#makeAgent(etomica.Atom)
	 */
	public Object makeAgent(Atom a) {
		return integratormd.makeAgent(a);
	}
	
	public MyMCMove[] mcMoves() {
		return new MyMCMove[] {mcMove1, mcMove2, mcMove3, mcMove4};
	}

    IntegratorIntervalEvent intervalEventMD, intervalEventMC;

}
