/*
 * Created on Apr 12, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package etomica.modules.dcvgcmd;
import etomica.Atom;
import etomica.Integrator;
import etomica.Modifier;
import etomica.PotentialMaster;
import etomica.Species;
import etomica.graphics.DisplayPhaseCanvas3DOpenGL;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.units.Dimension;


/**
 * @author ecc4
 *
 * To change the template for this generated type comment go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
public class IntegratorDCVGCMD extends Integrator{
	
	private int i = 0;
	IntegratorMC integratormc;
	IntegratorMD integratormd;
	double zFraction = 0.1;
	private MyMCMove mcMove1, mcMove2, mcMove3, mcMove4;
	private boolean doMD = true;
	private Species speciesA, speciesB;
	public DisplayPhaseCanvas3DOpenGL display;
	private double elapsedTime = 0.0;
	
	public IntegratorDCVGCMD(PotentialMaster parent, Species species1, Species species2) {
		super(parent);
		this.speciesA = species1;
		this.speciesB = species2;
	}
	
	public void doStep() {
		if(!doMD){
//		if(false) {
			mcMove1.setupActiveAtoms();
			mcMove2.setupActiveAtoms();
			mcMove3.setupActiveAtoms();
			mcMove4.setupActiveAtoms();
//			display.setSuspended(true);
			for(int i=0; i<50; i++) integratormc.doStep();
			integratormd.reset();
//			display.setSuspended(false);
	 	} else {	
	 		integratormd.doStep();
	 		elapsedTime += integratormd.getTimeStep();	
		} 
		doMD = !doMD;
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
		//two more mcmoves here
	}
	
	public void setMu(double mu1, double mu2) {
		mcMove1.setMu(mu1);
		mcMove2.setMu(mu2);
		mcMove3.setMu(mu2);
		mcMove4.setMu(mu1);
	}
		
	/**
	 * @see etomica.Integrator#doReset()
	 */
	public void reset() {
		integratormc.reset();
		integratormd.reset();
		elapsedTime = 0.0;
	}

	/**
	 * @see etomica.Integrator#makeAgent(etomica.Atom)
	 */
	public Object makeAgent(Atom a) {
		return integratormd.makeAgent(a);
	}
	
	public MyMCMove[] mcMoves() {
		return new MyMCMove[] {mcMove1, mcMove2, mcMove3, mcMove4};
	}

//	/**
//	 * @see etomica.Integrator#addPhase(etomica.Phase)
//	 */
//	public boolean addPhase(Phase p) {
//		return super.addPhase(p);
//		return true;
//		
//	}

}
