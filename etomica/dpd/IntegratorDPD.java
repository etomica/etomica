/*
 * Created on Apr 22, 2003
 */
package etomica.dpd;

import etomica.Atom;
import etomica.AtomIterator;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Simulation;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterTemperature;
import etomica.integrator.IntegratorMD;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.Vector;

/**
 * Dissipative particle dynamic integrator.
 * Based on the modified velocity Verlet algorithm proposed by Groot and Warren
 * (J. Chem Phys. 107 (11), 15 September 1997, pp.4423-4435)
 * lambdaV is the stepping parameter used to compute the velocity.
 * The timestep is defaulted in the constructor to 0.04;
 * 
 * @author cribbin
 *
 */

/*
 * History
 * 11/17/03 (DAK) added moveNotify call in doStep method
 */

public class IntegratorDPD extends IntegratorMD implements EtomicaElement {
	
	public String getVersion() {return "IntegratorDPD:03.04.22/"+etomica.integrator.VERSION;}
	AtomIterator atomIterator;
	private double t2;
	public final PotentialCalculationForceSum forceSum;
	private final IteratorDirective allAtoms = new IteratorDirective();
	private final MeterTemperature meterTemperature = new MeterTemperature(this);
    
    double lambdaV;	//paramter for modified velocity Verlet integration
    
	public IntegratorDPD() {
		this(Simulation.instance, 0.5);
	}
		
	public IntegratorDPD(double l){
		this(Simulation.instance, l);
	}
		
	public IntegratorDPD(SimulationElement parent){
		this(Simulation.instance, 0.5);
	}
		
	public IntegratorDPD(SimulationElement parent, double l) {
		super(parent);
		forceSum = new PotentialCalculationForceSum(space);
		setTimeStep(0.04);  
		lambdaV = 0.65;
	}
		
		
    
	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo("Molecular dynamics using modified velocity Verlet integration algorithm");
		return info;
	}
    
	public boolean addPhase(Phase p) {
		if(!super.addPhase(p)) return false;
		atomIterator = new AtomIteratorLeafAtoms(p);
		meterTemperature.setPhase(new Phase[] {p});
		return true;
	}
 
	public void setTimeStep(double t) {
		super.setTimeStep(t);
		t2 = 0.5*timeStep*timeStep;
  	}


          
//	  --------------------------------------------------------------
//	   steps all particles across time interval tStep

	public void doStep() {

		atomIterator.reset();              //reset iterator of atoms
		while(atomIterator.hasNext()) {    //loop over all atoms
			Atom a = atomIterator.nextAtom();  //  advancing positions full step
			MyAgent agent = (MyAgent)a.ia;     //  and momenta half step
			Vector r = a.coord.position();
			Vector p = a.coord.momentum();
			r.PEa1Tv1(timeStep*a.type.rm(),p);         // r += p*dt/m
			r.PEa1Tv1(t2*a.type.rm(),agent.force);  // r += f(old)*dt^2/2m
			p.PEa1Tv1(lambdaV*timeStep, agent.force);  // p += f(old)*dt*lambdaV
			agent.fOld.E(agent.force);
			agent.force.E(0.0);				//sets all components of the force vector to 0
		}
                
		//Compute forces on each atom
		potential.calculate(firstPhase, allAtoms, forceSum);
       
		//Finish integration step
		atomIterator.reset();
		while(atomIterator.hasNext()) {     //loop over atoms again
			Atom a = atomIterator.nextAtom();   //  finishing the momentum step
			Vector p = a.coord.momentum();
			MyAgent agent = (MyAgent)a.ia;
			p.PEa1Tv1((0.5-lambdaV)*timeStep, agent.fOld);//p += f(old)*(1/2 - lambdaV)*dt
			p.PEa1Tv1(0.5*timeStep, agent.force);//p += f(new)*dt/2
//			a.coord.momentum().PEa1Tv1(0.5*timeStep,((MyAgent)a.ia).force);  //p += f(new)*dt/2  We do not need to update by lambda.
		}
	}
    
//--------------------------------------------------------------

	public void reset() {
		atomIterator.reset();
		while(atomIterator.hasNext()) {
			Atom a = atomIterator.nextAtom();
			MyAgent agent = (MyAgent)a.ia;
			agent.force.E(0.0);
		}
		potential.calculate(firstPhase, allAtoms, forceSum);//assumes only one phase
        super.reset();
	}
             
//	  --------------------------------------------------------------

	public final Object makeAgent(Atom a) {
		return new MyAgent(simulation(),a);
	}
            
	public final static class MyAgent implements Integrator.Forcible {  //need public so to use with instanceof
		public Atom atom;
		public Vector force, fOld;

		public MyAgent(Simulation sim, Atom a) {
			atom = a;
			force = sim.space().makeVector();
			fOld = sim.space().makeVector();
		}
        
		public Vector force() {return force;}
	}//end of MyAgent
 
	/**
	 * @return
	 */
	public double getLambdaV() {
		return lambdaV;
	}
	/**
	 * @param d
	 */
	public void setLambdaV(double d) {
		lambdaV = d;
	}
		
}//end class-IntegratorDPD
