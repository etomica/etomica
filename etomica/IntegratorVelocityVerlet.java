package etomica;

/* History of changes
 * 08/29/02 (DAK) changed Andersen thermostat to velocity-scaling thermostat
 * 01/10/03 (DAK) reintroduced Andersen thermostat with flag to allow
 * selection of it or velocity-rescaling as the thermostat mechanism
 * */

public final class IntegratorVelocityVerlet extends IntegratorMD implements EtomicaElement {

    AtomIterator atomIterator;
    
    public final PotentialCalculationForceSum forceSum;
    private final Space space;
    private final IteratorDirective allAtoms = new IteratorDirective();
    private final MeterTemperature meterTemperature = new MeterTemperature();
    
    //Fields for Andersen thermostat
    double nu = 0.001;  //heat bath "collision" frequency
    private boolean andersenThermostat = false;
                
    public IntegratorVelocityVerlet(PotentialMaster potentialMaster, Space space) {
        super(potentialMaster);
        this.space = space;
        forceSum = new PotentialCalculationForceSum(space);
        
        setTimeStep(etomica.units.systems.LJ.SYSTEM.time().toSim(2.0));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using velocity Verlet integration algorithm");
        return info;
    }
    
    public boolean addPhase(Phase p) {
        if(!super.addPhase(p)) return false;
        atomIterator = new AtomIteratorLeafAtoms(p);
        meterTemperature.setPhase(phase);
        return true;
    }
    /**
     * Sets flag indicating which thermostat to use.  If true, Andersen
     * collision thermostat is used; if false (default) velocity rescaling is
     * used.  Whether to use thermostat or not is set independently; thermostat
     * is imposed only if a call to set(isothermal) is made with a "true"
     * argument.
     * @param b flag indicating type of thermostat to use.
     */
    public void setAndersenThermostat(boolean b) {
    	andersenThermostat = b;
    }
    /** 
     * @return boolean current value of andersenThermostat flag.  See mutator
     * method for more detail.
     */
    public boolean isAndersenThermostat() {
    	return andersenThermostat;
    }
     
  private double t2;
  public final void setTimeStep(double t) {
    super.setTimeStep(t);
    t2 = timeStep*timeStep;
  }
  

          
//--------------------------------------------------------------
// steps all particles across time interval tStep

    public void doStep() {

        atomIterator.reset();              //reset iterator of atoms
        while(atomIterator.hasNext()) {    //loop over all atoms
            Atom a = atomIterator.nextAtom();  //  advancing positions full step
            MyAgent agent = (MyAgent)a.ia;     //  and momenta half step
            Space.Vector r = a.coord.position();
            Space.Vector p = a.coord.momentum();
            p.PEa1Tv1(0.5*timeStep,agent.force);  // p += f(old)*dt/2
            r.PEa1Tv1(timeStep*a.coord.rm(),p);         // r += p*dt/m
            agent.force.E(0.0);
        }
                
        //Compute forces on each atom
        potential.calculate(firstPhase, allAtoms, forceSum);
        
        //Finish integration step
        atomIterator.reset();
        while(atomIterator.hasNext()) {     //loop over atoms again
            Atom a = atomIterator.nextAtom();   //  finishing the momentum step
            a.coord.momentum().PEa1Tv1(0.5*timeStep,((MyAgent)a.ia).force);  //p += f(new)*dt/2
        }
        if(isothermal) {
        	if(!andersenThermostat) {
	            //velocity-rescaling thermostat
	            double s = Math.sqrt(this.temperature/meterTemperature.getDataAsScalar(firstPhase));
	            atomIterator.reset();
	            while(atomIterator.hasNext()) {
	                Atom a = atomIterator.nextAtom();
	                a.coord.momentum().TE(s); //scale momentum
	            }
	        } else {
            //Andersen thermostat
            atomIterator.reset();
            double nut = nu*timeStep;
            while(atomIterator.hasNext()) {
                Atom a = atomIterator.nextAtom();
                if(Simulation.random.nextDouble() < nut) a.coord.randomizeMomentum(temperature);  //this method in Atom needs some work
            }
        	}
        }
        return;
    }
    
    public void setAndersenNu(double n) {nu = n;}
    public double getAndersenNu() {return nu;}
            
//--------------------------------------------------------------

    protected void doReset() {
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            MyAgent agent = (MyAgent)a.ia;
            agent.force.E(0.0);
        }
        potential.calculate(firstPhase, allAtoms, forceSum);//assumes only one phase
    }
              
//--------------------------------------------------------------

    public final Object makeAgent(Atom a) {
        return new MyAgent(space,a);
    }
            
    public final static class MyAgent implements Integrator.Forcible {  //need public so to use with instanceof
        public Atom atom;
        public Space.Vector force;

        public MyAgent(Space space, Atom a) {
            atom = a;
            force = space.makeVector();
        }
        
        public Space.Vector force() {return force;}
    }//end of MyAgent
    
}//end of IntegratorVelocityVerlet

