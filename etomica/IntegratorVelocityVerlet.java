package etomica;

/* History of changes
 * 08/29/02 (DAK) changed Andersen thermostat to velocity-scaling thermostat
 */

public final class IntegratorVelocityVerlet extends IntegratorMD implements EtomicaElement {

    public String getVersion() {return "IntegratorVelocityVerlet:01.07.05/"+IntegratorMD.VERSION;}

    AtomIterator atomIterator;
    
    public final PotentialCalculationForceSum forceSum;
    private final IteratorDirective allAtoms = new IteratorDirective();
    private final MeterTemperature meterTemperature = new MeterTemperature((Space)null);
    
    //Fields for Andersen thermostat
    double nu = 0.001;  //heat bath "collision" frequency
                
    public IntegratorVelocityVerlet() {
        this(Simulation.instance);
    }
    public IntegratorVelocityVerlet(Simulation sim) {
        super(sim);
        forceSum = new PotentialCalculationForceSum(sim.space());
        
        setTimeStep(etomica.units.LennardJones.Time.UNIT.toSim(2.0));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Molecular dynamics using velocity Verlet integration algorithm");
        return info;
    }
    
    public boolean addPhase(Phase p) {
        if(!super.addPhase(p)) return false;
        atomIterator = p.makeAtomIterator();
        meterTemperature.setPhase(p);
        return true;
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
            Atom a = atomIterator.next();  //  advancing positions full step
            Agent agent = (Agent)a.ia;     //  and momenta half step
            Space.Vector r = a.coord.position();
            Space.Vector p = a.coord.momentum();
            p.PEa1Tv1(0.5*timeStep,agent.force);  // p += f(old)*dt/2
            r.PEa1Tv1(timeStep*a.coord.rm(),p);         // r += p*dt/m
            agent.force.E(0.0);
        }
                
        //Compute forces on each atom
        potential.calculate(allAtoms, forceSum);
        
        //Finish integration step
        atomIterator.reset();
        while(atomIterator.hasNext()) {     //loop over atoms again
            Atom a = atomIterator.next();   //  finishing the momentum step
            a.coord.momentum().PEa1Tv1(0.5*timeStep,((Agent)a.ia).force);  //p += f(new)*dt/2
        }
        if(isothermal) {
            //velocity-rescaling thermostat
            double s = Math.sqrt(this.temperature/meterTemperature.currentValue(firstPhase.speciesMaster));
            atomIterator.reset();
            while(atomIterator.hasNext()) {
                Atom a = atomIterator.next();
                a.coord.momentum().TE(s); //scale momentum
            }
            
            //Andersen thermostat
            /* 
            atomIterator.reset();
            double nut = nu*timeStep;
            while(atomIterator.hasNext()) {
                Atom a = atomIterator.next();
                if(Simulation.random.nextDouble() < nut) a.coord.randomizeMomentum(temperature);  //this method in Atom needs some work
            }
            */
        }
        return;
    }
    
    public void setAndersenNu(double n) {nu = n;}
    public double getAndersenNu() {return nu;}
            
//--------------------------------------------------------------

    protected void doReset() {
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
            Agent agent = (Agent)a.ia;
            agent.force.E(0.0);
        }
        potential.set(firstPhase).calculate(allAtoms, forceSum);//assumes only one phase
    }
              
//--------------------------------------------------------------

    public final Integrator.Agent makeAgent(Atom a) {
        return new Agent(parentSimulation(),a);
    }
            
    public final static class Agent implements Integrator.Agent.Forcible {  //need public so to use with instanceof
        public Atom atom;
        public Space.Vector force;

        public Agent(Simulation sim, Atom a) {
            atom = a;
            force = sim.space().makeVector();
        }
        
        public Space.Vector force() {return force;}
    }//end of Agent
    
}//end of IntegratorVelocityVerlet

