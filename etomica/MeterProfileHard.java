package etomica;

/**
 * Profile meter for a collision-based property.
 * Tabulates values based on location of collisions.  May also include an atom-based contribution.
 *
 * @author Rob Riggleman
 */
public class MeterProfileHard extends MeterProfile implements IntegratorHardAbstract.CollisionListener, EtomicaElement {
    
    private MeterCollisional cMeter;
    private int nCollisions;
    /**
     * Atom-based contribution to profile property.  Computed using the superclass meter (if it is set to something).
     */
    private double[] z;
    /**
     * Collision-based contribution to profile property.  Accumulated in collisionAction method and averaged in
     * by currentValue method.
     */
    private double[] w;
    /**
     * Integrator time of last zeroing of collision-property accumulator (w).
     */
    private double t0;
    /** 
     * Convenience handle to integrator, to avoid repeated casting to IntegratorHard
     */
    private IntegratorHard integratorHard;
    
    
    public MeterProfileHard() {
        this(Simulation.instance);
    }
    public MeterProfileHard(Simulation sim) {
        super(sim);
        setX(0, 1, 100);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Breaks a hard meter's measurements into a profile along some direction in the phase");
        return info;
    }

    /**
     * Sets the meter used to tabulate the property.  
     * May be called with either a Meter.MeterCollisional or a MeterAbstract.  If both types of meter
     * are desired to contribute to average, should be called once with each.
     */
    public void setMeter(MeterAbstract m) {
        if (m instanceof MeterCollisional) {cMeter = (MeterCollisional)m;}
        if (m instanceof Atomic) {meter = (Atomic)m;}
    }
    
    public double[] getData() {
   //    double earlySum=0., middleSum=0., lateSum=0.;
        double t = integratorHard.elapsedTime();
        if (t > t0) {
            if (meter != null) {z = super.getData();}  //atomic contribution
            else {for (int i=0; i<nPoints; i++) { z[i] = 0.;}}
            double norm = 1.0/(deltaX*(t-t0));
            for (int i = 0; i < nPoints; i++) {
                w[i] *= norm;
                z[i] += w[i];
                w[i] = 0.0;
            }
            t0 = t;
        }
   //     for (int i = 5; i <= 15; i++) {earlySum += z[i];}
   //     for (int i = 16; i <= 35; i++) {middleSum += z[i];}
   //     for (int i = 36; i<= 45; i++) {lateSum += z[i];}
        return z;
    }
    
    /**
     * Tabulates the contribution of the last collision involving the given atom pair/potential to the profile.
     * Location of the collision is defined by the average of the atoms' coordinates.
     */
    public void collisionAction(IntegratorHardAbstract.Agent agent) {
        if(agent == null || agent.atom() == null || agent.collisionPartner() == null) return;
    //this way (adding once at midpoint of atoms) doesn't handle PBC well; 
    //causes peak in middle due to collisions at boundaries...
     //   double dot = 0.5*(agent.atom().coord.position().dot(profileVector) 
     //                       + agent.collisionPartner().coord.position().dot(profileVector));
     //...instead add twice, once at position of each atom
        double value = 0.5*cMeter.collisionValue(agent); //0.5 is to compensate for double addition
        
        double dot = agent.atom().coord.position().dot(profileVector);
        int j = getXIndex(dot);
        w[j] += value;
        dot = agent.collisionPartner().coord.position().dot(profileVector);
        j = getXIndex(dot);
        w[j] += value;
    }
    
    protected void resizeArrays() {
        super.resizeArrays();
        z = new double[nPoints];
        w = new double[nPoints];
    }
    
    /**
     * Invokes superclass method and registers meter as a collisionListener to integrator.
     * Performs only superclass method if integrator is not an instance of IntegratorHard.
     */
	protected void setPhaseIntegrator(Integrator newIntegrator) {
	    super.setPhaseIntegrator(newIntegrator);
	    if(newIntegrator instanceof IntegratorHard) {
	        integratorHard = (IntegratorHard)newIntegrator;
	        integratorHard.addCollisionListener(this);
    	    t0 = integratorHard.elapsedTime();
	    }
	    else throw new IllegalArgumentException("Error in integrator type in MeterPressureHard");
	}

}