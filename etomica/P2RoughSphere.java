package etomica;

/**
 * Basic hard-(rod/disk/sphere) potential, with surface roughness to couple rotation and translational motions.
 * Suitable for use in space of any dimension.
 *
 * @author David Kofke
 */
public class P2RoughSphere extends P2HardSphere implements EtomicaElement {

    private final Space3D.Vector omegaSum = new Space3D.Vector();
    private final Space.Vector v12Surface;
    private final Space.Vector v12Par;
    private final Space.Vector v12Perp;
    private final Space.Vector impulse;
    
    public P2RoughSphere() {
        this(Simulation.instance.hamiltonian.potential, Default.ATOM_SIZE);
    }
    public P2RoughSphere(double d) {
        this(Simulation.instance.hamiltonian.potential, d);
    }
    public P2RoughSphere(PotentialGroup parent, double d) {
        super(parent,d);
        v12Surface = parentSimulation().space.makeVector();
        v12Par = parentSimulation().space.makeVector();
        v12Perp = parentSimulation().space.makeVector();
        impulse = parentSimulation().space.makeVector();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard spheres with a roughness that allows collisions to transfer momentum to/from rotational motion");
        return info;
    }

    /**
     * Implements collision dynamics and updates lastCollisionVirial
     * Assumes atoms have same size and mass
     */
    public void bump(AtomPair pair) {
        Atom a1 = pair.atom1();
        Atom a2 = pair.atom2();
        double kappa = 4*((AtomType.Rotator)a1.type).momentOfInertia()[0]*a1.coord.rm()/(collisionDiameter*collisionDiameter);
        double r2 = pair.r2();
        Space.Vector p1 = a1.coord.momentum();
        Space.Vector p2 = a2.coord.momentum();
        dr.E(pair.dr());
        omegaSum.E(((Space.Coordinate.Angular)a1.coord).angularVelocity());
        omegaSum.PE(((Space.Coordinate.Angular)a2.coord).angularVelocity());
        // v12Surface should come to equal v2 - v1 - 1/2*(omega2+omega1) X (r2-r1)
        v12Surface.E(dr); // (r2 - r1)
        v12Surface.XE(omegaSum); //(r2-r1) X (omega2+omega1)
        v12Surface.TE(0.5); // +1/2 (r2-r1) X (omega2+omega1) [which equals -1/2*(omega2+omega1) X (r2-r1)]
        v12Surface.PEa1Tv1(+a2.coord.rm(),p2);// p2/m2 +1/2 (r2-r1) X (omega2+omega1)
        v12Surface.PEa1Tv1(-a1.coord.rm(),p1);// p2/m2 - p1/m1 +1/2 (r2-r1) X (omega2+omega1)
        //component of v12Surface parallel to r2-r1: v12Par = (v12Surface . dr) dr / |dr|^2
        v12Par.E(dr);
        v12Par.TE(v12Surface.dot(dr)/r2);
        //component of v12Surface perpendicular to r2-r1:  v12Perp = v12Surface - v12Par
        v12Perp.E(v12Surface);
        v12Perp.ME(v12Par);
        
        impulse.E(v12Par);
        impulse.PEa1Tv1(kappa/(1+kappa),v12Perp);
        impulse.TE(-a1.coord.mass());
        
        p2.PE(impulse);
        p1.ME(impulse);
        
        //here omegaSum is used to hold the angular impulse
        omegaSum.E(dr.cross(impulse));
        omegaSum.TE(-0.5);
        ((Space.Coordinate.Angular)a1.coord).angularAccelerateBy(omegaSum);
        ((Space.Coordinate.Angular)a2.coord).angularAccelerateBy(omegaSum);
        
        lastCollisionVirial = 2.0/(pair.atom1().coord.rm() + pair.atom2().coord.rm())*pair.vDotr();
        lastCollisionVirialr2 = lastCollisionVirial/r2;
    }
    
    //need to consider if hard-sphere virial is same as rough sphere virial
    public final double lastCollisionVirial() {
        return Double.NaN;
      //  return lastCollisionVirial;
    }
    
    //especially need to consider more carefully this method
    public final Space.Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.E(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;        
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
/*    public static void main(String[] args) {
	    IntegratorHard integratorHard1 = new IntegratorHard();
	    SpeciesSpheresRotating speciesSpheres1 = new SpeciesSpheresRotating();
	    Phase phase = new Phase();
	    P2RoughSphere potential = new P2RoughSphere();
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard1.new Timer(integratorHard1.chronoMeter());
	    timer.setUpdateInterval(10);
	    MeterEnergy meterEnergy = new MeterEnergy();
	    meterEnergy.setPhase(phase);
	    DisplayBox displayEnergy = new DisplayBox();
	    displayEnergy.setMeter(meterEnergy);
		Simulation.instance.panel().setBackground(java.awt.Color.yellow);
		Simulation.instance.elementCoordinator.go();
		                                            
        Simulation.makeAndDisplayFrame(Simulation.instance);
    }//end of main
    */
}