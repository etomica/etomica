package etomica;

/**
 * Basic hard-(rod/disk/sphere) potential, with surface roughness to couple rotation and translational motions.
 * Suitable for use in space of any dimension.
 */
public class PotentialRoughSphere extends PotentialHardDisk implements EtomicaElement
{

    private final Space3D.Vector omegaSum = new Space3D.Vector();
    private final Space.Vector v12Surface;
    private final Space.Vector v12Par;
    private final Space.Vector v12Perp;
    private final Space.Vector impulse;
    
    public PotentialRoughSphere() {
        this(Simulation.instance, Default.ATOM_SIZE);
    }
    public PotentialRoughSphere(double d) {
        this(Simulation.instance, d);
    }
    public PotentialRoughSphere(Simulation sim, double d) {
        super(sim,d);
        v12Surface = sim.space.makeVector();
        v12Par = sim.space.makeVector();
        v12Perp = sim.space.makeVector();
        impulse = sim.space.makeVector();
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
        double kappa = 4*((AtomType.Rotator)a1.type).momentOfInertia()[0]*a1.rm()/(collisionDiameter*collisionDiameter);
        double r2 = pair.r2();
        dr.E(pair.dr());
        omegaSum.E(((Space.Coordinate.Angular)a1.coordinate()).angularVelocity());
        omegaSum.PE(((Space.Coordinate.Angular)a2.coordinate()).angularVelocity());
        // v12Surface should come to equal v2 - v1 - 1/2*(omega2+omega1) X (r2-r1)
        v12Surface.E(dr); // (r2 - r1)
        v12Surface.XE(omegaSum); //(r2-r1) X (omega2+omega1)
        v12Surface.TE(0.5); // +1/2 (r2-r1) X (omega2+omega1) [which equals -1/2*(omega2+omega1) X (r2-r1)]
        v12Surface.PEa1Tv1(+a2.rm(),a2.p);// p2/m2 +1/2 (r2-r1) X (omega2+omega1)
        v12Surface.PEa1Tv1(-a1.rm(),a1.p);// p2/m2 - p1/m1 +1/2 (r2-r1) X (omega2+omega1)
        //component of v12Surface parallel to r2-r1: v12Par = (v12Surface . dr) dr / |dr|^2
        v12Par.E(dr);
        v12Par.TE(v12Surface.dot(dr)/r2);
        //component of v12Surface perpendicular to r2-r1:  v12Perp = v12Surface - v12Par
        v12Perp.E(v12Surface);
        v12Perp.ME(v12Par);
        
        impulse.E(v12Par);
        impulse.PEa1Tv1(kappa/(1+kappa),v12Perp);
        impulse.TE(-a1.mass());
        
        a2.p.PE(impulse);
        a1.p.ME(impulse);
        
        //here omegaSum is used to hold the angular impulse
        omegaSum.E(dr.cross(impulse));
        omegaSum.TE(-0.5);
        ((Space.Coordinate.Angular)a1.coordinate()).angularAccelerateBy(omegaSum);
        ((Space.Coordinate.Angular)a2.coordinate()).angularAccelerateBy(omegaSum);
        
        lastCollisionVirial = 2.0/(pair.atom1().rm() + pair.atom2().rm())*pair.vDotr();
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
    public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
	    IntegratorHard integratorHard1 = new IntegratorHard();
	    SpeciesSpheresRotating speciesDisks1 = new SpeciesSpheresRotating(3);
	    Phase phase1 = new Phase();
	    P2RoughSphere p2RoughSphere1 = new P2RoughSphere();
	    Controller controller1 = new Controller();
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard1.new Timer(integratorHard1.chronoMeter());
	    timer.setUpdateInterval(10);
	    MeterEnergy meterEnergy = new MeterEnergy();
	    meterEnergy.setPhase(phase1);
	    DisplayBox displayEnergy = new DisplayBox();
	    displayEnergy.setMeter(meterEnergy);
		Simulation.instance.setBackground(java.awt.Color.yellow);
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }//end of main
    
}