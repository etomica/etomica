package etomica.potential;

import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeLeaf;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateAngularKinetic;
import etomica.space.ICoordinateKinetic;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Basic hard-(rod/disk/sphere) potential, with surface roughness to couple rotation and translational motions.
 * Suitable for use in space of any dimension.
 *
 * @author David Kofke
 */
public class P2RoughSphere extends P2HardSphere {

    private final etomica.space3d.Vector3D omegaSum = new etomica.space3d.Vector3D();
    private final Vector v12Surface;
    private final Vector v12Par;
    private final Vector v12Perp;
    private final Vector impulse;
    
    public P2RoughSphere(Simulation sim) {
        this(sim.space, sim.getDefaults().atomSize, sim.getDefaults().ignoreOverlap);
    }
    public P2RoughSphere(Space space, double d, boolean ignoreOverlap) {
        super(space,d,ignoreOverlap);
        v12Surface = space.makeVector();
        v12Par = space.makeVector();
        v12Perp = space.makeVector();
        impulse = space.makeVector();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard spheres with a roughness that allows collisions to transfer momentum to/from rotational motion");
        return info;
    }

    /**
     * Implements collision dynamics and updates lastCollisionVirial
     * Assumes atoms have same size and mass
     */
    public void bump(AtomSet pair, double falseTime) {
        Atom a0 = ((AtomPair)pair).atom0;
        Atom a1 = ((AtomPair)pair).atom1;
		cPair.reset((AtomPair)pair);
        cPair.resetV();
        dr.E(cPair.dr());
        Vector dv = cPair.dv();
        dr.PEa1Tv1(falseTime,dv);
        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double rm0 = ((AtomTypeLeaf)a0.type).rm();
        double rm1 = ((AtomTypeLeaf)a1.type).rm();
        double kappa = 4*((AtomType.Rotator)a0.type).momentOfInertia()[0]*rm0/(collisionDiameter*collisionDiameter);
        Vector v1 = ((ICoordinateKinetic)a0.coord).velocity();
        Vector v2 = ((ICoordinateKinetic)a1.coord).velocity();
        dr.E(cPair.dr());
        omegaSum.E(((ICoordinateAngularKinetic)a0.coord).angularVelocity());
        omegaSum.PE(((ICoordinateAngularKinetic)a1.coord).angularVelocity());
        // v12Surface should come to equal v2 - v1 - 1/2*(omega2+omega1) X (r2-r1)
        v12Surface.E(dr); // (r2 - r1)
        v12Surface.XE(omegaSum); //(r2-r1) X (omega2+omega1)
        v12Surface.TE(0.5); // +1/2 (r2-r1) X (omega2+omega1) [which equals -1/2*(omega2+omega1) X (r2-r1)]
        v12Surface.PE(v2);// p2/m2 +1/2 (r2-r1) X (omega2+omega1)
        v12Surface.ME(v1);// p2/m2 - p1/m1 +1/2 (r2-r1) X (omega2+omega1)
        //component of v12Surface parallel to r2-r1: v12Par = (v12Surface . dr) dr / |dr|^2
        v12Par.E(dr);
        v12Par.TE(v12Surface.dot(dr)/r2);
        //component of v12Surface perpendicular to r2-r1:  v12Perp = v12Surface - v12Par
        v12Perp.E(v12Surface);
        v12Perp.ME(v12Par);
        
        impulse.E(v12Par);
        impulse.PEa1Tv1(kappa/(1+kappa),v12Perp);
        impulse.TE(((AtomTypeLeaf)a0.type).getMass());
        
        ((ICoordinateKinetic)a0.coord).velocity().PEa1Tv1( rm0,impulse);
        ((ICoordinateKinetic)a1.coord).velocity().PEa1Tv1(-rm1,impulse);
        a0.coord.position().PEa1Tv1(-falseTime*rm0,impulse);
        a1.coord.position().PEa1Tv1( falseTime*rm1,impulse);
        
        //here omegaSum is used to hold the angular impulse
        omegaSum.E(dr.cross(impulse));
        omegaSum.TE(-0.5);
        ((ICoordinateAngularKinetic)a0.coord).angularVelocity().PE(omegaSum);
        ((ICoordinateAngularKinetic)a1.coord).angularVelocity().PE(omegaSum);
        
        lastCollisionVirial = 2.0/(rm0 + rm1)*bij;
        lastCollisionVirialr2 = lastCollisionVirial/r2;
    }
    
    //need to consider if hard-sphere virial is same as rough sphere virial
    public final double lastCollisionVirial() {
        return Double.NaN;
      //  return lastCollisionVirial;
    }
    
    //especially need to consider more carefully this method
    public final Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.Ev1v2(dr, dr);
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
