package etomica;
import etomica.units.*;

/**
 * Meter that takes a (scalar) Meter and records its property as a 1-dimensional function of position in the simulation volume.
 * The measured property must be a quantity that can be associated with a single atom.
 * The position coordinate lies along an arbitrary direction vector.  The profile abscissa is a ratio of the position relative
 * to its maximum value along the chosen direction, and thus lies between zero and unity.
 * 
 * @author Rob Riggleman
 */
public class MeterProfile extends MeterFunction implements EtomicaElement {
    
    /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    final Space.Vector profileVector;
    /**
     * Meter that defines the property being profiled.
     */
    MeterScalar.Atomic meter;
    
    private double profileNorm = 1.0;
    private final AtomIteratorList ai1 = new AtomIteratorList();
    
    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterProfile() {
        this(Simulation.instance);
    }
    public MeterProfile(Simulation sim) {
        super(sim);
        setX(0, 1, 100);
        profileVector = sim.space().makeVector();
        profileVector.setX(0, 1.0);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Breaks a meter's measurements into a profile taken along some direction in phase");
        return info;
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        ai1.setBasis(p.speciesMaster.atomList);
    }

    /**
     * Declares that this meter uses the boundary of the phase, as it sizes the profile length 
     * according to the dimensions of the simulation cell.
     */
    public boolean usesPhaseBoundary() {return true;}
    
    /**
     * Updates the profile vector when informed that the phase boundary has changed.
     */
    public void setPhaseBoundary(Space.Boundary b) {
        setProfileVector(profileVector);
    }
    
    /**
     * Returns the ordinate label for the profile, obtained from the associated meter.
     */
    public String getLabel() {return ((MeterScalar)meter).getLabel();}
    
    /**
     * Indicates that the abscissa coordinate is dimensionless.
     * Abscissa is formed as the ratio of the profile position relative to its maximum value.
     */
    public Dimension getXDimension() {return Dimension.NULL;}
    
    /**
     * Returns the dimensions of the ordinate, obtained from the associated meter.
     */
    public Dimension getDimension() {return (meter==null) ? null : ((MeterScalar)meter).getDimension();}
        
    /**
     * The meter that defines the profiled quantity
     */
    public MeterAbstract getMeter() {return (MeterAbstract)meter;}
    
    /**
     * Accessor method for the meter that defines the profiled quantity.
     */
    public void setMeter(MeterScalar.Atomic m) {
        if (m instanceof MeterScalar.Atomic) meter = (MeterScalar.Atomic)m;
        else throw new IllegalArgumentException("Error in meter type in MeterProfile.  Must be Meter.Atomic");
    }
    
    
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     */
    public Space.Vector getProfileVector() {return profileVector;}
    
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     * The given vector is converted to a unit vector, if not already.
     */
    public void setProfileVector(Space.Vector v) {
        profileVector.E(v);
        profileVector.normalize();
        profileNorm = 1.0/phase.boundary().dimensions().dot(profileVector);
    }
    
    /**
     * Returns the profile for the current configuration.
     */
    public double[] currentValue() {
        for (int i = 0; i <nPoints; i++) {
            y[i] = 0.0;
        }
        ai1.reset();
        while(ai1.hasNext()) {
            Atom a = ai1.next();
            double value = meter.currentValue(a);
            int i = getXIndex(a.coord.position().dot(profileVector)*profileNorm);
            y[i] += value;
        }
        double norm = 1/((double)phase.atomCount()*deltaX);
        for (int i =0; i < nPoints; i++) {
            y[i] *= norm;
        }
        return y;
    }
}//end of MeterProfile