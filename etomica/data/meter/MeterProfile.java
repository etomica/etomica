package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.DataSourceAtomic;
import etomica.data.DataSourceUniform;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.phase.Phase;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Dimension;

/**
 * Meter that takes a (scalar) Meter and records its property as a 1-dimensional function of position in the simulation volume.
 * The measured property must be a quantity that can be associated with a single atom.
 * The position coordinate lies along an arbitrary direction vector.  The profile abscissa is a ratio of the position relative
 * to its maximum value along the chosen direction, and thus lies between zero and unity.
 * 
 * @author Rob Riggleman
 */
public class MeterProfile implements DataSource, Meter, java.io.Serializable {
    
    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterProfile(Space space) {
        xDataSource = new DataSourceUniform("x", Dimension.LENGTH);
        xData = (DataDoubleArray)xDataSource.getData();
        profileVector = space.makeVector();
        profileVector.setX(0, 1.0);
        position = space.makeVector();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Breaks a meter's measurements into a profile taken along some direction in phase");
        return info;
    }

    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }

    /**
     * The meter that defines the profiled quantity
     */
    public DataSourceAtomic getDataSource() {return meter;}
    
    /**
     * Accessor method for the meter that defines the profiled quantity.
     */
    public void setDataSource(DataSourceAtomic m) {
        if (!(m.getDataInfo().getDataClass().equals(DataDouble.class))) {
            throw new IllegalArgumentException("data source must return a DataDouble");
        }
        yData = new DataDoubleArray(m.getDataInfo().getLabel()+" Profile", m.getDataInfo().getDimension(),xDataSource.getNValues());
        data = new DataFunction(new DataDoubleArray[] {xData}, yData);
        meter = m;
    }
    
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     */
    public Vector getProfileVector() {return profileVector;}
    
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     * The given vector is converted to a unit vector, if not already.
     */
    public void setProfileVector(Vector v) {
        profileVector.E(v);
        profileVector.normalize();
    }
    
    /**
     * Returns the profile for the current configuration.
     */
    public Data getData() {
        Boundary boundary = phase.getBoundary();
        profileNorm = 1.0/boundary.getDimensions().dot(profileVector);
        yData.E(0);
        double[] y = yData.getData();
        ai1.reset();
        while(ai1.hasNext()) {
            Atom a = ai1.nextAtom();
            double value = ((DataDouble)meter.getData(a)).x;
            position.E(a.coord.position());
            position.PE(boundary.centralImage(position));
            int i = xDataSource.getIndex(position.dot(profileVector)*profileNorm);
            y[i] += value;
        }
        double dx = (xDataSource.getXMax() - xDataSource.getXMin())/y.length;
        double norm = 1.0/(phase.atomCount()*dx);
        yData.TE(norm);
        return data;
    }
    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
        ai1.setPhase(phase);
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    private Phase phase;
    private String name;
    private DataSourceUniform xDataSource;
    private DataFunction data;
    private DataDoubleArray xData;
    private DataDoubleArray yData;
    /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    final Vector profileVector;
    final Vector position;
    /**
     * Meter that defines the property being profiled.
     */
    DataSourceAtomic meter;
    
    private double profileNorm = 1.0;
    private final AtomIteratorPhaseDependent ai1 = new AtomIteratorLeafAtoms();
    
}//end of MeterProfile
