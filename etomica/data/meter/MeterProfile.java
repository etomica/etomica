package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataSourceAtomic;
import etomica.data.DataSourceIndependent;
import etomica.data.DataSourceUniform;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.phase.Phase;
import etomica.space.Boundary;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.units.Length;

/**
 * Meter that takes a (scalar) Meter and records its property as a 1-dimensional function of position in the simulation volume.
 * The measured property must be a quantity that can be associated with a single atom.
 * The position coordinate lies along an arbitrary direction vector.  The profile abscissa is a ratio of the position relative
 * to its maximum value along the chosen direction, and thus lies between zero and unity.
 * 
 * @author Rob Riggleman
 */
public class MeterProfile implements DataSource, DataSourceIndependent, java.io.Serializable {
    
    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterProfile(Space space) {
        xDataSource = new DataSourceUniform("x", Length.DIMENSION);
        profileVector = space.makeVector();
        profileVector.setX(0, 1.0);
        position = space.makeVector();
        tag = new DataTag();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Breaks a meter's measurements into a profile taken along some direction in phase");
        return info;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }

    /**
     * The meter that defines the profiled quantity
     */
    public DataSourceAtomic getDataSource() {return meter;}
    
    /**
     * Accessor method for the meter that defines the profiled quantity.
     */
    public void setDataSource(DataSourceAtomic m) {
        if (!(m.getDataInfo() instanceof DataInfoDouble)) {
            throw new IllegalArgumentException("data source must return a DataDouble");
        }
        data = new DataFunction(new int[] {xDataSource.getNValues()});
        dataInfo = new DataInfoFunction(m.getDataInfo().getLabel()+" Profile", m.getDataInfo().getDimension(), this);
        meter = m;
        dataInfo.addTag(meter.getTag());
        dataInfo.addTag(tag);
    }
    
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     */
    public IVector getProfileVector() {return profileVector;}
    
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     * The given vector is converted to a unit vector, if not already.
     */
    public void setProfileVector(IVector v) {
        profileVector.E(v);
        profileVector.normalize();
        double halfBox = 0.5*phase.getBoundary().getDimensions().dot(profileVector);
        xDataSource.setXMin(-halfBox);
        xDataSource.setXMax(halfBox);
    }
    
    /**
     * Returns the profile for the current configuration.
     */
    public Data getData() {
        Boundary boundary = phase.getBoundary();
        data.E(0);
        double[] y = data.getData();
        ai1.reset();
        while(ai1.hasNext()) {
            AtomLeaf a = (AtomLeaf)ai1.nextAtom();
            double value = ((DataDouble)meter.getData(a)).x;
            position.E(a.getPosition());
            position.PE(boundary.centralImage(position));
            int i = xDataSource.getIndex(position.dot(profileVector));
            y[i] += value;
        }
        double dx = (xDataSource.getXMax() - xDataSource.getXMin())/y.length;
        double norm = 1.0/(phase.atomCount()*dx);
        data.TE(norm);
        return data;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray)xDataSource.getData();
    }
    
    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray)xDataSource.getDataInfo();
    }
    
    public int getIndependentArrayDimension() {
        return 1;
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
        double halfBox = 0.5*phase.getBoundary().getDimensions().dot(profileVector);
        xDataSource.setXMin(-halfBox);
        xDataSource.setXMax(halfBox);
        ai1.setPhase(phase);
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    private static final long serialVersionUID = 1L;
    private Phase phase;
    private String name;
    private DataSourceUniform xDataSource;
    private DataFunction data;
    private IDataInfo dataInfo;
    /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    final IVector profileVector;
    final IVector position;
    /**
     * Meter that defines the property being profiled.
     */
    DataSourceAtomic meter;
    protected final DataTag tag;
    
    
    private final AtomIteratorPhaseDependent ai1 = new AtomIteratorLeafAtoms();
    
}//end of MeterProfile
