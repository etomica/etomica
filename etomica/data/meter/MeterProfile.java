package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.api.IAtomSet;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IVector;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataSourceAtomic;
import etomica.data.DataSourceIndependent;
import etomica.data.DataSourceUniform;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.space.ISpace;
import etomica.units.Length;

/**
 * Meter that takes a (scalar) Meter and records its property as a
 * 1-dimensional function of position in the simulation volume. The measured
 * property must be a quantity that can be associated with a single molecule.
 * The position coordinate lies along one dimension (x,y,z).
 * 
 * @author Rob Riggleman
 */
public class MeterProfile implements DataSource, DataSourceIndependent, java.io.Serializable {
    
    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterProfile(ISpace space) {
        xDataSource = new DataSourceUniform("x", Length.DIMENSION);
        position = space.makeVector();
        tag = new DataTag();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Breaks a meter's measurements into a profile taken along some direction in box");
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
    public int getProfileDim() {return profileDim;}
    
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     * The given vector is converted to a unit vector, if not already.
     */
    public void setProfileDim(int dim) {
        profileDim = dim;
        double halfBox = 0.5*box.getBoundary().getDimensions().x(dim);
        xDataSource.setXMin(-halfBox);
        xDataSource.setXMax(halfBox);
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
    }
    
    /**
     * Returns the profile for the current configuration.
     */
    public Data getData() {
        IBoundary boundary = box.getBoundary();
        data.E(0);
        double[] y = data.getData();
        IAtomSet moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.getAtomCount();
        for (int iMolecule=0; iMolecule<nMolecules; iMolecule++) {
            IMolecule a = (IMolecule)moleculeList.getAtom(iMolecule);
            double value = ((DataDouble)meter.getData(a)).x;
            position.E(a.getType().getPositionDefinition().position(a));
            position.PE(boundary.centralImage(position));
            int i = xDataSource.getIndex(position.x(profileDim));
            y[i] += value;
        }
        double dV = (xDataSource.getXMax() - xDataSource.getXMin())/y.length;
        for (int i=0; i<boundary.getDimensions().getD(); i++) {
            if (i==profileDim) continue;
            dV *= boundary.getDimensions().x(i);
        }
        data.TE(1.0/dV);
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
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(IBox box) {
        this.box = box;
        double halfBox = 0.5*box.getBoundary().getDimensions().x(profileDim);
        xDataSource.setXMin(-halfBox);
        xDataSource.setXMax(halfBox);
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    private static final long serialVersionUID = 1L;
    private IBox box;
    private String name;
    private DataSourceUniform xDataSource;
    private DataFunction data;
    private IDataInfo dataInfo;
    /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    protected int profileDim;
    protected final IVector position;
    /**
     * Meter that defines the property being profiled.
     */
    protected DataSourceAtomic meter;
    protected final DataTag tag;
}
