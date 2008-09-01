package etomica.data.meter;
import etomica.api.IAtomSet;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISpecies;
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
 * Data is averaged in each bin as
 * sum(values)/numAtoms
 * 
 * So, bins that don't have atoms will have an average of NaN.  The averaging
 * is done within this class, so use of an AccumulatorAverage is not needed.
 * The data returned by this class /is/ the average.  To reset the average,
 * this class' reset method should be called.
 * 
 * @author Rob Riggleman
 * @author Andrew Schultz
 */
public class MeterProfileByAtoms implements DataSource, DataSourceIndependent, java.io.Serializable {
    
    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterProfileByAtoms(ISpace space) {
        xDataSource = new DataSourceUniform("x", Length.DIMENSION);
        nAtoms = new int[xDataSource.getNValues()];
        y = new double[xDataSource.getNValues()];
        position = space.makeVector();
        tag = new DataTag();
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
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
        if (!(m.getAtomDataInfo() instanceof DataInfoDouble)) {
            throw new IllegalArgumentException("data source must return a DataDouble");
        }
        meter = m;
        reset();
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
        reset();
    }

    /**
     * Returns the profile for the current configuration.
     */
    public Data getData() {
        IBoundary boundary = box.getBoundary();
        IAtomSet moleculeList = box.getMoleculeList();
        if (species != null) {
            moleculeList = box.getMoleculeList(species);
        }
        int nMolecules = moleculeList.getAtomCount();
        for (int iMolecule=0; iMolecule<nMolecules; iMolecule++) {
            IMolecule a = (IMolecule)moleculeList.getAtom(iMolecule);
            double value = ((DataDouble)meter.getData(a)).x;
            position.E(a.getType().getPositionDefinition().position(a));
            position.PE(boundary.centralImage(position));
            int i = xDataSource.getIndex(position.x(profileDim));
            y[i] += value;
            nAtoms[i]++;
        }
        double[] yData = data.getData();
        for (int i=0; i<y.length; i++) {
            yData[i] = y[i] / nAtoms[i];
        }
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
    }

    public ISpecies getSpecies() {
        return species;
    }

    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }

    public void reset() {
        if (box == null) return;
        
        double halfBox = 0.5*box.getBoundary().getDimensions().x(profileDim);
        xDataSource.setXMin(-halfBox);
        xDataSource.setXMax(halfBox);
        nAtoms = new int[xDataSource.getNValues()];
        y = new double[xDataSource.getNValues()];
        
        if (meter != null) {
            data = new DataFunction(new int[] {xDataSource.getNValues()});
            dataInfo = new DataInfoFunction(meter.getAtomDataInfo().getLabel()+" Profile", meter.getAtomDataInfo().getDimension(), this);
            dataInfo.addTag(meter.getTag());
            dataInfo.addTag(tag);
        }
    }

    private static final long serialVersionUID = 1L;
    protected IBox box;
    protected DataSourceUniform xDataSource;
    protected DataFunction data;
    protected double[] y;
    protected int[] nAtoms;
    protected IDataInfo dataInfo;
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
    protected ISpecies species;
}
