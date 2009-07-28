package etomica.data.meter;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IData;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.api.IVectorMutable;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.IAtomPositionDefinition;
import etomica.data.DataSourceIndependent;
import etomica.data.DataSourceMolecular;
import etomica.data.DataSourceUniform;
import etomica.data.DataTag;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
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
 * Data is averaged as
 * sum(values)/dV
 * 
 * where dV is the volume of the bin.  As such, this meter is measuring the
 * "density" of the quantity returned by the atomic meter.
 * 
 * @author Rob Riggleman
 * @author Andrew Schultz
 */
public class MeterProfileByVolume implements IEtomicaDataSource, DataSourceIndependent, java.io.Serializable {
    
    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterProfileByVolume(ISpace space) {
        xDataSource = new DataSourceUniform("x", Length.DIMENSION);
        position = space.makeVector();
        tag = new DataTag();
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
        positionDefinition = new AtomPositionGeometricCenter(space);
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }

    /**
     * The meter that defines the profiled quantity
     */
    public DataSourceMolecular getDataSource() {return meter;}
    
    /**
     * Accessor method for the meter that defines the profiled quantity.
     */
    public void setDataSource(DataSourceMolecular m) {
        if (!(m.getMoleculeDataInfo() instanceof DataInfoDouble)) {
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
    public IData getData() {
        IBoundary boundary = box.getBoundary();
        data.E(0);
        double[] y = data.getData();
        IMoleculeList moleculeList = box.getMoleculeList();
        if (species != null) {
            moleculeList = box.getMoleculeList(species);
        }
        int nMolecules = moleculeList.getMoleculeCount();
        for (int iMolecule=0; iMolecule<nMolecules; iMolecule++) {
            IMolecule a = moleculeList.getMolecule(iMolecule);
            double value = ((DataDouble)meter.getData(a)).x;
            position.E(positionDefinition.position(a));
            position.PE(boundary.centralImage(position));
            int i = xDataSource.getIndex(position.getX(profileDim));
            y[i] += value;
        }
        double dV = (xDataSource.getXMax() - xDataSource.getXMin())/y.length;
        for (int i=0; i<boundary.getBoxSize().getD(); i++) {
            if (i==profileDim) continue;
            dV *= boundary.getBoxSize().getX(i);
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
    }
    
    public ISpecies getSpecies() {
        return species;
    }

    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }

    public void reset() {
        if (box == null) return;
        
        double halfBox = 0.5*box.getBoundary().getBoxSize().getX(profileDim);
        xDataSource.setXMin(-halfBox);
        xDataSource.setXMax(halfBox);
        
        if (meter != null) {
            data = new DataFunction(new int[] {xDataSource.getNValues()});
            dataInfo = new DataInfoFunction(meter.getMoleculeDataInfo().getLabel()+" Profile", meter.getMoleculeDataInfo().getDimension(), this);
            dataInfo.addTag(meter.getTag());
            dataInfo.addTag(tag);
        }
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

    public void setPositionDefinition(IAtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }

    public IAtomPositionDefinition getPositionDefinition() {
        return positionDefinition;
    }

    private static final long serialVersionUID = 1L;
    private IBox box;
    private DataSourceUniform xDataSource;
    private DataFunction data;
    private IEtomicaDataInfo dataInfo;
    /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    protected int profileDim;
    protected final IVectorMutable position;
    /**
     * Meter that defines the property being profiled.
     */
    protected DataSourceMolecular meter;
    protected final DataTag tag;
    protected ISpecies species;
    private IAtomPositionDefinition positionDefinition;
}
