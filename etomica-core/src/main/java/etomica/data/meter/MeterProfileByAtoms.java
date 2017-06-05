/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.space.Vector;
import etomica.atom.MoleculePositionGeometricCenter;
import etomica.atom.IMoleculePositionDefinition;
import etomica.data.DataSourceIndependent;
import etomica.data.DataSourceMolecular;
import etomica.data.DataSourceUniform;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.space.Space;
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
public class MeterProfileByAtoms implements IEtomicaDataSource, DataSourceIndependent, java.io.Serializable {
    
    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterProfileByAtoms(Space space) {
        xDataSource = new DataSourceUniform("x", Length.DIMENSION);
        nAtoms = new int[xDataSource.getNValues()];
        y = new double[xDataSource.getNValues()];
        position = space.makeVector();
        tag = new DataTag();
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
        positionDefinition = new MoleculePositionGeometricCenter(space);
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
        Boundary boundary = box.getBoundary();
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

    public DataTag getIndependentTag() {
        return xDataSource.getTag();
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
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
        nAtoms = new int[xDataSource.getNValues()];
        y = new double[xDataSource.getNValues()];
        
        if (meter != null) {
            data = new DataFunction(new int[] {xDataSource.getNValues()});
            dataInfo = new DataInfoFunction(meter.getMoleculeDataInfo().getLabel()+" Profile", meter.getMoleculeDataInfo().getDimension(), this);
            dataInfo.addTag(meter.getTag());
            dataInfo.addTag(tag);
        }
    }

    public void setPositionDefinition(IMoleculePositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }

    public IMoleculePositionDefinition getPositionDefinition() {
        return positionDefinition;
    }

    private static final long serialVersionUID = 1L;
    protected Box box;
    protected DataSourceUniform xDataSource;
    protected DataFunction data;
    private IMoleculePositionDefinition positionDefinition;
    protected double[] y;
    protected int[] nAtoms;
    protected IEtomicaDataInfo dataInfo;
    /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    protected int profileDim;
    protected final Vector position;
    /**
     * Meter that defines the property being profiled.
     */
    protected DataSourceMolecular meter;
    protected final DataTag tag;
    protected ISpecies species;
}
