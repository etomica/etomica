/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Length;

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
public class MeterProfileByVolume implements IDataSource, DataSourceIndependent, java.io.Serializable {
    
    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterProfileByVolume(Space space) {
        xDataSource = new DataSourceUniform("x", Length.DIMENSION);
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
            if (i>=0 && i<y.length) y[i] += value;
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
        
        Boundary boundary = box.getBoundary();
        double halfBox = 0.5*boundary.getBoxSize().getX(profileDim);
        xDataSource.setXMin(-halfBox);
        xDataSource.setXMax(halfBox);
        
        if (meter != null) {
            data = new DataFunction(new int[] {xDataSource.getNValues()});
            dataInfo = new DataInfoFunction(meter.getMoleculeDataInfo().getLabel()+" Profile", meter.getMoleculeDataInfo().getDimension(), this);
            dataInfo.addTag(meter.getTag());
            dataInfo.addTag(tag);

            dV = 2*halfBox/data.getLength();
            for (int i=0; i<boundary.getBoxSize().getD(); i++) {
                if (i==profileDim) continue;
                dV *= boundary.getBoxSize().getX(i);
            }
        }
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
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
    protected IMoleculePositionDefinition positionDefinition;
    protected double dV;
}
