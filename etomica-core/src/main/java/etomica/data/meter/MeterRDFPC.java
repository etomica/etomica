/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;
import etomica.action.IAction;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataSourceIndependent;
import etomica.data.DataSourceUniform;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.space.Space;
import etomica.units.Null;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).  The
 * meter takes data via actionPerformed and returns the average RDF via
 * getData.  This implementation uses a PotentialCalculation to perform the
 * actual work and so can take advantage of neighbor/cell lists so long as the
 * RDF range is not greater than the PotentialMaster range.
 *
 * @author David Kofke
 * @author Andrew Schultz
 */
public class MeterRDFPC implements IAction, IEtomicaDataSource, DataSourceIndependent {
	
	/**
	 * Creates meter with default to compute pair correlation for all
	 * leaf atoms in a box.
	 */
    public MeterRDFPC(Space space, PotentialMaster potentialMaster, Box box) {
	    this.space = space;
	    
	    this.potentialMaster = potentialMaster;
	    
	    pc = new PotentialCalculationRDF(space, box);

        xDataSource = pc.getXDataSource();
        
        rData = (DataDoubleArray)xDataSource.getData();
        data = new DataFunction(new int[] {rData.getLength()});
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);
        
        allAtoms = new IteratorDirective();

        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    /**
     * Zero's out the RDF sum tracked by this meter.
     */
    public void reset() {
        rData = (DataDoubleArray)xDataSource.getData();
        xMax = xDataSource.getXMax();
        data = new DataFunction(new int[] {rData.getLength()});
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        callCount = 0;
        pc.reset();
    }
    
    /**
     * Takes the RDF for the current configuration of the given box.
     */
    public void actionPerformed() {
        if (rData != xDataSource.getData() ||
                data.getLength() != rData.getLength() ||
                xDataSource.getXMax() != xMax) {
            reset();
        }
        potentialMaster.calculate(box, allAtoms, pc);
        
        callCount++;
    }
    
	/**
	 * Returns the RDF, averaged over the calls to actionPerformed since the
     * meter was reset or had some parameter changed (xMax or # of bins).
	 */
	public IData getData() {
        if (rData != xDataSource.getData() ||
            data.getLength() != rData.getLength() ||
            xDataSource.getXMax() != xMax) {
            reset();
            //that zeroed everything.  just return the zeros.
            return data;
        }
        
        final double[] y = data.getData();
        long numAtomPairs = 0;
        long numAtoms = box.getLeafList().getAtomCount();
        numAtomPairs = numAtoms*(numAtoms-1)/2;
	    double norm = numAtomPairs * callCount / box.getBoundary().volume();
	    double[] r = rData.getData();
	    double dx2 = 0.5*(xMax - xDataSource.getXMin())/r.length;
	    long[] gSum = pc.getGSum();
	    for(int i=0;i<r.length; i++) {
	        double vShell = space.sphereVolume(r[i]+dx2)-space.sphereVolume(r[i]-dx2);
	        y[i] = gSum[i] / (norm*vShell);
	    }
	    return data;
	}
    
    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }
	
    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray)xDataSource.getData();
    }
    
    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray)xDataSource.getDataInfo();
    }

    public DataTag getIndependentTag() {
        return xDataSource.getTag();
    }

    public int getIndependentArrayDimension() {
        return 1;
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
        pc.setBox(box);
        this.box = box;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    
    protected Box box;
    protected final Space space;
    protected DataFunction data;
    private IEtomicaDataInfo dataInfo;
    protected DataDoubleArray rData;
    protected IteratorDirective allAtoms;
    protected final DataSourceUniform xDataSource;
    protected double xMax;
    private String name;
    protected final DataTag tag;
    protected long callCount;
    protected final PotentialMaster potentialMaster;
    protected final PotentialCalculationRDF pc;
}
