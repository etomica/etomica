/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomsetIteratorBoxDependent;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Angle;
import etomica.units.dimensions.Null;

/**
 * Meter for tabulation of a dihedral angle distribution between nearest neighbors.  The
 * meter takes data via actionPerformed and returns the average dihedral angle via
 * getData.
 *
 * @author Michael Sellers, adapted from David Kofke's MeterRDF
 */
public class MeterDihedralAngle implements IAction, IDataSource, DataSourceIndependent, java.io.Serializable {

    /**
	 * Creates meter with default to compute dihedral angle for all
	 * leaf atoms in a box.
	 * @param space
	 */
    public MeterDihedralAngle(Space space) {
	    this.space = space;

        xDataSource = new DataSourceUniform("phi", Angle.DIMENSION);
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
        
        phiData = (DataDoubleArray)xDataSource.getData();
        data = new DataFunction(new int[] {phiData.getLength()});
        gSum = new long[phiData.getLength()];
        dataInfo = new DataInfoFunction("g(phi)", Null.DIMENSION, this);

        dr1 = space.makeVector();
        dr2 = space.makeVector();
        dr3 = space.makeVector();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public void setRmax(double rmax){
    	this.rMax = rmax;
    }
    /**
     * Zero's out the G(phi) sum tracked by this meter.
     */
    public void reset() {
    	xDataSource.setXMin(0);
    	xDataSource.setXMax(180);
        phiData = (DataDoubleArray)xDataSource.getData();
        data = new DataFunction(new int[] {phiData.getLength()});
        gSum = new long[phiData.getLength()];
        dataInfo = new DataInfoFunction("g(phi)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        tetraNum = 0;
    }
    
    /**
     * Measures the dihedral angle's for the current configuration of the given box.
     */
    public void actionPerformed() {
        if (phiData != xDataSource.getData() ||
            data.getLength() != phiData.getLength()) {
            reset();
        }
        int atomCount = box.getLeafList().size();
        
        double rMaxSquared = rMax*rMax;
        System.out.println(rMax);
        // iterate over all tetra's
        for (int i=0; i<atomCount; i++) {
        	IAtom atom0 = box.getLeafList().get(i);
        	
        	for (int j=0; j<atomCount; j++){
        		IAtom atom1 = box.getLeafList().get(j);
        		if(atom1==atom0){continue;}
        		dr1.Ev1Mv2(atom0.getPosition(),atom1.getPosition());
        		boundary.nearestImage(dr1);
        		if(dr1.squared()>rMaxSquared){continue;}
        		for (int k=0; k<atomCount; k++){

        			IAtom atom2 = box.getLeafList().get(k);
        			if(atom2==atom1 || atom2==atom0){continue;}
            		dr2.Ev1Mv2(atom1.getPosition(),atom2.getPosition());
            		boundary.nearestImage(dr2);
            		if(dr2.squared()>rMaxSquared){;continue;}
            		for (int l=0; l<atomCount; l++){
               			IAtom atom3 = box.getLeafList().get(l);
               			if(atom3==atom2 || atom3==atom1 || atom3==atom0){continue;}
                		dr3.Ev1Mv2(atom2.getPosition(),atom3.getPosition());
                		boundary.nearestImage(dr3);
            			if(dr3.squared()>rMaxSquared){continue;}
            			//compute dihedral angle
            			Vector tanY = space.makeVector();
            			tanY.Ea1Tv1(Math.sqrt(dr2.squared()),dr1);
            			Vector cross12 = space.makeVector();
            			Vector cross23 = space.makeVector();
            			cross12.E(dr1);
            			cross12.XE(dr2);
            			cross23.E(dr2);
            			cross23.XE(dr3);
            			
            			double angle1 = Math.atan2(tanY.dot(cross23),cross12.dot(cross23));
            			//determine histogram index, and add 1.
            			angle1 = Math.abs(angle1)*180/Math.PI;
            			int index = xDataSource.getIndex(angle1); 
            			gSum[index]++;
            			tetraNum++;
            		}
        		}
        	}
            
        }
    }
    
	/**
	 * Returns the Dihedral distribution, averaged over the calls to actionPerformed since the
     * meter was reset or had some parameter changed (# of bins).
	 */
	public IData getData() {
        if (phiData != xDataSource.getData() ||
            data.getLength() != phiData.getLength() || tetraNum == 0) {
            reset();
            //that zeroed everything.  just return the zeros.
            return data;
        }
        final double[] y = data.getData();
        for(int i=0;i<gSum.length; i++) {
	        y[i] = gSum[i] * phiData.getLength() / (double)tetraNum;
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
        boundary = box.getBoundary();
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    
    private static final long serialVersionUID = 1L;
    protected Box box;
    protected final Space space;
    protected long[] gSum;
    protected DataFunction data;
    private IDataInfo dataInfo;
    protected DataDoubleArray phiData;
    protected AtomsetIteratorBoxDependent iterator;
    private final Vector dr1, dr2, dr3;
    private Boundary boundary;
    protected final DataSourceUniform xDataSource;
    protected double rMax;
    private String name;
    protected final DataTag tag;
    protected long tetraNum;
}
