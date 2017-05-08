/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveGeneral;
import etomica.space.ISpace;
import etomica.units.Null;

/**
 * Meter for calculation of structure factor of atoms for all wave vectors less
 * than a cutoff.
 *
 * @author Michael Sellers
 * @author Andrew Schultz
 */
public class MeterStructureFactor implements IEtomicaDataSource, DataSourceIndependent {
	
	protected final ISpace space;
    protected IBox box;
    protected double[] struct;
    protected IVectorMutable [] waveVec;
    protected IAtomList atomList;
    protected DataFunction data;
    protected DataInfoFunction dataInfo;
    protected final DataTag tag, xTag;
    protected DataDoubleArray xData;
    protected DataInfoDoubleArray xDataInfo;

    /**
     * Creates meter with default to compute the structure factor for all atoms
     * in the box.  All wave vectors consistent with the box shape and with
     * magnitude less than cutoff are included.
     */
	public MeterStructureFactor(ISpace space, IBox aBox, double cutoff) {
	    this.space = space;
	    this.box = aBox;
        atomList = box.getLeafList();
        tag = new DataTag();
        xTag = new DataTag();
	    setCutoff(cutoff);
	}
	
	protected void resetData() {
        xData = new DataDoubleArray(waveVec.length);
        xDataInfo = new DataInfoDoubleArray("q", Null.DIMENSION, new int[]{waveVec.length});
        xDataInfo.addTag(xTag);

	    dataInfo = new DataInfoFunction("Structure Factor", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        data = new DataFunction(new int[]{waveVec.length}, struct);
	}

	protected int makeWaveVector(double cutoff) {
        int nVec = 0;
        double[] x = xData == null ? null : xData.getData();
        IVector[] edges = new IVector[3];
        edges[0] = box.getBoundary().getEdgeVector(0);
        edges[1] = box.getBoundary().getEdgeVector(1);
        edges[2] = box.getBoundary().getEdgeVector(2);
        Primitive primitiveBox = new PrimitiveGeneral(space, edges);
        Primitive recip = primitiveBox.makeReciprocal();
        IVector[] basis = recip.vectors();

        double cutoff2 = cutoff*cutoff;

        int[] iMax = new int[space.D()];
        // Be aggressive when look for wave vectors.  If the box is slanty,
        // we will need to go beyond cutoff/basis, but it's hard to know how
        // much.
        for (int i=0; i<space.D(); i++) {
            iMax[i] = 1+2*(int)(cutoff/Math.sqrt(basis[i].squared()));
        }

        int[] idx = new int[space.D()];
        idx[0] = 0;
        idx[1] = 0;
        idx[2] = 1;
        while (true) {
            IVectorMutable v = space.makeVector();
            boolean success = false;
            for  (int i=idx.length-1; i>=0; i--) {
                idx[i]++;
                if (idx[i] <= iMax[i]) {
                    success = true;
                    break;
                }
                
                idx[i] = -iMax[i];
            }
            if (!success) break;
            v.E(0);
            for (int i=0; i<idx.length; i++) {
                v.PEa1Tv1(idx[i], basis[i]);
            }
            double foo = Math.sqrt(v.squared());
            if (v.squared() > cutoff2) {
                continue;
            }
            if (waveVec != null) {
                waveVec[nVec] = v;
                x[nVec] = Math.sqrt(v.squared());
            }
            nVec++;
        }
        return nVec;
	}

    /**
     * Sets the wave vector cutoff.  All wave vectors consistent with the box
     * shape that have a magnitude less than the cutoff will be computed.
     * @param cutoff the cutoff for the wave vector magnitude
     */
	public void setCutoff(double cutoff) {
	    waveVec = null;
	    int nVec = makeWaveVector(cutoff);
        struct = new double[nVec];
	    waveVec = new IVectorMutable[nVec];
        resetData();
        makeWaveVector(cutoff);
	}
	
	/**
	 * @param waveVec Sets a custom wave vector array.
	 */
	public void setWaveVec(IVectorMutable [] waveVec){
	    this.waveVec = space.makeVectorArray(waveVec.length);
	    struct = new double[waveVec.length];
		for(int i=0; i<waveVec.length; i++){
			this.waveVec[i].E(waveVec[i]);
		}
		resetData();
	}
	
	/**
	 * @param atomList Sets the list of atoms for factor calculation.
	 */
	public void setAtoms(IAtomList atomList){
		this.atomList = atomList;
	}

    public IData getData() {
        double term1 = 0;
        double term2 = 0;
        double dotprod = 0;
        int numAtoms = atomList.getAtomCount();
        for(int k=0; k<waveVec.length; k++){
            term1 = 0;
            term2 = 0;
            dotprod = 0;
            for(int i=0; i<numAtoms; i++){
                dotprod = waveVec[k].dot(atomList.getAtom(i).getPosition());
                term1 += Math.cos(dotprod);
                term2 += Math.sin(dotprod);
            }
            struct[k] = ((term1*term1) + (term2*term2))/(numAtoms*numAtoms);
        }
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo;
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataTag getIndependentTag() {
        return xTag;
    }

}
