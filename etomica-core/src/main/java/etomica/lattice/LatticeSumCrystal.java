/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;

import etomica.space.Vector;
import etomica.data.FunctionData;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataGroup;

public class LatticeSumCrystal {

    public LatticeSumCrystal(BravaisLatticeCrystal lattice) {
        this.lattice = lattice;
        spaceDim = lattice.getSpace().D();
        IndexIteratorTriangularPermutations iteratorT = new IndexIteratorTriangularPermutations(spaceDim);
        setIterator(new IndexIteratorReflecting(iteratorT));
        coreIterator = iteratorT.getCoreIterator();
        setMaxLatticeShell(50);
        kVector = lattice.getSpace().makeVector();
        dr = lattice.getSpace().makeVector();
        siteIndex = new int[lattice.D()];//lattice.D() should be spaceDim+1
        
        //get coordinates of basis at the origin
        basisDim = lattice.getBasis().getScaledCoordinates().length;
        basis0 = new Vector[basisDim];
        for(int j=0; j<basisDim; j++) {
            siteIndex[spaceDim] = j;
            basis0[j] = lattice.getSpace().makeVector();
            basis0[j].E((Vector)lattice.site(siteIndex));
        }
        
    }

    public DataGroup calculateSum(FunctionData<Object> function) {
        IDataInfo dataInfo = function.getDataInfo();
        IData work = dataInfo.makeData();
        IData[][] sumR = new IData[basisDim][basisDim];
        IData[][] sumI = new IData[basisDim][basisDim];
        for(int jp=0; jp<basisDim; jp++) {
            for(int j=0; j<basisDim; j++) {
                sumR[jp][j] = dataInfo.makeData();
                sumI[jp][j] = dataInfo.makeData();
            }
        }
        
        //interactions among sites in origin cell
        //do all pairs twice, contributing once to each site of pair
        //cell-cell distance is zero, so exp(I k.r) = 1 + 0I
        for(int jp=0; jp<basisDim; jp++) {
            for(int j=0; j<basisDim; j++) {
                if(jp==j) continue;
                dr.Ev1Mv2(basis0[jp], basis0[j]);
                sumR[jp][j].PE(function.f(dr));//don't try to add to sumR[i0][i] for efficiency because we don't know what function does for -dr
            }
        }
        //loop over shells
        for(int m=1; m<=maxLatticeShell; m++) {
            //loop over cells in shell
            coreIterator.setMaxElement(m);
            coreIterator.setMaxElementMin(m);
            iterator.reset();
            while(iterator.hasNext()) {
                System.arraycopy(iterator.next(), 0, siteIndex, 0, spaceDim);
                double ckr = 0;
                double skr = 0;
                //loop over sites in lattice cell
                for(int jp=0; jp<basisDim; jp++) {
                    siteIndex[spaceDim] = jp;
                    Vector site = (Vector)lattice.site(siteIndex);
                    //loop over sites in origin cell
                    for(int j=0; j<basisDim; j++) {
                        dr.Ev1Mv2(site, basis0[j]);
                        if(jp == 0 && j == 0) {
                            //define cell distance as distance between 0th sites in cells
                            double kDotr = kVector.dot(dr);
                            ckr = Math.cos(kDotr);
                            skr = Math.sin(kDotr);
                        }
                        IData value = function.f(dr);
                        work.E(value);
                        work.TE(ckr);
                        sumR[jp][j].PE(work);
                        work.E(value);
                        work.TE(skr);
                        sumI[jp][j].PE(work);
                    }
                }
            }
        }
        
        return new DataGroupLSC(sumR, sumI);
    }
    
    public void setK(Vector k) {
        kVector.E(k);
    }
    
    public Vector getK() {
        return kVector;
    }
    
    public IndexIterator getIterator() {
        return iterator;
    }
    private void setIterator(IndexIterator iterator) {
        if(iterator.getD() != lattice.getSpace().D()) {
            throw new IllegalArgumentException("Given iterator does not produce index arrays of a dimension consistent with the lattice: iterator.getD() = "+iterator.getD()+"; lattice.D() = "+lattice.D());
        }
        this.iterator = iterator;
    }
        
    public SpaceLattice getLattice() {
        return lattice;
    }

    public int getMaxLatticeShell() {
        return maxLatticeShell;
    }

    /**
     * Specifies the largest index element that will be included in the lattice sum.
     * Default is 50. 
     */
    public void setMaxLatticeShell(int maxElement) {
        this.maxLatticeShell = maxElement;
    }

    private final BravaisLatticeCrystal lattice;
    private IndexIterator iterator;
    private IndexIteratorTriangular coreIterator;
    private final Vector kVector;
    private final Vector[] basis0;
    private final int[] siteIndex;
    private final Vector dr;
    private final int basisDim;
    private final int spaceDim;
    private int maxLatticeShell;
    
    /**
     * Helper class that encapsulates the complex basis-basis data in a manner that
     * permits easy access to real and imaginary components for a given pair of basis
     * elements.
     */
    public class DataGroupLSC extends DataGroup {
        private DataGroupLSC(IData[][] sumR, IData[][] sumI) {
            super(makeDataArray(sumR, sumI));
        }
        
        /**
         * Returns real part of complex data.
         * @param j index of basis element in origin cell
         * @param jp index of basis element in lattice cell
         */
        public IData getDataReal(int j, int jp) {
            return ((DataGroup)((DataGroup)this.getData(jp)).getData(j)).getData(0);
        }
        
        /**
         * Returns imaginary part of complex data.
         * @param j index of basis element in origin cell
         * @param jp index of basis element in lattice cell
         */
        public IData getDataImaginary(int j, int jp) {
            return ((DataGroup)((DataGroup)this.getData(jp)).getData(j)).getData(1);
        }

        private static final long serialVersionUID = 1L;
        
    }
    
    //used by DataGroupLSC constructor
    private static IData[] makeDataArray(IData[][] sumR, IData[][] sumI) {
        int basisDim = sumR.length;
        DataGroup[] dataArray = new DataGroup[basisDim];
        
        DataGroup[] data = new DataGroup[basisDim];
        for(int jp=0; jp<basisDim; jp++) {
            for(int j=0; j<basisDim; j++) {
                data[j] = new DataGroup(new IData[] {sumR[jp][j], sumI[jp][j]});
            }
            dataArray[jp] = new DataGroup(data);
        }
        return dataArray;
    }


}
