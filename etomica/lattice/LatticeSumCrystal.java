package etomica.lattice;

import etomica.data.Data;
import etomica.data.IDataInfo;
import etomica.data.types.DataGroup;
import etomica.space.IVector;
import etomica.util.FunctionGeneral;

public class LatticeSumCrystal {

    public LatticeSumCrystal(BravaisLatticeCrystal lattice) {
        this.lattice = lattice;
        spaceDim = lattice.getSpace().D();
        IndexIteratorTriangularPermutations iteratorT = new IndexIteratorTriangularPermutations(spaceDim);
        setIterator(new IndexIteratorReflecting(iteratorT));
        coreIterator = iteratorT.getCoreIterator();
        setMaxElement(50);
        kVector = lattice.getSpace().makeVector();
        dr = lattice.getSpace().makeVector();
        siteIndex = new int[lattice.D()];//lattice.D() should be spaceDim+1
        
        //get coordinates of basis at the origin
        basisDim = lattice.getBasis().getScaledCoordinates().length;
        basis0 = new IVector[basisDim];
        for(int i=0; i<basisDim; i++) {
            siteIndex[spaceDim] = i;
            basis0[i] = lattice.getSpace().makeVector();
            basis0[i].E((IVector)lattice.site(siteIndex));
        }
        
    }

    public DataGroup calculateSum(FunctionGeneral function) {
        IDataInfo dataInfo = function.getDataInfo();
        Data work = dataInfo.makeData();
        Data[][] sumR = new Data[basisDim][basisDim];
        Data[][] sumI = new Data[basisDim][basisDim];
        for(int i=0; i<basisDim; i++) {
            for(int j0=0; j0<basisDim; j0++) {
                sumR[i][j0] = dataInfo.makeData();
                sumI[i][j0] = dataInfo.makeData();
            }
        }
        
        //interactions among sites in origin cell
        //do all pairs twice, contributing once to each site of pair
        //cell-cell distance is zero, so exp(I k.r) = 1 + 0I
        for(int i=0; i<basisDim; i++) {
            for(int j0=0; j0<basisDim; j0++) {
                if(i==j0) continue;
                dr.Ev1Mv2(basis0[i], basis0[j0]);
                sumR[i][j0].PE(function.f(dr));//don't try to add to sumR[i0][i] for efficiency because we don't know what function does for -dr
            }
        }
        //loop over shells
        for(int m=1; m<=maxElement; m++) {
            //loop over cells in shell
            coreIterator.setMaxElement(m);
            coreIterator.setMaxElementMin(m);
            iterator.reset();
            while(iterator.hasNext()) {
                System.arraycopy(iterator.next(), 0, siteIndex, 0, spaceDim);
                double ckr = 0;
                double skr = 0;
                //loop over sites in lattice cell
                for(int i=0; i<basisDim; i++) {
                    siteIndex[spaceDim] = i;
                    IVector site = (IVector)lattice.site(siteIndex);
                    //loop over sites in origin cell
                    for(int j0=0; j0<basisDim; j0++) {
                        dr.Ev1Mv2(site, basis0[j0]);
                        if(i == 0 && j0 == 0) {
                            //define cell distance as distance between 0th sites in cells
                            double kDotr = kVector.dot(dr);
                            ckr = Math.cos(kDotr);
                            skr = Math.sin(kDotr);
                        }
                        Data value = function.f(dr);
                        work.E(value);
                        work.TE(ckr);
                        sumR[i][j0].PE(work);
                        work.E(value);
                        work.TE(skr);
                        sumI[i][j0].PE(work);
                    }
                }
            }
        }
        Data[] dataArr = new Data[2*basisDim*basisDim];
        int count = 0;
        for(int i=0; i<basisDim; i++) {
            for(int j0=0; j0<basisDim; j0++) {
                dataArr[count] = sumR[i][j0];
                dataArr[basisDim*basisDim+count] = sumI[i][j0];
                count++;
            }
        }
        
        //real in first half of group, imaginary in second half
        return new DataGroup(dataArr);
    }
    
    public void setK(IVector k) {
        kVector.E(k);
    }
    
    public IVector getK() {
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

    public int getMaxElement() {
        return maxElement;
    }

    /**
     * Specifies the largest index element that will be included in the lattice sum.
     * Default is 50. 
     */
    public void setMaxElement(int maxElement) {
        this.maxElement = maxElement;
    }

    private final BravaisLatticeCrystal lattice;
    private IndexIterator iterator;
    private IndexIteratorTriangular coreIterator;
    private final IVector kVector;
    private final IVector[] basis0;
    private final int[] siteIndex;
    private final IVector dr;
    private final int basisDim;
    private final int spaceDim;
    private int maxElement;
}
