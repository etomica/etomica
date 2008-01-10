package etomica.paracetamol;

import etomica.action.AtomGroupAction;
import etomica.atom.AtomArrayList;
import etomica.data.Data;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIterator;
import etomica.lattice.IndexIteratorReflecting;
import etomica.lattice.IndexIteratorTriangular;
import etomica.lattice.IndexIteratorTriangularPermutations;
import etomica.lattice.SpaceLattice;
import etomica.space.IVector;

public class LatticeSumCrystalParacetamol {

    public LatticeSumCrystalParacetamol(BravaisLatticeCrystal lattice) {
        this.lattice = lattice;
        spaceDim = lattice.getSpace().D();
        IndexIteratorTriangularPermutations iteratorT = new IndexIteratorTriangularPermutations(spaceDim);
        setIterator(new IndexIteratorReflecting(iteratorT));
        coreIterator = iteratorT.getCoreIterator();
        setMaxLatticeShell(50);
        kVector = lattice.getSpace().makeVector();
//        dr = lattice.getSpace().makeVector();
        siteIndex = new int[lattice.D()];//lattice.D() should be spaceDim+1

        atomGroupAction = new AtomGroupAction(new AtomActionTransformed(lattice.getSpace()));
        atomArrayList = new AtomArrayList();
        
        //get coordinates of basis at the origin
        System.out.println("At  LatticeSumCrystalParacetamol Constructor");
        basisDim = lattice.getBasis().getScaledCoordinates().length;
        basis0 = new IVector[basisDim];
        

        for(int j=0; j<basisDim; j++) {
            siteIndex[spaceDim] = j;
            basis0[j] = lattice.getSpace().makeVector();
            basis0[j].E((IVector)lattice.site(siteIndex));
        }
        
    }

    public DataGroup calculateSum(LatticeEnergyParacetamol lattice2ndDerivative) {
        IDataInfo dataInfo = lattice2ndDerivative.getDataInfo();
    	Data work = dataInfo.makeData();
        Data[][] sumR = new Data[basisDim][basisDim];
        Data[][] sumI = new Data[basisDim][basisDim];
        for(int jp=0; jp<basisDim; jp++) {
            for(int j=0; j<basisDim; j++) {
                sumR[jp][j] = dataInfo.makeData();
                sumI[jp][j] = dataInfo.makeData();
            }
        }
        
        int singleMolDim = coordinateDefinitionParacetamol.getCoordinateDim()
        				  /coordinateDefinitionParacetamol.getBasisCells()[0].molecules.getAtomCount();
        DataDoubleArray u = new DataDoubleArray(singleMolDim);
        
        /*
         *interactions among sites in origin cell
         *do all pairs twice, contributing once to each site of pair
         *cell-cell distance is zero, so exp(I k.r) = 1 + 0I
         *
         */
        for(int jp=0; jp<basisDim; jp++) {
            for(int j=0; j<basisDim; j++) {
                if(jp==j) continue;

                u.E(0);
                System.out.println ("In LatticeSumCrystalParacetmol calculateSum function!");
                
                //need to assign function to the particular molecules that I am looking?
	    		lattice2ndDerivative.setMolecule(j, jp);
	    		sumR[jp][j].PE(lattice2ndDerivative.getData());//don't try to add to sumR[i0][i] for efficiency because we don't know what function does for -dr
            }
        }
        
        //loop over shells
        for(int m=1; m<=maxLatticeShell; m++) {
            //loop over cells in shell
            coreIterator.setMaxElement(m);
            coreIterator.setMaxElementMin(m);
            iterator.reset();
            
//            while(iterator.hasNext()) {
//                System.arraycopy(iterator.next(), 0, siteIndex, 0, spaceDim);
//                double ckr = 0;
//                double skr = 0;
//                //loop over sites in lattice cell
//                for(int jp=0; jp<basisDim; jp++) {
//                    siteIndex[spaceDim] = jp;
//                    IVector site = (IVector)lattice.site(siteIndex);
//                    //loop over sites in origin cell
//                    for(int j=0; j<basisDim; j++) {
//                    	dr.Ev1Mv2(site, basis0[j]);
//                        if(jp == 0 && j == 0) {
//                            //define cell distance as distance between 0th sites in cells
//                            double kDotr = kVector.dot(dr);
//                            ckr = Math.cos(kDotr);
//                            skr = Math.sin(kDotr);
//                        }
//                        Data value = function.f(dr);
//                        work.E(value);
//                        work.TE(ckr);
//                        sumR[jp][j].PE(work);
//                        work.E(value);
//                        work.TE(skr);
//                        sumI[jp][j].PE(work);
//                    }
//                }
//            }
        }
        
        System.out.println ("In LatticeSumCrystalParacetmol calculateSum function!");
        
        return new DataGroupLSCParacetamol(sumR, sumI);
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
    
	public CoordinateDefinitionParacetamol getCoordinateDefinitionParacetamol() {
		return coordinateDefinitionParacetamol;
	}

	public void setCoordinateDefinitionParacetamol(
			CoordinateDefinitionParacetamol coordinateDefinitionParacetamol) {
		this.coordinateDefinitionParacetamol = coordinateDefinitionParacetamol;
	}
	
    private CoordinateDefinitionParacetamol coordinateDefinitionParacetamol;
    private final BravaisLatticeCrystal lattice;
    private IndexIterator iterator;
    private IndexIteratorTriangular coreIterator;
    private final IVector kVector;
    private final IVector[] basis0;
    private final int[] siteIndex;
//    private final IVector dr;
    private final int basisDim;
    private final int spaceDim;
    private int maxLatticeShell;
    protected double[][] basisOrientation;
    protected final AtomGroupAction atomGroupAction;
    protected final AtomArrayList atomArrayList;
    
    /**
     * Helper class that encapsulates the complex basis-basis data in a manner that
     * permits easy access to real and imaginary components for a given pair of basis
     * elements.
     */
    public class DataGroupLSCParacetamol extends DataGroup {
        private DataGroupLSCParacetamol(Data[][] sumR, Data[][] sumI) {
            super(makeDataArray(sumR, sumI));
        }
        
        /**
         * Returns real part of complex data.
         * @param j index of basis element in origin cell
         * @param jp index of basis element in lattice cell
         */
        public Data getDataReal(int j, int jp) {
            return ((DataGroup)((DataGroup)this.getData(jp)).getData(j)).getData(0);
        }
        
        /**
         * Returns imaginary part of complex data.
         * @param j index of basis element in origin cell
         * @param jp index of basis element in lattice cell
         */
        public Data getDataImaginary(int j, int jp) {
            return ((DataGroup)((DataGroup)this.getData(jp)).getData(j)).getData(1);
        }

        private static final long serialVersionUID = 1L;
        
    }
    
    //used by DataGroupLSCParacetamol constructor
    private static Data[] makeDataArray(Data[][] sumR, Data[][] sumI) {
        int basisDim = sumR.length;
        DataGroup[] dataArray = new DataGroup[basisDim];
        
        DataGroup[] data = new DataGroup[basisDim];
        for(int jp=0; jp<basisDim; jp++) {
            for(int j=0; j<basisDim; j++) {
                data[j] = new DataGroup(new Data[] {sumR[jp][j], sumI[jp][j]});
            }
            dataArray[jp] = new DataGroup(data);
        }
        return dataArray;
    }

}
