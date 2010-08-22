package etomica.models.nitrogen;

import etomica.action.MoleculeActionTranslateTo;
import etomica.action.MoleculeChildAtomAction;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.MoleculePair;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataGroup;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.IndexIterator;
import etomica.lattice.IndexIteratorReflecting;
import etomica.lattice.IndexIteratorTriangular;
import etomica.lattice.IndexIteratorTriangularPermutations;
import etomica.lattice.SpaceLattice;
import etomica.paracetamol.AtomActionTransformed;
import etomica.simulation.Simulation;
import etomica.space.Tensor;
import etomica.util.FunctionGeneral;

/**
 * Lattice sum crystal for molecular model
 * 
 * 
 * @author taitan
 *
 */

public class LatticeSumCrystalMolecular extends Simulation{

	public LatticeSumCrystalMolecular(BravaisLatticeCrystal lattice, CoordinateDefinitionNitrogen coordinateDef, IBox box) {
        
    	super(lattice.getSpace());
    	this.lattice = lattice;
    	this.coordinateDef = coordinateDef;
		this.ghostBox = box;
		
    	xzOrientationTensor = coordinateDef.getXzOrientationTensor();
    	yOrientationTensor = coordinateDef.getyOrientationTensor();

    	atomGroupAction = new MoleculeChildAtomAction(new AtomActionTransformed(lattice.getSpace()));
         
        spaceDim = lattice.getSpace().D();
        IndexIteratorTriangularPermutations iteratorT = new IndexIteratorTriangularPermutations(spaceDim);
        setIterator(new IndexIteratorReflecting(iteratorT));
        coreIterator = iteratorT.getCoreIterator();
        
        kVector = lattice.getSpace().makeVector();
        
        dr = lattice.getSpace().makeVector();
        siteIndex = new int[lattice.D()];//lattice.D() should be spaceDim+1
        
        offset = lattice.getSpace().makeVector();
        position = lattice.getSpace().makeVector();
        
        IVector[] primitiveVectors = lattice.getPrimitive().vectors();
        for (int i=0; i<primitiveVectors.length; i++) {
            offset.PEa1Tv1(1,primitiveVectors[i]);
        }
        
        offset.TE(-0.5);
        
        //get coordinates of basis at the origin
        basisDim = lattice.getBasis().getScaledCoordinates().length;
        basis0 = new IVectorMutable[basisDim];
        moleculeCell0 = new IMolecule[basisDim];
        
        for(int j=0; j<basisDim; j++) {
            siteIndex[spaceDim] = j;
            basis0[j] = lattice.getSpace().makeVector();
            basis0[j].E((IVectorMutable)lattice.site(siteIndex));
            moleculeCell0[j] = coordinateDef.getBasisCells()[0].molecules.getMolecule(j);
   
        }
        
        atomActionTranslateTo = new MoleculeActionTranslateTo(lattice.getSpace());
        
    }

    public DataGroup calculateSum(FunctionGeneral function) {
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
        
        MoleculePair pair;
        for(int jp=0; jp<basisDim; jp++) {
            for(int j=0; j<basisDim; j++) {
                if(jp==j) continue;
                pair = new MoleculePair(moleculeCell0[j], moleculeCell0[jp]);
                sumR[jp][j].PE(function.f(pair));//don't try to add to sumR[i0][i] for efficiency because we don't know what function does for -dr
            }
        }
    	System.out.println("maxLatticeShell: " + getMaxLatticeShell());
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
                    
                    IMolecule ghostMol = ghostBox.getMoleculeList(ghostBox.getMoleculeList().getMolecule(0).getType()).getMolecule(0);
                    ghostMol.getType().initializeConformation(ghostMol);
                    
         
                    int rotationNum = jp%4;
//                    System.out.println("rotationNum: " + rotationNum);
                    ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(yOrientationTensor[rotationNum]);
	                atomGroupAction.actionPerformed(ghostMol);
	                
	                ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(xzOrientationTensor[rotationNum]);
	                atomGroupAction.actionPerformed(ghostMol);
                    
	        
//	                for(int i=0; i<moleculeCell0[jp].getChildList().getAtomCount(); i++){
//                    	System.out.println("realMol(before)["+i+"]: " + moleculeCell0[jp].getChildList().getAtom(i).getPosition().toString());
//                    }
//                    for(int i=0; i<ghostMol.getChildList().getAtomCount(); i++){
//                    	System.out.println("ghostMol(before)["+i+"]: " + ghostMol.getChildList().getAtom(i).getPosition().toString());
//                    }
                    IVectorMutable site = (IVectorMutable)lattice.site(siteIndex);
                    position.E(site);
                    position.PE(offset);
//                    System.out.println("jp: "+jp+ " position: " + position.toString());
                
                    atomActionTranslateTo.setDestination(position);
                    atomActionTranslateTo.actionPerformed(ghostMol);  
//                    for(int i=0; i<moleculeCell0[jp].getChildList().getAtomCount(); i++){
//                    	System.out.println("realMol(after)["+i+"]: " + moleculeCell0[jp].getChildList().getAtom(i).getPosition().toString());
//                    }
//                    for(int i=0; i<ghostMol.getChildList().getAtomCount(); i++){
//                    	System.out.println(jp + " ghostMol(after)["+i+"]: " + ghostMol.getChildList().getAtom(i).getPosition().toString());
//                    }
//                    System.out.println("");
//                    //System.exit(1);
                    
                    
                  //  System.out.println("siteIndex: " + siteIndex[0] + " " + siteIndex[1] + " "+ siteIndex[2] + " "+ siteIndex[3] + " ");
                    //loop over sites in origin cell
                    for(int j=0; j<basisDim; j++) {
                        dr.Ev1Mv2(site, basis0[j]);
                        if(jp == 0 && j == 0) {
                            //define cell distance as distance between 0th sites in cells
                            double kDotr = kVector.dot(dr);
                            ckr = Math.cos(kDotr);
                            skr = Math.sin(kDotr);
                        }

                        MoleculePair molPair = new MoleculePair(moleculeCell0[j], ghostMol);
                        
//                        IVectorMutable leaf0Atom0 = molPair.atom0.getChildList().getAtom(0).getPosition();
//                        IVectorMutable leaf0Atom1 = molPair.atom0.getChildList().getAtom(1).getPosition();
//                        IVectorMutable leaf1Atom0 = molPair.atom1.getChildList().getAtom(0).getPosition();
//                        IVectorMutable leaf1Atom1 = molPair.atom1.getChildList().getAtom(1).getPosition();

//                      for(int i=0; i<moleculeCell0[jp].getChildList().getAtomCount(); i++){
//                    	  System.out.println(j+" realMol["+i+"]: " + molPair.atom0.getChildList().getAtom(i).getPosition().toString());
//                      }
//                      System.out.println("");
//                      for(int i=0; i<ghostMol.getChildList().getAtomCount(); i++){
//                    	  System.out.println("ghostMol["+i+"]: " + molPair.atom1.getChildList().getAtom(i).getPosition().toString());
//                      }
                        
//                        IVectorMutable delta = space.makeVector();
//                        IVectorMutable diff = space.makeVector();
//                        IVectorMutable comLeaf0 = space.makeVector();
//                        delta.Ev1Mv2(leaf0Atom1, leaf0Atom0);
//                        comLeaf0.E(leaf0Atom0);
//                        comLeaf0.PEa1Tv1(0.5, delta);
//                        
//                        IVectorMutable comLeaf1 = space.makeVector();
//                        delta.Ev1Mv2(leaf1Atom1, leaf1Atom0);
//                        comLeaf1.E(leaf1Atom0);
//                        comLeaf1.PEa1Tv1(0.5, delta);
//                        diff.Ev1Mv2(comLeaf1, comLeaf0);
//                        
//                        System.out.println("\ndiff: " + diff.toString());
//                        System.out.println("dr: " + dr.toString());
                    
                        IData value = function.f(molPair);
                        //System.out.println("energy: " + value.getValue(0));
                        work.E(value);
                        work.TE(ckr);
                        sumR[jp][j].PE(work);
                        work.E(value);
                        work.TE(skr);
                        sumI[jp][j].PE(work);
                        
                    }
                }
               // System.exit(1);
            }
            //System.exit(1);
        }
        
        return new DataGroupLSC(sumR, sumI);
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
    
	private static final long serialVersionUID = 1L;
    private final BravaisLatticeCrystal lattice;
    private IndexIterator iterator;
    private IndexIteratorTriangular coreIterator;
    private final IVectorMutable kVector;
    private final IVectorMutable[] basis0;
    private final int[] siteIndex;
    private final IVectorMutable dr;
    private final int basisDim;
    private final int spaceDim;
    private int maxLatticeShell;
    protected final IMolecule[] moleculeCell0;
    protected IVectorMutable offset, position;
    protected final MoleculeActionTranslateTo atomActionTranslateTo;
    protected CoordinateDefinitionNitrogen coordinateDef;
    protected IBox ghostBox;
    protected Tensor[] xzOrientationTensor, yOrientationTensor;
    protected MoleculeChildAtomAction atomGroupAction;
    
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
