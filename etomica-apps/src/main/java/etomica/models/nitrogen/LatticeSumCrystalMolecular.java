/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeActionTranslateTo;
import etomica.action.MoleculeChildAtomAction;
import etomica.box.Box;
import etomica.data.FunctionData;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataGroup;
import etomica.lattice.*;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePair;
import etomica.paracetamol.AtomActionTransformed;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Lattice sum crystal for molecular model
 * 
 * 
 * @author taitan
 *
 */

public class LatticeSumCrystalMolecular{

	public LatticeSumCrystalMolecular(BravaisLatticeCrystal lattice, CoordinateDefinitionNitrogen coordinateDef, Box box) {
       
    	this.lattice = lattice;
    	this.coordinateDef = coordinateDef;
		this.ghostBox = box;
		
    	xzOrientationTensor = coordinateDef.getXzOrientationTensor();
    	yOrientationTensor = coordinateDef.getyOrientationTensor();

    	if(coordinateDef.isBetaLatticeSum){
    		positionVector = coordinateDef.positionVector;
    	}
    	
    	atomGroupAction = new MoleculeChildAtomAction(new AtomActionTransformed(lattice.getSpace()));
         
        spaceDim = lattice.getSpace().D();
        IndexIteratorTriangularPermutations iteratorT = new IndexIteratorTriangularPermutations(spaceDim);
        setIterator(new IndexIteratorReflecting(iteratorT));
        coreIterator = iteratorT.getCoreIterator();
        
        kVector = lattice.getSpace().makeVector();
        
        dr = lattice.getSpace().makeVector();
        siteIndex = new int[lattice.D()];//lattice.D() should be spaceDim+1
        
        position = lattice.getSpace().makeVector();
        
        //get coordinates of basis at the origin
        basisDim = lattice.getBasis().getScaledCoordinates().length;
        basis0 = new Vector[basisDim];
        moleculeCell0 = new IMolecule[basisDim];
                
        for(int j=0; j<basisDim; j++) {
            siteIndex[spaceDim] = j;
            basis0[j] = lattice.getSpace().makeVector();
            basis0[j].E((Vector)lattice.site(siteIndex));
            moleculeCell0[j] = coordinateDef.getBasisCells()[0].molecules.get(j);
   
        }
        
        atomActionTranslateTo = new MoleculeActionTranslateTo(lattice.getSpace());
        
        translateBy = new AtomActionTranslateBy(lattice.getSpace());
        atomGroupActionTranslate = new MoleculeChildAtomAction(translateBy); 
        
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
        
        MoleculePair pair;
        for(int jp=0; jp<basisDim; jp++) {
            for(int j=0; j<basisDim; j++) {
                if(jp==j) continue;
                pair = new MoleculePair(moleculeCell0[j], moleculeCell0[jp]);
                sumR[jp][j].PE(function.f(pair));
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
                    
                    IMolecule ghostMol = ghostBox.getMoleculeList(ghostBox.getMoleculeList().get(0).getType()).get(0);
                    ghostMol.getType().initializeConformation(ghostMol);
                    
                    int rotationNum = jp%4;
                    ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(yOrientationTensor[rotationNum]);
	                atomGroupAction.actionPerformed(ghostMol);
	                
	                ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(xzOrientationTensor[rotationNum]);
	                atomGroupAction.actionPerformed(ghostMol);
                    
	                //Putting the molecule to its lattice site
                    Vector site = (Vector)lattice.site(siteIndex);
                    position.E(site);
                
                    atomActionTranslateTo.setDestination(position);
                    atomActionTranslateTo.actionPerformed(ghostMol);  
                    
                    if(coordinateDef.isBetaLatticeSum){
                    	translateBy.setTranslationVector(positionVector[rotationNum]);
                		atomGroupActionTranslate.actionPerformed(ghostMol);
                    }
                    
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
                    
                        IData value = function.f(molPair);
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
    
	private static final long serialVersionUID = 1L;
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
    protected final IMolecule[] moleculeCell0;
    protected Vector offset, position;
    protected final MoleculeActionTranslateTo atomActionTranslateTo;
    protected final AtomActionTranslateBy translateBy;
    protected MoleculeChildAtomAction atomGroupActionTranslate;
    protected CoordinateDefinitionNitrogen coordinateDef;
    protected Box ghostBox;
    protected Tensor[] xzOrientationTensor, yOrientationTensor;
    protected MoleculeChildAtomAction atomGroupAction;
    protected Vector[] positionVector;
    
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
