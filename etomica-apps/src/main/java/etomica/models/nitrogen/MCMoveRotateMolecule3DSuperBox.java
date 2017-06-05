/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;
import etomica.atom.*;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.potential.PotentialMaster;
import etomica.api.IPotentialMolecular;
import etomica.api.IRandom;
import etomica.space.Vector;
import etomica.atom.IMoleculePositionDefinition;
import etomica.atom.iterator.MoleculeIterator;
import etomica.atom.iterator.MoleculeIteratorArrayListSimple;
import etomica.integrator.mcmove.MCMoveMolecular;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.space.Space;
import etomica.space.RotationTensor;


/**
 * MC Rotate Move 3D for superbox
 * 
 * The potential master that being passed in this class is a "full-strength" potential without
 *  multiplying by 0.5 for the speciesA and speciesB interaction.
 *  
 * @author Tai Boon Tan
 *
 */
public class MCMoveRotateMolecule3DSuperBox extends MCMoveMolecule implements MCMoveMolecular{
    
    private static final long serialVersionUID = 2L;
    protected transient Vector r0;
    protected transient RotationTensor rotationTensor;
    protected IMoleculePositionDefinition positionDefinition;
    public int count;
    public int count1;
    public boolean flag = false;
    public boolean flag1 = false;
    protected int nA;
    protected int[][] molIndex;
    protected int molNum;
    protected IRandom random;
    protected BasisCell[] basisCell;
    protected IMolecule affectedMol;
    protected CoordinateDefinitionNitrogenSuperBox coordinateDef;
    protected MoleculeIteratorArrayListSimple affectedMoleculeIterator;
    protected MoleculeArrayList affectedMoleculeList;
    protected MoleculePair pairAB;
    protected IPotentialMolecular potentialAA;
    protected double uOldPE, uNewPE;
    protected double uCorrect;
    
    public MCMoveRotateMolecule3DSuperBox(PotentialMaster potentialMaster, IRandom random,
                                          Space _space, int nC, int basis, CoordinateDefinitionNitrogenSuperBox coordinateDef) {
        super(potentialMaster, random, _space, Math.PI/2, Math.PI);
        this.basisCell = coordinateDef.getBasisCells();
        this.random = random;
        this.coordinateDef = coordinateDef;
        box = coordinateDef.getBox();
        setBox(coordinateDef.getBox());
        
        affectedMoleculeIterator = new MoleculeIteratorArrayListSimple();
        affectedMoleculeList = new MoleculeArrayList();
        positionDefinition = new MoleculePositionGeometricCenter(space);
        
        rotationTensor = _space.makeRotationTensor();
        r0 = _space.makeVector();
        pairAB = new MoleculePair();
   
        nA = (nC*nC*nC)*basis;
        molIndex = new int[27][nA];
        
        int layerNum = (nC*basis*3)*(3*nC);
        int oneThirdLayerNum = (nC*basis*3)*(nC);
        
        int axisNum = nC*basis*3;
        int cellAxisNum = nC*basis;
        int ix= 0;
        int iy= 0;
        int iz= 0;
        
        for (int iCell=0; iCell<molIndex.length; iCell++){
        	
        	int counter = 0;
        	int cellConst = 0;
        	
        	if(iCell%27 > 17){
        		ix = 2;
        	} else if (iCell%27 > 8 && iCell%27 <=17){
        		ix = 1;
        	} else {
        		ix = 0;
        	}
        	
        	if(iCell%9 > 5){
        		iy = 2;
        	} else if (iCell%9 > 2 && iCell%9 <=5){
        		iy = 1;
        	} else {
        		iy = 0;
        	}
        	
        	iz = iCell%3;
        
        	cellConst = iz*(cellAxisNum) + iy*oneThirdLayerNum + ix*nC*layerNum;
        	
	        for(int xnC=0; xnC<nC; xnC++){
		        for(int ynC=0; ynC<nC; ynC++){
			        for(int i=0; i<nC*basis; i++){
			        	molIndex[iCell][counter] = (i+(ynC*axisNum)+(xnC*layerNum)+cellConst);
			        	++counter;
			        }
		        }
	        }
        }
       
    }
     
    public boolean doTrial() {
    	
        if(box.getMoleculeList().getMoleculeCount()==0) {molecule = null; return false;}
            
        molNum = random.nextInt(nA);
        molecule = basisCell[0].molecules.getMolecule(molIndex[13][molNum]);
      
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        
        uCorrect = 0.0;
        pairAB.atom0 = molecule;
        for(int i=0; i<molIndex.length; i++){
        	if(i==13) continue;
        	pairAB.atom1 = basisCell[0].molecules.getMolecule(molIndex[i][molNum]);
        	uCorrect += potentialAA.energy(pairAB);
     
        }
        
        uOld -= 0.5*uCorrect;
        
        if(Double.isInfinite(uOld)) {
            throw new RuntimeException("Overlap in initial state");
        }
        
        double dTheta = (2*random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);
        
        affectedMoleculeList.clear();
        
        for(int i=0; i<molIndex.length; i++){
        	molecule = basisCell[0].molecules.getMolecule(molIndex[i][molNum]);
        	affectedMoleculeList.add(affectedMol);
        	r0.E(positionDefinition.position(molecule));
        	doTransform(molecule, r0);
        }
        
        /*
         * After the rotational move
         */
      
        molecule = basisCell[0].molecules.getMolecule(molIndex[13][molNum]);
        energyMeter.setTarget(molecule);
        uNew = energyMeter.getDataAsScalar();
        
        /*
         * uCorrect is account for the double-counting of the molecular energy
         */
        uCorrect = 0.0;
        pairAB.atom0 = molecule;
        for(int i=0; i<molIndex.length; i++){
        	if(i==13) continue;
        	pairAB.atom1 = basisCell[0].molecules.getMolecule(molIndex[i][molNum]);
        	uCorrect += potentialAA.energy(pairAB);

        }
        uNew -= 0.5*uCorrect;
        
        return true;
    }

    public double getB(){
    	return -(uNew - uOld);
    }
    
	public MoleculeIterator affectedMolecules(Box aBox) {
	   if (box == aBox) {
		   affectedMoleculeIterator.setList(affectedMoleculeList);
		   return affectedMoleculeIterator;
        }
	    
	   return null;
	   
	}
    
    protected void doTransform(IMolecule molecule, Vector r0) {
        IAtomList childList = molecule.getChildList();
        for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {
            IAtom a = childList.getAtom(iChild);
            Vector r = a.getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
        }
    }
    
    public void rejectNotify() {
    	rotationTensor.invert();
    	for(int i=0; i<molIndex.length; i++){
        	molecule = basisCell[0].molecules.getMolecule(molIndex[i][molNum]);
        	r0.E(positionDefinition.position(molecule));
	        doTransform(molecule, r0);
    	}
    }

    public void setPotential(IPotentialMolecular newPotential){
        potentialAA = newPotential;
    }
    
}
