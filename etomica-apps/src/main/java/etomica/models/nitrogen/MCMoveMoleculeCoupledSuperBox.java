/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.api.*;
import etomica.box.Box;
import etomica.atom.MoleculeArrayList;
import etomica.atom.MoleculePair;
import etomica.atom.MoleculeSource;
import etomica.atom.MoleculeSourceRandomMolecule;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.MoleculeIterator;
import etomica.atom.iterator.MoleculeIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.integrator.mcmove.MCMoveMolecular;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Standard Monte Carlo molecule-displacement trial move for superbox. Two molecules are 
 * moved at a time in such a way that the geometric center of the system is not changed.
 *
 * When move one molecule in the center cell; the same molecule in the other 26 unit cells 
 * will move too!
 *
 * The potential master that being passed in this class is a "full-strength" potential without
 *  multiplying by 0.5 for the speciesA and speciesB interaction.
 * 
 * @author Tai Boon Tan
 */
public class MCMoveMoleculeCoupledSuperBox extends MCMoveBoxStep implements MCMoveMolecular{

    
    public MCMoveMoleculeCoupledSuperBox(PotentialMaster potentialMaster, IRandom nRandom,
                                         Space _space, Box box, int nC, int basis, CoordinateDefinitionNitrogenSuperBox coordinateDef){
        super(potentialMaster);
        this.random = nRandom;
        this.box = box;
        this.basisCell = coordinateDef.getBasisCells();
        
        super.setBox(box);
        
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule)moleculeSource).setRandomNumberGenerator(random);
        moleculeSource.setBox(box);
        
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        energyMeter.setBox(box);
        
        affectedMoleculeIterator = new MoleculeIteratorArrayListSimple();
        affectedMoleculeList = new MoleculeArrayList();
        
        singleAction = new AtomActionTranslateBy(_space);
        groupTransVect = singleAction.getTranslationVector();
        
        moveMoleculeAction = new MoleculeChildAtomAction(singleAction);
        
        pair = new MoleculePair();
        pairAB = new MoleculePair();
        
        perParticleFrequency = true;
        
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

    public void setBox(Box newBox) {
 
    }
    
    public void setPotential(IPotentialMolecular newPotential){
        potential = newPotential;
    }
    
    public AtomIterator affectedAtoms() {
    	/*
    	 * I'm NOT USING THIS METHOD!!!
    	 */
    	System.err.println("<MCMoveMoleculeCoupledSuperBox> I'm not using AtomIterator but MoleculIterator!!");
        return null;
    }

	public MoleculeIterator affectedMolecules(Box aBox) {
		   if (box == aBox) {
			   affectedMoleculeIterator.setList(affectedMoleculeList);
			   return affectedMoleculeIterator;
	        }
		    
		   return null;
	}
    
    public double energyChange() {return uNew - uOld;}

    public void acceptNotify() {
        // I do believe nothing needs to happen here.
    }

    public boolean doTrial() {
    	
        randomMol0 = random.nextInt(nA);
        randomMol1 = random.nextInt(nA);

        molecule0 = basisCell[0].molecules.getMolecule(molIndex[13][randomMol0]);
        molecule1 = basisCell[0].molecules.getMolecule(molIndex[13][randomMol1]);

        affectedMoleculeList.clear();
        affectedMoleculeList.add(molecule0);
        affectedMoleculeList.add(molecule1);
        
        if(molecule0==null || molecule1==null || molecule0==molecule1) return false;
        
        energyMeter.setTarget(molecule0);
        uOld = energyMeter.getDataAsScalar();
        
        //molSpeciesA
        /*
         * uCorrect is account for the double-counting of the molecular energy
         */
        uCorrect = 0.0;
        pairAB.atom0 = molecule0;
        for(int i=0; i<molIndex.length; i++){
        	if(i==13) continue;
        	pairAB.atom1 = basisCell[0].molecules.getMolecule(molIndex[i][randomMol0]);
        	uCorrect += potential.energy(pairAB);
        	
        	pairAB.atom1 = basisCell[0].molecules.getMolecule(molIndex[i][randomMol1]);
        	uCorrect += potential.energy(pairAB);
        }
        pairAB.atom0 = molecule1;
        for(int i=0; i<molIndex.length; i++){
        	if(i==13) continue;
        	pairAB.atom1 = basisCell[0].molecules.getMolecule(molIndex[i][randomMol0]);
        	uCorrect += potential.energy(pairAB);
        	
        	pairAB.atom1 = basisCell[0].molecules.getMolecule(molIndex[i][randomMol1]);
        	uCorrect += potential.energy(pairAB);
        }
        
        energyMeter.setTarget(molecule1);
        uOld += energyMeter.getDataAsScalar();
        uOld -= 0.5*uCorrect;
        pair.atom0 = molecule0;
        pair.atom1 = molecule1;
        uOld -= potential.energy(pair);
        
        
        if(uOld > 1e10){
            throw new ConfigurationOverlapException(box);
        }
        
        groupTransVect.setRandomCube(random);
        groupTransVect.TE(stepSize);
        for(int i=0; i<molIndex.length; i++){
        	moveMoleculeAction.actionPerformed(basisCell[0].molecules.getMolecule(molIndex[i][randomMol0]));
        }
        groupTransVect.TE(-1.0);
        for(int i=0; i<molIndex.length; i++){
        	moveMoleculeAction.actionPerformed(basisCell[0].molecules.getMolecule(molIndex[i][randomMol1]));
        }
  
        molecule0 = basisCell[0].molecules.getMolecule(molIndex[13][randomMol0]);
        molecule1 = basisCell[0].molecules.getMolecule(molIndex[13][randomMol1]);

        uCorrect = 0.0;
        pairAB.atom0 = molecule0;
        for(int i=0; i<molIndex.length; i++){
        	if(i==13) continue;
        	pairAB.atom1 = basisCell[0].molecules.getMolecule(molIndex[i][randomMol0]);
        	uCorrect += potential.energy(pairAB);
        	
        	pairAB.atom1 = basisCell[0].molecules.getMolecule(molIndex[i][randomMol1]);
        	uCorrect += potential.energy(pairAB);
        }
        pairAB.atom0 = molecule1;
        for(int i=0; i<molIndex.length; i++){
        	if(i==13) continue;
        	pairAB.atom1 = basisCell[0].molecules.getMolecule(molIndex[i][randomMol0]);
        	uCorrect += potential.energy(pairAB);
        	
        	pairAB.atom1 = basisCell[0].molecules.getMolecule(molIndex[i][randomMol1]);
        	uCorrect += potential.energy(pairAB);
        }
        
        energyMeter.setTarget(molecule0);
        uNew = energyMeter.getDataAsScalar();
        uNew -= 0.5*uCorrect;
        energyMeter.setTarget(molecule1);
        uNew += energyMeter.getDataAsScalar();
        uNew -= potential.energy(pair);
        /*
         * Because we have uNew is infinity, and we don't want to have to 
         * worry about the system subtracting infinity from infinity, and
         * setting uNew equal to zero, and accepting the move.
         */
       // if(Double.isInfinite(uNew)) {return true;}  
        
        return true;
    }

    public double getA() {
        return 1.0;
    }

    public double getB() {
        return -(uNew - uOld);
    }

    public void rejectNotify() {
        for(int i=0; i<molIndex.length; i++){
        	moveMoleculeAction.actionPerformed(basisCell[0].molecules.getMolecule(molIndex[i][randomMol0]));
        }
        groupTransVect.TE(-1.0);

        for(int i=0; i<molIndex.length; i++){
        	moveMoleculeAction.actionPerformed(basisCell[0].molecules.getMolecule(molIndex[i][randomMol1]));
        }
    }
    
    private static final long serialVersionUID = 1L;
    protected final MoleculeChildAtomAction moveMoleculeAction;
    protected final Vector groupTransVect;
    protected IMolecule molecule0, molecule1;
    protected final MeterPotentialEnergy energyMeter;
    protected MoleculeSource moleculeSource;
    protected double uOld, uNew;
    protected final IRandom random;
    protected MoleculeIteratorArrayListSimple affectedMoleculeIterator;
    protected MoleculeArrayList affectedMoleculeList;
    protected final AtomActionTranslateBy singleAction;
    protected final MoleculePair pair, pairAB;
    protected IPotentialMolecular potential;
    protected int[][] molIndex;
    protected int randomMol0, randomMol1;
    protected IMolecule molSpeciesA, molSpeciesB;
    protected int nA;
    protected BasisCell[] basisCell;
    protected double uCorrect;
    
}
