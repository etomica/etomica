/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

 package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.api.IRandom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomPair;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.Potential2;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Standard Monte Carlo atom-displacement trial move. Two atoms are moved at a
 * time in such a way that all the atoms' "images" are being moved at the same time too.
 * 
 * This class is written specifically for 32 system size FCC packing
 * 
 * 
 * @author Tai Boon Tan
 */
public class MCMoveAtomSuperBox extends MCMoveBoxStep {


    public MCMoveAtomSuperBox(PotentialMaster potentialMaster, IRandom random,
                              Space _space, CoordinateDefinitionLeafSuperBox coordinateDefinition) {
        super(potentialMaster);
        this.random = random;
        this.coordinateDefinition = coordinateDefinition;
 
        cells = coordinateDefinition.getBasisCells();
        
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        translationVector = _space.makeVector();
        setStepSizeMax(0.5);
        setStepSizeMin(0.0);
        setStepSize(0.1);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        affectedAtomList = new AtomArrayList(2);
        affectedAtomIterator = new AtomIteratorArrayListSimple(affectedAtomList);
        pair = new AtomPair();
        pairAB = new AtomPair();
        
        /*
         * determine the 'imaginary' box index
         * 	these boxes are not the actual simulation box!!!
         * 
         * boxIndex[a] = x
         * 	a is the cell # [total number of cells is 3x3x3x8 =216]
         * 	x is the box # [total number of boxes is 3x3x3 = 27]
         */
        if (coordinateDefinition.is864()){ //assignment for 864 atoms
        	//boxIndex = new int[216];
        	boxCells = new BasisCell[27][8];
        	int[] n = new int[27];
        	for (int j=0; j<6; j++){
	
        		for(int i=0; i<cells.length/6; i++){
    		
        			int thisBoxIndex =  (i/12)*3 + ((i%6)/2) + ((j/2)*9);
        			//boxIndex[i+(j*36)] = thisBoxIndex; //thisBoxIndex is the index for the 27 boxes
        			boxCells[thisBoxIndex][n[thisBoxIndex]] = cells[i+(j*36)];
        			n[thisBoxIndex]++;
        		}
    	
        	}
		
        } else { //assignment for 256 atoms
        	/*
        	 * boxCells[a][b]
        	 * 
        	 * a is the total number of equivalent cells
        	 * b is the total number of distinctive(or individual) cells
        	 */
        	boxCells = new BasisCell[8][8];
        	int[] n = new int[8];
        	
        	for (int j=0; j<4; j++){                    //loop over the layers of cells  
            	int plus;
        		for (int i=0; i<16; i++){   //loop over the cells in each layer of 16 cells
        			if (((i/8)%2)==0){
        				plus = 0;
        			} else {
        				plus = 1;
        			}
        			int thisBoxIndex = ((i%16)/8) + ((i%4)/2) + plus + ((j/2)*4);
        			
        			boxCells[thisBoxIndex][n[thisBoxIndex]] = cells[i+(j*16)];
        			n[thisBoxIndex]++;
        			
        		}
        		
        	}
        	        	
        }
        
    }
    
    public void setPotential(Potential2 newPotential) {
        potential = newPotential;
    }
    
    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
    	
    	randomNumber0 = random.nextInt(8);
    	randomAtom0 = random.nextInt(4);
    	
    	if(coordinateDefinition.is864()){
    		atom0 = boxCells[13][randomNumber0].molecules.getMolecule(randomAtom0).getChildList().getAtom(0);
    	} else {
    		atom0 = boxCells[7 - randomNumber0][randomNumber0].molecules.getMolecule(randomAtom0).getChildList().getAtom(0);
    	}
    	
    	randomNumber1 = random.nextInt(8);
    	randomAtom1 = random.nextInt(4);
    	
    	if(coordinateDefinition.is864()){
    		atom1 = boxCells[13][randomNumber1].molecules.getMolecule(randomAtom1).getChildList().getAtom(0);
    	} else{
    		atom1 = boxCells[7 - randomNumber1][randomNumber1].molecules.getMolecule(randomAtom1).getChildList().getAtom(0);
    	}
    	
        if (atom0 == null || atom1 == null || atom0 == atom1) return false;
        
        energyMeter.setTarget(atom0);
        uOld = energyMeter.getDataAsScalar();
        
        /*
         * pair atomSpeciesA atomSpeciesB
         */
        
        BasisCell cellA, cellB;
        double uCorrect = 0;
        	
    	for (int nInnerCells=0; nInnerCells < boxCells[0].length; nInnerCells++){
    		
    		if (coordinateDefinition.is864()){
    			cellA = boxCells[13][nInnerCells];
    		} else {
    			cellA = boxCells[7-nInnerCells][nInnerCells];
    		}
    		
    		for (int i=0; i<4; i++){
    			atomSpeciesA = cellA.molecules.getMolecule(i).getChildList().getAtom(0);
    			pairAB.atom0 = atomSpeciesA;
    			
    			if (atomSpeciesA == atom0 || atomSpeciesA == atom1) continue;
    			    			    			
				for (int nCellsBox=0; nCellsBox<boxCells.length; nCellsBox++){
					boolean isCenterCell=true;
					
					if(coordinateDefinition.is864()){
						if(nCellsBox == 13){
							isCenterCell = true;
						} else {
							isCenterCell = false;
						}
						
					} else {
						if(nCellsBox==(7-nInnerCells) && nInnerCells==(7-nCellsBox)){ //need to check
							isCenterCell = true;
						} else {
							isCenterCell = false;
						}
						
					}
					
					if (!isCenterCell){
						cellB = boxCells[nCellsBox][randomNumber0];
						atomSpeciesB = cellB.molecules.getMolecule(randomAtom0).getChildList().getAtom(0);
						pairAB.atom1 = atomSpeciesB;
						uCorrect += potential.energy(pairAB)/2;
					
						cellB = boxCells[nCellsBox][randomNumber1];
						atomSpeciesB = cellB.molecules.getMolecule(randomAtom1).getChildList().getAtom(0);
						pairAB.atom1 = atomSpeciesB;
						uCorrect += potential.energy(pairAB)/2;
					
					}
				
				}
    			
    		}
    	}
    	       
        energyMeter.setTarget(atom1);
        uOld += energyMeter.getDataAsScalar();
        uOld += uCorrect;
        pair.atom0 = atom0;
        pair.atom1 = atom1;
        uOld -= potential.energy(pair); //to take care of the double counting
        
        if(uOld > 1e10) {
            throw new ConfigurationOverlapException(box);
        }
        
        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        
        /*
         * Moving all the corresponding atoms in other boxes
         * 	at the same time with the atom that was randomly 
         * 	generated 
         * 
         * there are total of 27 boxes
         */
   
        for (int i=0; i<boxCells.length;i++){
        	boxCells[i][randomNumber0].molecules.getMolecule(randomAtom0).getChildList().getAtom(0).getPosition().PE(translationVector);
        	boxCells[i][randomNumber1].molecules.getMolecule(randomAtom1).getChildList().getAtom(0).getPosition().ME(translationVector);
        	
        }
        
        uCorrect = 0;
        for (int nInnerCells=0; nInnerCells < boxCells[0].length; nInnerCells++){
    		
    		if (coordinateDefinition.is864()){
    			cellA = boxCells[13][nInnerCells];
    		} else {
    			cellA = boxCells[7-nInnerCells][nInnerCells];
    		}
    		
    		for (int i=0; i<4; i++){
    			atomSpeciesA = cellA.molecules.getMolecule(i).getChildList().getAtom(0);
    			pairAB.atom0 = atomSpeciesA;
    			
    			if (atomSpeciesA == atom0 || atomSpeciesA == atom1) continue;
    			    			    			
				for (int nCellsBox=0; nCellsBox<boxCells.length; nCellsBox++){
					boolean isCenterCell=true;
					
					if(coordinateDefinition.is864()){
						if(nCellsBox == 13){
							isCenterCell = true;
						} else {
							isCenterCell = false;
						}
						
					} else {
						if(nCellsBox==(7-nInnerCells) && nInnerCells==(7-nCellsBox)){ //need to check
							isCenterCell = true;
						} else {
							isCenterCell = false;
						}
						
					}
					
					if (!isCenterCell){
						cellB = boxCells[nCellsBox][randomNumber0];
						atomSpeciesB = cellB.molecules.getMolecule(randomAtom0).getChildList().getAtom(0);
						pairAB.atom1 = atomSpeciesB;
						uCorrect += potential.energy(pairAB)/2;
					
						cellB = boxCells[nCellsBox][randomNumber1];
						atomSpeciesB = cellB.molecules.getMolecule(randomAtom1).getChildList().getAtom(0);
						pairAB.atom1 = atomSpeciesB;
						uCorrect += potential.energy(pairAB)/2;
					
					}
				
				}
    			
    		}
    	}
     
        energyMeter.setTarget(atom0);
        uNew = energyMeter.getDataAsScalar();
        uNew += uCorrect;
        energyMeter.setTarget(atom1);
        uNew += energyMeter.getDataAsScalar(); 
        uNew -= potential.energy(pair);
        
        return true; 
    }//end of doTrial
    
    
    /**
     * Returns log of the ratio of the trial probabilities, ln(Tij/Tji) for the
     * states encountered before (i) and after (j) the most recent call to doTrial(). 
     * Tij is the probability that this move would generate state j from state i, and
     * Tji is the probability that a subsequent call to doTrial would return to state i
     * from state j.
     */
    public double getA() {return 1.0;}
    
    /**
     * Returns the log of the limiting-distribution probabilities of states, ln(Pj/Pi), 
     * for the states encountered before (i) and after (j) the most recent call to 
     * doTrial.
     */
    public double getB() {
        return -(uNew - uOld);
    }
    
    public double energyChange() {return uNew - uOld;}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
    }
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.  This method should cause the system to be restored to the condition
     * before the most recent call to doTrial.
     */
    public void rejectNotify() {
    	
    	//System.err.println("Rejected atom is in cell["+randomNumber0+"]["+randomAtom0+"]");
        for (int i=0; i<boxCells.length;i++){
        	boxCells[i][randomNumber0].molecules.getMolecule(randomAtom0).getChildList().getAtom(0).getPosition().ME(translationVector);
        	boxCells[i][randomNumber1].molecules.getMolecule(randomAtom1).getChildList().getAtom(0).getPosition().PE(translationVector);
 
        }
    }
        
    
    public AtomIterator affectedAtoms() {
        affectedAtomList.clear();
        affectedAtomList.add(atom0);
        affectedAtomList.add(atom1);
        return affectedAtomIterator;
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        potential.setBox(p);
    }
    
    
    private static final long serialVersionUID = 2L;
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final AtomArrayList affectedAtomList;
    protected final MeterPotentialEnergy energyMeter;
    protected final Vector translationVector;
    protected IAtom atom0, atom1, atomSpeciesA, atomSpeciesB;
    protected double uOld;
    protected double uNew;
    protected final IRandom random;
    protected Potential2 potential;
    protected final AtomPair pair, pairAB;
    protected CoordinateDefinitionLeafSuperBox coordinateDefinition;
    protected BasisCell[] cells;
    //protected int[] boxIndex;
    protected BasisCell[][] boxCells; 
    protected int randomNumber0, randomNumber1, randomAtom0, randomAtom1;
}
