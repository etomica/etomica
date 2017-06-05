/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.api.*;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.atom.MoleculeArrayList;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.util.random.IRandom;

/**
 * Monte Carlo molecule-displacement trial move for a molecule in smer
 * 
 * @author Hye Min Kim
 */
public class MCMoveMoleculeSmer extends MCMoveMolecule {
	protected AssociationManagerMolecule associationManager;
	protected final MoleculeArrayList bondList, smerList;
	public static boolean dodebug = true;;
	protected final Vector dr;
	protected int maxLength = Integer.MAX_VALUE;
	protected IAssociationHelperMolecule associationHelper;
	

	public MCMoveMoleculeSmer(Simulation sim, PotentialMaster potentialMaster,
                              Space _space) {
		this(potentialMaster, sim.getRandom(), _space, 1.0, 15.0);
	}


	public MCMoveMoleculeSmer(PotentialMaster potentialMaster, IRandom random,
                              Space _space, double stepSize, double stepSizeMax) {
		super(potentialMaster, random, _space, stepSize, stepSizeMax);
		this.smerList = new MoleculeArrayList();
		this.dr = _space.makeVector();
		bondList = new MoleculeArrayList();
	}
	public void setAssociationManager(AssociationManagerMolecule associationManager, IAssociationHelperMolecule associationHelper){
		this.associationManager = associationManager;
		this.associationHelper = associationHelper;
		MoleculeSourceRandomDimer moleculeSourceRandomDimer = new MoleculeSourceRandomDimer();
		moleculeSourceRandomDimer.setAssociationManager(associationManager);
		if (box != null) {
			moleculeSourceRandomDimer.setBox(box);
		}
		moleculeSourceRandomDimer.setRandomNumberGenerator(random);
		setMoleculeSource(moleculeSourceRandomDimer);
	}
	public boolean doTrial() {
        molecule = moleculeSource.getMolecule();
        if (molecule == null) return false;
        bondList.clear();
        bondList.addAll(associationManager.getAssociatedMolecules(molecule));//making a copy of the list of bonded atoms
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        if(uOld > 1e8) {
        	//PotentialCalculationEnergySum.dodebug = true;
        	System.out.println("molecule "+molecule+" bondList "+bondList.getMolecule(0));
        	energyMeter.getDataAsScalar();
            throw new RuntimeException("molecule "+molecule+" in box "+box+" has an overlap");
        }
        groupTranslationVector.setRandomCube(random);
        groupTranslationVector.TE(stepSize);
        moveMoleculeAction.actionPerformed(molecule);
        return true;
    }
	
	public double getA(){
		if (associationHelper.populateList(smerList,molecule,true)){
    		return 0;
    	}
		IMoleculeList newBondList = associationManager.getAssociatedMolecules(molecule);
		if (bondList.getMoleculeCount() != newBondList.getMoleculeCount()){
			return 0;
		}
		for (int i = 0; i < bondList.getMoleculeCount(); i+=1){
			IMolecule b = bondList.getMolecule(i);
			boolean success = false;
			for (int j = 0; j<newBondList.getMoleculeCount(); j+=1){
				if (b == newBondList.getMolecule(j)){
					success = true;//to check all the atoms in the bondList are still associating
					break;
				}
			}
			if (!success){
				return 0;
			}
		}
		return 1.0;
	}
}
