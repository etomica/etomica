/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.molecule.MoleculeArrayList;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Monte Carlo molecule-displacement trial move for monomer molecules
 *
 * @author Hye Min Kim
 */
public class MCMoveMoleculeMonomer extends MCMoveMolecule {
	protected AssociationManagerMolecule associationManager;
	protected final MoleculeArrayList smerList;
	protected final Vector dr;
	protected int maxLength = Integer.MAX_VALUE;
	protected IAssociationHelperMolecule associationHelper;


	public MCMoveMoleculeMonomer(IRandom random, PotentialMaster potentialMaster,
								 Space _space) {
		this(potentialMaster, random, _space, 1.0, 15.0);
	}


	public MCMoveMoleculeMonomer(PotentialMaster potentialMaster, IRandom random,
								 Space _space, double stepSize, double stepSizeMax) {
		super(potentialMaster, random, _space, stepSize, stepSizeMax);
		this.smerList = new MoleculeArrayList();
		this.dr = _space.makeVector();
	}
	public void setAssociationManager(AssociationManagerMolecule associationManager, IAssociationHelperMolecule associationHelper){
		this.associationManager = associationManager;
		this.associationHelper = associationHelper;
		MoleculeSourceRandomMonomer moleculeSourceRandomMonomer = new MoleculeSourceRandomMonomer();
		moleculeSourceRandomMonomer.setAssociationManager(associationManager);
		if (box != null) {
			moleculeSourceRandomMonomer.setBox(box);
		}
		moleculeSourceRandomMonomer.setRandomNumberGenerator(random);
		setMoleculeSource(moleculeSourceRandomMonomer);
	}

    public double getChi(double temperature) {
        if (associationHelper.populateList(smerList,molecule,true)){
    		return 0;
    	}
		if (associationManager.getAssociatedMolecules(molecule).size() > 0) {
        	return 0.0;
        } 
		if (smerList.size() > maxLength) {
    		return 0.0;
		}
        return super.getChi(temperature);
    }

	
}
