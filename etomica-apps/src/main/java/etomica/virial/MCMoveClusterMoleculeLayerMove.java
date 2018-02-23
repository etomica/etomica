/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * Monte Carlo molecule-displacement from 1 layer to the other layer trial move for cluster integrals.
 * 
 * copied and modified from MCMoveClusterMolecule class
 * for virial coefficients calculation of phenanthrene and anthracene 
 * for n = 4 , the molecules tend to form either 2 layers or 3 layers structures 
 * add this monte carlo move, force one molecule to move out the existing layer to the other layer to explore more configurations
 * @author shu
 * date : april 27 2011
 */
public class MCMoveClusterMoleculeLayerMove extends MCMoveClusterMolecule {
    private static final long serialVersionUID = 1L;
    protected final Vector vector1 , vector2 , crossVector; //crossVector = vector1 * vector2
    public MCMoveClusterMoleculeLayerMove(Simulation sim, Space _space) {
    	super (sim.getRandom(), _space, 1.0); //superclass parameter stepsize
    	
    	vector1     =  space.makeVector(); // Initialize 
    	vector2     =  space.makeVector();// Initialize 
    	crossVector =  space.makeVector(); // Initialize 
    }
      
    public boolean doTrial() {
    	if(box.getMoleculeList().getMoleculeCount()==1) return false;
        
        molecule = moleculeSource.getMolecule();
        while (molecule.getIndex() == 0) {
            molecule = moleculeSource.getMolecule();// keep on picking up molecules randomly if picking up the molecule in the origin
        }
        
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box); // calculate integrand (-1/2) * f 12 , i.e. B2 
       
         // generate a random integer, 0 or 1 call nextInt in IRandom interfce
        int direction = random.nextInt(2);/////tai :cannot be the same 
       // now the  _random value is either 0 or 1
       // if is 0, then the chosen molecule is moved toward the adjacent plane 
        IAtomList atoms = molecule.getChildList(); // get atoms
     // public IAtom getAtom(int i);
        IAtom atom0 = atoms.get(0);
        IAtom atom1 = atoms.get(1);
        IAtom atom2 = atoms.get(2);
      //public IVectorMutable getPosition();
        Vector position0 = atom0.getPosition();
        Vector position1 = atom1.getPosition();
        Vector position2 = atom2.getPosition();
        vector1.Ev1Mv2(position1, position0);
        vector2.Ev1Mv2(position2, position1);
        crossVector.E(vector1);
        crossVector.XE(vector2);
        crossVector.normalize();
        crossVector.TE(3.85 * (1 + random.nextInt(2)));
        groupTranslationVector.E(crossVector);
        if (direction == 0 ){
        	groupTranslationVector.TE(-1);
        }
        moveMoleculeAction.actionPerformed(molecule);
        //groupTranslationVector.setRandomCube(random); //  random unit vector
      //  groupTranslationVector.TE(stepSize);  // multiply by a number, this will be tuned to get a 50 % trail move acceptance
      //  moveMoleculeAction.actionPerformed(molecule); // execute???????????
       
        
        ((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }
    
    
        
}
