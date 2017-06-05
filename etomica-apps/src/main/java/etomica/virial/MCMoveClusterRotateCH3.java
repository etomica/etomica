/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.*;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Vector3D;
import etomica.util.random.IRandom;

/**
 * Monte Carlo CH3 rotation for cluster integrals(Alkane TraPPE-EH).
 * 1.pick up a CH3 group randomly(2 choices)
 * 2.rotate HHH plane (perpendicular to C0-C1 or C[n-2]-C[n-1]) to a random angle from -60 to 60 degree
 * 
 * @author shu
 * Feb 2013
 */
public class MCMoveClusterRotateCH3 extends MCMoveMolecule {

    public MCMoveClusterRotateCH3(Simulation sim, PotentialMaster potentialMaster, int nAtoms, Space _space) {
    	this(potentialMaster,sim.getRandom(), 1.0, nAtoms, _space);
    }
    
    public MCMoveClusterRotateCH3(PotentialMaster potentialMaster,
                                  IRandom random, double stepSize, int nAtoms, Space _space) {
        super(potentialMaster,random,_space, stepSize,Double.POSITIVE_INFINITY);
        this.space = _space;
        setStepSizeMax(Math.PI/3);
        energyMeter = new MeterPotentialEnergy(potential);
        axis = space.makeVector();
        rotateTensor = (RotationTensor3D)(space.makeRotationTensor());

    }

    public void setBox(Box p) {
    	super.setBox(p);
        selectedAtoms = new IAtom[box.getMoleculeList().getMoleculeCount()];
        translationVectors = new Vector3D[box.getMoleculeList().getMoleculeCount()];
        for (int i=0; i<translationVectors.length; i++) {
            translationVectors[i] = space.makeVector();
        }
        energyMeter.setBox(p);
    }
    
    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }

    //note that total energy is calculated
    public boolean doTrial() {
        uOld = energyMeter.getDataAsScalar();
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);

        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i=0; i<moleculeList.getMoleculeCount(); i++) {
            if (species != null && moleculeList.getMolecule(i).getType() != species) {
                continue;
            }
            IAtomList childList = moleculeList.getMolecule(i).getChildList();
            Vector position = space.makeVector();
            Vector positionNeighbor = space.makeVector();
            int numChildren = childList.getAtomCount();// total atoms in the i-th molecule
            int numCarbons = (numChildren-2)/3;// number of carbons in the i-th molecule
            int j = random.nextInt(2);// 0 or 1
            double theta = (random.nextDouble()-0.5)* 2.0 * stepSize;// [-60,60)

            selectedAtoms[i] = childList.getAtom(j);
//            System.out.println(selectedAtoms[i]+" "+j+" before "+selectedAtoms[i].coord.position());
            if (j==0){
            	position = childList.getAtom(0).getPosition();
            	positionNeighbor = childList.getAtom(1).getPosition();
            }
            else if (j==1){
            	position = childList.getAtom(numCarbons-1).getPosition();
            	positionNeighbor = childList.getAtom(numCarbons-2).getPosition();
            }
            else {
            	throw new RuntimeException("wrong random number in CH3 MC move!");
            }
            translationVectors[i].Ea1Tv1(-1,position);
  
            ///// ######################################################################### ///////////////////
            ///// #########  GET the new position of j molecule ########################### ///////////////////
            ///// ######################################################################### ///////////////////
            axis.Ev1Mv2(position, positionNeighbor);
            axis.normalize();
            rotateTensor.setRotationAxis(axis, theta);

            IAtom[] hydrogen = new IAtom[3];
            if (j==0){
            	///// ######################################################################### ///////////////////
            	///// #########  GET the new positions of 3 H on C0  ########################### ///////////////////
            	///// ######################################################################### ///////////////////
            	for (int s=0;s<3; s++){
            		hydrogen[s] =  childList.getAtom(numCarbons * (s+1));// [n], [2n], [3n]
            		Vector r = hydrogen[s].getPosition();
            		r.ME(position);//position is position of C0
            		rotateTensor.transform(r);
            		r.PE(position);
            	}
            }
            else if (j==1){
            	///// ######################################################################### ///////////////////
            	///// #########  GET the new positions of 3 H on C[n-1]      ###################### ///////////////////
            	///// ######################################################################### ///////////////////
            	hydrogen[0]=childList.getAtom(numCarbons*2-1);//[n-1+n]
            	hydrogen[1]=childList.getAtom(numCarbons*3-1);//[n-1+2n]
            	hydrogen[2]=childList.getAtom(numCarbons*3+1);//[3n+1]
            	for (int s=0;s<3; s++){
            		Vector r = hydrogen[s].getPosition();
            		r.ME(position);//position is position of C[n-1]
            		rotateTensor.transform(r);
            		r.PE(position);
            	}
            }
            
            // shift the whole molecule so that the center of mass doesn't change
            translationVectors[i].PE(position);
            axis.E(translationVectors[i]);
            axis.TE(1.0/childList.getAtomCount());
            for (int k=0; k<childList.getAtomCount(); k++) {
                childList.getAtom(k).getPosition().ME(axis);
            }
        
        }
        ((BoxCluster)box).trialNotify();
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        uNew = energyMeter.getDataAsScalar();
        return true;
    }
    
    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i=0; i<selectedAtoms.length; i++) {
            if (species != null && moleculeList.getMolecule(i).getType() != species) continue;
            IAtomList childList = moleculeList.getMolecule(i).getChildList();
            axis.E(translationVectors[i]);
            axis.TE(1.0/childList.getAtomCount());
            for (int k=0; k<childList.getAtomCount(); k++) {
                childList.getAtom(k).getPosition().PE(axis);
            }
            selectedAtoms[i].getPosition().ME(translationVectors[i]);
        }
        ((BoxCluster)box).rejectNotify();
    }

    public double getB() {
        return -(uNew - uOld);
    }
    
    public double getA() {
        return wNew/wOld;
    }
	
    private static final long serialVersionUID = 1L;
    protected final MeterPotentialEnergy energyMeter;
    protected IAtom[] selectedAtoms;
    protected final Vector axis;
    protected Vector[] translationVectors;
    protected double wOld, wNew;
    protected final Space space;
    protected ISpecies species;
    protected RotationTensor3D rotateTensor;

}
