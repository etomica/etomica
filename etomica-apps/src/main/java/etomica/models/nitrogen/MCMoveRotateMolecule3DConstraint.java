/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.MoleculeChildAtomAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.paracetamol.AtomActionTransformed;
import etomica.potential.PotentialMaster;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Tensor3D;
import etomica.units.Degree;
import etomica.util.random.IRandom;


/**
 * MCMoveRotate that moves the molecules according to pre-set constraint angle
 *  which is given by the parameter "angle"
 *  
 * This class does not consider rotational energy, so the variable "energyChange"
 *  always returns zero
 *  
 * getB() always returns positive one is to ensure the proposed move is always 
 *  accepted
 * 
 * @author Tai Boon Tan
 *
 */
public class MCMoveRotateMolecule3DConstraint extends MCMoveMolecule {
    
    private static final long serialVersionUID = 2L;
    protected transient Vector r0;
    protected transient RotationTensor rotationTensor;
    protected IMoleculePositionDefinition positionDefinition;
    protected double constraintAngle;
    protected CoordinateDefinitionNitrogen coordinateDef;
    protected Vector[][] initMolecOrientation;
    protected Vector molecOrientation, rotationAxis, workVector;
	protected RotationTensor3D rotation;
	protected Tensor3D tensor;
    protected final MoleculeChildAtomAction atomGroupAction;
    
    public MCMoveRotateMolecule3DConstraint(PotentialMaster potentialMaster, IRandom random,
                                            Space _space, double angle, CoordinateDefinitionNitrogen coordinateDef, Box box) {
        super(potentialMaster, random, _space, 0.5*Math.PI, Math.PI);
        rotationTensor = _space.makeRotationTensor();
        r0 = _space.makeVector();
        positionDefinition = new MoleculePositionGeometricCenter(space);
        constraintAngle = angle;
        this.coordinateDef = coordinateDef;
        
        int numMolec = box.getMoleculeList().getMoleculeCount();
     	initMolecOrientation = new Vector[numMolec][3];
     	molecOrientation = space.makeVector();
     	rotationAxis = space.makeVector();
     	workVector = space.makeVector();
     	
    	tensor = new Tensor3D(new double[][]{{1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{0.0, 0.0, 1.0}});
		rotation = new RotationTensor3D();
		rotation.E(tensor);
    	/*
		 * initializing the initial orientation of the molecule
		 */
		for (int i=0; i<numMolec; i++){
			initMolecOrientation[i] = space.makeVectorArray(3);
			initMolecOrientation[i] = coordinateDef.getMoleculeOrientation(box.getMoleculeList().getMolecule(i));
		}
		atomGroupAction = new MoleculeChildAtomAction(new AtomActionTransformed(space));
	        
    }
     
    public boolean doTrial() {

        if(box.getMoleculeList().getMoleculeCount()==0) {molecule = null; return false;}
        int iMolecule = random.nextInt(box.getMoleculeList().getMoleculeCount());
        
        molecule = coordinateDef.getBox().getMoleculeList().getMolecule(iMolecule);
        r0.E(positionDefinition.position(molecule));
        
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        if(Double.isInfinite(uOld)) {
            throw new RuntimeException("Overlap in initial state");
        }
        
        Vector leafPos0 = molecule.getChildList().get(0).getPosition();
		Vector leaftPos1 = molecule.getChildList().get(1).getPosition();
		
		molecOrientation.Ev1Mv2(leaftPos1, leafPos0);
		molecOrientation.normalize();
		
		workVector.Ev1Mv2(molecOrientation, initMolecOrientation[iMolecule][0]);
		
		double check = Math.sqrt(workVector.squared());
		double angleMol;
		if(check < 1e-7){
			angleMol = 0.0;
		} else {
			angleMol = Math.acos(molecOrientation.dot(initMolecOrientation[iMolecule][0]));
		}
		
		if(check !=0.0){
            doTransformToInitial(iMolecule, angleMol);
	        
		}

        double randomValue = (2*random.nextDouble() - 1.0);
        double constant = (constraintAngle/180.0) *Math.PI;
        double dTheta = randomValue*constant;
        //rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

        r0.E(positionDefinition.position(molecule));
        //doTransform();
		doTransformMolec(iMolecule, dTheta);
        
		molecOrientation.Ev1Mv2(leaftPos1, leafPos0);
		molecOrientation.normalize();
		workVector.Ev1Mv2(molecOrientation, initMolecOrientation[iMolecule][0]);
		angleMol = Math.acos(molecOrientation.dot(initMolecOrientation[iMolecule][0]));
		
        energyMeter.setTarget(molecule);
        
        if(Double.isNaN(energyMeter.getDataAsScalar())){
        	double energy = energyMeter.getDataAsScalar();
        	System.out.println("energy"+ energy);
        	System.out.println("dTheta: " + dTheta + " angleMol: "+ Degree.UNIT.fromSim(angleMol));
        	throw new RuntimeException("<MCMoveRotate> energy is NaN!!!!!!!!!!!!");
        }
        return true;
    }
    
    protected void doTransform() {
        IAtomList childList = molecule.getChildList();
        for (int iChild = 0; iChild<childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
        }
    }
    
    protected void doTransformMolec(int iMolecule, double dTheta) {
        IAtomList childList = molecule.getChildList();
        for (int iChild = 0; iChild<childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.ME(r0);
        }
        double s = (2*random.nextDouble() - 1.0);
        double t = (2*random.nextDouble() - 1.0);
        
        rotationAxis.E(new double[]{0.0, 0.0, 0.0});
		rotationAxis.PEa1Tv1(s, initMolecOrientation[iMolecule][1]);
		rotationAxis.PEa1Tv1(t, initMolecOrientation[iMolecule][2]);
		rotationAxis.normalize();
		
		Vector vec = space.makeVector();
		vec.E(initMolecOrientation[iMolecule][0]);

		double check = vec.dot(rotationAxis);
		if(check > 1e-7){
			System.out.println("check: " + check);
			throw new RuntimeException("<MCMoveRotateMolecule3DConstraint> The rotation is NOT perpendicular to" +
					" the nominal molecular axis");
		}
		
		rotation.setRotationAxis(rotationAxis, dTheta);
        ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(rotation);
        atomGroupAction.actionPerformed(molecule);
        
        for (int iChild = 0; iChild<childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.PE(r0);
        }
        
    }
    
    protected void doTransformToInitial(int iMolecule, double angleMol) {
        IAtomList childList = molecule.getChildList();
        for (int iChild = 0; iChild<childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.ME(r0);
        }
        
		rotationAxis.E(molecOrientation);
		rotationAxis.XE(initMolecOrientation[iMolecule][0]);
		rotationAxis.normalize();
      
		rotation.setRotationAxis(rotationAxis, angleMol);
        ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(rotation);
        atomGroupAction.actionPerformed(molecule);
        
        for (int iChild = 0; iChild<childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.PE(r0);
        }
        
    }
    
    public void rejectNotify() {
        rotationTensor.invert();
        doTransform();
    }

    public double getChi(double temperature) {
//    	double energy = energyMeter.getDataAsScalar();
//    	if(Double.isInfinite(energy)){
//    		return -1.0;
//    	}
        return 1.0; //always accept the rotational move, See the acceptance criteria in IntegratorMC
    }
    
    public double energyChange() { return 0.0;}
}
