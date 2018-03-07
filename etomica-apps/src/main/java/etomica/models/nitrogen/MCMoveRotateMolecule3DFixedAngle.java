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
 * the proposed angle criteria is based on u3 and u4 distribution
 * 
 * @author Tai Boon Tan
 *
 */
public class MCMoveRotateMolecule3DFixedAngle extends MCMoveMolecule {
    
    private static final long serialVersionUID = 2L;
    protected transient Vector r0;
    protected transient RotationTensor rotationTensor;
    protected IMoleculePositionDefinition positionDefinition;
    protected double constraintAngle;
    protected CoordinateDefinitionNitrogen coordinateDef;
    protected Vector[][] initMolecOrientation;
    protected Vector molecOrientation, rotationAxis;
	protected RotationTensor3D rotation;
	protected Tensor3D tensor;
    protected final MoleculeChildAtomAction atomGroupAction;
    
    public MCMoveRotateMolecule3DFixedAngle(PotentialMaster potentialMaster, IRandom random,
                                            Space _space, double angle, CoordinateDefinitionNitrogen coordinateDef, Box box) {
        super(potentialMaster, random, _space, 0.5*Math.PI, Math.PI);
        rotationTensor = _space.makeRotationTensor();
        r0 = _space.makeVector();
        positionDefinition = new MoleculePositionGeometricCenter(space);
        constraintAngle = angle;
        this.coordinateDef = coordinateDef;
        
        int numMolec = box.getMoleculeList().size();
     	initMolecOrientation = new Vector[numMolec][3];
     	molecOrientation = space.makeVector();
     	rotationAxis = space.makeVector();
     	
    	tensor = new Tensor3D(new double[][]{{1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{0.0, 0.0, 1.0}});
		rotation = new RotationTensor3D();
		rotation.E(tensor);
    	/*
		 * initializing the initial orientation of the molecule
		 */
		for (int i=0; i<numMolec; i++){
			initMolecOrientation[i] = space.makeVectorArray(3);
			initMolecOrientation[i] = coordinateDef.getMoleculeOrientation(box.getMoleculeList().get(i));
		}
		atomGroupAction = new MoleculeChildAtomAction(new AtomActionTransformed(space));
	        
    }
     
    public boolean doTrial() {

        if(box.getMoleculeList().size()==0) {molecule = null; return false;}
        int iMolecule = random.nextInt(box.getMoleculeList().size());
        
        molecule = coordinateDef.getBox().getMoleculeList().get(iMolecule);
        r0.E(positionDefinition.position(molecule));
        
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        if(Double.isInfinite(uOld)) {
            throw new RuntimeException("Overlap in initial state");
        }
        
        double u3 = 2.0; // set u3 and u4 to large values
        double u4 = 2.0;
        double theta = 361;
        while (theta > constraintAngle){
        	u3 = (2*random.nextDouble() - 1.0)*(constraintAngle/180.0) *Math.PI;
        	u4 = (2*random.nextDouble() - 1.0)*(constraintAngle/180.0) *Math.PI;
        	theta = Degree.UNIT.fromSim(Math.acos(1.0000000000000004 - (u3*u3 + u4*u4)*0.5));
        } 
        
        r0.E(positionDefinition.position(molecule));
        setToU(iMolecule, u3, u4);
		
	    Vector leafPos0 = molecule.getChildList().get(0).getPosition();
		Vector leaftPos1 = molecule.getChildList().get(1).getPosition();
		molecOrientation.Ev1Mv2(leaftPos1, leafPos0);
		molecOrientation.normalize();
		
		double angleMol = Math.acos(molecOrientation.dot(initMolecOrientation[iMolecule][0]));
		
        energyMeter.setTarget(molecule);
        
        if(Double.isNaN(energyMeter.getDataAsScalar())){
        	System.out.println("theta: " + theta+ " angleMol: "+ Degree.UNIT.fromSim(angleMol));
        	System.out.println("molecOrientation: " + molecOrientation.toString());
        	System.out.println("initMolecOrientation: " + initMolecOrientation[iMolecule][0].toString());
        	System.out.println("dot: " + molecOrientation.dot(initMolecOrientation[iMolecule][0]));
        	double energy = energyMeter.getDataAsScalar();
        	System.out.println("energy"+ energy);
        	throw new RuntimeException("<MCMoveRotate3DFixedAngle> energy is NaN!!!!!!!!!!!!");
        }
        return true;
    }

    protected void setToU(int iMolecule, double u3, double u4){
    	IAtomList childList = molecule.getChildList();
    	for (int iChild = 0; iChild<childList.size(); iChild++) {
    		IAtom a = childList.get(iChild);
    		Vector r = a.getPosition();
    		r.ME(r0);
    	}
    	
    	Vector rotationAxis = space.makeVector();
    	RotationTensor3D rotation = new RotationTensor3D();
    	rotation.E(tensor);
	         
    	Vector leafPos0 = molecule.getChildList().get(0).getPosition();
    	Vector leafPos1 = molecule.getChildList().get(1).getPosition();
		    	
    	/*
         * a.
         */
	   	Vector axis = space.makeVector();
	   	axis.Ev1Mv2(leafPos1, leafPos0);
	   	axis.normalize();
		    		    	
    	double angle = Math.acos(axis.dot(initMolecOrientation[iMolecule][0]));

    	if (Math.abs(angle) > 5e-8){ // make sure we DO NOT cross-product vectors with very small angle
    		rotationAxis.E(axis);
	    	rotationAxis.XE(initMolecOrientation[iMolecule][0]);
	    	rotationAxis.normalize();
			    	
	    	/*
	    	 * 	 c. rotating clockwise.
		 	 */
	    	rotation.setRotationAxis(rotationAxis, angle);
			    	
	    	if(rotation.isNaN()){	
	    		System.out.println("Step 1 Rotation tensor is BAD!");
				System.out.println("Rotation Angle is too small, angle: "+ angle);
	    		System.out.println("Rotation is not necessary");
	    		System.out.println(rotation);
	    		throw new RuntimeException();
	    	}
	        ((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(rotation);
	        atomGroupAction.actionPerformed(molecule);
    	}

        if (Math.abs(u3)>1e-7 || Math.abs(u4)>1e-7){
		       
       	/*
         * a.	
         */
       	axis.E(0);
    	axis.Ea1Tv1(u3, initMolecOrientation[iMolecule][1]);
    	axis.PEa1Tv1(u4, initMolecOrientation[iMolecule][2]);
    	axis.normalize();
			    	
    	/*
    	 * b.
    	 */
    	angle = Math.acos(1.0000000000000004 - (u3*u3 + u4*u4)*0.5);
    	if(Math.abs(angle) > 1e-7){
	    	rotationAxis.E(0);
	    	rotationAxis.E(axis);
	    	rotationAxis.XE(initMolecOrientation[iMolecule][0]);
	    	rotationAxis.normalize();
				    	
	    	rotation.setRotationAxis(rotationAxis, -angle);
				    	
	    	if(rotation.isNaN()){
	    		System.out.println("Step 2 Rotation tensor is BAD!");
	    		System.out.println("Rotation Angle is too small, angle: "+ angle);
	    		System.out.println("Rotation is not necessary");
	    		System.out.println(rotation);
	    		throw new RuntimeException();
	    	}
				    	
	    	((AtomActionTransformed)atomGroupAction.getAtomAction()).setTransformationTensor(rotation);
	        atomGroupAction.actionPerformed(molecule);
	   		}
        }
    
        for (int iChild = 0; iChild<childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.PE(r0);
        }
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
