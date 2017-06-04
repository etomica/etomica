/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;
import etomica.api.*;
import etomica.atom.MoleculePositionGeometricCenter;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IMoleculePositionDefinition;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.RotationTensor;


/**
 * Impose constraint to prevent the alpha-N2 phase from becoming orientationally-disordered
 *  the constraint imposed is the average cosine theta, where
 *  theta is the orientation angle off the molecular nominal orientation
 * 
 * @author taitan
 *
 */
public class MCMoveRotateMolecule3DN2AveCosThetaConstraint extends MCMoveMolecule {
    
    private static final long serialVersionUID = 2L;
    protected transient Vector r0, molAxis;
    protected transient RotationTensor rotationTensor;
    protected IMoleculePositionDefinition positionDefinition;
    protected CoordinateDefinitionNitrogen coordinateDef;
    protected Vector[] initMolecOrientation;
    protected int numMolecule;
    protected double aveCosTheta = 1.0;
    protected double newTotalCosTheta, oldTotalCosTheta;
    protected double cosThetaConstraint, workAveCosTheta;
    
    public MCMoveRotateMolecule3DN2AveCosThetaConstraint(PotentialMaster potentialMaster, IRandom random,
                                                         Space _space, CoordinateDefinitionNitrogen coordinateDef, double cosThetaConstraint) {
        super(potentialMaster, random, _space, Math.PI/2, Math.PI);
        rotationTensor = _space.makeRotationTensor();
        r0 = _space.makeVector();
        this.coordinateDef = coordinateDef;
        this.cosThetaConstraint = cosThetaConstraint;
        
        positionDefinition = new MoleculePositionGeometricCenter(space);
        
        IMoleculeList moleculeList = coordinateDef.getBox().getMoleculeList();
        numMolecule = moleculeList.getMoleculeCount();
        initMolecOrientation = new Vector[numMolecule];
        molAxis = space.makeVector();
        
        for (int i=0; i<numMolecule; i++){
			initMolecOrientation[i] = space.makeVector();
			initMolecOrientation[i] = coordinateDef.getMoleculeOrientation(moleculeList.getMolecule(i))[0];
		}
        
    }
     
    public boolean doTrial() {
//        System.out.println("doTrial MCMoveRotateMolecule called");
        
        if(box.getMoleculeList().getMoleculeCount()==0) {molecule = null; return false;}
           
        int iMol = random.nextInt(numMolecule);
        molecule = coordinateDef.getBox().getMoleculeList().getMolecule(iMol);
        
        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();
        
        Vector leafPos0 = molecule.getChildList().getAtom(0).getPosition();
    	Vector leafPos1 = molecule.getChildList().getAtom(1).getPosition();

    	molAxis.Ev1Mv2(leafPos1, leafPos0);
       	molAxis.normalize();
       	
        double l2 = molAxis.Mv1Squared(initMolecOrientation[iMol]);
        double oldWorkTotalCosTheta;
       
       	oldWorkTotalCosTheta = aveCosTheta*numMolecule - (1.0 - 0.5*l2);
        
        if(Double.isInfinite(uOld)) {
            throw new RuntimeException("Overlap in initial state");
        }
        double dTheta = (2*random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

        r0.E(positionDefinition.position(molecule));
        doTransform();
        energyMeter.setTarget(molecule);
                
        molAxis.Ev1Mv2(leafPos1, leafPos0);
       	molAxis.normalize();
        
       	l2 = molAxis.Mv1Squared(initMolecOrientation[iMol]);
       	
   		workAveCosTheta = (oldWorkTotalCosTheta + (1.0 - 0.5*l2))/numMolecule;

//       	System.out.println("workAveCosTheta: " + workAveCosTheta);
        // hitting the constraint
        if(workAveCosTheta<cosThetaConstraint){
       		uNew = Double.POSITIVE_INFINITY;
//       		System.out.println("hit constraint");
       	
        } else {
       		uNew = energyMeter.getDataAsScalar();
    
       	}

        return true;
    }
    
    public double getB() {
        return -(uNew - uOld);
        
    }
    
    protected void doTransform() {
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
        doTransform();
    }
    
    public void acceptNotify() { 
    	aveCosTheta = workAveCosTheta;
    }
    
    public void calcAveCosThetaInitial(){
        double totalCosTheta = 0.0;
        for (int i=0; i<numMolecule; i++){
		    IMolecule molec = coordinateDef.getBox().getMoleculeList().getMolecule(i);
		        
		    Vector leafPos0 = molec.getChildList().getAtom(0).getPosition();
		    Vector leafPos1 = molec.getChildList().getAtom(1).getPosition();

		    molAxis.Ev1Mv2(leafPos1, leafPos0);
		    molAxis.normalize();
		       	
		    double l2 = molAxis.Mv1Squared(initMolecOrientation[i]);
		    totalCosTheta += (1.0 - 0.5*l2);
        }
        
        aveCosTheta = totalCosTheta/numMolecule;
    }
}
