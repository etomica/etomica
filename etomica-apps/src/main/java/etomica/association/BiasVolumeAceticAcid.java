/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.models.OPLS.SpeciesAceticAcid;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePair;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.util.random.IRandom;

public class BiasVolumeAceticAcid extends BiasVolumeMolecule {
    
    /**insert acetic acid with certain angles and distances to satisfy dimerization criteria
     * 
	 * @author Hye Min Kim
	 */
	private static final long serialVersionUID = 1L;
	private double radius;
    private double innerRadius;
    private final IRandom random;
    private Boundary boundary;
    private double maxCosTheta,maxCosPhi;
    protected final MoleculePair pair;
    protected final Vector H1O2, H2O1, C1CH31, C2CH32, OO1, OO2, C2C1, secondAxis, thirdAxis, newPositionC, work1, work2, C1SBO1, C1DBO1, C2SBO2, C2DBO2,dv;
    protected Vector groupTranslationVector;
    protected MoleculeChildAtomAction moveMoleculeAction;
    protected final RotationTensor3D rotationTensor;
    
    public BiasVolumeAceticAcid(Space space, IRandom random, Box box){
        super(space);
        this.random = random;
        pair = new MoleculePair();
		H1O2 = space.makeVector();
		H2O1 = space.makeVector();
		C1CH31 = space.makeVector();
		C2CH32 = space.makeVector();
		C1SBO1 = space.makeVector();
		C1DBO1 = space.makeVector();
		C2SBO2 = space.makeVector();
		C2DBO2 = space.makeVector();
		OO1 = space.makeVector();
		OO2 = space.makeVector();
		C2C1 = space.makeVector();
		secondAxis = space.makeVector();
		thirdAxis = space.makeVector();
		newPositionC = space.makeVector();
		work1 = space.makeVector();
		work2 = space.makeVector();
		dv = space.makeVector();
        rotationTensor = (RotationTensor3D)(space.makeRotationTensor());
    	radius= 4.2;
    	innerRadius = 3.8;
    	maxCosTheta = -0.9;
    	maxCosPhi = -0.9;
    	AtomActionTranslateBy translator = new AtomActionTranslateBy(space);
        groupTranslationVector = translator.getTranslationVector();
        moveMoleculeAction = new MoleculeChildAtomAction(translator);
    }
    
    public void setBox(Box box) {
    	boundary = box.getBoundary();
    }
    
    public void setBiasSphereRadius(double radius) {
        this.radius = radius;
    }
    
    public void setBiasSphereInnerRadius(double innerRadius) {
        this.innerRadius = innerRadius;
    }
    
    public void setmaxCosPhi(double maxCosPhi) {
        this.maxCosPhi = maxCosPhi;
    }
    
    public void setmaxCosTheta(double maxCostheta) {
        this.maxCosTheta = maxCosTheta;
    }
    
    public double getBiasSphereRadius() {return radius;}
    
    public double biasVolume() {
    	//bonding volume = (volume of shell) * (2 cone fractions from 2 theta) * (2*phi/360)
    	double partialVolume = Math.PI*(radius*radius*radius-innerRadius*innerRadius*innerRadius)*(2.0/3.0*(1+maxCosTheta))*(2.0/3.0*(1+maxCosTheta))*Math.acos(maxCosPhi)/Math.PI;
        double totalVolume = 4.0/3.0*Math.PI*(radius*radius*radius-innerRadius*innerRadius*innerRadius);//volume of shell
        return partialVolume*partialVolume/totalVolume;
    }
    
    // Insert molecule1 in to the Bonding region of molecule2 , randomly rotate/move molecule1 and check the molecule2 and then check the molecule1
    public void biasInsert(IMolecule molecule1, IMolecule molecule2) {
    	double randomX = random.nextDouble();
    	double r = Math.pow(innerRadius*innerRadius*innerRadius+(radius*radius*radius-innerRadius*innerRadius*innerRadius)*randomX, 1.0/3.0);
    	double a = 1.0/(1+maxCosTheta);
    	randomX = random.nextDouble();
    	double cosTheta = (a*maxCosTheta-randomX)/a;//uniform theta distribution(X is linear)
		C2CH32.Ev1Mv2(molecule2.getChildList().get(SpeciesAceticAcid.indexCH3).getPosition(), molecule2.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
		C2CH32.normalize();
		secondAxis.E(0.0);
		if (Math.abs(C2CH32.getX(0))> 0.5){
			secondAxis.setX(1, 1);//y direction
		}
		else {
			secondAxis.setX(0, 1);//x direction
		}
		secondAxis.PEa1Tv1(-secondAxis.dot(C2CH32), C2CH32);
		secondAxis.normalize();
    	thirdAxis.E(secondAxis);
    	thirdAxis.XE(C2CH32);
    	double randomAlpha = random.nextDouble()*2*Math.PI;//from 0 to 2pi
    	
    	//newPositionC=R*[a1*cosTheta+(a2*cosalpha+a3*sinalpha)*sintheta], a1, a2, a3: axis for rotation
    	newPositionC.Ea1Tv1(Math.cos(randomAlpha), secondAxis);
    	newPositionC.PEa1Tv1(Math.sin(randomAlpha), thirdAxis);
    	newPositionC.TE(Math.sqrt(1-cosTheta*cosTheta));
    	newPositionC.PEa1Tv1(cosTheta, C2CH32);
    	newPositionC.TE(r);
    	newPositionC.PE(molecule2.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
        groupTranslationVector.Ev1Mv2(newPositionC, molecule1.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
        moveMoleculeAction.actionPerformed(molecule1);//molecule1 is translated, only C is perfect.The location of C is between 3.8 and 4.2. We should check other atoms. 
        //pick a phi for the perfect case
        //phi: angles between OO of 1 and oo of 2
        
        C1CH31.Ev1Mv2(molecule1.getChildList().get(SpeciesAceticAcid.indexCH3).getPosition(), molecule1.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
		C1CH31.normalize();
        C2C1.Ev1Mv2(molecule2.getChildList().get(SpeciesAceticAcid.indexC).getPosition(), molecule1.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
        boundary.nearestImage(C2C1);
        C2C1.normalize();
        
    	work1.Ea1Tv1(-1.0, C2CH32);
    	double psi = calcAngle(C1CH31, work1);
    	work1.XE(C1CH31);
    	if (work1.isZero() && psi > Math.PI/2.0){
    		doFlip(molecule1);
    	}
    	else if (!work1.isZero()){
        	work1.normalize();
        	doTransform(molecule1,newPositionC,work1,-psi);//Now C1CH31 and C2CH32 are parallel.
    	}
    	
        C1CH31.Ev1Mv2(molecule1.getChildList().get(SpeciesAceticAcid.indexCH3).getPosition(), molecule1.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
		C1CH31.normalize();

    	//perpendicular vector of C1CH31 on the place C1dBO1SBO1
    	C1SBO1.Ev1Mv2(molecule1.getChildList().get(SpeciesAceticAcid.indexSBO).getPosition(), molecule1.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
    	C1DBO1.Ev1Mv2(molecule1.getChildList().get(SpeciesAceticAcid.indexDBO).getPosition(), molecule1.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
    	work1.E(C1SBO1);
    	work1.XE(C1DBO1);
    	work1.XE(C1CH31);
    	work1.normalize();
    	
    	//perpendicular vector of C2CH32 on the place C2dBO2SBO2
    	C2SBO2.Ev1Mv2(molecule2.getChildList().get(SpeciesAceticAcid.indexSBO).getPosition(), molecule2.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
    	C2DBO2.Ev1Mv2(molecule2.getChildList().get(SpeciesAceticAcid.indexDBO).getPosition(), molecule2.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
    	work2.E(C2SBO2);
    	work2.XE(C2DBO2);
    	work2.XE(C2CH32);
    	work2.normalize();
    	secondAxis.E(work1);
    	secondAxis.XE(work2);
    	secondAxis.normalize();

    	double phi = calcAngle(work1, work2);
    	
    	//System.out.println("phi "+phi+ " work1 "+work1);
    	if (!secondAxis.isNaN() && secondAxis.dot(C1CH31) < 0) phi = -phi;//to distinguish negative phi and positive phi
    	randomX = random.nextDouble();
    	double minPhi = Math.acos(maxCosPhi);
    	double newPhi = minPhi+ randomX*2*(Math.PI-minPhi);//range of phi = 2*maxPhi, this is phi that we want
    	doTransform(molecule1,newPositionC,C1CH31,phi-newPhi);//check if the rotation is fine
        rotationTensor.transform(work1);

        randomX = random.nextDouble();
        double cosTheta2 = (a*maxCosTheta-randomX)/a;
        randomAlpha = random.nextDouble()*2*Math.PI;
		secondAxis.E(0.0);
		if (Math.abs(C2C1.getX(0))> 0.5){
			secondAxis.setX(1, 1);//y direction
		}
		else {
			secondAxis.setX(0, 1);//x direction
		}
		secondAxis.PEa1Tv1(-secondAxis.dot(C2C1), C2C1);
		secondAxis.normalize();
    	thirdAxis.E(secondAxis);
    	thirdAxis.XE(C2C1);
        work1.Ea1Tv1(Math.cos(randomAlpha), secondAxis);
    	work1.PEa1Tv1(Math.sin(randomAlpha), thirdAxis);
    	work1.TE(Math.sqrt(1-cosTheta2*cosTheta2));
    	work1.PEa1Tv1(cosTheta2, C2C1);
    	
    	double psi2 = calcAngle(C1CH31, work1);
    	work1.XE(C1CH31);
    	if (!work1.isZero()){
        	work1.normalize();
        	doTransform(molecule1,newPositionC,work1,-psi2);
    	}

    	C1CH31.Ev1Mv2(molecule1.getChildList().get(SpeciesAceticAcid.indexCH3).getPosition(), molecule1.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
		C1CH31.normalize();
        
    	if (!isAssociated(molecule1,molecule2)){
    		System.out.println("r "+r+" costheta "+cosTheta+" cosTheta2 "+cosTheta2+" cosPhi "+Math.cos(newPhi));
    		isAssociated(molecule1,molecule2);
    		throw new RuntimeException();
    	}   	
    }
    
    protected void doFlip(IMolecule molecule){
    	IAtomList childList = molecule.getChildList();
    	Vector rC = childList.get(SpeciesAceticAcid.indexC).getPosition();
    	for (int i = 0; i<childList.size(); i+=1){
    		if (i == SpeciesAceticAcid.indexC)continue;
    		childList.get(i).getPosition().TE(-1);
    		childList.get(i).getPosition().PEa1Tv1(2, rC);
    	}
    }
    
    protected double calcAngle(Vector v1, Vector v2){
    	dv.Ev1Mv2(v1, v2);
    	double l = 0.5*Math.sqrt(dv.squared());
    	return Math.asin(l)*2;
    }
    
    protected void doTransform(IMolecule molecule, Vector r0, Vector axis, double theta) {
        IAtomList childList = molecule.getChildList();
        rotationTensor.setRotationAxis(axis, theta);
        for (int iChild = 0; iChild<childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.ME(r0);
            rotationTensor.transform(r);
            r.PE(r0);
        }
    }

    /**
     *Function to check for bonding
     *calculate the distance between atom1 and atom2
     */
    
    public boolean isAssociated(IMolecule molecule1, IMolecule molecule2){
    	
		C2CH32.Ev1Mv2(molecule2.getChildList().get(SpeciesAceticAcid.indexCH3).getPosition(), molecule2.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
		C2CH32.normalize();
        C2C1.Ev1Mv2(molecule1.getChildList().get(SpeciesAceticAcid.indexC).getPosition(), molecule2.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
        boundary.nearestImage(C2C1);
        double distance = C2C1.squared();
        if (distance < innerRadius*innerRadius || distance > radius*radius){
        	return false;
        }
        C2C1.normalize();
        double cosTheta = C2CH32.dot(C2C1);
        if (cosTheta > maxCosTheta){
        	return false;
        }
        
        C1CH31.Ev1Mv2(molecule1.getChildList().get(SpeciesAceticAcid.indexCH3).getPosition(), molecule1.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
		C1CH31.normalize();
        double cosTheta2 = -C1CH31.dot(C2C1);
        if (cosTheta2 > maxCosTheta){
        	return false;
        }
    	//perpendicular vector of C1CH31 on the place C1dBO1SBO1
    	C1SBO1.Ev1Mv2(molecule1.getChildList().get(SpeciesAceticAcid.indexSBO).getPosition(), molecule1.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
    	C1DBO1.Ev1Mv2(molecule1.getChildList().get(SpeciesAceticAcid.indexDBO).getPosition(), molecule1.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
    	work1.E(C1SBO1);
    	work1.XE(C1DBO1);
    	work1.XE(C1CH31);
    	work1.normalize();
    	
    	double psi = Math.PI-calcAngle(C1CH31, C2CH32);
    	secondAxis.E(C2CH32);
    	secondAxis.XE(C1CH31);
    	//double cosPsi = C1CH31.dot(C2CH32);//C1CH31 and C2CH32 are pointing opposite
    	if (!secondAxis.isZero()){
        	//new axis to rotate molecule1 to make CH3 happy
        	secondAxis.normalize();
	    	rotationTensor.setRotationAxis(secondAxis, psi);
	        rotationTensor.transform(work1);
	        rotationTensor.transform(C1CH31);
    	}
    	
    	//perpendicular vector of C2CH32 on the place C2dBO2SBO2
    	C2SBO2.Ev1Mv2(molecule2.getChildList().get(SpeciesAceticAcid.indexSBO).getPosition(), molecule2.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
    	C2DBO2.Ev1Mv2(molecule2.getChildList().get(SpeciesAceticAcid.indexDBO).getPosition(), molecule2.getChildList().get(SpeciesAceticAcid.indexC).getPosition());
    	work2.E(C2SBO2);
    	work2.XE(C2DBO2);
    	work2.XE(C2CH32);
    	work2.normalize();
    	double cosPhi = work1.dot(work2);
    	work1.XE(work2);
    	work1.normalize();
		if ( cosPhi <maxCosPhi){
    		return true;
		}
    	return false;
    }
}
