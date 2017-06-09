/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.util.random.IRandom;
import etomica.space.Vector;
import etomica.atom.MoleculePair;
import etomica.models.OPLS.SpeciesAceticAcid;
import etomica.space.Space;
import etomica.space3d.RotationTensor3D;

public class BiasVolume2SiteAceticAcid implements AssociationDefinitionMolecule {
    
    /**insert acetic acid with certain angles and distances to satisfy dimerization criteria
     * 
	 * @author Hye Min Kim
	 */
	private static final long serialVersionUID = 1L;
	private double radius;
    private double innerRadius;
    private double innerOHLength;
    private double outerOHLength;
    private final IRandom random;
    private Boundary boundary;
    private double maxCosTheta,maxCosPhi, maxCosOHO, maxCosHOC, minCosHOC;
    protected final MoleculePair pair;
    protected final Vector H1O2, H2O1, C1CH31, C2CH32, OO1, OO2, C2C1, secondAxis, thirdAxis, newPositionC, work1, work2;
    protected final Vector C1SBO1, C1DBO1, C2SBO2, C2DBO2,dv, H1DBO2, H1SBO1, DBO2H1, DBO2C2;
    protected Vector groupTranslationVector;
    protected MoleculeChildAtomAction moveMoleculeAction;
    protected final RotationTensor3D rotationTensor;

    
    public BiasVolume2SiteAceticAcid(Space space, IRandom random){
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
		H1DBO2 = space.makeVector();
		H1SBO1 = space.makeVector();
		DBO2H1 = space.makeVector();
		DBO2C2 = space.makeVector();
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
    	innerOHLength = 1.8;
    	outerOHLength = 2.3;
    	maxCosTheta = -0.9;
    	maxCosPhi = -0.9;
    	maxCosOHO = -0.7;
    	maxCosHOC = -0.4;
    	minCosHOC = -0.9;
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
    
    public void setmaxCosTheta(double maxCosTheta) {
        this.maxCosTheta = maxCosTheta;
    }
    
    public double getBiasSphereRadius() {return radius;}


    
    protected void doFlip(IMolecule molecule){
    	IAtomList childList = molecule.getChildList();
    	Vector rC = childList.getAtom(SpeciesAceticAcid.indexC).getPosition();
    	for (int i = 0;i<childList.getAtomCount();i+=1){
    		if (i == SpeciesAceticAcid.indexC)continue;
    		childList.getAtom(i).getPosition().TE(-1);
    		childList.getAtom(i).getPosition().PEa1Tv1(2, rC);
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
        for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {
            IAtom a = childList.getAtom(iChild);
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
    	
        H1DBO2.Ev1Mv2(molecule2.getChildList().getAtom(SpeciesAceticAcid.indexDBO).getPosition(), molecule1.getChildList().getAtom(SpeciesAceticAcid.indexH).getPosition());
        boundary.nearestImage(H1DBO2);
        double distance1 = H1DBO2.squared();
        H1DBO2.normalize();
        H1SBO1.Ev1Mv2(molecule1.getChildList().getAtom(SpeciesAceticAcid.indexSBO).getPosition(), molecule1.getChildList().getAtom(SpeciesAceticAcid.indexH).getPosition());
        H1SBO1.normalize();
        double cosOHO1 = H1SBO1.dot(H1DBO2);
        DBO2H1.Ev1Mv2(molecule1.getChildList().getAtom(SpeciesAceticAcid.indexH).getPosition(), molecule2.getChildList().getAtom(SpeciesAceticAcid.indexDBO).getPosition());
        DBO2C2.Ev1Mv2(molecule2.getChildList().getAtom(SpeciesAceticAcid.indexC).getPosition(), molecule2.getChildList().getAtom(SpeciesAceticAcid.indexDBO).getPosition());
        DBO2H1.normalize();
        DBO2C2.normalize();
        double cosHOC1 = DBO2H1.dot(DBO2C2);
        if (distance1 > innerOHLength*innerOHLength && distance1 < outerOHLength*outerOHLength && cosOHO1 < maxCosOHO && cosHOC1 < maxCosHOC && cosHOC1 > minCosHOC){
        	return true;
        }
    	return false;
    }
}
