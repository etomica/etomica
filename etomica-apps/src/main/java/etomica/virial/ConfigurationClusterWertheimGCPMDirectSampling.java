/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.action.MoleculeActionTranslateTo;
import etomica.api.*;
import etomica.atom.*;
import etomica.box.Box;
import etomica.models.water.PNWaterGCPMThreeSite;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space.RotationTensor;
import etomica.util.random.IRandom;

public class ConfigurationClusterWertheimGCPMDirectSampling extends ConfigurationCluster {
	
	protected IRandom random;
	protected PNWaterGCPMThreeSite associationPotential;
	protected PNWaterGCPMThreeSite associationPotential2;
	protected PNWaterGCPMThreeSite associationPotential3;
	protected PNWaterGCPMThreeSite nonAssociationPotential;
	protected double diagramIndex;

	public ConfigurationClusterWertheimGCPMDirectSampling(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential) {
		super(_space);
		this.associationPotential = associationPotential;
		this.random = random;
	}

	public ConfigurationClusterWertheimGCPMDirectSampling(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential, PNWaterGCPMThreeSite associationPotential2 ) {
		super(_space);
		this.associationPotential = associationPotential;
		this.associationPotential2 = associationPotential2;
		this.random = random;
	}
	public ConfigurationClusterWertheimGCPMDirectSampling(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential, PNWaterGCPMThreeSite associationPotential2, PNWaterGCPMThreeSite associationPotential3 ) {
		super(_space);
		this.associationPotential = associationPotential;
		this.associationPotential2 = associationPotential2;
		this.associationPotential3 = associationPotential3;
		this.random = random;
	}
	
	public ConfigurationClusterWertheimGCPMDirectSampling(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential, PNWaterGCPMThreeSite associationPotential2, double diagramIndex ) {
		super(_space);
		this.associationPotential = associationPotential;
		this.associationPotential2 = associationPotential2;
		this.random = random;
		this.diagramIndex = diagramIndex;
	}

	public void initializeCoordinatesTrimer(Box box) {

		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		associationPotential2.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair1 = new MoleculePair();
		MoleculePair pair2 = new MoleculePair();
		pair1.atom0 = list.getMolecule(0);
		pair1.atom1 = list.getMolecule(1);
		pair2.atom0 = list.getMolecule(1);
		pair2.atom1 = list.getMolecule(2);
		
		double[] d = new double[] {3.5,0,0};
		double[] e = new double[] {7,0.0,0};
		translationTrimer(d,e,box);
		
		association(pair1,box);
		association2(pair2,box);

		 clusterBox.trialNotify();
		 clusterBox.acceptNotify();
		 System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	public void initializeCoordinatesTetramer(Box box) {

		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		associationPotential2.setBox(box);
		associationPotential3.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair1 = new MoleculePair();
		MoleculePair pair2 = new MoleculePair();
		MoleculePair pair3 = new MoleculePair();
		pair1.atom0 = list.getMolecule(0);
		pair1.atom1 = list.getMolecule(1);
		pair2.atom0 = list.getMolecule(1);
		pair2.atom1 = list.getMolecule(2);
		pair3.atom0 = list.getMolecule(2);
		pair3.atom1 = list.getMolecule(3);
		
		double[] d = new double[] {3.5,0,0};
		double[] e = new double[] {3.5,3.5,0};
		double[] f = new double[] {0,3.5,0};
		translationTetramer(d,e,f,box);
		
		association(pair1,box);
		association2(pair2,box);
		association3(pair3,box);

		 clusterBox.trialNotify();
		 clusterBox.acceptNotify();
		 System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	public void initializeCoordinatesBranch(Box box) {

		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		associationPotential2.setBox(box);
		associationPotential3.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair1 = new MoleculePair();
		MoleculePair pair2 = new MoleculePair();
		MoleculePair pair3 = new MoleculePair();
		pair1.atom0 = list.getMolecule(0);
		pair1.atom1 = list.getMolecule(1);
		pair2.atom0 = list.getMolecule(1);
		pair2.atom1 = list.getMolecule(2);
		pair3.atom0 = list.getMolecule(1);
		pair3.atom1 = list.getMolecule(3);
		
		double[] d = new double[] {3.5,0,0};
		double[] e = new double[] {3.5,3.5,0};
		double[] f = new double[] {0,0,0};
		translationTetramer(d,e,f,box);
		
		association(pair1,box);
		association2(pair2,box);
		association4(pair3,box);

		 clusterBox.trialNotify();
		 clusterBox.acceptNotify();
		 System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	public void translationTrimer(double[] d,double[] e, Box box){//place molecule1,2 at some position
		IMoleculeList list = box.getMoleculeList();
        MoleculeActionTranslateTo translationA = new MoleculeActionTranslateTo(space);
        MoleculeActionTranslateTo translationB = new MoleculeActionTranslateTo(space);
        Vector a = space.makeVector();
        Vector b = space.makeVector();
        a.E(d);
        b.E(e);
        translationA.setDestination(a);
        translationB.setDestination(b);
		IMolecule mol1 = list.getMolecule(1);
		IMolecule mol2 = list.getMolecule(2);
        translationA.actionPerformed(mol1);
        translationB.actionPerformed(mol2);
	}
	public void translationTetramer(double[] d,double[] e,double[] f, Box box){//place molecule1,2 at some position
		IMoleculeList list = box.getMoleculeList();
        MoleculeActionTranslateTo translationA = new MoleculeActionTranslateTo(space);
        MoleculeActionTranslateTo translationB = new MoleculeActionTranslateTo(space);
        MoleculeActionTranslateTo translationC = new MoleculeActionTranslateTo(space);
        Vector a = space.makeVector();
        Vector b = space.makeVector();
        Vector c = space.makeVector();
        a.E(d);
        b.E(e);
        c.E(f);
        translationA.setDestination(a);
        translationB.setDestination(b);
        translationC.setDestination(c);
		IMolecule mol1 = list.getMolecule(1);
		IMolecule mol2 = list.getMolecule(2);
		IMolecule mol3 = list.getMolecule(3);
        translationA.actionPerformed(mol1);
        translationB.actionPerformed(mol2);
        translationC.actionPerformed(mol3);
	}
	public void association(MoleculePair pair, Box box){
		RotationTensor rotationTensor = space.makeRotationTensor();
		Vector r0 = space.makeVector();
		IMoleculePositionDefinition positionDefinition = new MoleculePositionGeometricCenter(space);
		IMoleculeList list = box.getMoleculeList();
		pair.atom0 = list.getMolecule(0);
		IMolecule water = list.getMolecule(1);
		pair.atom1 = water;

		
        while (true){
	        Vector positionWater = space.makeVector();
	        positionWater.setRandomInSphere(random);
	        positionWater.TE(4.0);//place water molecule within a sphere with r = 4A
	        MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
	        translation.setDestination(positionWater);
	        translation.actionPerformed(water);
	       
	        if (associationPotential.energy(pair) != 0.0){//when there is an association, fix the position of water
        		break;
        	}
	        double dTheta = (2*random.nextDouble() - 1.0)*Math.PI;
	        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

	        r0.E(positionDefinition.position(water));
		    IAtomList childList = water.getChildList();
		    for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {//free rotation until finding association
		        IAtom a = childList.getAtom(iChild);
		        Vector r = a.getPosition();
		        r.ME(r0);
		        box.getBoundary().nearestImage(r);
		        rotationTensor.transform(r);
		        r.PE(r0);
	    	}
	        if (associationPotential.energy(pair) != 0.0){
        		break;
        	}
		}
	}
		
	public void association2(MoleculePair pair, Box box){
		RotationTensor rotationTensor = space.makeRotationTensor();
		Vector r0 = space.makeVector();
		IMoleculePositionDefinition positionDefinition = new MoleculePositionGeometricCenter(space);
		IMoleculeList list = box.getMoleculeList();
		pair.atom0 = list.getMolecule(1);
		IMolecule water = list.getMolecule(2);
		pair.atom1 = water;
	
        while (true){
        	Vector positionWater = space.makeVector();
	        positionWater.setRandomInSphere(random);
	        positionWater.TE(4.0);//place water molecule within a sphere with r = 8A
	        positionWater.PE(pair.atom0.getChildList().getAtom(0).getPosition());
	        MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
	        translation.setDestination(positionWater);
	        translation.actionPerformed(water);

	       
	        if (associationPotential2.energy(pair) != 0.0){//when there is an association, fix the position of water
        		break;
        	}
	        double dTheta = (2*random.nextDouble() - 1.0)*Math.PI;
	        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

	        r0.E(positionDefinition.position(water));
		    IAtomList childList = water.getChildList();
		    for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {//free rotation until finding association
		        IAtom a = childList.getAtom(iChild);
		        Vector r = a.getPosition();
		        r.ME(r0);
		        box.getBoundary().nearestImage(r);
		        rotationTensor.transform(r);
		        r.PE(r0);
	    	}
	        if (associationPotential2.energy(pair) != 0.0){
        		break;
        	}
		}
	}
	public void association3(MoleculePair pair, Box box){
		RotationTensor rotationTensor = space.makeRotationTensor();
		Vector r0 = space.makeVector();
		IMoleculePositionDefinition positionDefinition = new MoleculePositionGeometricCenter(space);
		IMoleculeList list = box.getMoleculeList();
		pair.atom0 = list.getMolecule(2);
		IMolecule water = list.getMolecule(3);
		pair.atom1 = water;
	
        while (true){
        	Vector positionWater = space.makeVector();
	        positionWater.setRandomInSphere(random);
	        positionWater.TE(4.0);//place water molecule within a sphere with r = 8A
	        positionWater.PE(pair.atom0.getChildList().getAtom(0).getPosition());
	        MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
	        translation.setDestination(positionWater);
	        translation.actionPerformed(water);

	       
	        if (associationPotential3.energy(pair) != 0.0){//when there is an association, fix the position of water
        		break;
        	}
	        double dTheta = (2*random.nextDouble() - 1.0)*Math.PI;
	        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

	        r0.E(positionDefinition.position(water));
		    IAtomList childList = water.getChildList();
		    for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {//free rotation until finding association
		        IAtom a = childList.getAtom(iChild);
		        Vector r = a.getPosition();
		        r.ME(r0);
		        box.getBoundary().nearestImage(r);
		        rotationTensor.transform(r);
		        r.PE(r0);
	    	}
	        if (associationPotential3.energy(pair) != 0.0){
        		break;
        	}
		}
	}
	public void association4(MoleculePair pair, Box box){
		RotationTensor rotationTensor = space.makeRotationTensor();
		Vector r0 = space.makeVector();
		IMoleculePositionDefinition positionDefinition = new MoleculePositionGeometricCenter(space);
		IMoleculeList list = box.getMoleculeList();
		pair.atom0 = list.getMolecule(1);
		IMolecule water = list.getMolecule(3);
		pair.atom1 = water;
	
        while (true){
        	Vector positionWater = space.makeVector();
	        positionWater.setRandomInSphere(random);
	        positionWater.TE(4.0);//place water molecule within a sphere with r = 8A
	        positionWater.PE(pair.atom0.getChildList().getAtom(0).getPosition());
	        MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
	        translation.setDestination(positionWater);
	        translation.actionPerformed(water);

	       
	        if (associationPotential3.energy(pair) != 0.0){//when there is an association, fix the position of water
        		break;
        	}
	        double dTheta = (2*random.nextDouble() - 1.0)*Math.PI;
	        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

	        r0.E(positionDefinition.position(water));
		    IAtomList childList = water.getChildList();
		    for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {//free rotation until finding association
		        IAtom a = childList.getAtom(iChild);
		        Vector r = a.getPosition();
		        r.ME(r0);
		        box.getBoundary().nearestImage(r);
		        rotationTensor.transform(r);
		        r.PE(r0);
	    	}
	        if (associationPotential3.energy(pair) != 0.0){
        		break;
        	}
		}
	}

}
