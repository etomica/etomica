/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.action.MoleculeActionTranslateTo;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.models.water.PNWaterGCPMThreeSite;
import etomica.molecule.*;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class ConfigurationClusterWertheimGCPM4Pt extends ConfigurationCluster {
	
	protected IRandom random;
	protected PNWaterGCPMThreeSite associationPotential;
	protected PNWaterGCPMThreeSite associationPotential2;
	protected PNWaterGCPMThreeSite associationPotential3;
	protected PNWaterGCPMThreeSite nonAssociationPotential;
	protected double diagramIndex;

	public ConfigurationClusterWertheimGCPM4Pt(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential) {
		super(_space);
		this.associationPotential = associationPotential;
		this.random = random;
	}

	public ConfigurationClusterWertheimGCPM4Pt(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential, PNWaterGCPMThreeSite associationPotential2 ) {
		super(_space);
		this.associationPotential = associationPotential;
		this.associationPotential2 = associationPotential2;
		this.random = random;
	}
	public ConfigurationClusterWertheimGCPM4Pt(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential, PNWaterGCPMThreeSite associationPotential2, PNWaterGCPMThreeSite associationPotential3 ) {
		super(_space);
		this.associationPotential = associationPotential;
		this.associationPotential2 = associationPotential2;
		this.associationPotential3 = associationPotential3;
		this.random = random;
	}
	
	public ConfigurationClusterWertheimGCPM4Pt(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential, PNWaterGCPMThreeSite associationPotential2, double diagramIndex ) {
		super(_space);
		this.associationPotential = associationPotential;
		this.associationPotential2 = associationPotential2;
		this.random = random;
		this.diagramIndex = diagramIndex;
	}
	public void initializeCoordinates(Box box) {
		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair = new MoleculePair();
		pair.mol0 = list.get(0);
		pair.mol1 = list.get(1);
		double[] d = new double[] {3.5,0,0};
		double[] e = new double[] {5.0,5.0,0};
		double[] f = new double[] {0,5.0,0};
		translation(d,e,f,box);
		association(pair,box);
		clusterBox.trialNotify();
		clusterBox.acceptNotify();
		System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	
	public void initializeCoordinatesER(Box box) {
		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		double[] d = new double[] {3.5,0,0};
		double[] e = new double[] {5.0,5.0,0};
		double[] f = new double[] {0,5.0,0};
		translation(d,e,f,box);
		clusterBox.trialNotify();
		clusterBox.acceptNotify();
		System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	
	public void initializeCoordinates2(Box box) {
		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		associationPotential2.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair1 = new MoleculePair();
		MoleculePair pair2 = new MoleculePair();
		pair1.mol0 = list.get(0);
		pair1.mol1 = list.get(1);
		pair2.mol0 = list.get(1);
		pair2.mol1 = list.get(2);
		double[] d = new double[] {3.5,0,0};
		double[] e = new double[] {3.5,3.5,0};
		double[] f = new double[] {0,10.0,0};
		translation(d,e,f,box);
		association(pair1,box);
		association2(pair2,box);
		clusterBox.trialNotify();
		clusterBox.acceptNotify();
		System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	public void initializeCoordinates3(Box box) {
		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		associationPotential2.setBox(box);
		associationPotential3.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair1 = new MoleculePair();
		MoleculePair pair2 = new MoleculePair();
		MoleculePair pair3 = new MoleculePair();
		pair1.mol0 = list.get(0);
		pair1.mol1 = list.get(1);
		pair2.mol0 = list.get(1);
		pair2.mol1 = list.get(2);
		pair3.mol0 = list.get(2);
		pair3.mol1 = list.get(3);
		double[] d = new double[] {3.5,0,0};
		double[] e = new double[] {7,0,0};
		double[] f = new double[] {10.5,0,0};
		translation(d,e,f,box);
		association(pair1,box);
		association2(pair2,box);
		association3(pair3,box);
		clusterBox.trialNotify();
		clusterBox.acceptNotify();
		System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	public void initializeCoordinatesRef(Box box) {
		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		associationPotential2.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		for (int i = 1; i<list.size(); i++){
			list.get(i).getChildList().get(0).getPosition().setX(0, 0.9*i);
		 }
		clusterBox.trialNotify();
		clusterBox.acceptNotify();
		System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	public void initializeCoordinates4(Box box) {
		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		associationPotential2.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair1 = new MoleculePair();
		MoleculePair pair2 = new MoleculePair();
		pair1.mol0 = list.get(0);
		pair1.mol1 = list.get(1);
		pair2.mol0 = list.get(2);
		pair2.mol1 = list.get(3);
		double[] d = new double[] {3.5,0,0};
		double[] e = new double[] {3.5,10.0,0};
		double[] f = new double[] {0,10.0,0};
		translation(d,e,f,box);
		association(pair1,box);
		association3(pair2,box);
		clusterBox.trialNotify();
		clusterBox.acceptNotify();
		System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	
	public void initializeCoordinates5(Box box) {
		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		associationPotential2.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair1 = new MoleculePair();
		MoleculePair pair2 = new MoleculePair();
		pair1.mol0 = list.get(0);
		pair1.mol1 = list.get(1);
		pair2.mol0 = list.get(2);
		pair2.mol1 = list.get(3);
		double[] d = new double[] {3.5,0,0};
		double[] e = new double[] {3.5,10.0,0};
		double[] f = new double[] {0,10.0,0};
		translation(d,e,f,box);
		association(pair1,box);
		association3(pair2,box);
		clusterBox.trialNotify();
		clusterBox.acceptNotify();
		System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	
	public void translation(double[] d,double[] e,double[] f, Box box){//place molecule1,2,3 at some position
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
		IMolecule mol1 = list.get(1);
		IMolecule mol2 = list.get(2);
		IMolecule mol3 = list.get(3);
        translationA.actionPerformed(mol1);
        translationB.actionPerformed(mol2);
        translationC.actionPerformed(mol3);
	}
	public void association(MoleculePair pair, Box box){
		RotationTensor rotationTensor = space.makeRotationTensor();
		Vector r0 = space.makeVector();
		IMoleculePositionDefinition positionDefinition = new MoleculePositionGeometricCenter(space);
		IMoleculeList list = box.getMoleculeList();
		pair.mol0 = list.get(0);
		IMolecule water = list.get(1);
		pair.mol1 = water;
		
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
		    for (int iChild = 0; iChild<childList.size(); iChild++) {//free rotation until finding association
		        IAtom a = childList.get(iChild);
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
		pair.mol0 = list.get(1);
		IMolecule water = list.get(2);
		pair.mol1 = water;
	
        while (true){
        	Vector positionWater = space.makeVector();
	        positionWater.setRandomInSphere(random);
	        positionWater.TE(4.0);//place water molecule within a sphere with r = 4A
	        positionWater.PE(pair.mol0.getChildList().get(0).getPosition());
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
		    for (int iChild = 0; iChild<childList.size(); iChild++) {//free rotation until finding association
		        IAtom a = childList.get(iChild);
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
		pair.mol0 = list.get(2);
		IMolecule water = list.get(3);
		pair.mol1 = water;
	
        while (true){
        	Vector positionWater = space.makeVector();
	        positionWater.setRandomInSphere(random);
	        positionWater.TE(4.0);//place water molecule within a sphere with r = 8A
	        positionWater.PE(pair.mol0.getChildList().get(0).getPosition());
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
		    for (int iChild = 0; iChild<childList.size(); iChild++) {//free rotation until finding association
		        IAtom a = childList.get(iChild);
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
	
	public void association4(MoleculePair pair,MoleculePair pair2, Box box){
		RotationTensor rotationTensor = space.makeRotationTensor();
		Vector r0 = space.makeVector();
		IMoleculePositionDefinition positionDefinition = new MoleculePositionGeometricCenter(space);
		IMoleculeList list = box.getMoleculeList();
		pair.mol0 = list.get(1);
		IMolecule water = list.get(2);
		pair.mol1 = water;
		pair2.mol0 =list.get(0);
		pair2.mol1 = water;
		
	
        while (true){
        	Vector positionWater = space.makeVector();
	        positionWater.setRandomInSphere(random);
	        positionWater.TE(4.0);//place water molecule within a sphere with r = 8A
	        positionWater.PE(pair.mol0.getChildList().get(0).getPosition());
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
		    for (int iChild = 0; iChild<childList.size(); iChild++) {//free rotation until finding association
		        IAtom a = childList.get(iChild);
		        Vector r = a.getPosition();
		        r.ME(r0);
		        box.getBoundary().nearestImage(r);
		        rotationTensor.transform(r);
		        r.PE(r0);
	    	}
	        if (associationPotential2.energy(pair) != 0.0 &&nonAssociationPotential.energy(pair2) != 0.0){
        		break;
        	}
		}
	}
	public void association5(MoleculePair pair,MoleculePair pair2, Box box){
		RotationTensor rotationTensor = space.makeRotationTensor();
		Vector r0 = space.makeVector();
		IMoleculePositionDefinition positionDefinition = new MoleculePositionGeometricCenter(space);
		IMoleculeList list = box.getMoleculeList();
		pair.mol0 = list.get(1);
		IMolecule water = list.get(2);
		pair.mol1 = water;
		pair2.mol0 =list.get(0);
		pair2.mol1 = water;
	
        while (true){
        	Vector positionWater = space.makeVector();
	        positionWater.setRandomInSphere(random);
	        positionWater.TE(3.5);//place water molecule within a sphere with r = 4A
	        //positionWater.PE(pair.atom0.getChildList().getAtom(0).getPosition());
	        MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
	        translation.setDestination(positionWater);
	        translation.actionPerformed(water);

	       
	        if (associationPotential2.energy(pair) != 0.0 && associationPotential2.energy(pair2) != 0.0 ){//when there is an association, fix the position of water
        		break;
        	}
	        double dTheta = (2*random.nextDouble() - 1.0)*Math.PI;
	        rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

	        r0.E(positionDefinition.position(water));
		    IAtomList childList = water.getChildList();
		    for (int iChild = 0; iChild<childList.size(); iChild++) {//free rotation until finding association
		        IAtom a = childList.get(iChild);
		        Vector r = a.getPosition();
		        r.ME(r0);
		        box.getBoundary().nearestImage(r);
		        rotationTensor.transform(r);
		        r.PE(r0);
	    	}
	        if (associationPotential2.energy(pair) != 0.0&& associationPotential2.energy(pair2) != 0.0 ){
        		break;
        	}
		}
	}
}
