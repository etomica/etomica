/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.action.MoleculeActionTranslateTo;
import etomica.api.*;
import etomica.atom.*;
import etomica.box.Box;
import etomica.models.water.PNWaterGCPMThreeSite;
import etomica.models.water.SpeciesWater4P;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space.RotationTensor;
import etomica.util.random.IRandom;

public class ConfigurationClusterWertheimGCPM extends ConfigurationCluster {
	
	protected IRandom random;
	protected PNWaterGCPMThreeSite associationPotential;
	protected PNWaterGCPMThreeSite associationPotential2;
	protected PNWaterGCPMThreeSite associationPotential3;
	protected PNWaterGCPMThreeSite nonAssociationPotential;
	protected double diagramIndex;

	public ConfigurationClusterWertheimGCPM(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential) {
		super(_space);
		this.associationPotential = associationPotential;
		this.random = random;
	}

	public ConfigurationClusterWertheimGCPM(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential, PNWaterGCPMThreeSite associationPotential2 ) {
		super(_space);
		this.associationPotential = associationPotential;
		this.associationPotential2 = associationPotential2;
		this.random = random;
	}
	public ConfigurationClusterWertheimGCPM(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential, PNWaterGCPMThreeSite associationPotential2, PNWaterGCPMThreeSite associationPotential3 ) {
		super(_space);
		this.associationPotential = associationPotential;
		this.associationPotential2 = associationPotential2;
		this.associationPotential3 = associationPotential3;
		this.random = random;
	}
//	public ConfigurationClusterWertheimGCPM(ISpace _space, IRandom random, PNWaterGCPMThreeSite associationPotential,PNWaterGCPMThreeSite associationPotential2, PNWaterGCPMThreeSite nonAssociationPotential ) {
//		super(_space);
//		this.associationPotential = associationPotential;
//		this.associationPotential2 = associationPotential2;
//		this.nonAssociationPotential = nonAssociationPotential;
//		this.random = random;
//	}
	
	public ConfigurationClusterWertheimGCPM(Space _space, IRandom random, PNWaterGCPMThreeSite associationPotential, PNWaterGCPMThreeSite associationPotential2, double diagramIndex ) {
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
		pair.atom0 = list.getMolecule(0);
		pair.atom1 = list.getMolecule(1);
		if (list.getMoleculeCount()==3){
			IMolecule mol3 = list.getMolecule(2);
			Vector v1 = space.makeVector();
			v1.E(10.0);
	        MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
	        translation.setDestination(v1);
	        translation.actionPerformed(mol3);
		}
		if (list.getMoleculeCount()==4 && diagramIndex == 0){
			IMolecule mol3 = list.getMolecule(2);
			IMolecule mol4 = list.getMolecule(3);
			Vector v1 = space.makeVector();
			v1.E(10.0);
	        MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
	        translation.setDestination(v1);
	        translation.actionPerformed(mol3);
	        translation.actionPerformed(mol4);
		} else if (list.getMoleculeCount()==4 && diagramIndex == 1){
			IMolecule mol3 = list.getMolecule(2);
			IMolecule mol4 = list.getMolecule(3);
			Vector v1 = space.makeVector();
			v1.E(10.0);
	        MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
	        translation.setDestination(v1);
	        translation.actionPerformed(mol3);
	        translation.actionPerformed(mol4);
		}

		association(pair,box);

//		for (int i=1;i<list.getMoleculeCount();i++){
//			pair.atom0 = list.getMolecule(i-1);
//	        IMolecule water = list.getMolecule(i);
//	        pair.atom1 = water;
//
//		 }
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
		translation(d,e,box);
		clusterBox.trialNotify();
		clusterBox.acceptNotify();
		System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	
	public void initializeCoordinates2(Box box) {

		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair1 = new MoleculePair();
		MoleculePair pair2 = new MoleculePair();
		pair1.atom0 = list.getMolecule(0);
		pair1.atom1 = list.getMolecule(1);
		pair2.atom0 = list.getMolecule(1);
		pair2.atom1 = list.getMolecule(2);
		
		association(pair1,box);
		association2(pair2,box);
		
		Vector v1 = space.makeVector();
		Vector v2 = space.makeVector();
		Vector v3 = space.makeVector();
		
		v1.E(box.getMoleculeList().getMolecule(1).getChildList().getAtom(SpeciesWater4P.indexO).getPosition());
		v1.ME(box.getMoleculeList().getMolecule(0).getChildList().getAtom(SpeciesWater4P.indexO).getPosition());
		double distance12 = Math.sqrt(v1.squared());
		
		v2.E(box.getMoleculeList().getMolecule(2).getChildList().getAtom(SpeciesWater4P.indexO).getPosition());
		v2.ME(box.getMoleculeList().getMolecule(1).getChildList().getAtom(SpeciesWater4P.indexO).getPosition());
		double distance23 = Math.sqrt(v2.squared());
		
		v3.E(box.getMoleculeList().getMolecule(2).getChildList().getAtom(SpeciesWater4P.indexO).getPosition());
		v3.ME(box.getMoleculeList().getMolecule(0).getChildList().getAtom(SpeciesWater4P.indexO).getPosition());
		double distance13 = Math.sqrt(v3.squared());
		
		System.out.println("distance12 "+distance12+" distance23 "+distance23+" distance13 "+distance13);

		 clusterBox.trialNotify();
		 clusterBox.acceptNotify();
		 System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	public void initializeCoordinates3(Box box) {

		super.initializeCoordinates(box);
		associationPotential.setBox(box);
		associationPotential2.setBox(box);
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
		double[] e = new double[] {7,0.0,0};
		translation(d,e,box);
		
		association(pair1,box);
		association3(pair2,box);
		
		IMolecule mol4 = list.getMolecule(3);
		
		if (list.getMoleculeCount()==4&& diagramIndex == 1){
			Vector v4 = space.makeVector();
			v4.E(10.0);
	        MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
	        translation.setDestination(v4);
	        translation.actionPerformed(mol4);
		}
		
		if (list.getMoleculeCount()==4 && diagramIndex == 2){
			association3(pair3,box);
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
		MoleculePair pair3 = new MoleculePair();
		pair1.atom0 = list.getMolecule(0);
		pair1.atom1 = list.getMolecule(1);
		pair2.atom0 = list.getMolecule(1);
		pair2.atom1 = list.getMolecule(2);
		pair3.atom0 = list.getMolecule(0);
		pair3.atom1 = list.getMolecule(2);
		
		association(pair1,box);
		association5(pair2,pair3,box);

		 clusterBox.trialNotify();
		 clusterBox.acceptNotify();
		 System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	public void initializeCoordinates5(Box box) {

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
		pair3.atom0 = list.getMolecule(0);
		pair3.atom1 = list.getMolecule(2);
		
		double[] d = new double[] {1.5,2.5,0};
		double[] e = new double[] {3.0,0.0,0};
		translation(d,e,box);
		
		association(pair1,box);
		association3(pair2,box);
		association6(pair2,pair3,box);

		 clusterBox.trialNotify();
		 clusterBox.acceptNotify();
		 System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
	public void translation(double[] d,double[] e, Box box){//place molecule1,2 at some position
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
	        positionWater.TE(4.0);//place water molecule within a sphere with r = 4A
	        positionWater.PE(pair.atom0.getChildList().getAtom(0).getPosition());
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
	
	public void association3(MoleculePair pair, Box box){
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
	
	public void association4(MoleculePair pair,MoleculePair pair2, Box box){
		RotationTensor rotationTensor = space.makeRotationTensor();
		Vector r0 = space.makeVector();
		IMoleculePositionDefinition positionDefinition = new MoleculePositionGeometricCenter(space);
		IMoleculeList list = box.getMoleculeList();
		pair.atom0 = list.getMolecule(1);
		IMolecule water = list.getMolecule(2);
		pair.atom1 = water;
		pair2.atom0 =list.getMolecule(0);
		pair2.atom1 = water;
		
	
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
		pair.atom0 = list.getMolecule(1);
		IMolecule water = list.getMolecule(2);
		pair.atom1 = water;
		pair2.atom0 =list.getMolecule(0);
		pair2.atom1 = water;
	
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
		    for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {//free rotation until finding association
		        IAtom a = childList.getAtom(iChild);
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
	public void association6(MoleculePair pair,MoleculePair pair2, Box box){
		RotationTensor rotationTensor = space.makeRotationTensor();
		Vector r0 = space.makeVector();
		IMoleculePositionDefinition positionDefinition = new MoleculePositionGeometricCenter(space);
		IMoleculeList list = box.getMoleculeList();
		pair.atom0 = list.getMolecule(1);
		IMolecule water = list.getMolecule(2);
		pair.atom1 = water;
		pair2.atom0 =list.getMolecule(0);
		pair2.atom1 = water;
	
        while (true){
        	Vector positionWater = space.makeVector();
	        positionWater.setRandomInSphere(random);
	        positionWater.TE(3.5);//place water molecule within a sphere with r = 4A
	        //positionWater.PE(pair.atom0.getChildList().getAtom(0).getPosition());
	        MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
	        translation.setDestination(positionWater);
	        translation.actionPerformed(water);

	       
	        if (associationPotential2.energy(pair) != 0.0 && associationPotential3.energy(pair2) != 0.0 ){//when there is an association, fix the position of water
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
	        if (associationPotential2.energy(pair) != 0.0&& associationPotential3.energy(pair2) != 0.0 ){
        		break;
        	}
		}
	}
}
