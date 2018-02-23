/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.action.MoleculeActionTranslateTo;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.*;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class ConfigurationClusterAceticAcid extends ConfigurationCluster {
	 
	protected IRandom random;
	protected MayerFunction f;
	protected MayerFunction f2;
	protected MayerFunction f3;
	
	public ConfigurationClusterAceticAcid(Space _space, IRandom random, MayerFunction f) {
		super(_space);
		this.f = f;
		this.random = random;
	}
	public ConfigurationClusterAceticAcid(Space _space, IRandom random, MayerFunction f, MayerFunction f2) {
		super(_space);
		this.f = f;
		this.f2 = f2;
		this.random = random;
	}
	public ConfigurationClusterAceticAcid(Space _space, IRandom random, MayerFunction f, MayerFunction f2, MayerFunction f3) {
		super(_space);
		this.f = f;
		this.f2 = f2;
		this.f3 = f3;
		this.random = random;
	}

	public void initializeCoordinates(Box box, boolean a, boolean b) {
		super.initializeCoordinates(box);
		f.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair1 = new MoleculePair();
		MoleculePair pair2 = new MoleculePair();
		MoleculePair pair3 = new MoleculePair();
		pair1.atom0 = list.get(0);
		pair1.atom1 = list.get(1);
		association(f,pair1,box);
		if (list.size()==3){
			pair2.atom0 = list.get(1);
			pair2.atom1 = list.get(2);
			double[] e = new double[] {1.0,1.0,0};
			translation2Mol(2, e, box);
			if(a){
				association(f,pair2,box);
			}
		}
		if (list.size()==4){
			pair2.atom0 = list.get(1);
			pair2.atom1 = list.get(2);
			pair3.atom0 = list.get(2);
			pair3.atom1 = list.get(3);
			double[] e = new double[] {1.0,1.0,0};
			double[] g = new double[] {0,1.0,0};
			translation3Mol(e,g,box);
			if(a){
				association(f,pair2,box);
			}
			if (b){
				association(f,pair3,box); 
			}
		}
		clusterBox.trialNotify();
		clusterBox.acceptNotify();
		System.out.println(clusterBox+" "+clusterBox.getSampleCluster().value(clusterBox));
	}
	 
	public void initializeCoordinates4Mol(Box box, char a, char b, char c, char d, char e, int molNumber, int molNumber2) {
		super.initializeCoordinates(box);
		f.setBox(box);
		f2.setBox(box);
		f3.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair1 = new MoleculePair();
		MoleculePair pair2 = new MoleculePair();
		MoleculePair pair3 = new MoleculePair();
		MoleculePair pair4 = new MoleculePair();
		MoleculePair pair5 = new MoleculePair();
		pair1.atom0 = list.get(0);
		pair1.atom1 = list.get(1);
		pair2.atom0 = list.get(1);
		pair2.atom1 = list.get(2);
		pair4.atom0 = list.get(0);
		pair4.atom1 = list.get(2);
		  
		if (list.size() >3){
			pair3.atom0 = list.get(0);
			pair3.atom1 = list.get(3);
			pair5.atom0 = list.get(2);
			pair5.atom1 = list.get(3);
		}
		if (molNumber2 >-1){
			double[] amount = new double[] {5.0,5.0,0};
			translation2Mol(molNumber2, amount, box);
		}
		
		if (a == 'E'){
			association(f,pair1,box);
		}
		if (a == 'C'){
			association(f2,pair1,box);
		}
		if (a == 'H'){
			association(f3,pair1,box);
		}
		if (b == 'E'){
			association(f,pair2,box);
		}
		if (b == 'C'){
			association(f2,pair2,box);
		}
		if (b == 'H'){
			association(f3,pair2,box);
		}
		if (c == 'E'){
			association(f,pair3,box);
		}
		if (c == 'C'){
			association(f2,pair3,box);
		}
		if (c == 'H'){
			association(f3,pair3,box);
		}
		if (d == 'E'){
			association(f,pair4,box);
		}
		if (d == 'C'){
			association(f2,pair4,box);
		}
		if (d == 'H'){
			association(f3,pair4,box);
		}
		if (e == 'E'){
			association(f,pair5,box);
		}
		if (e == 'C'){
			association(f2,pair5,box);
		}
		if (e == 'H'){
			association(f3,pair5,box);
		}
		if (molNumber > -1){
			translationRandom(list.get(molNumber), box);
		}
		clusterBox.trialNotify();
		clusterBox.acceptNotify();
		System.out.println(clusterBox+" "+clusterBox.getSampleCluster().value(clusterBox));
	}
	 
	public void initializeCoordinates2(Box box, boolean a, boolean b) {
		super.initializeCoordinates(box);
		f.setBox(box);
		f2.setBox(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IMoleculeList list = box.getMoleculeList();
		MoleculePair pair1 = new MoleculePair();
		MoleculePair pair2 = new MoleculePair();
		MoleculePair pair3 = new MoleculePair();
		pair1.atom0 = list.get(0);
		pair1.atom1 = list.get(1);
		pair2.atom0 = list.get(1);
		pair2.atom1 = list.get(2);
		pair3.atom0 = list.get(2);
		pair3.atom1 = list.get(3);
		double[] e = new double[] {1.0,1.0,0};
		double[] g = new double[] {0,1.0,0};
		translation3Mol(e,g,box);
		association(f,pair1,box);
		association(f2,pair3,box); 
		clusterBox.trialNotify();
		clusterBox.acceptNotify();
		System.out.println(clusterBox+" "+clusterBox.getSampleCluster().value(clusterBox));
	}

	public void translation2Mol(int a, double[] e, Box box){//place molecule a at some position
		IMoleculeList list = box.getMoleculeList();
		MoleculeActionTranslateTo translationB = new MoleculeActionTranslateTo(space);
		Vector b = space.makeVector();
		b.E(e);
		translationB.setDestination(b);
		IMolecule mol = list.get(a);
		translationB.actionPerformed(mol);
	}
	 
	public void translation3Mol(double[] e,double[] f, Box box){//place molecule2,3 at some position
		IMoleculeList list = box.getMoleculeList(); 
		MoleculeActionTranslateTo translationB = new MoleculeActionTranslateTo(space);
		MoleculeActionTranslateTo translationC = new MoleculeActionTranslateTo(space);
		Vector b = space.makeVector();
		Vector c = space.makeVector();
		b.E(e);
		c.E(f);
		translationB.setDestination(b);
		translationC.setDestination(c);
		IMolecule mol2 = list.get(2);
		IMolecule mol3 = list.get(3);
		translationB.actionPerformed(mol2);
		translationC.actionPerformed(mol3);
	}
	
	public void translationRandom(IMolecule mol, Box box){
		BoxCluster clusterBox =(BoxCluster) box; 
		while (true){
			Vector positionAceticAcid = space.makeVector();
			positionAceticAcid.setRandomInSphere(random);
			positionAceticAcid.TE(8.0);//place acetic acid molecule within a sphere with r = 8A
			MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
			translation.setDestination(positionAceticAcid);
			translation.actionPerformed(mol);
			clusterBox.trialNotify();
			clusterBox.acceptNotify();
			if (clusterBox.getSampleCluster().value(clusterBox) != 0.0){//when there is an association, fix position of the molecule
				break;
			}
		 }
	}
	
	public void association(MayerFunction f, MoleculePair pair, Box box){
		RotationTensor rotationTensor = space.makeRotationTensor();
		Vector r0 = space.makeVector();
		IMoleculePositionDefinition positionDefinition = new MoleculePositionGeometricCenter(space);
		while (true){
			Vector positionAceticAcid = space.makeVector();
			positionAceticAcid.setRandomInSphere(random);
			positionAceticAcid.TE(8.0);//place acetic acid molecule within a sphere with r = 8A
			positionAceticAcid.PE(positionDefinition.position(pair.atom0));
			MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
			translation.setDestination(positionAceticAcid);
			translation.actionPerformed(pair.atom1);
			
			if (f.f(pair, 0, 0.001) != 0.0){//when there is an association, fix position of the molecule
				f.f(pair, 0, 0.001);
				break;
			}
	         double dTheta = (2*random.nextDouble() - 1.0)*Math.PI;
	         rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

	         r0.E(positionDefinition.position(pair.atom1));
	         IAtomList childList = pair.atom1.getChildList();
	         for (int iChild = 0; iChild<childList.size(); iChild++) {//free rotation until finding association
	        	 IAtom a = childList.get(iChild);
	        	 Vector r = a.getPosition();
	        	 r.ME(r0);
	        	 box.getBoundary().nearestImage(r);
	        	 rotationTensor.transform(r);
	        	 r.PE(r0);
	         }
	         if (f.f(pair, 0, 0.001) != 0.0){
	        	 f.f(pair, 0, 0.001);
	        	 break;
	         }
		 }
	 }
}
