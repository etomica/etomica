package etomica.virial;

import etomica.action.MoleculeActionTranslateTo;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.IAtomPositionDefinition;
import etomica.atom.MoleculePair;
import etomica.models.OPLS.P2AceticAcidOneSite;
import etomica.models.water.PNWaterGCPMThreeSite;
import etomica.models.water.SpeciesWater4P;
import etomica.potential.PotentialGroup;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;
import etomica.space.RotationTensor;

public class ConfigurationClusterAceticAcid extends ConfigurationCluster {
	 
	 protected IRandom random;
	 protected MayerFunction f;
	 protected MayerFunction f2;

	 public ConfigurationClusterAceticAcid(ISpace _space, IRandom random, MayerFunction f) {
	  super(_space);
	  this.f = f;
	  this.random = random;
	 }
	 public ConfigurationClusterAceticAcid(ISpace _space, IRandom random, MayerFunction f, MayerFunction f2) {
	  super(_space);
	  this.f = f;
	  this.f2 = f2;
	  this.random = random;
	 }

	 public void initializeCoordinates(IBox box,boolean a,boolean b) {
	  super.initializeCoordinates(box);
	  f.setBox(box);
	  BoxCluster clusterBox =(BoxCluster) box;
	  IMoleculeList list = box.getMoleculeList();
	  MoleculePair pair1 = new MoleculePair();
	  MoleculePair pair2 = new MoleculePair();
	  MoleculePair pair3 = new MoleculePair();
	  pair1.atom0 = list.getMolecule(0);
	  pair1.atom1 = list.getMolecule(1);
	  association(f,pair1,box);
	  if (list.getMoleculeCount()==3){
	   pair2.atom0 = list.getMolecule(1);
	   pair2.atom1 = list.getMolecule(2);
	   double[] e = new double[] {1.0,1.0,0};
	   translation2Mol(e,box);
	         if(a){
	          association(f,pair2,box);
	         }
	  }
	  if (list.getMoleculeCount()==4){
	   pair2.atom0 = list.getMolecule(1);
	   pair2.atom1 = list.getMolecule(2);
	   pair3.atom0 = list.getMolecule(2);
	   pair3.atom1 = list.getMolecule(3);
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
	 
	 public void initializeCoordinates2(IBox box,boolean a,boolean b) {
	  super.initializeCoordinates(box);
	  f.setBox(box);
	  f2.setBox(box);
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
	  double[] e = new double[] {1.0,1.0,0};
	  double[] g = new double[] {0,1.0,0};
	  translation3Mol(e,g,box);
	  association(f,pair1,box);
	  association(f2,pair3,box); 
	   clusterBox.trialNotify();
	   clusterBox.acceptNotify();
	   System.out.println(clusterBox+" "+clusterBox.getSampleCluster().value(clusterBox));
	 }


	 public void translation2Mol(double[] e, IBox box){//place molecule1,2 at some position
	  IMoleculeList list = box.getMoleculeList();
	        MoleculeActionTranslateTo translationB = new MoleculeActionTranslateTo(space);
	        IVectorMutable b = space.makeVector();
	        b.E(e);
	        translationB.setDestination(b);
	  IMolecule mol2 = list.getMolecule(2);
	        translationB.actionPerformed(mol2);
	 }
	 
	 public void translation3Mol(double[] e,double[] f, IBox box){//place molecule1,2 at some position
	  IMoleculeList list = box.getMoleculeList();
	       
	        MoleculeActionTranslateTo translationB = new MoleculeActionTranslateTo(space);
	        MoleculeActionTranslateTo translationC = new MoleculeActionTranslateTo(space);
	        IVectorMutable b = space.makeVector();
	        IVectorMutable c = space.makeVector();
	        b.E(e);
	        c.E(f);
	        translationB.setDestination(b);
	        translationC.setDestination(c);
	  IMolecule mol2 = list.getMolecule(2);
	  IMolecule mol3 = list.getMolecule(3);
	        translationB.actionPerformed(mol2);
	        translationC.actionPerformed(mol3);
	 }
	 public void association(MayerFunction f, MoleculePair pair, IBox box){
	  RotationTensor rotationTensor = space.makeRotationTensor();
	  IVectorMutable r0 = space.makeVector();
	  IAtomPositionDefinition positionDefinition = new AtomPositionGeometricCenter(space); 
	        while (true){
	         IVectorRandom positionWater = (IVectorRandom)space.makeVector();
	         positionWater.setRandomInSphere(random);
	         positionWater.TE(8.0);//place water molecule within a sphere with r = 8A
	         positionWater.PE(positionDefinition.position(pair.atom0));
	         MoleculeActionTranslateTo translation = new MoleculeActionTranslateTo(space);
	         translation.setDestination(positionWater);
	         translation.actionPerformed(pair.atom1);
	       
	         if (f.f(pair, 0, 0.001) != 0.0){//when there is an association, fix the position of water
	          break;
	         }
	         double dTheta = (2*random.nextDouble() - 1.0)*Math.PI;
	         rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2,dTheta);

	         r0.E(positionDefinition.position(pair.atom1));
	      IAtomList childList = pair.atom1.getChildList();
	      for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {//free rotation until finding association
	          IAtom a = childList.getAtom(iChild);
	          IVectorMutable r = a.getPosition();
	          r.ME(r0);
	          box.getBoundary().nearestImage(r);
	          rotationTensor.transform(r);
	          r.PE(r0);
	      }
	         if (f.f(pair, 0, 0.001) != 0.0){
	          break;
	         }
	  }
	 }
	}