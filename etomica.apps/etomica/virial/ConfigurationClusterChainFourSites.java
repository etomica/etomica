package etomica.virial;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.atom.IAtomOriented;
import etomica.space.ISpace;
import etomica.space3d.IOrientationFull3D;

public class ConfigurationClusterChainFourSites extends ConfigurationCluster {

	public ConfigurationClusterChainFourSites(ISpace _space) {
		super(_space);
	}

	public void initializeCoordinates(IBox box) {
		double angle1 = Math.acos(-1.0/3.0);//theta
    	double cosAngle1 = -1.0/3.0;
    	double sinAngle1 = Math.sin(angle1);
    	//double angle2 = (Math.PI-angle1)/2.0;//(pi-theta)/2
    	//double cosAngle2 = Math.cos(angle2);
    	//double length = 2*cosAngle2;//length of tetrahedron
    	double coordX = cosAngle1;
        double coordY = cosAngle1*(1-cosAngle1)/sinAngle1;
        //double coordZ = Math.sqrt(1-coordX*coordX-coordY*coordY);
        double rotateAngle1 = 2*Math.PI-Math.acos(Math.abs(coordY)/Math.sqrt((1-coordX*coordX)));
        System.out.println("rotateAngle1= "+rotateAngle1);
        double rotateAngle2 = 2*Math.PI-Math.acos(Math.abs(coordX));
        System.out.println("rotateAngle2= "+rotateAngle2);
        //double rotateAngle2 = Math.PI-angle1;
		super.initializeCoordinates(box);
		BoxCluster clusterBox =(BoxCluster) box;
		IAtomList list = box.getLeafList();
		for (int i=1;i<list.getAtomCount();i++){
			((IVectorMutable)list.getAtom(i).getPosition()).setX(0, 0.9*i);
			if (list.getAtom(i) instanceof IAtomOriented){
				IVectorMutable xDirection = space.makeVector();
				IVectorMutable zDirection = space.makeVector();
				xDirection.setX(0,1);//X axis=0, Y=1, Z=2
				zDirection.setX(2,1);//Z
				 ((IOrientationFull3D)((IAtomOriented)list.getAtom(i)).getOrientation()).rotateBy(rotateAngle1*i, xDirection);
				 ((IOrientationFull3D)((IAtomOriented)list.getAtom(i)).getOrientation()).rotateBy(rotateAngle2*i, zDirection);
			}
		 }
		 clusterBox.trialNotify();
		 clusterBox.acceptNotify();
		 System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
}
