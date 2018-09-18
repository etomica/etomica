/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.atom.IAtomOriented;
import etomica.space.Space;
import etomica.space3d.IOrientationFull3D;

public class ConfigurationClusterChainFourSites extends ConfigurationCluster {

	public ConfigurationClusterChainFourSites(Space _space) {
		super(_space);
	}

	public void initializeCoordinates(Box box) {
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
        Vector direction = space.makeVector();
        Vector secondaryDirection = space.makeVector();
		for (int i = 1; i<list.size(); i++){
			list.get(i).getPosition().setX(0, 0.9*i);
			if (list.get(i) instanceof IAtomOriented){
				double x = 1.0/3.0;
				double y = -2.0/(3.0*Math.sqrt(2.0));
				double z = Math.sqrt(2.0/3.0);
				direction.E(new double[]{x,y,z});
				double norm = Math.sqrt(16/9.0*(x*x+y*y+0.25*z*z));
				secondaryDirection.E(new double[]{4.0*x/3.0, 4.0*y/3.0, -2.0*z/3.0});
				secondaryDirection.TE(1.0/norm);
				((IOrientationFull3D)((IAtomOriented)list.get(i)).getOrientation()).setDirections(direction, secondaryDirection);
			}
		 }
		 clusterBox.trialNotify();
		 clusterBox.acceptNotify();
		 System.out.println("box "+clusterBox.getSampleCluster().value(clusterBox));
	}
}
