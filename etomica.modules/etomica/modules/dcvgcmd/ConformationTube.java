package etomica.modules.dcvgcmd;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.config.Conformation;
import etomica.space.Space;

/**
 * @author nsives, msellers
 * 
 * Sets carbon atom coordinates in tubular form.
 *  
 */

public class ConformationTube extends Conformation {

	int atomsPerRing;

	double tubeRadius;

	double A, B, Q;

	double dtheta;

	public ConformationTube(Space space, int atomsPerRing) {
		
		super(space);	
		this.atomsPerRing = atomsPerRing;

		dtheta = Math.PI * 2.0 / atomsPerRing;
		
		
		A= 1.42;

		B= 0.705;
		
		//Distance between ORTHO carbons in hexagonal ring, x-z plane
		Q = 2.45;

		//Trig operation on right triangle formed, x-y plane
		tubeRadius = Q / 2.0/Math.sin(dtheta / 2);
		
	}	

	public void initializePositions(AtomArrayList atomList) {

		int size = atomList.size();
		if (size == 0)
			return;

		AtomIteratorArrayListSimple atomIterator = new AtomIteratorArrayListSimple();
		atomIterator.setList(atomList);
		atomIterator.reset();

		int N = 0;
		double theta, theta0, x, y, z, dz;
		int ctr;
		
//		tubeRadius = atomsPerRing*.5;
		x = tubeRadius;
		y = 0;
		z = 0;
		dz = 0;
		theta = 0;
		theta0 = 0;
		ctr = 0;

		while (atomIterator.hasNext()) {
			AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();

			a.coord.position().setX(0, x);
			a.coord.position().setX(1, y);
			a.coord.position().setX(2, z);

			N = N + 1;

			if (N == atomsPerRing) {
				N = 0;
				switch (ctr) {
				case (0):
					dz = A;
					theta0 = 0;
					ctr = 1;
					break;

				case (1):
					dz = B;
					theta0 = dtheta / 2;
					ctr = 2;
					break;

				case (2):
					dz = A;
					theta0 = dtheta / 2;
					ctr = 3;
					break;

				case (3):
					dz = B;
					theta0 = 0;
					ctr = 0;
					break;
				}
				z = z + dz;
			}

			theta = theta0 + N * dtheta;

			x = tubeRadius * Math.cos(theta);
			y = tubeRadius * Math.sin(theta);
		}
	}

}