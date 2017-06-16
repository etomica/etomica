/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.AtomHydrogen;
import etomica.atom.AtomTypeOriented;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.chem.elements.Hydrogen;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.math.Quaternion;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.IOrientation3D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.random.IRandom;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * MCMove that fully regrows the beads of a ring polymer by rotating the images, accepting or
 * rejecting the move based on the sampling weight.
 *
 * @author Ram
 */
public class MCMoveClusterRingRegrowOrientation extends MCMoveBox {

	// Private instance variables
	private final int P;
	private final Space space;
	private final IRandom random;
	private final AtomIteratorLeafAtoms leafIterator;
	private final Vector utilityVec1, utilityVec2, utilityVec3;
	private IOrientation3D[][] oldOrientations;
	private double weightOld, weightNew,kHarmonic,pacc;
	private int acc = 0, molIndexUntouched = -1;
	private boolean firstMove = true;
	private boolean[] doExchange;
	private double [] newAlpha;

	public MCMoveClusterRingRegrowOrientation(IRandom random, Space _space, int P) {
		super(null);
		this.space = _space;
		this.P = P;
		this.random = random;
		leafIterator = new AtomIteratorLeafAtoms();
		utilityVec1 = space.makeVector();
		utilityVec2 = space.makeVector();
		utilityVec3 = space.makeVector();
	}

	public static void main(String[] args) {
		Space space = Space3D.getInstance();
		ClusterWeight cluster = new ClusterWeight() {

			@Override
			public double value(BoxCluster box) {
				// TODO Auto-generated method stub
				return 1;
			}

			@Override
			public void setTemperature(double temperature) {
				// TODO Auto-generated method stub

			}

			@Override
			public int pointCount() {
				// TODO Auto-generated method stub
				return 1;
			}

			@Override
			public ClusterAbstract makeCopy() {
				// TODO Auto-generated method stub
				return null;
			}
		};
		BoxCluster box = new BoxCluster(cluster, space);
		Simulation sim = new Simulation(space);
		sim.addBox(box);
		AtomTypeOriented atype = new AtomTypeOriented(Hydrogen.INSTANCE, space);
		SpeciesSpheresHetero species = new SpeciesSpheresHetero(space, new AtomTypeOriented[]{atype});
		sim.addSpecies(species);
		File file1 = new File("acceptance.dat");
		if (file1.exists()) {
			file1.delete();
		}
		for (int p = 2; p <= 512; p *= 2) {
			box.setNMolecules(species, 0);
			species.setChildCount(new int[]{p});
			box.setNMolecules(species, 1);
			IntegratorMC integrator = new IntegratorMC(sim, null);
			integrator.setBox(box);
			MCMoveClusterRingRegrowOrientation move = new MCMoveClusterRingRegrowOrientation(sim.getRandom(), space, p);

			for (int iTemp = 40; iTemp <= 40; iTemp += 2) {
				move.acc = 0;
				move.setStiffness(Kelvin.UNIT.toSim(iTemp), species.getAtomType(0).getMass());
				integrator.getMoveManager().addMCMove(move);
				integrator.reset();
				int total = 100;
				for (int i = 0; i < total; i++) {
					integrator.doStep();
				}
				try {
					FileWriter Temp = new FileWriter("acceptance.dat", true);
					Temp.write(iTemp + " " + p + " " + move.getStiffness() + " " + ((double) move.acc) / total + "\n");
					Temp.close();
				} catch (IOException ex1) {
					throw new RuntimeException(ex1);
				}
				System.out.println("p = " + p + " ,Temp = " + iTemp + " ,acceptance ratio = " + ((double) move.acc) / total);
			}

		}
	}

	@Override
	public void setBox(Box p) {
		super.setBox(p);
		int nMolecules = box.getMoleculeList().getMoleculeCount();
		doExchange = new boolean[nMolecules];
		oldOrientations = new IOrientation3D[nMolecules][];
		for (int i=0; i<nMolecules; i++) {
			doExchange[i] = false;
			int nAtoms = box.getMoleculeList().getMolecule(i).getChildList().getAtomCount();
			oldOrientations[i] = new IOrientation3D[nAtoms+1];
			for (int j=0; j<nAtoms+1; j++) {
				oldOrientations[i][j] = (IOrientation3D) space.makeOrientation();
			}
		}
		leafIterator.setBox(p);
	}

	/**
	 * Set PI harmonic spring stiffness for given temperature and atomic mass.
	 * Dimer is assumed to be composed of two atoms of the given mass.
	 */
	public void setStiffness(double t, double mass) {
		double lambda = Constants.PLANCK_H/(Math.sqrt(2*Math.PI*mass*t));
		kHarmonic = Math.PI*P/(lambda*lambda);
	}

	public double getStiffness() {
		return kHarmonic;
	}

	@Override
	public boolean doTrial() {
		weightOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
		Vector axis = space.makeVector();
		Vector oldCenter = space.makeVector();
		Vector newCenter = space.makeVector();
		IMoleculeList molecules = box.getMoleculeList();
		IOrientation3D [][] newOrientations = new IOrientation3D[molecules.getMoleculeCount()][P+1];
		double [] oldAlpha = new double [P];
		newAlpha = new double [P];
		double [] theta = new double [P];
		int fromImage = 0;
		int toImage = 0;
		double uOld = 0;
		double uNew = 0;
		double pGenRatio = 1.00;
		Vector pVecOld = space.makeVector();
		Vector pVecNew = space.makeVector();
		int nMolecules = molecules.getMoleculeCount();
		molIndexUntouched = random.nextInt(nMolecules);
		for (int i=0; i<nMolecules; i++) {
			if (molIndexUntouched == i) continue;
			IMolecule molecule = molecules.getMolecule(i);
			IAtomList atoms = molecule.getChildList();
			for (int j=0; j<P; j++) {
				int prev = j-1;
				if (prev < 0) prev = P-1;
				AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
				AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
				uOld += kHarmonic*dist(jAtom,jPrev, i);
			}
			oldOrientations[i][0].setDirection(((IAtomOriented) atoms.getAtom(0)).getOrientation().getDirection());
			Vector rV1 = space.makeVector();
			rV1.setRandomSphere(random);
			newOrientations[i][0] = (IOrientation3D)((IAtomOriented) atoms.getAtom(0)).getOrientation();
			newOrientations[i][0].setDirection(rV1);
			oldOrientations[i][P].setDirection(oldOrientations[i][0].getDirection());
			newOrientations[i][P] = (IOrientation3D) space.makeOrientation();
			newOrientations[i][P].setDirection(newOrientations[i][0].getDirection());
			pVecOld.E(oldOrientations[i][0].getDirection());
			pVecNew.E(newOrientations[i][0].getDirection());
			for (int dr = 2; dr<= P; dr*=2){
				double kEff = 4*kHarmonic*dr/P;
				double y0 = 0;
				double kEff1Old = 0;
				double kEff1New = 0;
				double pGenNew = 0;
				double y1 = 0;
				for (int nr = 1; nr<dr; nr+=2){
					boolean piFlagNew = false;
					boolean piFlagOld = false;
					int imageIndex = nr*P/dr;
					//                    System.out.println("image # = "+imageIndex);
					IAtomOriented jAtom = ((IAtomOriented)atoms.getAtom(imageIndex));
					oldOrientations[i][imageIndex].setDirection(jAtom.getOrientation().getDirection());
					newOrientations[i][imageIndex] = (IOrientation3D) jAtom.getOrientation();
					fromImage = (nr-1)*P/dr;
					toImage = (nr+1)*P/dr;
					double bl = 0;
					for (int k = fromImage; k<= toImage; k++) {
						if (k == P) {
							bl += ((AtomHydrogen)atoms.getAtom(0)).getBondLength();
						}
						else {
							bl += ((AtomHydrogen)atoms.getAtom(k)).getBondLength();
						}

					}
					bl *= dr/(2.0*P + dr);
					// same as r /= (toImage - fromImage + 1)


					if (imageIndex == P/2 && doExchange[i]) {
						pVecOld.TE(-1);
						pVecNew.TE(-1);
						oldOrientations[i][P].setDirection(pVecOld);
						newOrientations[i][P].setDirection(pVecNew);
					}
					else {
						oldCenter.Ev1Pv2(oldOrientations[i][fromImage].getDirection(), oldOrientations[i][toImage].getDirection());
						oldCenter.normalize();
						if (oldCenter.isNaN()) throw new RuntimeException("Fix for solving this problem: Set starting orientations in molecule "+i);

						y0 = oldOrientations[i][fromImage].getDirection().dot(oldOrientations[i][toImage].getDirection());
						if (y0 > 1.0) y0 = 1.0;
						if (y0 < -1.0) y0 = -1.0;
						double oldY0 = y0;
						if (y0 == -1.0) piFlagOld = true;
						if (!piFlagOld) kEff1Old = kEff*bl*bl*Math.sqrt((1+y0)/2.0);
						y0 = newOrientations[i][fromImage].getDirection().dot(newOrientations[i][toImage].getDirection());
						if (y0 > 1.0) y0 = 1.0;
						if (y0 < -1.0) y0 = -1.0;
						double newY0 = y0;
						if (y0 == -1.0) piFlagNew = true;
						if (!piFlagNew) kEff1New = kEff*bl*bl*Math.sqrt((1+y0)/2.0);
						y0 = oldCenter.dot(oldOrientations[i][imageIndex].getDirection());
						if (y0 > 1.0) y0 = 1.0;
						if (y0 < -1.0) y0 = -1.0;
						double oldCenterY0 = y0;
						oldAlpha[imageIndex] = Math.acos(y0);
						double x = random.nextDouble();

						if (!piFlagNew) {
							pGenNew = Math.log(1 - x)/kEff1New + Math.log(1 + x*Math.exp(-2*kEff1New)/(1-x))/kEff1New;
							if (pGenNew > 0) {
								pGenNew = 0;
							}
							if (pGenNew < -2.0) {
								pGenNew = -2.0;
							}
							y1 = 1 + pGenNew;

							newAlpha[imageIndex] = Math.acos(y1);

							if (newAlpha[imageIndex] != newAlpha[imageIndex] || y1 != y1) {
								System.out.println("x = "+x);
								System.out.println("pGenNew = "+pGenNew);
								System.out.println("y1 = "+y1);
								System.out.println("imageIndex = "+imageIndex);
								System.out.println("newAlpha = "+newAlpha[imageIndex]);
								System.out.println("oldY0 = "+oldY0);
								System.out.println("newY0 = "+newY0);
								System.out.println("oldCenterY0 = "+oldCenterY0);
								System.out.println("kEff = "+kEff);
								System.out.println("kEffOld = "+kEff1Old);
								System.out.println("kEffNew = "+kEff1New);
								System.out.println("r = "+bl);
								throw new RuntimeException("kEff1New = " +kEff1New);
							}
						}
						else {
							newAlpha[imageIndex] = 2*Math.PI*x;
						}
					}
					if (imageIndex == P/2) {
						newOrientations[i][imageIndex].setDirection(rV1);
						Vector rV2 = space.makeVector();
						if (Math.abs(rV1.getX(0)) > 0.5) {
							rV2.setX(1,1);
						}
						else {
							rV2.setX(0, 1);
						}
						rV2.PEa1Tv1(-rV2.dot(rV1), rV1);
						rV2.normalize();
						double dummyAlpha = 2*Math.PI*random.nextDouble();
						rotateVectorV(dummyAlpha, rV1, rV2);

						if (!doExchange[i]) {
							rotateVectorV(newAlpha[imageIndex],rV2, newOrientations[i][imageIndex].getDirection());
						}
						else {
							double angle = 2*Math.PI*random.nextDouble();
							rotateVectorV(angle,rV2, newOrientations[i][imageIndex].getDirection());
						}
						theta[imageIndex] = 0; // without the loss of generality
					}
					else {
						newCenter.Ev1Pv2(newOrientations[i][fromImage].getDirection(), newOrientations[i][toImage].getDirection());
						newCenter.normalize();
						newOrientations[i][imageIndex].setDirection(newCenter);

						axis.E(0);
						if (Math.abs(newCenter.getX(0)) > 0.5) {
							axis.setX(1,1);
						}
						else {
							axis.setX(0, 1);
						}
						axis.PEa1Tv1(-axis.dot(newCenter), newCenter);
						axis.normalize();

						rotateVectorV(newAlpha[imageIndex], axis, newOrientations[i][imageIndex].getDirection());
						theta[imageIndex] = random.nextDouble()*2*Math.PI;
						newOrientations[i][imageIndex].rotateBy(theta[imageIndex], newCenter);
					}
					if (newOrientations[i][imageIndex].getDirection().isNaN()) throw new RuntimeException("bead "+imageIndex+" orientation is NaN");
					if (!doExchange[i] || imageIndex != P/2) {
						oldCenter.ME(oldOrientations[i][imageIndex].getDirection());
						double v = oldCenter.squared();
						double v1 = oldCenter.dot(oldOrientations[i][imageIndex].getDirection());
						double s1 = v - v1*v1;
						double y0m1 = -s1/(1+y0);
						//                  	if (xOld[i][imageIndex] == 0) xOld[i][imageIndex] = xNew[i][imageIndex];
						//                  	pGenOld *= Math.exp(kEff1Old*y0)*kEff1Old/Math.sinh(kEff1Old);
						//                  	pGenNew *= Math.exp(kEff1New*y1)*kEff1New/Math.sinh(kEff1New);
						//                  	pGenOld *= 2*Math.exp(kEff1Old*y0m1)*kEff1Old/(1 - Math.exp(-2*kEff1Old));
						//                  	pGenNew *= 2*Math.exp(kEff1New*a)*kEff1New/(1 - Math.exp(-2*kEff1New));
						//                  	System.out.println(kEff1Old+" "+kEff1New+" "+y0m1+" "+(y0-1)+" "+a+" "+(y1_new-1));
						//                  	pGenRatio *= Math.exp(kEff1New*a - kEff1Old*y0m1)*kEff1New*(1-Math.exp(-2*kEff1Old))/(kEff1Old*(1-Math.exp(-2*kEff1New)));
						//                  	double aNew = kEff1New*(-xNew[i][imageIndex] + 1/(1 - Math.exp(-2*kEff1New)))/(kEff1Old*(-xOld[i][imageIndex] + 1/(1 - Math.exp(-2*kEff1Old))));
						double pRatio = 0;
						if (!piFlagOld && !piFlagNew) {
							pRatio = Math.exp(kEff1New*pGenNew - kEff1Old*y0m1)*kEff1New*(1-Math.exp(-2*kEff1Old))/(kEff1Old*(1-Math.exp(-2*kEff1New)));
						}
						else if (!piFlagOld && piFlagNew){
							pRatio = Math.exp(-kEff1Old*y0m1)*(1-Math.exp(-2*kEff1Old))/(2*kEff1Old);
						}
						else if (!piFlagNew && piFlagOld) {
							pRatio = Math.exp(kEff1New*pGenNew)*kEff1New*2/(1-Math.exp(-2*kEff1New));
						}
						else {
							pRatio = 1.0;
						}

						pGenRatio *= pRatio;
						if (Double.isNaN(pGenRatio)) {
							System.out.println("a = "+pGenNew);
							System.out.println("y1 = "+y1);
							System.out.println("imageIndex = "+imageIndex);
							System.out.println("newAlpha = "+newAlpha[imageIndex]);
							System.out.println("kEff = "+kEff);
							System.out.println("r = "+bl);
							System.out.println("v = "+v);
							System.out.println("v1 = "+v1);
							System.out.println("s1 = "+s1);
							System.out.println("y0m1 = "+y0m1);
							throw new RuntimeException("pGenRatio is not a number!");
						}
					}
				}
			}
			for (int j=0; j<P; j++) {
				int prev = j-1;
				if (prev < 0) prev = P-1;
				AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
				AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
				uNew += kHarmonic*dist(jAtom,jPrev, i);

			}
		}
		double pActRatio = Math.exp(uOld-uNew);
		pacc = pActRatio/pGenRatio;
		((BoxCluster)box).trialNotify();
		weightNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
		return true;
	}

	@Override
	public double getA() {
		if (firstMove) {
			return weightNew/weightOld;
		}
		else {
			return pacc*weightNew/weightOld;
		}

	}

	@Override
	public double getB() {
		return 0.0;
	}

	@Override
	public void rejectNotify() {
		IMoleculeList molecules = box.getMoleculeList();
		for (int i=0; i<molecules.getMoleculeCount(); i++) {
			if (molIndexUntouched == i) continue;
			IAtomList atoms = molecules.getMolecule(i).getChildList();
			for (int k=0; k<P; k++) {
				IAtomOriented kAtom = ((IAtomOriented)atoms.getAtom(k));
				kAtom.getOrientation().setDirection(oldOrientations[i][k].getDirection());
			}
		}
		((BoxCluster)box).rejectNotify();
	}

	@Override
	public void acceptNotify() {
		acc++;
		firstMove = false;
		((BoxCluster)box).acceptNotify();
	}

	@Override
	public double energyChange() {
		return 0;
	}

	@Override
	public AtomIterator affectedAtoms() {
		return leafIterator;
	}

	public void setDoExchange(boolean[] b) {
		doExchange = b;
	}

	private double dist(AtomHydrogen a0, AtomHydrogen a1, int moleculeIndex) {
		if (moleculeIndex < 0) throw new ArrayIndexOutOfBoundsException("Molecule Index: "+moleculeIndex);
		utilityVec3.E(a0.getOrientation().getDirection());
		if (a0.getIndex() == 0 && doExchange[moleculeIndex]) {
			utilityVec3.TE(-1);
		}
		utilityVec1.Ea1Tv1(a0.getBondLength(), utilityVec3);
		utilityVec2.Ea1Tv1(a1.getBondLength(), a1.getOrientation().getDirection());
		double r2 = utilityVec1.Mv1Squared(utilityVec2);
		return r2;
	}

	private void rotateVectorV(double angle, Vector axis, Vector v) {
		double q0 = Math.cos(angle/2.0);
		double sth2 = Math.sin(angle/2.0);
		Vector a1 = space.makeVector();
		a1.E(axis);
		a1.TE(sth2);
		double q1 = a1.getX(0);
		double q2 = a1.getX(1);
		double q3 = a1.getX(2);
		Quaternion q = new Quaternion(q0,q1,q2,q3);
		Quaternion vec = new Quaternion(0,v.getX(0),v.getX(1),v.getX(2));
		Quaternion w = q.preMultiply(vec).preMultiply(q.conjugate());
		if (Math.abs(w.getScalar()) > 1E-10 ) throw new RuntimeException("Quaternion product is not a vector!");
		v.E(w.getVector());
	}
}


