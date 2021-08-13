/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveMoleculeFasterer;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.potential.P4BondTorsion;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * An MC Move for cluster simulations that performs torsion moves on a chain
 * molecule (of length N>=4).  The move is performed on all molecules in the
 * Box.  The move needs the torsion potential in order to choose appropriate
 * torsion angles.  The torsion angle is chosen from a Boltzmann distribution,
 * determined numerically from the potential.  The angles are divided into bins
 * of equal probability.  For a move, a bin is chosen at random and accepted
 * based on the ratio of the
 * <p>
 * This class assumes that the torsion potential is symmetric about theta=0 --
 * that u(theta) = u(-theta)
 *
 * @author Andrew Schultz
 */
public class MCMoveTorsionFasterer extends MCMoveMoleculeFasterer {

    protected final P4BondTorsion torsionPotential;
    protected IMoleculePositionDefinition positionDefinition;
    protected final double[] probabilityBins;
    protected final double[] binSize;
    protected final int[] probabilityReverseMap;
    protected final Vector work1, work2, work3;
    protected final Vector dr21, dr23, dr34;
    protected final Vector oldCenter;
    protected double bias;
    protected ISpecies species;

    public MCMoveTorsionFasterer(PotentialCompute pc, Box box, IRandom random,
                                 P4BondTorsion torsionPotential) {
        this(pc, box, random, torsionPotential, 20);
    }

    public MCMoveTorsionFasterer(PotentialCompute pc, Box box,
                                 IRandom random, P4BondTorsion torsionPotential, int nBins) {
        super(random, pc, box);
        ((MCMoveStepTracker)getTracker()).setTunable(false);
        positionDefinition = new MoleculePositionGeometricCenter(space);
        probabilityBins = new double[nBins+1];
        binSize = new double[nBins];
        probabilityReverseMap = new int[nBins+1];
        this.torsionPotential = torsionPotential;
        work1 = space.makeVector();
        work2 = space.makeVector();
        work3 = space.makeVector();
        dr21 = space.makeVector();
        dr23 = space.makeVector();
        dr34 = space.makeVector();
        oldCenter = space.makeVector();
    }

    public void setBox(Box p) {
        super.setBox(p);
        ((MCMoveStepTracker)getTracker()).setTunable(false);
    }

    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }
    public void setTemperature(double temperature) {
        int nBins = probabilityBins.length - 1;
        int nSubBins = 100;
        double beta = 1.0/temperature;
        // numerically integrate P = exp(-beta U) from cosphi=1 to cosphi=-1 (0 to pi radians)
        double totP = 0.5 * Math.exp(-beta * torsionPotential.u(1));
        for (int i=1; i<nSubBins * nBins; i++) {
            double cosphi = Math.cos((Math.PI * i) / (nSubBins * nBins));
            double u = torsionPotential.u(cosphi);
            totP +=  Math.exp(-beta*u);
        }
        totP += 0.5 * Math.exp(-beta * torsionPotential.u(-1));
        // pPerBin is the amount of probability allocated to each bin.
        double pPerBin = totP / nBins;
        
        double thisBinP = 0;
        double previousP = 1;
        int iBin = 1;
        // last bin always ends at phi=PI
        int lastMapBin = 0;
        // probabilityBins is an array of the boundaries of the bins (in terms of cos phi)
        // bin i goes from probabilityBins[i] to probabilityBins[i+1]
        // the first bin is bin 0 and the last bin is bin nBins-1
        probabilityBins[0] = 0;
        probabilityBins[nBins] = Math.PI;
        
        // probabilityReverseMap maps an angle to a bin
        // probabilityReverseMap[phi / PI * nBins] will be less than or
        //   equal to the appropriate bin.
        // probabilityReverseMap[(phi / PI * nBins) + 1] will be greater than or
        //   equal to the appropriate bin.
        probabilityReverseMap[0] = 0;
        // sum up probability (=exp(-beta U)) until we exceed pPerBin.  Then
        // interpolate to find approximate when the sum is equal to pPerBin.
        // Call that the bin boundary and then begin summing for the next bin.
        for (int i=1; i<nSubBins * nBins + 1; i++) {
            double cosphi = Math.cos((Math.PI * i) / (nSubBins * nBins));
            double u = torsionPotential.u(cosphi);
            double thisP = Math.exp(-beta *u);
            double newP = 0.5 * (previousP + thisP);
//            System.out.println(i+" "+cosphi+" "+u+" "+newTot+" "+thisBinP+" "+pPerBin);
            // if i is too much probability, then ix would be the value
            // of i (based on interpolation) where we hit the limit
            double ix = (i-1) + (pPerBin-thisBinP)/newP;
            while (ix < i) {
                // we have too much probability for this bin.  back up (interpolate)
                // to find where we hit the limit for this bin.  dump the rest into the
                // next bin
                probabilityBins[iBin] = (ix*Math.PI)/(nSubBins*nBins);
                // we'll use this to correct acceptance (which is not correct due
                // to energy inhomogeniety within the bin)
                binSize[iBin-1] = (probabilityBins[iBin]-probabilityBins[iBin-1])/(Math.PI/nBins);

                // we found the lower bound of iBin.  Now mark 
                int nextMapBin = (int)((probabilityBins[iBin])/Math.PI * nBins);
                for (int j=lastMapBin+1; j<nextMapBin+1; j++) {
                    probabilityReverseMap[j] = iBin-1;
//                    System.out.println(j+" "+Math.acos(-1.0+(2.0*j)/nBins)+" "+(iBin-1));
                }
                lastMapBin = nextMapBin;
//                System.out.println(iBin+" "+probabilityBins[iBin]+" "+Math.acos(probabilityBins[iBin]) +" "+ torsionPotential.energyAtAngle(probabilityBins[iBin]));
//                double phi = Math.acos((probabilityBins[iBin-1]+probabilityBins[iBin])*0.5);
//                System.out.println(phi+" "+ torsionPotential.energyAtAngle(Math.cos(phi)) + " "+ (probabilityBins[iBin-1]-probabilityBins[iBin]));
                iBin++;
                if (iBin == nBins) {
                    // we just found the end of the next-to-the-last bin
                    // we know the last bin ends at PI (and already marked it as such)
//                    phi = Math.acos((probabilityBins[iBin-1]+probabilityBins[iBin])*0.5);
//                    System.out.println(phi+" "+ torsionPotential.energyAtAngle(Math.cos(phi))+" "+(probabilityBins[iBin-1]-probabilityBins[iBin]));
                    // we'll use this to correct acceptance (which is not correct due
                    // to energy inhomogeniety within the bin)
                    binSize[iBin-1] = (probabilityBins[iBin]-probabilityBins[iBin-1])/(Math.PI/nBins);
                    for (int j=lastMapBin+1; j<nBins+1; j++) {
                        probabilityReverseMap[j] = iBin-1;
//                        System.out.println(j+" "+Math.acos(-1.0+(2.0*j)/nBins)+" "+(iBin-1));
                    }
                    return;
                }
                thisBinP -= pPerBin;
                // if we still have too much probability (for the next bin)
                // ix will be the the new value where we hit the limit
                ix += pPerBin/newP;
            }
            thisBinP += newP;
            previousP = thisP;
        }
        throw new RuntimeException("oops");
    }
    
    //note that total energy is calculated
    public boolean doTrial() {
        molecule = moleculeSource.getMolecule();
        while (oldPositions.size() < molecule.getChildList().size()) {
            oldPositions.add(space.makeVector());
        }
        uOld = potentialCompute.computeOneOldMolecule(molecule);
        bias = 1;

        oldCenter.E(CenterOfMass.position(box, molecule));
        IAtomList childList = molecule.getChildList();
        int numChildren = childList.size();

        int j = random.nextInt(numChildren-3);  // j=0 ==> first torsion bond (atoms 0,1,2,3)
//            System.out.println("rotating about "+(j+1)+" and "+(j+2));
        // atoms j+1 and j+2 stay fixed.  atoms 0 to j and j+3 to N-1 move
        IAtom atom0 = childList.get(j+0);
        IAtom atom1 = childList.get(j+1);
        IAtom atom2 = childList.get(j+2);
        IAtom atom3 = childList.get(j+3);
        dr21.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
        box.getBoundary().nearestImage(dr21);
        dr23.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
        box.getBoundary().nearestImage(dr23);
        dr34.Ev1Mv2(atom3.getPosition(), atom2.getPosition());
        box.getBoundary().nearestImage(dr34);

//            System.out.println("|dr21| "+Math.sqrt(dr21.squared()));
//            if (Math.abs(Math.sqrt(dr21.squared())-1.54) > 0.0001) {
//                throw new RuntimeException("oops dr21 "+Math.sqrt(dr21.squared())+" "+i);
//            }
//            System.out.println("|dr23| "+Math.sqrt(dr23.squared()));
//            if (Math.abs(Math.sqrt(dr23.squared())-1.54) > 0.0001) {
//                throw new RuntimeException("oops dr23 "+Math.sqrt(dr23.squared()));
//            }
//            System.out.println("|dr34| "+Math.sqrt(dr34.squared()));
//            if (Math.abs(Math.sqrt(dr34.squared())-1.54) > 0.0001) {
//                throw new RuntimeException("oops dr34 "+Math.sqrt(dr34.squared()));
//            }

        double dr23Sq = dr23.squared();
        dr23.TE(1.0/Math.sqrt(dr23Sq));
        dr21.PEa1Tv1(-dr21.dot(dr23), dr23);
        dr34.PEa1Tv1(-dr34.dot(dr23), dr23);

        // current angle
        double oldcosphi = dr21.dot(dr34)/Math.sqrt(dr21.squared()*dr34.squared());
        if (oldcosphi > 0.99999 || oldcosphi < -0.99999) {
            // torsion angle is almost exactly 0 or PI.  skip this molecule
            return false;
        }
        double oldphi = Math.acos(oldcosphi);
        int oldTorsionBin = probabilityReverseMap[(int)((oldphi/Math.PI)*(probabilityBins.length-1))];
//            System.out.println("oldTorsionBin0 = "+oldTorsionBin);
        if (probabilityBins[oldTorsionBin] > oldphi) {
            throw new RuntimeException("reversemap overshot "+oldphi+" "+oldTorsionBin+" "+probabilityBins[oldTorsionBin]);
        }
        while (oldTorsionBin < probabilityReverseMap.length && probabilityBins[oldTorsionBin+1] < oldphi) {
            oldTorsionBin++;
        }

        work1.E(dr21);
        work1.XE(dr34);
        // do we need to rotate the first part of the molecule clockwise
        // or counter-clockwise?
        boolean flippit = work1.dot(dr23) > 0;

        // pick a new bin
        int torsionBin = random.nextInt(probabilityBins.length-1);
        // pick some angle within that bin
        double newphi = probabilityBins[torsionBin] + random.nextDouble()*
                    (probabilityBins[torsionBin+1]-probabilityBins[torsionBin]);
        double deltaphi = 0.5*(newphi - oldphi);

        // if we move from a high-U (wide) bin to a low-U (narrow) bin, the nominal
        // acceptance probability will be too high.  correct for that here.
        bias *= binSize[torsionBin]/binSize[oldTorsionBin];

        if (flippit) deltaphi = -deltaphi;

        for (int k=0; k<numChildren; k++) {
            IAtom atomk = childList.get(k);
            oldPositions.get(k).E(atomk.getPosition());
            if (k == j+3) {
                deltaphi = -deltaphi;
            }
            else if (k > j && k < j+3) {
                continue;
            }

//                System.out.println("rotating "+atomk);

            // we can use atom1 as a reference atom to rotate the first and
            // last part of the molecule.  anything on the atom2-atom3 axis
            // would work
            work1.Ev1Mv2(atomk.getPosition(), atom1.getPosition());
            box.getBoundary().nearestImage(work1);
//                System.out.println("|work1| "+Math.sqrt(work1.squared()));
            work2.E(dr23);   // dr23 is the axis of rotation
            // v1 = v1overAxis * axis
            double v1overAxis = work2.dot(work1);

//                if (Math.abs(Math.abs(v1overAxis)-Math.sqrt(work1.squared())) < 1E-10) {
//                    // axis is almost exactly parallel or anti-parallel to direction,
//                    // so just don't rotate.
//                    continue;
//                }

            work2.TE(-v1overAxis);
            work2.PE(work1);
            dr21.E(work2);
            dr21.PEa1Tv1(v1overAxis, dr23);
//                System.out.println("|v1+v2| "+Math.sqrt(dr21.squared()));

            // now temp = v2
            double v2Sq = work2.squared();
            work3.E(dr23);
            work3.XE(work1);
            work3.TE(Math.sqrt(v2Sq/work3.squared()));
            //System.out.println(Math.sqrt(work2.squared())+" "+Math.sqrt(work3.squared()));

            dr21.E(work3);
            dr21.PEa1Tv1(v1overAxis, dr23);
//                System.out.println("|v1+v3| "+Math.sqrt(dr21.squared()));

            // now temp = v3
//                System.out.println("dot "+dr23.dot(work2)+" "+work2.dot(work3)+" "+dr23.dot(work3));

            work1.Ea1Tv1(Math.cos(deltaphi), work2);
            work1.PEa1Tv1(Math.sin(deltaphi), work3);
            work1.PEa1Tv1(v1overAxis, dr23);
//                System.out.println("|work1| "+Math.sqrt(work1.squared()));
            work1.PE(atom1.getPosition());

            atomk.getPosition().E(work1);

            if (k==j || k == j+2) {
                dr21.Ev1Mv2(atom0.getPosition(), atom1.getPosition());
                box.getBoundary().nearestImage(dr21);
//                    System.out.println("|dr21| "+Math.sqrt(dr21.squared()));
//                    if (Math.abs(Math.sqrt(dr21.squared())-1.54) > 0.0001) {
//                        throw new RuntimeException("oops dr21 "+Math.sqrt(dr21.squared()));
//                    }
                dr23.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                box.getBoundary().nearestImage(dr23);
//                    System.out.println("|dr23| "+Math.sqrt(dr23.squared()));
//                    if (Math.abs(Math.sqrt(dr23.squared())-1.54) > 0.0001) {
//                        throw new RuntimeException("oops dr23 "+Math.sqrt(dr23.squared()));
//                    }
                dr34.Ev1Mv2(atom3.getPosition(), atom2.getPosition());
                box.getBoundary().nearestImage(dr34);
//                    System.out.println("|dr34| "+Math.sqrt(dr34.squared()));
//                    if (Math.abs(Math.sqrt(dr34.squared())-1.54) > 0.0001) {
//                        throw new RuntimeException("oops dr34 "+Math.sqrt(dr34.squared()));
//                    }

                dr23Sq = dr23.squared();
                dr23.TE(1.0/Math.sqrt(dr23Sq));
                dr21.PEa1Tv1(-dr21.dot(dr23), dr23);
                dr34.PEa1Tv1(-dr34.dot(dr23), dr23);

                // current angle
                double newcosphicheck = dr21.dot(dr34)/Math.sqrt(dr21.squared()*dr34.squared());
//                    System.out.println(Math.acos(oldcosphi)+" "+Math.acos(newcosphi)+" "+Math.acos(newcosphicheck));
//                    System.out.println(oldcosphi+" "+newcosphi+" "+newcosphicheck);
                if (k<j+1) {
//                        System.out.println("error: "+((Math.acos(oldcosphi) + Math.acos(newcosphi))*0.5-Math.acos(newcosphicheck)));
                    if (Math.abs((oldphi + newphi)*0.5-Math.acos(newcosphicheck)) > 1.e-8) {
                        throw new RuntimeException("shouldn't need to flip");
                    }
                }
                else {
//                        System.out.println("error: "+(Math.acos(newcosphi)-Math.acos(newcosphicheck)));
                    if (Math.abs(newphi-Math.acos(newcosphicheck)) > 1.e-8) {
                        throw new RuntimeException("shouldn't need to flip");
                    }
                }
            }
        }
        oldCenter.ME(CenterOfMass.position(box, molecule));
        for (int k=0; k<numChildren; k++) {
            // shift the whole molecule so that the center of mass doesn't change
            IAtom atomk = childList.get(k);
            Vector r = atomk.getPosition();
            r.PE(oldCenter);
            r.PE(box.getBoundary().centralImage(r));
            potentialCompute.updateAtom(atomk);
        }
        return true;
    }

    public void acceptNotify() {
//        System.out.println("accepted");
        potentialCompute.processAtomU(1);
        Vector[] newPositions = new Vector[oldPositions.size()];
        // put it back, then compute old contributions to energy
        molecule.getChildList().forEach(atom -> {
            newPositions[atom.getIndex()] = box.getSpace().makeVector();
            newPositions[atom.getIndex()].E(atom.getPosition());
            atom.getPosition().E(oldPositions.get(atom.getIndex()));
            potentialCompute.updateAtom(atom);
        });
        potentialCompute.computeOneMolecule(molecule);
        potentialCompute.processAtomU(-1);
        molecule.getChildList().forEach(atom -> {
            atom.getPosition().E(newPositions[atom.getIndex()]);
            potentialCompute.updateAtom(atom);
        });
    }
}
