/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;


/**
 * Class which handles orientation of linear molecules during volume change
 * trials.  During compression, the orientation is scaled back toward the
 * nominal orientation according to latticeScale as calculated by the
 * superClass.
 * 
 * @author Andrew Schultz
 */
public class MCMoveVolumeSolidNPTMolecularOriented extends
        MCMoveVolumeSolidNPTMolecular {

    protected final Vector orientation, dr, com;
    protected final Vector[] drSum, drSumSave;
    protected double maxPhi;
    protected double cosNominalTheta;
    protected transient RotationTensor rotationTensor;
    protected int maxMolecule;
    protected double thetaFrac = 1;
    
    public MCMoveVolumeSolidNPTMolecularOriented(PotentialMaster potentialMaster, IRandom random,
                                                 Space space, double pressure, Vector[] drSum) {
        super(potentialMaster, random, space, pressure, 5);
        rotationTensor = space.makeRotationTensor();
        this.drSum = drSum;
        orientation = space.makeVector();
        dr = space.makeVector();
        com = space.makeVector();
        drSumSave = new Vector[drSum.length];
        for (int i=0; i<drSumSave.length; i++) {
            drSumSave[i] = space.makeVector();
        }
        maxMolecule = -1;
    }
    
    public void setNominalTheta(double newNominalTheta) {
        cosNominalTheta = Math.cos(newNominalTheta);
    }
    
    public void setMaxPhi(double newMaxPhi) {
        maxPhi = newMaxPhi;
    }
    
    public void setThetaFrac(double newThetaFrac) {
        thetaFrac = newThetaFrac;
    }
    
    public double getThetaFrac() {
        return thetaFrac;
    }
    
    public void doTransform(double vScaleLocal) {
        super.doTransform(vScaleLocal);

        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.getMoleculeCount();
        if (maxMolecule == -1) {
            maxMolecule = nMolecules;
        }

        for (int i=0; i<maxMolecule; i++) {
            IMolecule molecule = moleculeList.getMolecule(i);
            
            IAtomList atoms = molecule.getChildList();

            dr.E(atoms.getAtom(1).getPosition());
            dr.ME(atoms.getAtom(0).getPosition());
            dr.normalize();
            double sintheta = Math.sqrt(dr.getX(0)*dr.getX(0) + dr.getX(1)*dr.getX(1));
            double costheta = dr.getX(2);
            double u3 = costheta - cosNominalTheta;
            
            double x = dr.getX(0);
            double y = dr.getX(1);
            double phi = Math.atan2(y, x);
            double u4 = phi * sintheta;
            
            u3 *= Math.pow(latticeScale, thetaFrac);
            u4 *= Math.pow(latticeScale, 2.0-thetaFrac);
            
            double newCosTheta = u3 + cosNominalTheta;
            double newSinTheta = Math.sqrt(1-newCosTheta*newCosTheta);
            double newPhi = u4 / newSinTheta;

            if (newCosTheta < 1e-8 || newCosTheta > 1 - 1e-8 || Math.abs(newPhi) > Math.PI - 1e-8) {
                // we don't rotate theta past vertical or horizontal.  also can't rotate phi past backwards
                // bail now (don't rotate).  we'll force move rejection.
                // when we go to unrotate, we'll stop at the previous molecule
                maxMolecule = i;
//                if (newCosTheta < 0 || newCosTheta > 1) {
//                    System.out.println("bailing because cosTheta "+costheta+" => "+newCosTheta+" for "+molecule);
//                }
//                else {
//                    System.out.println("bailing because phi "+phi+" => "+newPhi+" for "+molecule);
//                }
                return;
            }
            
            double a = newSinTheta/sintheta;
            double b = newCosTheta - newSinTheta/sintheta*costheta;
            doTransformTheta(molecule, a, b);

            rotationTensor.setAxial(2, newPhi - phi);
            doTransformPhi(molecule);
        }
        maxMolecule = -1;
    }

    protected void doTransformPhi(IMolecule molecule) {
        com.E(moleculeCenter.position(molecule));
        IAtomList childList = molecule.getChildList();
        for (int iChild = 0; iChild<childList.getAtomCount(); iChild++) {
            IAtom a = childList.getAtom(iChild);
            Vector r = a.getPosition();
            r.ME(com);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(com);
        }
    }

    protected void doTransformTheta(IMolecule molecule, double a, double b) {
        com.E(moleculeCenter.position(molecule));
        IAtomList childList = molecule.getChildList();
        for (int i=0; i<childList.getAtomCount(); i++) {
            Vector p1 = childList.getAtom(i).getPosition();
            p1.ME(com);
            double bl = Math.sqrt(p1.squared());
            if (p1.getX(2) < 0) bl = -bl;
            p1.TE(a);
            p1.setX(2, p1.getX(2) + b*bl);
            p1.PE(com);
        }
    }

    public double getChi(double temperature) {
        if (maxMolecule > -1) {
            // we over-rotated a molecule (past theta=0 or theta=90)
            // reject
            return 0;
        }
        IMoleculeList moleculeList = box.getMoleculeList();
        int nMolecules = moleculeList.getMoleculeCount();
        for (int i=0; i<drSum.length; i++) {
            drSumSave[i].E(drSum[i]);
            drSum[i].E(0);
        }
        // now recalculate drSum
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = moleculeList.getMolecule(i);
            
            int nPlanes = drSum.length;
            int iPlane = (molecule.getIndex()/2)%nPlanes;
            Vector drSumi = drSum[iPlane];

            IAtomList atoms = molecule.getChildList();
            // add the trial orientation
            drSumi.PE(atoms.getAtom(1).getPosition());
            drSumi.ME(atoms.getAtom(0).getPosition());
        }

        for (int i=0; i<drSum.length; i++) {
            double phi = Math.atan2(drSum[i].getX(1), drSum[i].getX(0));
            if (Math.abs(phi) > maxPhi) {
//                System.out.println("rejecting because phi "+i+" = "+phi+" (>"+maxPhi+")");
                return 0;
            }
        }
        return super.getChi(temperature);
    }

    public void rejectNotify() {
        boolean bailed = maxMolecule > -1;
        super.rejectNotify();
        if (!bailed) {
            // if the trial was successful then we calculated new drSum in getA
            // we rejected the trial, so restore drSum
            for (int i=0; i<drSum.length; i++) {
                drSum[i].E(drSumSave[i]);
            }
        }
        if (maxMolecule != -1) {
            throw new RuntimeException("oops");
        }
    }
    
    public void acceptNotify() {
        maxMolecule = -1;
    }
}
