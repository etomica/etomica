/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveMoleculeRotateFasterer;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionCOM;
import etomica.molecule.MoleculePositionCOMPBC;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveRotateMoleculePhiThetaFasterer extends MCMoveMoleculeRotateFasterer {

    protected final PotentialCompute potentialCompute;
    protected final Vector[] drSum;
    protected double maxAngle;
    protected boolean rotatePhi;
    protected final Vector dr, com;
    protected IMoleculePositionDefinition pos;
    protected double delta;

    public MCMoveRotateMoleculePhiThetaFasterer(PotentialCompute potentialCompute, Box box,
                                                IRandom random, Vector[] drSum, boolean doPhi) {
        super(random, potentialCompute, box);
        this.potentialCompute = potentialCompute;
        this.drSum = drSum;
        dr = box.getSpace().makeVector();
        rotatePhi = doPhi;
        com = box.getSpace().makeVector();
        pos = new MoleculePositionCOM(space);

        IMoleculeList molecules = box.getMoleculeList();
        int nPlanes = drSum.length;
        for (int iPlane=0; iPlane<nPlanes; iPlane++) {
            drSum[iPlane].E(0);
        }
        for (int i = 0; i<molecules.size(); i++) {
            IAtomList atoms = molecules.get(i).getChildList();
            int iPlane = (i/2)%nPlanes;
            drSum[iPlane].PE(atoms.get(1).getPosition());
            drSum[iPlane].ME(atoms.get(0).getPosition());
        }
        setBox(box);
    }
    
    public void setMaxPhi(double newMaxAngle) {
        maxAngle = newMaxAngle;
    }
    
    public double getConstraint() {
        return maxAngle;
    }
    
    public boolean doTrial() {
//      System.out.println("doTrial MCMoveRotateMolecule called");
      
        if(box.getMoleculeList().size()==0) {molecule = null; return false;}
          
        molecule = moleculeSource.getMolecule();
        uOld = potentialCompute.computeOneOldMolecule(molecule);
        while (oldPositions.size() < molecule.getChildList().size()) {
            oldPositions.add(space.makeVector());
        }
        for (IAtom a : molecule.getChildList()) {
            oldPositions.get(a.getIndex()).E(a.getPosition());
        }

        if(Double.isInfinite(uOld)) {
            throw new RuntimeException("Overlap in initial state");
        }
        
        int nPlanes = drSum.length;
        int iPlane = (molecule.getIndex()/2)%nPlanes;
        Vector drSumi = drSum[iPlane];
        IAtomList atoms = molecule.getChildList();
        // subtract the original orientation
        double planePhi = Math.atan2(drSumi.getX(1), drSumi.getX(0));
        if (Math.abs(planePhi) > maxAngle) {
            throw new RuntimeException("oops "+molecule.getIndex()+" "+planePhi);
        }
        drSumi.ME(atoms.get(1).getPosition());
        drSumi.PE(atoms.get(0).getPosition());
        Vector mydr = Vector.d(3);
        mydr.ME(atoms.get(1).getPosition());
        mydr.PE(atoms.get(0).getPosition());

        if (rotatePhi) {
            delta = (2*random.nextDouble() - 1.0)*stepSize;
            rotationTensor.setAxial(2,delta);

            r0.E(MoleculePositionCOMPBC.com(box.getBoundary(), molecule));
            doTransform();

            mydr.E(0);
            mydr.ME(atoms.get(1).getPosition());
            mydr.PE(atoms.get(0).getPosition());
        }
        else {
            dr.E(atoms.get(1).getPosition());
            dr.ME(atoms.get(0).getPosition());
            dr.normalize();
            double sintheta = Math.sqrt(dr.getX(0)*dr.getX(0) + dr.getX(1)*dr.getX(1));
            double costheta = dr.getX(2);
            delta = (2*random.nextDouble() - 1.0)*stepSize;

            double x = dr.getX(0);
            double y = dr.getX(1);
            double phi = Math.atan2(y, x);
            
            double newCosTheta = costheta + delta;
            double newSinTheta = Math.sqrt(1-newCosTheta*newCosTheta);

            double newPhi = phi * sintheta / newSinTheta;

            if (newCosTheta < 1e-8 || newCosTheta > 1 - 1e-8 || Math.abs(newPhi) > Math.PI - 1e-8) {
                // we don't rotate past vertical or horizontal, also can't rotate backwards
                drSumi.PE(atoms.get(1).getPosition());
                drSumi.ME(atoms.get(0).getPosition());
                return false;
            }
//            System.out.println(sintheta+" "+costheta+"  =>  "+newSinTheta+" "+newCosTheta);

            double a = newSinTheta/sintheta;
            double b = newCosTheta - newSinTheta/sintheta*costheta;

            doRotationTransform(a, b);
            rotationTensor.setAxial(2,newPhi-phi);

            r0.E(MoleculePositionCOMPBC.com(box.getBoundary(), molecule));
            doTransform();
        }

      
        return true;
    }

    protected void doRotationTransform(double a, double b) {
        com.E(MoleculePositionCOMPBC.com(box.getBoundary(), molecule));
        IAtomList childList = molecule.getChildList();
        for (int i = 0; i<childList.size(); i++) {
            Vector p = childList.get(i).getPosition();
            p.ME(com);
            double l = Math.sqrt(p.squared());
            if (p.getX(2) < 0) l = -l;
            p.TE(a);
            p.setX(2, p.getX(2) + b*l);
            p.PE(com);
        }
    }


    public double getChi(double temperature) {
        int nPlanes = drSum.length;
        int iPlane = (molecule.getIndex()/2)%nPlanes;
        Vector drSumi = drSum[iPlane];
        IAtomList atoms = molecule.getChildList();
        // add the trial orientation
        drSumi.PE(atoms.get(1).getPosition());
        drSumi.ME(atoms.get(0).getPosition());
        double phi = Math.atan2(drSumi.getX(1), drSumi.getX(0));

//        System.out.println("rotating "+molecule.getIndex()+" "+phi);
        if (Math.abs(phi) > maxAngle) {
            return 0;
        }
        return super.getChi(temperature);
    }

    public void acceptNotify() {
//        System.out.println("accepted");
        potentialCompute.processAtomU(1);
        // put it back, then compute old contributions to energy
        molecule.getChildList().forEach(atom -> {
            atom.getPosition().E(oldPositions.get(atom.getIndex()));
            potentialCompute.updateAtom(atom);
        });
        potentialCompute.computeOneMolecule(molecule);
        potentialCompute.processAtomU(-1);

        if (rotatePhi) {
            doTransform();
        }
        else {
            IAtomList atoms = molecule.getChildList();
            dr.E(atoms.get(1).getPosition());
            dr.ME(atoms.get(0).getPosition());
            dr.normalize();
            double sintheta = Math.sqrt(dr.getX(0)*dr.getX(0) + dr.getX(1)*dr.getX(1));
            double costheta = dr.getX(2);

            double x = dr.getX(0);
            double y = dr.getX(1);
            double phi = Math.atan2(y, x);

            double newCosTheta = costheta + delta;
            double newSinTheta = Math.sqrt(1-newCosTheta*newCosTheta);

            double newPhi = phi * sintheta / newSinTheta;

            double a = newSinTheta/sintheta;
            double b = newCosTheta - newSinTheta/sintheta*costheta;

            doRotationTransform(a, b);
            doTransform();
        }

//        double utot = potentialCompute.computeAll(false);
//        if (utot==Double.POSITIVE_INFINITY || Double.isNaN(utot)) {
//            throw new RuntimeException("oops "+utot);
//        }


    }

    public void rejectNotify() {
        int nPlanes = drSum.length;
        int iPlane = (molecule.getIndex()/2)%nPlanes;
        Vector drSumi = drSum[iPlane];
        IAtomList atoms = molecule.getChildList();
        // subtract the new orientation
        drSumi.ME(atoms.get(1).getPosition());
        drSumi.PE(atoms.get(0).getPosition());

        molecule.getChildList().forEach(atom -> {
            atom.getPosition().E(oldPositions.get(atom.getIndex()));
            potentialCompute.updateAtom(atom);
        });

        // subtract the old orientation
        drSumi.PE(atoms.get(1).getPosition());
        drSumi.ME(atoms.get(0).getPosition());
    }
    
    public String toString() {
        return super.toString()+" "+(rotatePhi ? "phi" : "theta");
    }
}
