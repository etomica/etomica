/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.potential.PotentialMaster;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.util.random.IRandom;


public class MCMoveRotateMolecule3D extends MCMoveMolecule {

    private static final long serialVersionUID = 2L;
    protected transient Vector r0;
    //TODO I changed RotationTensor to RotationTensor3D here, is it ok???
    protected transient RotationTensor3D rotationTensor;
    protected IMoleculePositionDefinition positionDefinition;

    public MCMoveRotateMolecule3D(PotentialMaster potentialMaster, IRandom random,
                                  Space _space) {
        super(potentialMaster, random, _space, Math.PI / 2, Math.PI);
        rotationTensor = (RotationTensor3D) _space.makeRotationTensor();
        r0 = _space.makeVector();
        positionDefinition = new MoleculePositionGeometricCenter(space);
    }

    public IMoleculePositionDefinition getPositionDefinition() {
        return positionDefinition;
    }

    public void setPositionDefinition(IMoleculePositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }

    public boolean doTrial() {
//        System.out.println("doTrial MCMoveRotateMolecule called");

        if (box.getMoleculeList().getMoleculeCount() == 0) {
            molecule = null;
            return false;
        }


        molecule = moleculeSource.getMolecule();

        //TODO  only move first molecule
        boolean onlyOneRotate = false;
        if (onlyOneRotate) {
            IMoleculeList molecules = box.getMoleculeList();
            molecule = molecules.getMolecule(0);
        }

        energyMeter.setTarget(molecule);
        uOld = energyMeter.getDataAsScalar();

        if (Double.isInfinite(uOld)) {
            throw new RuntimeException("Overlap in initial state");
        }
        double dTheta = (2 * random.nextDouble() - 1.0) * stepSize;


        boolean onlyKappa3 = false;
        if (onlyKappa3) {
            //        //TODO only rotate around om so only kappa3
//        Vector h1 = molecule.getChildList().getAtom(0).getPosition();
//        Vector h2 = molecule.getChildList().getAtom(1).getPosition();
            Vector o = molecule.getChildList().getAtom(2).getPosition();
            Vector m = molecule.getChildList().getAtom(3).getPosition();
            Vector dr = space.makeVector();
            dr.Ev1Mv2(m, o);
            dr.normalize();
            rotationTensor.setRotationAxis(dr, dTheta);

        } else {
            rotationTensor.setAxial(r0.getD() == 3 ? random.nextInt(3) : 2, dTheta);
        }

        r0.E(positionDefinition.position(molecule));
        doTransform();

        energyMeter.setTarget(molecule);
        return true;
    }

    protected void doTransform() {
        IAtomList childList = molecule.getChildList();
        for (int iChild = 0; iChild < childList.getAtomCount(); iChild++) {
            IAtom a = childList.getAtom(iChild);
            Vector r = a.getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
        }
    }

    public void rejectNotify() {
        rotationTensor.invert();
        doTransform();
    }
}
