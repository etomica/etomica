/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.potential.IPotential;
import etomica.potential.IPotentialField;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

public class P1Sinusoidal implements IPotential, IPotentialField {

    public P1Sinusoidal(Space space) {
        this.space = space;
        setB(1);
        this.offset = space.makeVector();
        r = space.makeVector();
        waveVectors = new Vector[3];
        setCellSize(1, 1);
        gradient = new Vector[1];
        gradient[0] = space.makeVector();
    }
    
    public void setOffset(Vector newOffset) {
        offset.E(newOffset);
    }
    
    public Vector getOffset() {
        return offset;
    }
    
    public void setB(double newB) {
        b45 = newB/4.5;
    }
    
    public double getB() {
        return 4.5*b45;
    }
    
    public void setCellSize(double xSize, double zSize) {
        waveVectors[0] = Vector.of(new double[]{1 / xSize, 0, -1 / zSize});
        waveVectors[1] = Vector.of(new double[]{-1 / xSize, 0, -1 / zSize});
        waveVectors[2] = Vector.of(new double[]{0, 0, 2.0 / zSize});
        waveVectors[0].TE(2.0*Math.PI);
        waveVectors[1].TE(2.0*Math.PI);
        waveVectors[2].TE(2.0*Math.PI);
    }
    
    public double energy(IAtomList atoms) {
        IAtom a = atoms.get(0);
        r.Ev1Mv2(a.getPosition(), offset);
        double sum = 0;
        for (int i=0; i<3; i++) {
            sum += Math.cos(r.dot(waveVectors[i]));
        }
        return b45 * (3.0 - sum);
    }

    public Vector[] gradient(IAtomList atoms) {
        IAtom a = atoms.get(0);
        r.Ev1Mv2(a.getPosition(), offset);
        gradient[0].E(0);
        for (int i=0; i<3; i++) {
            gradient[0].PEa1Tv1(Math.sin(r.dot(waveVectors[i])), waveVectors[i]);
        }
        gradient[0].TE(b45);
        return gradient;
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    public double getRange() {
        return 0;
    }

    protected final Space space;
    protected double b45;
    protected final Vector offset;
    protected final Vector r;
    protected final Vector[] waveVectors;
    protected final Vector[] gradient;
}
