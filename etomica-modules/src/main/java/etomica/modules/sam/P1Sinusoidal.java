/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.atom.IAtom;
import etomica.potential.IPotential1;
import etomica.space.Space;
import etomica.space.Vector;

public class P1Sinusoidal implements IPotential1 {

    public P1Sinusoidal(Space space) {
        this.space = space;
        setB(1);
        this.offset = space.makeVector();
        r = space.makeVector();
        waveVectors = new Vector[3];
        setCellSize(1, 1);
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

    @Override
    public double u(IAtom atom) {
        r.Ev1Mv2(atom.getPosition(), offset);
        double sum = 0;
        for (int i=0; i<3; i++) {
            sum += Math.cos(r.dot(waveVectors[i]));
        }
        return b45 * (3.0 - sum);
    }

    @Override
    public double udu(IAtom atom, Vector f) {
        r.Ev1Mv2(atom.getPosition(), offset);
        double sum = 0;
        for (int i=0; i<3; i++) {
            sum += Math.cos(r.dot(waveVectors[i]));
            f.PEa1Tv1(-b45*Math.sin(r.dot(waveVectors[i])), waveVectors[i]);
        }
        return b45 * (3.0 - sum);
    }

    protected final Space space;
    protected double b45;
    protected final Vector offset;
    protected final Vector r;
    protected final Vector[] waveVectors;
}
