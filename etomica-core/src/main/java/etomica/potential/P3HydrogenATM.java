/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.api.IPotentialAtomic;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.IAtomOriented;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;
import etomica.units.Kelvin;

public class P3HydrogenATM implements IPotential{
    public static void main(String[] args) {       
        try {
            FileWriter file = new FileWriter("GarberoglioATM.dat");
            int maxi = 10000;
            for (int i=0; i<maxi; i++) {
                double r = 1.5 + i*10.00/maxi;
                double r3 = r*r*r;
                double E = E0*(1 + 3/8)/(r3*r3*r3); 
                file.write(r+" "+E+"\n");            
            }
            file.close();
        } catch(IOException e) {
            throw new RuntimeException("caught IOException: " + e.getMessage());
        }
    }        
    protected IBoundary boundary;
    protected static final double E0 = 49400; //(K Angstorm^9)
    protected IVectorMutable dr,r1,r2,r3;

    public P3HydrogenATM(ISpace space) {
        dr = space.makeVector();
        r1 = space.makeVector();
        r2 = space.makeVector();
        r3 = space.makeVector();
    }    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(IBox box) {    
        boundary = box.getBoundary();
    }
    
    public int nBody() {
        return 3;
    }
    public static class P3HydrogenAtomic extends P3HydrogenATM implements IPotentialAtomic {
        public P3HydrogenAtomic(ISpace space) {
            super(space);     
        }
        
        public double energy(IAtomList atoms) {
            IAtom a0 = atoms.getAtom(0);
            IAtom a1 = atoms.getAtom(1);        
            IAtom a2 = atoms.getAtom(2);                    
            IVector com0 = a0.getPosition();               
            IVector com1 = a1.getPosition();
            IVector com2 = a2.getPosition();
            IVector hh0 = ((IAtomOriented)a0).getOrientation().getDirection();
            IVector hh1 = ((IAtomOriented)a1).getOrientation().getDirection();
            IVector hh2 = ((IAtomOriented)a2).getOrientation().getDirection();
            dr.Ev1Mv2(com0, com1);
            double r12 = dr.squared()*Math.sqrt(dr.squared());            
            dr.Ev1Mv2(com1, com2);
            double r23 = dr.squared()*Math.sqrt(dr.squared());
            dr.Ev1Mv2(com0, com2);
            double r13 = dr.squared()*Math.sqrt(dr.squared());
            double c12 = hh0.dot(hh1);
            double c23 = hh1.dot(hh2);
            double c13 = hh0.dot(hh2);
            double U = E0*(1 + 3*c12*c23*c13)/(r12*r23*r13);
            return U;
        }
    }

}
