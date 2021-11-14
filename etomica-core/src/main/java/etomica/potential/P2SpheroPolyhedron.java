/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomOrientedQuaternion;
import etomica.atom.AtomTypeSpheroPolyhedron;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.spaceNd.VectorND;

import java.util.List;

public class P2SpheroPolyhedron implements Potential2Soft {

    protected final Space space;
    protected final Vector dr;
    protected final Vector v, w;
    protected final Vector W0, W1, W2;
    
    protected final Vector sv, sqv, qsqv;
    
    protected final Vector dir0, dir1, dir2;
    protected final Vector norm0, norm1, norm2;
    
    public P2SpheroPolyhedron(Space space) {
        this.space = space;
        dr = space.makeVector();
        v = space.makeVector();
        W0 = space.makeVector();
        W1 = space.makeVector();
        W2 = space.makeVector();
        
        w = space.makeVector();
        sv = space.makeVector();
        sqv = space.makeVector();
        qsqv = space.makeVector();
        
        dir0 = space.makeVector();
        dir1 = space.makeVector();
        dir2 = space.makeVector();
        norm0 = space.makeVector();
        norm1 = space.makeVector();
        norm2 = space.makeVector();
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public double energy(IAtomList atoms) {
        AtomOrientedQuaternion atom0 = (AtomOrientedQuaternion)atoms.get(0);
        AtomOrientedQuaternion atom1 = (AtomOrientedQuaternion)atoms.get(1);
        AtomTypeSpheroPolyhedron atomType0 = (AtomTypeSpheroPolyhedron)atom0.getType();
        AtomTypeSpheroPolyhedron atomType1 = (AtomTypeSpheroPolyhedron)atom1.getType();
        Vector position0 = atom0.getPosition();
        Vector position1 = atom1.getPosition();
        
        dr.Ev1Mv2(position0, position1);
        double r2 = dr.squared();
        double din = atomType0.getInnerRadius() + atomType1.getInnerRadius();
        if (r2 < din*din) return Double.POSITIVE_INFINITY;

        double dout = atomType0.getOuterRadius() + atomType1.getOuterRadius();
        if (r2 > dout*dout) return 0;

        return gjke(atom0, atom1);
    }

    protected void inverseProduct4D(Vector q0, Vector q1, Vector q2) {
        // take inverse of q1, multiply by v and assign to v1
        double a, b, c;
        a = q1.getX(0) * q2.getX(0) + q1.getX(1) * q2.getX(1) + q1.getX(2) * q2.getX(2) + q1.getX(3) * q2.getX(3);
        b = q1.getX(0) * q2.getX(1) - q1.getX(1) * q2.getX(0) - q1.getX(2) * q2.getX(3) + q1.getX(3) * q2.getX(2);
        c = q1.getX(0) * q2.getX(2) + q1.getX(1) * q2.getX(3) - q1.getX(2) * q2.getX(0) - q1.getX(3) * q2.getX(1);
        q0.setX(3, q1.getX(0) * q2.getX(3) - q1.getX(1) * q2.getX(2) + q1.getX(2) * q2.getX(1) - q1.getX(3) * q2.getX(0));
        q0.setX(0, a);
        q0.setX(1, b);
        q0.setX(2, c);
    }
    
    protected void product3D(Vector v1, Vector q1, Vector vv) {
        // take inverse of q1, multiply by v and assign to v1
        double x1, y1, z1;
        x1 = q1.getX(0) * vv.getX(0) + q1.getX(2) * vv.getX(2) - q1.getX(3) * vv.getX(1);
        y1 = q1.getX(0) * vv.getX(1) + q1.getX(3) * vv.getX(0) - q1.getX(1) * vv.getX(2);
        z1 = q1.getX(0) * vv.getX(2) + q1.getX(1) * vv.getX(1) - q1.getX(2) * vv.getX(0);
        double vx, vy, vz;
        vx = vv.getX(0) + 2.0 * (q1.getX(2) * z1 - q1.getX(3) * y1);
        vy = vv.getX(1) + 2.0 * (q1.getX(3) * x1 - q1.getX(1) * z1);
        vz = vv.getX(2) + 2.0 * (q1.getX(1) * y1 - q1.getX(2) * x1);
        v1.setX(0, vx);
        v1.setX(1, vy);
        v1.setX(2, vz);
    }

    protected void inverseProduct3D(Vector v1, Vector q1, Vector vv) {
        // take inverse of q1, multiply by v and assign to v1
        double x1, y1, z1;
        /*System.out.println("q1 "+q1);
        System.out.println("vv "+vv); */
        x1 = q1.getX(0) * vv.getX(0) - q1.getX(2) * vv.getX(2) + q1.getX(3) * vv.getX(1);
        y1 = q1.getX(0) * vv.getX(1) - q1.getX(3) * vv.getX(0) + q1.getX(1) * vv.getX(2);
        z1 = q1.getX(0) * vv.getX(2) - q1.getX(1) * vv.getX(1) + q1.getX(2) * vv.getX(0);
        //System.out.println("xyz1 "+x1+" "+y1+" "+z1);
        double vx, vy, vz;
        vx = vv.getX(0) + 2.0 * (-q1.getX(2) * z1 + q1.getX(3) * y1);
        vy = vv.getX(1) + 2.0 * (-q1.getX(3) * x1 + q1.getX(1) * z1);
        vz = vv.getX(2) + 2.0 * (-q1.getX(1) * y1 + q1.getX(2) * x1);
        v1.setX(0, vx);
        v1.setX(1, vy);
        v1.setX(2, vz);
    }

    // Calculates the support in direction "direction"
    // The support is the vertex with the lowest distance from the origin
    //XXX farthest?
    // direction is replaced with the distance
    protected void calcSupport(AtomTypeSpheroPolyhedron atomType, Vector direction) {
        List<Vector> vertices = atomType.getVertices();
        if (vertices.size() == 0) {
            direction.E(0);
            return;
        }
        int index = 0;
        double distance = 0.0;
        for (int i = 0; i < vertices.size(); i++) {
            double d = vertices.get(i).dot(direction);
            if (d > distance) {
                index = i;
                distance = d;
            }
        }
        direction.E(vertices.get(index));
        return;
    }
    
    public double gjke(AtomOrientedQuaternion atom0, AtomOrientedQuaternion atom1) {
      //System.out.println("we're here");
        VectorND q = new VectorND(4);
        inverseProduct4D(q, atom0.getQuaternion(), atom1.getQuaternion());
        //System.out.println("q "+q);
        Vector t = space.makeVector();
        inverseProduct3D(t, atom0.getQuaternion(), dr);
        //System.out.println("t "+t);
        AtomTypeSpheroPolyhedron atomType0 = (AtomTypeSpheroPolyhedron)atom0.getType();
        AtomTypeSpheroPolyhedron atomType1 = (AtomTypeSpheroPolyhedron)atom1.getType();
        double sr = atomType0.getSweepRadius() + atomType1.getSweepRadius();
        //System.out.println("sr "+sr);
        
        v.E(t);
//        callcount++;
        
        for (int iteration=0; iteration<999; iteration++) {
            // Find a new vertex in direction v and store it in w
            // Vector3d w = support(-v) - q * shapeB.support(q.inverse() * v) + t;
            sv.Ea1Tv1(-1, v);
            calcSupport(atomType0, sv);
            inverseProduct3D(sqv, q, v);
            calcSupport(atomType1, sqv);
            w.E(sv);
            product3D(qsqv, q, sqv);
            w.ME(qsqv);
            w.PE(t);
            // Take care of rounding
            if (sr > 0.0) {
                w.PEa1Tv1(-sr/Math.sqrt(v.squared()), v);
            }
            /*System.out.println("iteration "+iteration);
            System.out.println(" v "+v);
            System.out.println(" w "+w);*/
            if (v.dot(w) >= 0) {
//                isum+=iteration;
                return 0;
            }
            
            /***** P O I N T *****/
            if (iteration == 0) {
                // Found: [0]
                v.E(w);
                W0.E(w);
            }
            
            /***** L I N E *****/
            else if (iteration == 1) {
                // Found: [0,1]
                dir0.Ev1Mv2(W0, w);
                v.Ea1Tv1(dir0.squared(), w);
                v.PEa1Tv1(-w.dot(dir0), dir0);
                W1.E(w);
            }
            
            /***** TRIANGLE *****/
            else if (iteration == 2) {
                // Found: [0,1,2]
                dir0.Ev1Mv2(W0, w);
                dir1.Ev1Mv2(W1, w);
                v.E(dir0);
                v.XE(dir1);
                if (v.dot(w) < 0) v.TE(-1);
                W2.E(w);
            }
            
            /***** T E T R A H E D R O N  *****/
            else {
                dir0.Ev1Mv2(W0, w);
                dir1.Ev1Mv2(W1, w);
                dir2.Ev1Mv2(W2, w);
                /*System.out.println("d0 "+dir0);
                System.out.println("d1 "+dir1);
                System.out.println("d2 "+dir2); */
                norm0.E(dir1);
                norm0.XE(dir2);
                boolean flipNorms = norm0.dot(dir0) < 0;
                if (flipNorms) {
                    norm0.TE(-1);
                }
                //System.out.println("n0 "+norm0);
                
                if (norm0.dot(w) > 0) {
                    // Found: [1,2,3]
                    v.E(norm0);
                    W0.E(w);
                    continue;
                }

                norm1.E(dir2);
                norm1.XE(dir0);
                if (flipNorms) {
                    norm1.TE(-1);
                }
                if (norm1.dot(w) > 0) {
                    // Found: [0,2,3]
                    v.E(norm1);
                    W1.E(w);
                    continue;
                }

                norm2.E(dir0);
                norm2.XE(dir1);
                if (flipNorms) {
                    norm2.TE(-1);
                }
                if (norm2.dot(w) > 0) {
                    // Found: [0,1,3]
                    v.E(norm2);
                    W2.E(w);
                    continue;
                }
                // vertex is inside: there is an overlap and we can return
//                isum+=iteration;
                return Double.POSITIVE_INFINITY;
            }
        }
        throw new RuntimeException("Infinite Loop Detected!");
    }

    @Override
    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        AtomOrientedQuaternion atom1q = (AtomOrientedQuaternion)atom1;
        AtomOrientedQuaternion atom2q = (AtomOrientedQuaternion)atom2;
        AtomTypeSpheroPolyhedron atomType1 = (AtomTypeSpheroPolyhedron)atom1q.getType();
        AtomTypeSpheroPolyhedron atomType2 = (AtomTypeSpheroPolyhedron)atom2q.getType();
        Vector position0 = atom1q.getPosition();
        Vector position1 = atom2q.getPosition();

        dr.Ev1Mv2(position0, position1);
        double r2 = dr.squared();
        double din = atomType1.getInnerRadius() + atomType2.getInnerRadius();
        if (r2 < din*din) return Double.POSITIVE_INFINITY;

        double dout = atomType1.getOuterRadius() + atomType2.getOuterRadius();
        if (r2 > dout*dout) return 0;

        return gjke(atom1q, atom2q);
    }

    @Override
    public double virial(IAtomList atoms) {
        return 0;
    }

    @Override
    public Vector[] gradient(IAtomList atoms) {
        return new Vector[0];
    }

    @Override
    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return new Vector[0];
    }
}
