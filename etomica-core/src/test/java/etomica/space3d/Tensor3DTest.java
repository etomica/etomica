/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space3d;

import etomica.space.Vector;
import etomica.space.Tensor;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;

import static junit.framework.TestCase.assertTrue;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

/**
 * Created by kofke on 5/7/17.
 */
public class Tensor3DTest {

    private Tensor zero, t1, t2, t3;
    private int dim;
    private double tolerance = 1e-10;

    @Before
    public void setUp() throws Exception {
        zero = new Tensor3D();
        zero.E(0.0);
        t1 = new Tensor3D(new double[][] {{1,2,3},{4,-5,6},{7,8,9}});
        t2 = new Tensor3D();
        t3 = new Tensor3D();
        dim = t1.D();
    }

    @Test
    public void testEquals() throws Exception {
        t2 = (Tensor)t1.clone();
        assertTrue(t2.equals(t1));
        assertFalse(t2.equals(zero));
        t2.setComponent(1,2,-50.);//check that clone made a deep copy
        assertFalse(t2.equals(t1));
    }

    //tests component, setComponent, component PE
    @Test
    public void testSetComponent() throws Exception {
        t2.E(0.0);
        t3.E(0.0);
        double k = 0.;
        for(int i=0; i<dim; i++) {
            for(int j=0; j<dim; j++) {
                t2.setComponent(i,j,4.1);
                assertEquals(t2.component(i,j),4.1,tolerance);
                t2.E(0.0);
                k++;
                t3.PE(i,j,k);
            }
        }
        t3.PE(1,1,-10.);
        assertTrue(t3.equals(t1));
    }

    // tests E(Vector[]) and assignTo(Vector[])
    @Test
    public void testE() throws Exception {
        Vector3D v0 = new Vector3D(1.0,4.0,7.0);
        Vector3D v1 = new Vector3D(2.0,-5.0,8.0);
        Vector3D v2 = new Vector3D(3.0,6.0,9.0);
        t2.E(new Vector[] {v0, v1, v2});
        assertTrue(t2.equals(t1));

        Vector3D[] vecs = new Vector3D[] {new Vector3D(),new Vector3D(), new Vector3D()};
        t2.assignTo(vecs);
        assertTrue(vecs[0].equals(v0));
        assertTrue(vecs[1].equals(v1));
        assertTrue(vecs[2].equals(v2));
    }

    //tests Ev1v2, E, ME, PEv1v2
    @Test
    public void testEv1v2() throws Exception {
        Vector3D v1 = new Vector3D(new double[] {1.0, 1.5, 2.0});
        Vector3D v2 = new Vector3D(new double[] {-3.0, 4.0, 1.7});
        t3.Ev1v2(v1,v2);
        double[][] d = {{-3.0,4.0,1.7},{-4.5,6.0,2.55},{-6,8.0,3.4}};
        t2.E(d);
        assertTrue(t3.equals(t2));
        t3.ME(t3);
        t3.PEv1v2(v1,v2);
        assertTrue(t3.equals(t2));
    }

    @Test
    public void testTrace() throws Exception {
        assertEquals(t1.trace(),5.0,tolerance);
    }

    @Test
    public void testTranspose() throws Exception {
        t3.E(t1);
        t1.transpose();
        double[][] d = {{1,4,7},{2,-5,8},{3,6,9}};
        t2.E(d);
        assertTrue(t1.equals(t2));
        t1.transpose();
        assertTrue(t1.equals(t3));
    }

    @Test
    public void testDeterminant() throws Exception {
        assertEquals(t1.determinant(),120.,tolerance);
    }

    @Test
    public void testTransform() throws Exception {
        Vector3D v1 = new Vector3D(new double[] {1.0, 1.5, 2.0});
        Vector3D v2 = new Vector3D(new double[] {10.,8.5,37.0});
        t1.transform(v1);
        assertTrue(v1.equals(v2));
    }

    @Test
    public void testInvert() throws Exception {
        t2 = new Tensor3D(new double[][] {{1.,1.,1.},{4.,5.,6.},{2.,6.,4.}});
        assertEquals(t2.determinant(),-6.0,tolerance);
        t2.invert();
        t2.TE(6.0);
        t3 = new Tensor3D(new double[][] {{16.,-2.,-1.},{4.,-2.,2.},{-14.,4.,-1.}});
        assertTrue(t2.equals(t3));

        t2.E(t1);
        t2.invert();
        t2.invert();
        t2.map(x -> Math.round(x));
        assertTrue(t2.equals(t1));
    }

    @Test
    public void testMap() throws Exception {
        t2.E(0);
        t2.PEa1Tt1(2.0,t1);
        t3.E(t1);
        t3.map(x -> 2.0*x);
        assertTrue(t2.equals(t3));
    }

    @Test
    public void testMEv1v2() throws Exception {
        Vector3D v1 = new Vector3D(new double[] {1.0, 1.5, 2.0});
        Vector3D v2 = new Vector3D(new double[] {10.,8.5,37.0});
        t2.E(0.0);
        t2.PEv1v2(v1,v2);
        t2.MEv1v2(v1,v2);
        assertTrue(t2.equals(zero));
    }

    //tests E, PE, TE
    @Test
    public void testTE() throws Exception {
        t2.E(t1);
        t2.PE(t1);
        t3.E(t1);
        t3.TE(2.0);
        assertTrue(t2.equals(t3));

        t2.E(0.0);
        t2.PE(5.0);
        t3.E(5.0);
        assertTrue(t2.equals(t3));

        t2.E(t1);
        t2.transpose();
        t2.PE(-1.0);
        t2.TE(t1);//matrix multiplication
        double[][] d0 = new double[][] {{54, 33, 72}, {26, 88, 30}, {78, 43, 108}};
        t3.E(d0);
        assertTrue(t2.equals(t3));

        double[] d1 = new double[] {54, 33, 72, 26, 88, 30, 78, 43, 108};
        double[] d2 = new double[dim*dim];
        ((Tensor3D)t3).E(d1);
        assertTrue(t2.equals(t3));
        t3.assignTo(d2);
        assertTrue(Arrays.equals(d1,d2));

        double[][] d3 = new double[dim][dim];
        t3.assignTo(d3);
        for(int i=0; i<dim; i++) {
            assertTrue(Arrays.equals(d0[i], d3[i]));
        }

        d2 = t3.toArray();
        assertTrue(Arrays.equals(d1,d2));
    }

    @Test
    public void testDiagE() throws Exception {
        t2 = new Tensor3D(new double[][] {{2.0,0.0,0.0},{0.0,-1.5,0.0},{0.0,0.0,3.0}});
        t3.diagE(new Vector3D(2.0,-1.5,3.0));
        assertTrue(t2.equals(t3));
    }

}
