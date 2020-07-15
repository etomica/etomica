package etomica.spin.heisenberg;

import etomica.integrator.Integrator;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space1d.Tensor1D;
import etomica.space1d.Vector1D;
import etomica.space2d.Vector2D;

/**
 * Holds fields needed by the meter that calculates dielectric constant using mapped averaging.
 */
public class MoleculeAgent implements Integrator.Torquable, Integrator.Forcible {
    public final Vector torque;
    public final Vector force;
    public final Vector vE;
    public final Vector vEE;
    public final Tensor phi;

    //dvEx means dvEx/dtheta and d2vEx means d2vEx/dtheta^2 etc.
    public final Vector vEx, vEEx, dvEx, dvEEx, d2vEx;
    public final Vector vEy, vEEy, dvEy, dvEEy, d2vEy;

    double[] Axc0, Axs0, dAxc0, dAxs0, Axc1, Axs1, dAxc1, dAxs1;
    double[] d2Axc0, d2Axs0, d3Axc0, d3Axs0, d2Axc1, d2Axs1;
    double[] Ayc0, Ays0, dAyc0, dAys0, Ayc1, Ays1, dAyc1, dAys1;
    double[] d2Ayc0, d2Ays0, d3Ayc0, d3Ays0, d2Ayc1, d2Ays1;
    double[] sinntheta, cosntheta;


    public MoleculeAgent(int nMax) {
        phi = new Tensor1D();
        vEE = new Vector2D();
        force = new Vector2D();
        torque = new Vector1D();
        vE = new Vector2D();

        vEx = new Vector1D();
        vEEx = new Vector1D();
        dvEx = new Vector1D();
        dvEEx = new Vector1D();
        d2vEx = new Vector1D();
        vEy = new Vector1D();
        vEEy = new Vector1D();
        dvEy = new Vector1D();
        dvEEy = new Vector1D();
        d2vEy = new Vector1D();

        sinntheta = new double[nMax + 1];
        cosntheta = new double[nMax + 1];


        Axc0 = new double[nMax + 1];
        Axs0 = new double[nMax + 1];
        dAxc0 = new double[nMax + 1];
        dAxs0 = new double[nMax + 1];
        Axc1 = new double[nMax + 1];
        Axs1 = new double[nMax + 1];
        dAxc1 = new double[nMax + 1];
        dAxs1 = new double[nMax + 1];
        d2Axc0 = new double[nMax + 1];
        d2Axs0 = new double[nMax + 1];
        d3Axc0 = new double[nMax + 1];
        d3Axs0 = new double[nMax + 1];
        d2Axc1 = new double[nMax + 1];
        d2Axs1 = new double[nMax + 1];
        Ayc0 = new double[nMax + 1];
        Ays0 = new double[nMax + 1];
        dAyc0 = new double[nMax + 1];
        dAys0 = new double[nMax + 1];
        Ayc1 = new double[nMax + 1];
        Ays1 = new double[nMax + 1];
        dAyc1 = new double[nMax + 1];
        dAys1 = new double[nMax + 1];
        d2Ayc0 = new double[nMax + 1];
        d2Ays0 = new double[nMax + 1];
        d3Ayc0 = new double[nMax + 1];
        d3Ays0 = new double[nMax + 1];
        d2Ayc1 = new double[nMax + 1];
        d2Ays1 = new double[nMax + 1];


    }

    public Vector torque() {
        return torque;
    }


    public Tensor phi() {
        return phi;
    }


    public Vector force() {
        return force;
    }


    public Vector vEx() {
        return vEx;
    }

    public Vector vEEx() {
        return vEEx;
    }

    public Vector dvEx() {
        return dvEx;
    }

    public Vector dvEEx() {
        return dvEEx;
    }

    public Vector d2vEx() {
        return d2vEx;
    }

    public Vector vEy() {
        return vEy;
    }

    public Vector vEEy() {
        return vEEy;
    }

    public Vector dvEy() {
        return dvEy;
    }

    public Vector dvEEy() {
        return dvEEy;
    }

    public Vector d2vEy() {
        return d2vEy;
    }

    public void zeroSum() {

        vEx().E(0);
        vEy().E(0);
        vEEx().E(0);
        vEEy().E(0);
        dvEx().E(0);
        dvEy().E(0);
        dvEEx().E(0);
        dvEEy().E(0);
        d2vEx().E(0);
        d2vEy().E(0);
    }


}
