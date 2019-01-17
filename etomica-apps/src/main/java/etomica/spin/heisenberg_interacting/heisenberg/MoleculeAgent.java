package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.integrator.Integrator;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space1d.Tensor1D;
import etomica.space1d.Vector1D;
import etomica.space2d.Tensor2D;
import etomica.space2d.Vector2D;
import etomica.spaceNd.TensorND;

/**
 * Holds fields needed by the meter that calculates dielectric constant using mapped averaging.
 */
public class MoleculeAgent implements Integrator.Torquable, Integrator.Forcible {
    public final Vector torque;
    public final Vector btPhi;
    public final Vector force;
    public final Vector vE;
    public final Vector vEE;
    public final Tensor phi;

    //dvEx means dvEx/dtheta and d2vEx means d2vEx/dtheta^2 etc.
    public final Vector vEx, vEEx, dvEx, dvEEx, d2vEx;
    public final Vector vEy, vEEy, dvEy, dvEEy, d2vEy;
    public final Vector fEx, fEy, btForce;


    public MoleculeAgent() {
        phi = new Tensor1D();
        vEE = new Vector2D();
        force = new Vector2D();
        torque = new Vector1D();
        vE = new Vector2D();

        btPhi = new Vector1D();
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
        fEx = new Vector1D();
        fEy = new Vector1D();
        btForce = new Vector1D();
    }

    public Vector torque() {
        return torque;
    }

    public Tensor phi() {
        return phi;
    }

    public Vector btPhi() {
        return btPhi;
    }

    public Vector force() {
        return force;
    }

    public Vector btForce() {
        return btForce;
    }

    public Vector fEx() {
        return fEx;
    }

    public Vector fEy() {
        return fEy;
    }

    public Vector vE() {
        return vE;
    }

    public Vector vEE() {
        return vEE;
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


}
