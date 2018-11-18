package etomica.spin.heisenberg_interacting.heisenberg;

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
    public final Tensor phi;
    public final Vector force;
    public final Vector vE;
    public final Vector vEE;

    public MoleculeAgent() {
        phi = new Tensor1D();
        vEE = new Vector2D();
        force = new Vector2D();
        torque = new Vector1D();
        vE = new Vector2D();
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

    public Vector vE() {
        return vE;
    }

    public Vector vEE() {
        return vEE;
    }
}
