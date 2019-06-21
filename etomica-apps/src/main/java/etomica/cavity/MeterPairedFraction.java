package etomica.cavity;

import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorHard;
import etomica.units.dimensions.Null;

public class MeterPairedFraction extends DataSourceScalar implements IntegratorHard.CollisionListener {

    protected final IntegratorHard integratorHard;
    protected double lastTime;
    protected double lastSwitchTime;
    protected boolean internal;
    protected double tInternal, tExternal;

    public MeterPairedFraction(IntegratorHard integrator) {
        super("Contact ratio", Null.DIMENSION);
        integrator.addCollisionListener(this);
        integratorHard = integrator;
    }

    public void reset() {
        lastTime = integratorHard.getCurrentTime();
        lastSwitchTime = integratorHard.getCurrentTime();
        tInternal = 0;
        tExternal = 0;

    }

    public void collisionAction(IntegratorHard.Agent agent) {
        P2HardSphereCavity p2 = (P2HardSphereCavity) agent.collisionPotential;
        P2HardSphereCavity.CollisionType cType = p2.getLastCollisionType();
        if (cType == P2HardSphereCavity.CollisionType.CAPTURE) {
            internal = true;
            double t = integratorHard.getCurrentTime() + agent.collisionTime();
            tExternal += t - lastSwitchTime;
            lastSwitchTime = t;
        } else if (cType == P2HardSphereCavity.CollisionType.ESCAPE) {
            internal = false;
            double t = integratorHard.getCurrentTime() + agent.collisionTime();
            tInternal += t - lastSwitchTime;
            lastSwitchTime = t;
        }
    }

    @Override
    public double getDataAsScalar() {
        double ti = tInternal, te = tExternal;
        if (internal) ti += integratorHard.getCurrentTime() - lastSwitchTime;
        else te += integratorHard.getCurrentTime() - lastSwitchTime;
        double frac = ti / (te + ti);
        reset();
        return frac;
    }
}
