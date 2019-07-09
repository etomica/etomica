package etomica.cavity;

import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorHard;
import etomica.units.dimensions.Null;

public class MeterContactRatio extends DataSourceScalar implements IntegratorHard.CollisionListener {

    protected long totalCollisions, internalCollisions;
    protected double sigma;

    public MeterContactRatio(IntegratorHard integrator) {
        super("Contact ratio", Null.DIMENSION);
        integrator.addCollisionListener(this);
    }

    public void reset() {
        totalCollisions = internalCollisions = 0;
    }

    public void collisionAction(IntegratorHard.Agent agent) {
        totalCollisions++;
        P2HardSphereCavity p2 = (P2HardSphereCavity) agent.collisionPotential;
        P2HardSphereCavity.CollisionType cType = p2.getLastCollisionType();
        sigma = p2.getCollisionDiameter();
        if (cType == P2HardSphereCavity.CollisionType.ESCAPE) {
            internalCollisions++;
        } else if (cType == P2HardSphereCavity.CollisionType.INTERNAL_BOUNCE) {
            internalCollisions++;
        }
    }

    @Override
    public double getDataAsScalar() {
        double fac = 1;
        if (internalCollisions > 0) {
            long externalCollision = totalCollisions - internalCollisions;
            fac = externalCollision / (double) internalCollisions;
        }
        reset();
        return fac;
    }
}
