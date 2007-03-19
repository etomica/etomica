package etomica.data.meter;

import etomica.integrator.IntegratorHard;

/**
 * Interface for a meter that can return a value based on a hard-potential collision.
 */
public interface DataSourceCollisional extends IntegratorHard.CollisionListener {
    public double collisionValue(IntegratorHard.Agent agent);
}