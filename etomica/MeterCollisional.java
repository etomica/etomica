/*
 * History
 * Created on Jul 25, 2004 from MeterAbstract.Collisional
 */
package etomica;

/**
 * Interface for a meter that can return a value based on a hard-potential collision.
 */
public interface MeterCollisional extends IntegratorHard.CollisionListener {
    public double collisionValue(IntegratorHard.Agent agent);
}