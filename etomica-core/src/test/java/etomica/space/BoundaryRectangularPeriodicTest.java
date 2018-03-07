package etomica.space;

import etomica.space3d.Space3D;
import org.junit.Before;
import org.junit.Test;

import static etomica.UnitTestUtil.DELTA;
import static org.junit.Assert.*;

public class BoundaryRectangularPeriodicTest {
    Boundary b;

    @Before
    public void setUp() {
        b = new BoundaryRectangularPeriodic(Space3D.getInstance(), 30);
    }

    @Test
    public void testNearestImage() {
        Vector v = Vector.of(3.5, 46, -12.4);
        b.nearestImage(v);
        assertArrayEquals(new double[]{3.5, -14.0, -12.4}, v.toArray(), DELTA);
    }

    @Test
    public void testCentralImage() {
        Vector v = Vector.of(3.5, 46, -12.4);
        Vector centralImage = b.centralImage(v);
        assertArrayEquals(new double[]{0, -60.0, 0}, centralImage.toArray(), DELTA);
    }
}
