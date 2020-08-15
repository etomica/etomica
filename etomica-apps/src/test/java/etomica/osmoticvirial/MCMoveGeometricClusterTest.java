package etomica.osmoticvirial;

import etomica.action.activity.ActivityIntegrate2;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.liquidLJ.LjMC3D;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class MCMoveGeometricClusterTest {
    private LjMC3D ljMC3D;
    private MCMoveGeometricCluster mcMoveGeometricCluster;

    @BeforeEach
    public void setUp() throws Exception {
        ljMC3D = new LjMC3D(200,2,0.7,3);
        ljMC3D.getController2().runActivityBlocking(new ActivityIntegrate2(ljMC3D.integrator), 10000);
    }

    @Test
    public void energyChange() throws Exception {
        mcMoveGeometricCluster = new MCMoveGeometricCluster(ljMC3D.potentialMasterCell, ljMC3D.getSpace(),
                ljMC3D.getRandom(), ljMC3D.integrator, null);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(ljMC3D.potentialMasterCell);
        meterPE.setBox(ljMC3D.box);
        double oldE = meterPE.getDataAsScalar();
        mcMoveGeometricCluster.setBox(ljMC3D.box);
        mcMoveGeometricCluster.doTrial();
        double newE = meterPE.getDataAsScalar();
        double expected = newE-oldE;
        double actual = mcMoveGeometricCluster.energyChange();
        assertEquals(expected, actual, 2e-12);
    }
}