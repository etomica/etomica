package etomica.modules.rosmosis;

import etomica.api.IAtom;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IIntegratorNonintervalListener;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.MoleculeAgentManager;
import etomica.atom.MoleculeAgentManager.MoleculeAgentSource;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.space.ISpace;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Pressure;
import etomica.units.Quantity;
import etomica.units.Time;

/**
 * Meter to measure flux across a boundary or boundaries.  If an atom is on one
 * side of the boundary at one time and on the other side of the boundary the
 * next time, then the meter counts that as a crossing.  Each boundary can have
 * a coefficient, such that flow into one region through boundaries on either
 * side can be considered "positive" flux.
 * 
 * If an IntegratorMD is used, flux will be given in terms of crossings per
 * area per time.  Otherwise, flux will be given in terms of crossings per
 * area per step.
 *
 * @author Andrew Schultz
 */
public class MeterOsmoticPressure implements DataSource {

    public MeterOsmoticPressure(PotentialCalculationTorqueSumWallForce pc, IBox box) {
        this.pc = pc;
        this.box = box;
        data = new DataDouble();
        dataInfo = new DataInfoDouble("osmotic pressure", Pressure.DIMENSION);
        tag = new DataTag();
    }
    
    public Data getData() {
        IVector dimensions = box.getBoundary().getDimensions();
        data.x = pc.getWallForce() / (dimensions.x(1) * dimensions.x(2));
        return data;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected final PotentialCalculationTorqueSumWallForce pc;
    protected final IBox box;
    protected DataDouble data;
    protected DataInfoDouble dataInfo;
    protected DataTag tag;
}
