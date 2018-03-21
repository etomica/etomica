package etomica.models.nitrogen;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorRigidIterative;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.normalmode.MoleculeSiteSource;
import etomica.normalmode.MoleculeSiteSourceNitrogen;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;
import etomica.units.dimensions.Null;

/**
 * @author Weisong Lin
 */
//TODO I replace IEtomicaDataSource to IDATASource here
public class MeterDADBNitrogen implements IDataSource, AgentSource<MyAgent> {

    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final MoleculeAgentManager latticeCoordinates;
    protected final DataSourceScalar meterPE;
    protected final PotentialCalculationTorqueSum pcTorqueSum;
    protected final PotentialMaster potentialMaster;
    protected final MoleculeAgentManager forceManager;
    protected final IteratorDirective id;
    protected final Vector dr;
    protected double latticeEnergy;
    protected final double temperature;
    public static boolean justDADB = true;
    public static boolean justU = false;
    protected final Vector ph1h2;
    protected final Vector q;
    protected final Vector totalForce;
    private final Vector centerMass;
    public boolean doTranslation;
    public boolean doRotation;

    protected final Vector torque;

    public MeterDADBNitrogen(Simulation sim, DataSourceScalar meterPE, PotentialMaster potentialMaster, double temperature, MoleculeAgentManager latticeCoordinates) {
        int nData = justDADB ? 1 : 9;
        //default is 1?
        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.space = sim.getSpace();
        this.latticeCoordinates = latticeCoordinates;
        this.meterPE = meterPE;
        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pcTorqueSum = new PotentialCalculationTorqueSum();
        MoleculeAgentManager.MoleculeAgentSource molAgentSource = new MoleculeAgentManager.MoleculeAgentSource() {

            public void releaseAgent(Object agent, IMolecule molecule) {
            }

            public Object makeAgent(IMolecule mol) {
                return new IntegratorRigidIterative.MoleculeAgent(space);
            }

            public Class getMoleculeAgentClass() {
                return IntegratorRigidIterative.MoleculeAgent.class;
            }
        };
        forceManager = new MoleculeAgentManager(sim, latticeCoordinates.getBox(), molAgentSource);
        pcTorqueSum.setMoleculeAgentManager(forceManager);
        dr = space.makeVector();
        MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(potentialMaster);
        meterPE2.setBox(latticeCoordinates.getBox());
        latticeEnergy = meterPE2.getDataAsScalar();
        this.temperature = temperature;
        ph1h2 = space.makeVector();
        q = space.makeVector();
        totalForce = space.makeVector();
        centerMass = space.makeVector();
        torque = space.makeVector();
    }

    public void setLatticeEnergy(double newLatticeEnergy) {
        latticeEnergy = newLatticeEnergy;
    }

    public IData getData() {
        Box box = latticeCoordinates.getBox();
        double[] x = data.getData();
        double x0 = meterPE.getDataAsScalar() - latticeEnergy;
        pcTorqueSum.reset();
        potentialMaster.calculate(box, id, pcTorqueSum);
        IMoleculeList molecules = box.getMoleculeList();
        double ForceSum = 0;
        double orientationSum = 0;


        for (int i = 0; i < molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);

            totalForce.E(((Integrator.Forcible) forceManager.getAgent(molecule)).force());
            Vector n1 = molecule.getChildList().getAtom(0).getPosition();
            Vector n2 = molecule.getChildList().getAtom(1).getPosition();

            centerMass.E(n1);
            centerMass.PE(n2);
            centerMass.TE(0.5);


            Vector axis = space.makeVector();
            Vector nn = space.makeVector();
            nn.E(n1);
            nn.ME(n2);
            nn.normalize();

            if (doTranslation) {
                q.E(((Integrator.Torquable) forceManager.getAgent(molecule)).torque());
                Vector lPos = ((MoleculeSiteSourceNitrogen.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).position;
                dr.Ev1Mv2(centerMass, lPos);
                ForceSum += totalForce.dot(dr);
//                System.out.println("doTranslation");
                //get the forcesum!!
            }
            if (doRotation) {
                Orientation3D or = ((MoleculeSiteSourceNitrogen.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).orientation;
                Vector a0 = or.getDirection();//NN
                double theta = Math.acos(Math.abs(nn.dot(a0)));
                axis.E(a0);
                axis.XE(nn);
                axis.normalize();
                double DUDT = q.dot(axis);
                orientationSum += (1 - Math.cos(theta)) / Math.sin(theta) * DUDT;//OrientationSum
//                System.out.println("doRotation");
            }
        }


        if (!doTranslation) ForceSum = 0;
        if (!doRotation) orientationSum = 0;
        if (justDADB) {
            if (justU) {
                int N = molecules.getMoleculeCount();
                double fac = (doTranslation ? 1.5 : 0) * (N - 1) + (doRotation ? 1.5 : 0) * N;
                x[0] = (x0 + latticeEnergy) + (fac * temperature) + 0.5 * ForceSum + orientationSum;
            } else {
                x[0] = x0 + 0.5 * ForceSum + orientationSum;
            }

            if (data.isNaN()) {
                throw new RuntimeException();
            }
            return data;
        }
        x[0] = x0;
        x[1] = ForceSum;
        x[2] = x[0] + 0.5 * x[1];
        // so this is the energy
        x[3] = (2 * temperature - x[0]) * x[0];   // x[0]*x[0] = 2*T*x[0] - x[3]
        x[4] = (2 * temperature - x[0]) * x[2];
//        x[5] = (x[0]*x[0]+2*temperature*temperature)*x[2];
//        x[5] = (-6*temperature*(temperature - x[0]) - x[0]*x[0])*x[2];
        x[5] = x[0] * x[0];
        x[6] = x[0] * x[0] * x[2];
        //what does x[3-7] used for?

//        x[7] = x[0]*x[1];
        //dAc/dT = -(uAvg+0.5*fdrAvg))/ T^2
        //       = -<x2>/T^2
        //d2Ac/dT2 = ( (2*temperature*(uAvg + 0.5*fdrAvg) - u2Avg + uAvg*uAvg + 0.5*(uAvg*fdrAvg - ufdrAvg) ) / T^4
        //         = (<x4> + <x0>*<x2>)/T^4;
        //d3Ac/dT3 = ( -6*T*T*(uAvg+0.5fdrAvg) + 6*T*(u2Avg - uAvg + 0.5*(ufdrAvg - uAvg*fdrAvg)) - uuu# - 0.5*uuf#)
        //  xyz# = <xyz> - <x><yz> - <y><xz> - <z><xy> + 2<x><y><z>
        //d3Ac/dT3 = (<x5> - 6T*<x0><x2> + 3*<u2>*<x2> + 2*<x0>*(-<x0>*<x2> + 0.5*<uF>) ) / T^6
        //d3Ac/dT3 = (<x5> - 6T*<x0><x2> + 3*<x6>*<x2> + 2*<x0>*(-<x0>*<x2> + 0.5*<x7>) ) / T^6
        return data;
    }


    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public final MyAgent makeAgent(IAtom a, Box box) {
        return new MyAgent(space);
    }

    public void releaseAgent(MyAgent agent, IAtom atom, Box box) {
    }
}
