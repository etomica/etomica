package etomica.models.nitrogen;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.normalmode.MoleculeSiteSource;
import etomica.normalmode.MoleculeSiteSourceNitrogen;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;
import etomica.units.dimensions.Null;

import java.util.Arrays;

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
    protected final PotentialCalculationForceSum pcForceSum;
    protected final PotentialMaster potentialMaster;
    protected final AtomLeafAgentManager<MyAgent> forceManager;
    protected final IteratorDirective id;
    protected final Vector dr;
    protected double latticeEnergy;
    protected final double temperature;
    public static boolean justDADB = true;
    public static boolean justU = false;
    protected final Vector ph1h2;
    protected final Vector q;
    protected final Vector totalforce;
    private final Vector centerMass;
    public boolean doTranslation;
    public boolean doRotation;

    protected final Vector torque;

    public MeterDADBNitrogen(Space space, DataSourceScalar meterPE, PotentialMaster potentialMaster, double temperature, MoleculeAgentManager latticeCoordinates) {
        int nData = justDADB ? 1 : 9;
        //default is 1?
        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.space = space;
        this.latticeCoordinates = latticeCoordinates;
        this.meterPE = meterPE;
        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pcForceSum = new PotentialCalculationForceSum();
        forceManager = new AtomLeafAgentManager<MyAgent>(this, latticeCoordinates.getBox(), MyAgent.class);
        pcForceSum.setAgentManager(forceManager);
        dr = space.makeVector();
        MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(potentialMaster);
        meterPE2.setBox(latticeCoordinates.getBox());
        latticeEnergy = meterPE2.getDataAsScalar();
        this.temperature = temperature;
        ph1h2 = space.makeVector();
        q = space.makeVector();
        totalforce = space.makeVector();
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
        pcForceSum.reset();
        potentialMaster.calculate(box, id, pcForceSum);
        IMoleculeList molecules = box.getMoleculeList();
        double ForceSum = 0;
        double orientationSum = 0;


        for (int i = 0; i < molecules.getMoleculeCount(); i++) {
            IMolecule molecule = molecules.getMolecule(i);
            IAtomList leafList = molecule.getChildList();
            Vector n1Force = forceManager.getAgent(leafList.getAtom(0)).force();
            Vector n2Force = forceManager.getAgent(leafList.getAtom(1)).force();
            Vector p1lForce = forceManager.getAgent(leafList.getAtom(2)).force();
            Vector p2lForce = forceManager.getAgent(leafList.getAtom(3)).force();
            Vector p1rForce = forceManager.getAgent(leafList.getAtom(4)).force();
            Vector p2rForce = forceManager.getAgent(leafList.getAtom(5)).force();
            totalforce.E(n1Force);
            totalforce.PE(n2Force);
            totalforce.PE(p1lForce);
            totalforce.PE(p2lForce);
            totalforce.PE(p1rForce);
            totalforce.PE(p2rForce);
            Vector n1 = leafList.getAtom(0).getPosition();
            Vector n2 = leafList.getAtom(1).getPosition();
            Vector p1l = leafList.getAtom(2).getPosition();
            Vector p2l = leafList.getAtom(3).getPosition();
            Vector p1r = leafList.getAtom(4).getPosition();
            Vector p2r = leafList.getAtom(5).getPosition();

            Vector axis = space.makeVector();
            Vector nn = space.makeVector();
            nn.E(n1);
            nn.ME(n2);

            centerMass.E(n1);
            centerMass.PE(n2);
            centerMass.TE(0.5);

            dr.Ev1Mv2(n1, centerMass);
            torque.E(n1Force);
            torque.XE(dr);
            q.E(torque);
            dr.Ev1Mv2(n2, centerMass);
            torque.E(n2Force);
            torque.XE(dr);
            q.PE(torque);
            dr.Ev1Mv2(p1l, centerMass);
            torque.E(p1lForce);
            torque.XE(dr);
            q.PE(torque);
            dr.Ev1Mv2(p2l, centerMass);
            torque.E(p2lForce);
            torque.XE(dr);
            q.PE(torque);
            dr.Ev1Mv2(p1r, centerMass);
            torque.E(p1rForce);
            torque.XE(dr);
            q.PE(torque);
            dr.Ev1Mv2(p2r, centerMass);
            torque.E(p2rForce);
            torque.XE(dr);
            q.PE(torque);
            //for the total torque q

            Vector lPos = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).position;
            dr.Ev1Mv2(centerMass, lPos);
            ForceSum += totalforce.dot(dr);
            //get the forcesum!!


            Orientation3D or = ((MoleculeSiteSourceNitrogen.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).orientation;
            Vector a0 = or.getDirection();//NN

            //TODO define vetor NN for melucle orientation


            double theta = Math.acos(Math.abs(nn.dot(a0)));

            axis.E(a0);
            axis.ME(nn);
            double DUDT = q.dot(axis);

            orientationSum += Math.sin(theta) * DUDT;//TODO careful about the sign

        }


        if (!doTranslation) ForceSum = 0;
        if (!doRotation) orientationSum = 0;
        if (justDADB) {
            if (justU) {
                int N = molecules.getMoleculeCount();
                double fac = (doTranslation ? 1.5 : 0) * (N - 1) + (doRotation ? 1.5 : 0) * N;
                x[0] = (x0 + latticeEnergy) + (fac * temperature) + 0.5 * ForceSum + orientationSum;
            } else {
//                System.out.println(x0+" "+(0.5*sum)+" "+(x0+0.5*sum)/atoms.getAtomCount());
                x[0] = x0 + 0.5 * ForceSum + orientationSum; //translation and rotation
//                x[0] = x0 + 0.5*ForceSum ;//only translation
//                x[0] = x0 + orientationSum;//only rotation
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

    private void Ev1Mv2(Vector a0, double d) {
        // TODO Auto-generated method stub

    }

    public void debug() {
        Box box = latticeCoordinates.getBox();
        IMoleculeList molecules = box.getMoleculeList();
        for (int j = 99; j > 0; j--) {
            for (int i = 0; i < molecules.getMoleculeCount(); i++) {
                IMolecule molecule = molecules.getMolecule(i);
                IAtomList leafList = molecule.getChildList();
                Vector h1 = leafList.getAtom(0).getPosition();
                Vector h2 = leafList.getAtom(1).getPosition();
                Vector o = leafList.getAtom(2).getPosition();
                Vector m = leafList.getAtom(3).getPosition();
                OrientationFull3D or = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).orientation;
                Vector a0 = or.getDirection();//om
                Vector a1 = or.getSecondaryDirection();//h1h2
//                dr.Ev1Mv2(m, o);
//                dr.normalize();
//                System.out.println("om = " + dr);
//                System.out.println("a0 =  " + a0);
//                System.out.println("a1 = " + a1);
//                if(i == 3){System.exit(2);}

                Vector om = space.makeVector();
                Vector hh = space.makeVector();
                Vector axis = space.makeVector();
                Vector a2 = space.makeVector();

                a2.E(a0);
                a2.XE(a1);
                a2.normalize();
                double[][] array = new double[3][3];
                a0.assignTo(array[0]);
                a1.assignTo(array[1]);
                a2.assignTo(array[2]);
                Matrix a = new Matrix(array).transpose();
                om.Ev1Mv2(m, o);
                om.normalize();
                hh.Ev1Mv2(h2, h1);
                //try to test shake tolerance effect the purpendicular or not
//                Vector p = space.makeVector();
//                p.E(om);
//                p.TE(hh.dot(om));
//                hh.ME(p);


                Vector p = space.makeVector();
                p.E(om);
                p.TE(hh.dot(om));
                hh.ME(p);


                hh.normalize();
                a2.E(om);
                a2.XE(hh);
                a2.normalize();

                double[][] array1 = new double[3][3];
                om.assignTo(array1[0]);
                hh.assignTo(array1[1]);
                a2.assignTo(array1[2]);
                Matrix newa = new Matrix(array1).transpose();
                a = a.inverse();
                Matrix matrix = newa.times(a);
                double beta = 0;

                EigenvalueDecomposition eigenvalueDecomposition = matrix.eig();
                Matrix eigenvectormatrix = eigenvalueDecomposition.getV();
                double[] eigenValueArray = eigenvalueDecomposition.getRealEigenvalues();
                double[][] eigenVectors = eigenvectormatrix.transpose().getArrayCopy();

                int best = 0;
                double value = Math.abs(eigenValueArray[0] - 1);
                if (value > Math.abs(eigenValueArray[1] - 1)) {
                    best = 1;
                    value = Math.abs(eigenValueArray[1] - 1);
                }
                if (value > Math.abs(eigenValueArray[2] - 1)) {
                    best = 2;
                    value = Math.abs(eigenValueArray[2] - 1);
                }


                axis.E(eigenVectors[best]);


                om.Ev1Mv2(m, o);
                hh.Ev1Mv2(h2, h1);
                if (om.dot(a0) > hh.dot(a0)) {
                    dr.E(om);
                    dr.XE(a0);
                } else {
                    dr.E(hh);
                    dr.XE(a1);
                }

                beta = Math.signum(axis.dot(dr)) * Math.acos((matrix.trace() - 1) / 2.0);//rotation angle


//                dr.Ev1Mv2(m, o);
//    			ph1h2.Ev1Mv2(h2, h1);
//    			System.out.println("old  =" + dr.dot(ph1h2));
                double hmass = leafList.getAtom(0).getType().getMass();
                double omass = leafList.getAtom(2).getType().getMass();
                centerMass.Ea1Tv1(hmass, h1);
                centerMass.PEa1Tv1(hmass, h2);
                centerMass.PEa1Tv1(omass, o);
                centerMass.TE(1 / (2 * hmass + omass));

//              h1.ME(o);
//    			double lenth = Math.sqrt(h1.squared());
//    			h1.normalize();
//    			Orientation3D orientation = new Orientation3D(space);
//    			orientation.setDirection(h1);
//    			orientation.rotateBy(-beta/(j+1), a0);
//    			h1.Ea1Tv1(lenth, orientation.getDirection());
//    			h1.PE(o);
//    			h2.ME(o);
//    		    h2.normalize();
//    		    orientation.setDirection(h2);
//    			orientation.rotateBy(-beta/(j+1), a0);
//    			h2.Ea1Tv1(lenth, orientation.getDirection());
//    			h2.PE(o);
                h1.ME(centerMass);
                h2.ME(centerMass);
                o.ME(centerMass);
                m.ME(centerMass);
                double hlength = Math.sqrt(h1.squared());
                double olength = Math.sqrt(o.squared());
                double mlength = Math.sqrt(m.squared());
                h1.normalize();
                h2.normalize();
                o.normalize();
                m.normalize();
                Orientation3D orientation = new Orientation3D(space);
                orientation.setDirection(h1);
                orientation.rotateBy(beta / (j + 1), axis);
                h1.Ea1Tv1(hlength, orientation.getDirection());
                h1.PE(centerMass);
                orientation.setDirection(h2);
                orientation.rotateBy(beta / (j + 1), axis);
                h2.Ea1Tv1(hlength, orientation.getDirection());
                h2.PE(centerMass);
                orientation.setDirection(o);
                orientation.rotateBy(beta / (j + 1), axis);
                o.Ea1Tv1(olength, orientation.getDirection());
                o.PE(centerMass);
                orientation.setDirection(m);
                orientation.rotateBy(beta / (j + 1), axis);
                m.Ea1Tv1(mlength, orientation.getDirection());
                m.PE(centerMass);
//    			
//    			dr.Ev1Mv2(m, o);
//    			ph1h2.Ev1Mv2(h2, h1);
//    			System.out.println("new =" +dr.dot(ph1h2));


//    			if(i == 0){
//    				System.out.println(beta);
//    			}
            }
            System.out.print(j * 0.01 + " ");
            double y = getData().getValue(0);
//        	System.out.println( j*0.01 + " " + y );  

        }

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
