package etomica.models.clathrates;
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
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.normalmode.AtomSiteSource;
import etomica.normalmode.LatticeSumMolecularCrystal;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.RotationTensor3D;
import etomica.units.dimensions.Null;

import java.util.Arrays;
public class MeterDADBWaterTIP4P  implements IDataSource {


    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final MoleculeAgentManager latticeCoordinates;
    protected final DataSourceScalar meterPE;
    protected final PotentialCalculationForceSum pcForceSum;
    protected final PotentialMaster potentialMaster;
    protected final AtomLeafAgentManager<Vector> forceManager;
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

    public MeterDADBWaterTIP4P(Space space, DataSourceScalar meterPE, PotentialMaster potentialMaster, double temperature, MoleculeAgentManager latticeCoordinates) {
        int nData = justDADB ? 1 : 9;
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
        forceManager = new AtomLeafAgentManager<>(a -> space.makeVector(), latticeCoordinates.getBox());
        pcForceSum.setAgentManager(forceManager);
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
        pcForceSum.reset();
        potentialMaster.calculate(box, id, pcForceSum);
        IMoleculeList molecules = box.getMoleculeList();
        double ForceSum = 0;
        double orientationSum = 0;
        AtomLeafAgentManager<Vector> atomLatticeCoordinates = new AtomLeafAgentManager<>(new AtomSiteSource(space), box, Vector.class);

        for (int i = 0; i < molecules.size(); i++) {    //#########GET GETATOM/MOLECULE WORKING
            IMolecule molecule = molecules.get(i);
            IAtomList leafList = molecule.getChildList();

            Vector h1Force = agentManager.getAgent(leafList.get(0));
            Vector h2Force = agentManager.getAgent(leafList.get(1));
            Vector oForce = agentManager.getAgent(leafList.get(2));
            Vector mForce = agentManager.getAgent(leafList.get(3));

            totalForce.E(h1Force); //########################WHAT IS E,PE,TE?
            totalForce.PE(h2Force);
            totalForce.PE(mForce);
            totalForce.PE(oForce);
            Vector h1 = leafList.get(0).getPosition();
            Vector h2 = leafList.get(1).getPosition();
            Vector o = leafList.get(2).getPosition();
            Vector m = leafList.get(3).getPosition();
            double hMass = leafList.get(0).getType().getMass();
            double oMass = leafList.get(2).getType().getMass();

            centerMass.Ea1Tv1(hMass, h1);
            centerMass.PEa1Tv1(hMass, h2);
            centerMass.PEa1Tv1(oMass, o);   //WHAT IS THIS PE AND TE
            centerMass.TE(1 / (2 * hMass + oMass));   //##################WHAT IS TE?

            dr.Ev1Mv2(m, centerMass);
            torque.E(mForce);
            torque.XE(dr);  //#############################WHAT DOES THIS XE DO?
            q.E(torque);
            dr.Ev1Mv2(h1, centerMass);
            torque.E(h1Force);
            torque.XE(dr);
            q.PE(torque);
            dr.Ev1Mv2(h2, centerMass);
            torque.E(h2Force);
            torque.XE(dr);
            q.PE(torque);
            dr.Ev1Mv2(o, centerMass);
            torque.E(oForce);
            torque.XE(dr);
            q.PE(torque);
            //for the total torque q
//#####################MOLECULE SITE SOURCE
            Vector lPos = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).position;
            dr.Ev1Mv2(centerMass, lPos);
            ForceSum += totalForce.dot(dr);

//bring class
            OrientationFull3D or = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).orientation;
            Vector a0 = or.getDirection();//om
            Vector a1 = or.getSecondaryDirection();//h1h2

            dr.Ev1Mv2(h2, h1);
            dr.normalize();

            Vector axis = space.makeVector();
            Vector om = space.makeVector();
            Vector hh = space.makeVector();
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
            hh.normalize();
            a2.E(om);
            a2.XE(hh);
            a2.normalize();
            double[][] array1 = new double[3][3];
            om.assignTo(array1[0]);
            hh.assignTo(array1[1]);
            a2.assignTo(array1[2]);
            Matrix newa = new Matrix(array1).transpose();

            // newa = rotationMatrix * a so rotationMatrix = newa * a^-1
            a = a.inverse();
            Matrix matrix = newa.times(a);

            //https://en.wikipedia.org/wiki/Rotation_matrix#Determining_the_axis
            double d, b, c, g, h, f;
            b = matrix.get(0, 1);
            c = matrix.get(0, 2);
            d = matrix.get(1, 0);
            f = matrix.get(1, 2);
            g = matrix.get(2, 0);
            h = matrix.get(2, 1);
            axis.setX(0, h - f);
            axis.setX(1, c - g);
            axis.setX(2, d - b);
            double beta = Math.asin(0.5 * Math.sqrt(axis.squared()));
            if (beta == 0) continue;
            axis.normalize();


            double DUDT = -q.dot(axis);
            double denominator = 1 - Math.cos(beta);
            if (denominator == 0) continue;
            orientationSum += 1.5 * (beta - Math.sin(beta)) / denominator * DUDT;  //TRY ANOTHER 3/2 THAT U WERE GETTING FOR TORQUE ..DERIVATION
//            orientationSum += 0.5 * betaNew * DUDT;  //Only kappa3


        }

        if (!doTranslation) ForceSum = 0;
        if (!doRotation) orientationSum = 0;
        if (justDADB) {
            if (justU) {
                int N = molecules.size();
                double fac = (doTranslation ? 1.5 : 0) * (N - 1) + (doRotation ? 1.5 : 0) * N;
                x[0] = (x0 + latticeEnergy) + (fac * temperature) + 0.5 * ForceSum + orientationSum;

                for (int i = 0; i < molecules.size(); i++) {    //#########GET GETATOM/MOLECULE WORKING
                    for (int j = 0; j > i; j++) {  //12, 13...

                        LatticeSumMolecularCrystal.atomicToMolecularD(aTensor, i, j);
                    }
                }

                for (int i = 0; i < molecules.size(); i++) {    //#########GET GETATOM/MOLECULE WORKING
                    LatticeSumMolecularCrystal.atomicToMolecularD(aTensor, i, i);  //11, 22, etc. (call same molecules)
                }

            } else {
                x[0] = x0 + 0.5 * ForceSum + orientationSum; //translation and rotation
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
        //##########################3WHAT IS THE FOLLOWING STUFF???/
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


    public void debug() {
        Box box = latticeCoordinates.getBox();
        IMoleculeList molecules = box.getMoleculeList();
        //test make k3 rotate back to its nominal orientation
//###############################WHAT IS DOKAPPA3?
        System.out.println("You're doing the bebug in MeterDADBWaterTIP4P");
        boolean doKappa3 = false;
        boolean onlyTranslation = false;
        if (doKappa3) System.out.println("You're doing kappa3");
        for (int j = 99; j >= 0 && doKappa3; j--) {
            for (int i = 0; i < molecules.size(); i++) {
                IMolecule molecule = molecules.get(i);
                IAtomList leafList = molecule.getChildList();
                Vector h1 = leafList.get(0).getPosition();
                Vector h2 = leafList.get(1).getPosition();
                Vector o = leafList.get(2).getPosition();
                Vector m = leafList.get(3).getPosition();
                OrientationFull3D or = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).orientation;
                Vector a0 = or.getDirection();//om
                Vector a1 = or.getSecondaryDirection();//h1h2
                dr.Ev1Mv2(h2, h1);
                dr.normalize();
                double acos = a1.dot(dr);
                if (acos > 1) acos = 1;
                if (acos < -1) acos = -1;
                dr.XE(a1);
                double beta = -Math.signum(dr.dot(a0)) * Math.acos(acos);

                h1.ME(o);
                double length = Math.sqrt(h1.squared());
                h1.normalize();
                Orientation3D orientation = new Orientation3D(space);
                orientation.setDirection(h1);
                orientation.rotateBy(-beta / (j + 1), a0);
                h1.Ea1Tv1(length, orientation.getDirection());
                h1.PE(o);
                h2.ME(o);
                h2.normalize();
                orientation.setDirection(h2);
                orientation.rotateBy(-beta / (j + 1), a0);
                h2.Ea1Tv1(length, orientation.getDirection());
                h2.PE(o);
//    			if(i == 0){
//    				System.out.println(kappa);
//    			}
            }
            double u = meterPE.getDataAsScalar() - latticeEnergy;
            System.out.println(j * 0.01 + " " + u);
        }
        if (doKappa3) System.exit(2);

//############################WHAT IS THIS TEST? WHY DOING IT
        //test make configuration translate back to its nominal position
        if (onlyTranslation) System.out.println("You're doing only translation");
        for (int j = 99; j > 0 && onlyTranslation; j--) {
            for (int i = 0; i < molecules.size(); i++) {
                IMolecule molecule = molecules.get(i);
                IAtomList leafList = molecule.getChildList();
                Vector h1 = leafList.get(0).getPosition();
                Vector h2 = leafList.get(1).getPosition();
                Vector o = leafList.get(2).getPosition();
                Vector m = leafList.get(3).getPosition();
                double hmass = leafList.get(0).getType().getMass();
                double omass = leafList.get(2).getType().getMass();
                centerMass.Ea1Tv1(hmass, h1);
                centerMass.PEa1Tv1(hmass, h2);
                centerMass.PEa1Tv1(omass, o);
                centerMass.TE(1 / (2 * hmass + omass));
                Vector nominalCenterMass = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).position;
                dr.Ev1Mv2(nominalCenterMass, centerMass);
                dr.TE(1.0 / (j + 1));
                o.PE(dr);
                h1.PE(dr);
                h2.PE(dr);
                m.PE(dr);

//    			if(i == 0){
//    				System.out.println(latticeEnergy);
//    			}
            }

            double u = meterPE.getDataAsScalar() - latticeEnergy;

            pcForceSum.reset();
            potentialMaster.calculate(box, id, pcForceSum);
            double fdr = 0, fdr2 = 0;
            for (int i = 0; i < molecules.size(); i++) {

                IMolecule molecule = molecules.get(i);
                IAtomList leafList = molecule.getChildList();
                Vector h1 = leafList.get(0).getPosition();
                Vector h2 = leafList.get(1).getPosition();
                Vector o = leafList.get(2).getPosition();
                Vector m = leafList.get(3).getPosition();
                double hmass = leafList.get(0).getType().getMass();
                double omass = leafList.get(2).getType().getMass();
                centerMass.Ea1Tv1(hmass, h1);
                centerMass.PEa1Tv1(hmass, h2);
                centerMass.PEa1Tv1(omass, o);
                centerMass.TE(1 / (2 * hmass + omass));
                Vector nominalCenterMass = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).position;
                dr.Ev1Mv2(nominalCenterMass, centerMass);

                Vector h1Force = agentManager.getAgent(leafList.get(0));
                Vector h2Force = agentManager.getAgent(leafList.get(1));
                Vector oForce = agentManager.getAgent(leafList.get(2));
                Vector mForce = agentManager.getAgent(leafList.get(3));
                totalForce.E(h1force);
                totalForce.PE(h2force);
                totalForce.PE(mforce);
                totalForce.PE(oforce);
                fdr += totalForce.dot(dr);
                dr.normalize();
                fdr2 += totalForce.dot(dr);
            }
            System.out.println(j * 0.01 + " " + u + " " + fdr + " " + fdr2);

        }

        if (onlyTranslation) System.exit(2);

//debug for rotation angle and axis
        boolean allKappa = false;
        if (allKappa) System.out.println("You're doing all kappa");
        double[] beta0 = new double[46];
        double forceSum = 0;
        for (int j = 99; j > 0 && allKappa; j--) {
            for (int i = 0; i < molecules.size(); i++) {
                IMolecule molecule = molecules.get(i);
                IAtomList leafList = molecule.getChildList();
                Vector h1 = leafList.get(0).getPosition();
                Vector h2 = leafList.get(1).getPosition();
                Vector o = leafList.get(2).getPosition();
                Vector m = leafList.get(3).getPosition();
                OrientationFull3D or = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).orientation;
                Vector a0 = or.getDirection();//om
                Vector a1 = or.getSecondaryDirection();//h1h2
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

                double d, b, c, g, h, f;
                b = matrix.get(0, 1);
                c = matrix.get(0, 2);
                d = matrix.get(1, 0);
                f = matrix.get(1, 2);
                g = matrix.get(2, 0);
                h = matrix.get(2, 1);
                axis.setX(0, h - f);
                axis.setX(1, c - g);
                axis.setX(2, d - b);
                double beta = Math.asin(0.5 * Math.sqrt(axis.squared()));
                if (beta == 0) continue;
                axis.normalize();

                if (j == 99) beta0[i] = beta;
                double hMass = leafList.get(0).getType().getMass();
                double oMass = leafList.get(2).getType().getMass();
                centerMass.Ea1Tv1(hMass, h1);
                centerMass.PEa1Tv1(hMass, h2);
                centerMass.PEa1Tv1(oMass, o);
                centerMass.TE(1 / (2 * hMass + oMass));

                h1.ME(centerMass);
                h2.ME(centerMass);
                o.ME(centerMass);
                m.ME(centerMass);
                double hLength = Math.sqrt(h1.squared());
                double oLength = Math.sqrt(o.squared());
                double mLength = Math.sqrt(m.squared());
                h1.normalize();
                h2.normalize();
                o.normalize();
                m.normalize();
                Orientation3D orientation = new Orientation3D(space);
                orientation.setDirection(h1);
                orientation.rotateBy(beta / (j + 1), axis);
                h1.Ea1Tv1(hLength, orientation.getDirection());
                h1.PE(centerMass);
                orientation.setDirection(h2);
                orientation.rotateBy(beta / (j + 1), axis);
                h2.Ea1Tv1(hLength, orientation.getDirection());
                h2.PE(centerMass);
                orientation.setDirection(o);
                orientation.rotateBy(beta / (j + 1), axis);
                o.Ea1Tv1(oLength, orientation.getDirection());
                o.PE(centerMass);
                orientation.setDirection(m);
                orientation.rotateBy(beta / (j + 1), axis);
                m.Ea1Tv1(mLength, orientation.getDirection());
                m.PE(centerMass);

//    			dr.Ev1Mv2(m, o);
//    			ph1h2.Ev1Mv2(h2, h1);
//    			System.out.println("new =" +dr.dot(ph1h2));


                if (i == 0) {
                    System.out.println("i" + i + " " + beta);
                }
            }


            double u = meterPE.getDataAsScalar() - latticeEnergy;

            System.out.println(j * 0.01 + " " + u);

            pcForceSum.reset();
            potentialMaster.calculate(box, id, pcForceSum);

            double orientationSum = 0;
            double DUDTsum = 0, fooSum = 0;
            for (int i = 0; i < molecules.size(); i++) {

                IMolecule molecule = molecules.get(i);
                IAtomList leafList = molecule.getChildList();
                Vector h1Force = agentManager.getAgent(leafList.get(0));
                Vector h2Force = agentManager.getAgent(leafList.get(1));
                Vector oForce = agentManager.getAgent(leafList.get(2));
                Vector mForce = agentManager.getAgent(leafList.get(3));

                totalForce.E(h1Force);
                totalForce.PE(h2Force);
                totalForce.PE(mForce);
                totalForce.PE(oForce);
                Vector h1 = leafList.get(0).getPosition();
                Vector h2 = leafList.get(1).getPosition();
                Vector o = leafList.get(2).getPosition();
                Vector m = leafList.get(3).getPosition();
                double hMass = leafList.get(0).getType().getMass();
                double oMass = leafList.get(2).getType().getMass();

                centerMass.Ea1Tv1(hMass, h1);
                centerMass.PEa1Tv1(hMass, h2);
                centerMass.PEa1Tv1(oMass, o);
                centerMass.TE(1 / (2 * hMass + oMass));

                dr.Ev1Mv2(m, centerMass);
                torque.E(mForce);
                torque.XE(dr);
                q.E(torque);
                dr.Ev1Mv2(h1, centerMass);
                torque.E(h1Force);
                torque.XE(dr);
                q.PE(torque);
                dr.Ev1Mv2(h2, centerMass);
                torque.E(h2Force);
                torque.XE(dr);
                q.PE(torque);
                dr.Ev1Mv2(o, centerMass);
                torque.E(oForce);
                torque.XE(dr);
                q.PE(torque);
                //for the total torque q

                Vector lPos = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).position;
                dr.Ev1Mv2(centerMass, lPos);
                forceSum += totalForce.dot(dr);


                OrientationFull3D or = ((MoleculeSiteSource.LatticeCoordinate) latticeCoordinates.getAgent(molecule)).orientation;
                Vector a0 = or.getDirection();//om
                Vector a1 = or.getSecondaryDirection();//h1h2

                dr.Ev1Mv2(h2, h1);
                dr.normalize();

                Vector axis = space.makeVector();
                Vector om = space.makeVector();
                Vector hh = space.makeVector();
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
                hh.normalize();
                a2.E(om);
                a2.XE(hh);
                a2.normalize();
                double[][] array1 = new double[3][3];
                om.assignTo(array1[0]);
                hh.assignTo(array1[1]);
                a2.assignTo(array1[2]);
                Matrix newa = new Matrix(array1).transpose();

                // newa = rotationMatrix * a so rotationMatrix = newa * a^-1
                a = a.inverse();
                Matrix matrix = newa.times(a);

                //https://en.wikipedia.org/wiki/Rotation_matrix#Determining_the_axis
                //get axis in a newer and cleaner way.
                double d, b, c, g, h, f;
                b = matrix.get(0, 1);
                c = matrix.get(0, 2);
                d = matrix.get(1, 0);
                f = matrix.get(1, 2);
                g = matrix.get(2, 0);
                h = matrix.get(2, 1);
                axis.setX(0, h - f);
                axis.setX(1, c - g);
                axis.setX(2, d - b);
                double beta = Math.asin(0.5 * Math.sqrt(axis.squared()));
                if (beta == 0) continue;
                axis.normalize();


//                System.out.println(beta);
//                System.exit(2);

                double DUDT = -q.dot(axis);
                double denominator = 1 - Math.cos(beta);
                if (denominator == 0) continue;
                orientationSum += 1.5 * (beta - Math.sin(beta)) / denominator * DUDT;
                fooSum += (beta - Math.sin(beta)) / (1 - Math.cos(beta));
                DUDTsum += DUDT * beta0[i];
            }
//            System.out.println(j * 0.01 + " " + u + " " + orientationSum + " " + DUDTsum);
        }
        if (allKappa) System.exit(2);


    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }










}
