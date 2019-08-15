/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity.crystal;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.dimensions.*;

public class MeterDensityAnisotropicHMA implements IDataSource, DataSourceIndependent, AtomLeafAgentManager.AgentSource<Vector> {

    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective id;
    protected DataSourceUniform xDataSourcer;
    protected DataSourceUniform xDataSourcetheta;
    protected DataSourceUniform xDataSourcephi;

    protected final PotentialCalculationForceSum pc;
    protected final AtomLeafAgentManager<Vector> agentManager;
    protected final Box box;
    protected DataFunction data;
    protected IDataInfo dataInfo;
    protected double Rmax;
    protected Vector rivector;
    protected Vector rvector;
    protected Vector perpendtorvectorinriplane;
    protected Vector crossprodrandri;
    protected Vector temporary;
    protected Vector rivectortransformed;
    protected Vector fivectortransformed;
    protected Vector rvectortransformed;
    protected int profileDim;
    protected final OrientationFull3D orientation;
    protected final RotationTensor3D rotationTensor;
    protected final Vector vel;
    protected final DataTag tag;
    protected double temperature;
    protected CoordinateDefinition latticesite;
    protected double[] arraymsd;
    protected int rnumberofbins;
    protected int thetaphinumberofbins;

    public MeterDensityAnisotropicHMA(double[] arraymsd, int rnumberofbins, int thetaphinumberofbins, Box box, PotentialMaster potentialMaster, double temperature, CoordinateDefinition latticesite) {
        this.box = box;
        this.temperature = temperature;
        this.arraymsd = arraymsd;
        this.rnumberofbins = rnumberofbins;
        this.thetaphinumberofbins = thetaphinumberofbins;

        this.Rmax = Math.sqrt(arraymsd[0])*4;
        this.latticesite = latticesite;
        this.rivector = box.getSpace().makeVector();
        this.rvector = box.getSpace().makeVector();
        this.perpendtorvectorinriplane = box.getSpace().makeVector();
        this.crossprodrandri = box.getSpace().makeVector();
        this.temporary = box.getSpace().makeVector();
        this.rivectortransformed = box.getSpace().makeVector();
        this.rvectortransformed = box.getSpace().makeVector();
        this.fivectortransformed = box.getSpace().makeVector();

        xDataSourcer = new DataSourceUniform("r", Length.DIMENSION);
        tag = new DataTag();
        xDataSourcer.setTypeMax(LimitType.HALF_STEP);
        xDataSourcer.setTypeMin(LimitType.HALF_STEP);
        xDataSourcer.setXMin(0);
        xDataSourcer.setXMax(Rmax);

        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pc = new PotentialCalculationForceSum();
        agentManager = new AtomLeafAgentManager<Vector>(this, box);
        pc.setAgentManager(agentManager);

        xDataSourcetheta = new DataSourceUniform("theta", Length.DIMENSION);
        xDataSourcetheta.setTypeMax(LimitType.HALF_STEP);
        xDataSourcetheta.setTypeMin(LimitType.HALF_STEP);

        xDataSourcephi = new DataSourceUniform("phi", Length.DIMENSION);
        xDataSourcephi.setTypeMax(LimitType.HALF_STEP);
        xDataSourcephi.setTypeMin(LimitType.HALF_STEP);

        xDataSourcetheta.setXMin(0);
        xDataSourcetheta.setXMax(Math.PI);
        xDataSourcephi.setXMin(0);
        xDataSourcephi.setXMax(2*Math.PI);
        xDataSourcer.setNValues(rnumberofbins);
        xDataSourcetheta.setNValues(thetaphinumberofbins);
        xDataSourcephi.setNValues(thetaphinumberofbins);

        rotationTensor = new RotationTensor3D();
        orientation = new OrientationFull3D(Space3D.getInstance());
        vel = new Vector3D();
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }
//sim.coordinateDefinition - latticesites

    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     */
    public int getProfileDim() {
        return profileDim;
    }

    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     * The given vector is converted to a unit vector, if not already.
     */
    public void setProfileDim(int dim) {
        profileDim = dim;
        reset();
    }

    /**
     * Returns the profile for the current configuration.
     */
    public IData getData() {
        pc.reset();
        potentialMaster.calculate(box, id, pc);
        data.E(0);
        double[] y = data.getData();

        IAtomList atoms = box.getLeafList();
        double L = box.getBoundary().getBoxSize().getX(profileDim);
        double dz = Rmax / xDataSourcer.getNValues();
         double [][][] ytemporary=new double[rnumberofbins][thetaphinumberofbins][thetaphinumberofbins];

        for (IAtom atom : atoms) {
            rivector.Ev1Mv2(atom.getPosition(), latticesite.getLatticePosition(atom));
            box.getBoundary().nearestImage(rivector);
            double ri = Math.sqrt(rivector.squared());
            //           double fr = agentManager.getAgent(atom).dot(rivector) / ri;
            //           double thetai = Math.acos(rivector.getX(2) / ri);
            //           double phii = Math.atan2(rivector.getX(1), rivector.getX(0));
            //           if (phii < 0) {
            //               phii = phii + 2 * Math.PI;
            //           }
            //
            for (int i = 0; i < rnumberofbins; i++) {
                double r = (i + 0.5) * dz;

                for (int j = 0; j < thetaphinumberofbins; j++) {
                    double theta = (j + 0.5) * Math.PI / thetaphinumberofbins;

                    for (int k = 0; k < thetaphinumberofbins; k++) {
                        double phi = (k + 0.5) * 2 * Math.PI / thetaphinumberofbins;

                        int mm = j * thetaphinumberofbins + k;
                        double pr = (Math.exp(-3 * r * r / (2 * arraymsd[mm])));
                        double q = Math.pow(2 * Math.PI * arraymsd[mm] / 3, 1.5);
                        int nmax = 5;
                        //    System.out.println(j + k + "r" + r + "theta" + theta + "phi" + phi + " " + "ri" + ri + "thetai" + thetai + "phii" + phii);
//point where dens measured from lat site
                        rvector.setX(0, Math.sin(theta) * Math.cos(phi));
                        rvector.setX(1, Math.sin(theta) * Math.sin(phi));
                        rvector.setX(2, Math.cos(theta));

                        //projection of rivector on rvector
                        double dotproduct = rvector.dot(rivector);
                        temporary.Ea1Tv1(dotproduct, rvector);
                        perpendtorvectorinriplane.Ev1Mv2(rivector, temporary);
                        perpendtorvectorinriplane.normalize();
                        crossprodrandri.E(rvector);
                        crossprodrandri.XE(perpendtorvectorinriplane);

                        orientation.setDirections(perpendtorvectorinriplane, crossprodrandri);
                        rotationTensor.setOrientation(orientation);
                        rotationTensor.invert();

                        double ti = Math.acos(dotproduct / ri);
                        double[] mappingvelocitytransformed = Singlet3Dmapping.xyzDot(nmax, ri, ti, 0, r, Math.sqrt(arraymsd[mm] / 3.0));
                        vel.E(mappingvelocitytransformed);
                        rotationTensor.transform(vel);

                        double rdot = vel.dot(rivector) / ri;
                        double dlnpridr = -3 * ri / arraymsd[mm];

                        ytemporary[i][j][k] += (pr / q) - agentManager.getAgent(atom).dot(vel) + (temperature * rdot * dlnpridr);
                    }
                }
            }
        }
        int n = 0;
        for (int i = 0; i < rnumberofbins; i++) {
            for (int j = 0; j < thetaphinumberofbins; j++) {
                for (int k = 0; k < thetaphinumberofbins; k++) {
                    y[n] = ytemporary[i][j][k] / box.getLeafList().size();
                    n = n + 1;
                }
            }
        }
         return data;
    }

    public DataDoubleArray getIndependentData(int i) {
        if(i==0){return (DataDoubleArray) xDataSourcer.getData();}
        if(i==1){return (DataDoubleArray) xDataSourcetheta.getData();}
        if(i==2){return (DataDoubleArray) xDataSourcephi.getData();}
        throw new RuntimeException();
    }
    public DataInfoDoubleArray getIndependentDataInfo(int i) {

        if(i==0){return (DataInfoDoubleArray) xDataSourcer.getDataInfo();}
        if(i==1){return (DataInfoDoubleArray) xDataSourcetheta.getDataInfo();}
        if(i==2){return (DataInfoDoubleArray) xDataSourcephi.getDataInfo();}
        throw new RuntimeException();
    }
    public int getIndependentArrayDimension() {
        return 3;
    }

    public DataTag getIndependentTag(){throw new RuntimeException();}
    public DataTag getIndependentTag(int i) {

        if(i==0){return xDataSourcer.getTag();}
        if(i==1){return xDataSourcetheta.getTag();}
        if(i==2){return xDataSourcephi.getTag();}
        throw new RuntimeException();
    }
    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    public void reset() {
        if (box == null) return;

        Boundary boundary = box.getBoundary();

        xDataSourcer.setXMin(0);
        xDataSourcer.setXMax(Rmax);
        xDataSourcetheta.setXMin(0);
        xDataSourcetheta.setXMax(Math.PI);
        xDataSourcephi.setXMin(0);
        xDataSourcephi.setXMax(2*Math.PI);

        data = new DataFunction(new int[]{xDataSourcer.getNValues(),xDataSourcetheta.getNValues(),xDataSourcephi.getNValues()});
        dataInfo = new DataInfoFunction("Mapped Average Profile", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{1, -1}), this);
        dataInfo.addTag(tag);
    }

    public DataSourceUniform getXDataSource(int i) {
        if(i==0)return xDataSourcer;
        if(i==1)return xDataSourcetheta;
        if(i==2)return xDataSourcephi;
        throw new RuntimeException();

    }

    @Override
    public Vector makeAgent(IAtom a, Box agentBox) {
        return box.getSpace().makeVector();
    }

    @Override
    public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {

    }
}
