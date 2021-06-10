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
    protected final Vector rivector;
    protected final Vector deltaVec;
    protected final Vector rVec;
    protected final Vector RVec;
    protected int profileDim;
    protected final OrientationFull3D orientation;
    protected final RotationTensor3D rotationTensor;
    protected final Vector vel;
    protected final DataTag tag;
    protected double temperature;
    protected CoordinateDefinition latticesite;
    protected double msd;
    protected int rnumberofbins;
    protected int costhetaphinumberofbins;

    public MeterDensityAnisotropicHMA(double msd, int rnumberofbins, int costhetaphinumberofbins, Box box, PotentialMaster potentialMaster, double temperature, CoordinateDefinition latticesite) {
        this.box = box;
        this.temperature = temperature;
        this.msd = msd;
        this.rnumberofbins = rnumberofbins;
        this.costhetaphinumberofbins = costhetaphinumberofbins;

        this.Rmax = Math.sqrt(msd)*4;
        this.latticesite = latticesite;
        this.rivector = box.getSpace().makeVector();
        this.deltaVec = box.getSpace().makeVector();
        this.rVec = box.getSpace().makeVector();
        this.RVec = box.getSpace().makeVector();

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

        xDataSourcetheta.setXMin(-1);
        xDataSourcetheta.setXMax(+1);
        xDataSourcephi.setXMin(0);
        xDataSourcephi.setXMax(2*Math.PI);
        xDataSourcer.setNValues(rnumberofbins);
        xDataSourcetheta.setNValues(costhetaphinumberofbins);
        xDataSourcephi.setNValues(costhetaphinumberofbins);

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
        double dz = Rmax / xDataSourcer.getNValues();
        double [][][] ytemporary=new double[rnumberofbins][costhetaphinumberofbins][costhetaphinumberofbins];

        double sigma2 = msd / 3;
        double sigma = Math.sqrt(sigma2);
        double q = Math.pow(2 * Math.PI * sigma2, 1.5);

        for (IAtom atom : atoms) {
            rivector.Ev1Mv2(atom.getPosition(), latticesite.getLatticePosition(atom));
            box.getBoundary().nearestImage(rivector);
            double ri = Math.sqrt(rivector.squared());
            double dlnpridr = - ri / sigma2;

            for (int i = 0; i < rnumberofbins; i++) {
                double r = (i + 0.5) * dz;

                for (int j = 0; j < costhetaphinumberofbins; j++) {
                    double costheta = -1 + (j+0.5)*2.0 / costhetaphinumberofbins;
                    double sintheta = Math.sqrt(1.0 - costheta*costheta);

                    for (int k = 0; k < costhetaphinumberofbins; k++) {
                        double phi = (k + 0.5) * 2 * Math.PI / costhetaphinumberofbins;

                        int mm = j * costhetaphinumberofbins + k;
                        double pr = (Math.exp(-r * r / (2 * sigma2)));

                        //point where dens measured, relative to lat site
                        deltaVec.setX(0, r * sintheta * Math.cos(phi));
                        deltaVec.setX(1, r * sintheta * Math.sin(phi));
                        deltaVec.setX(2, r * costheta);

                        rVec.Ev1Mv2(rivector, deltaVec);
                        RVec.Ea1Tv1(-1, deltaVec);
                        vel.E(Singlet3DmappingDelta0.xyzDot(rVec, RVec, sigma));

                        double rdot = vel.dot(rivector) / ri;
                        double a1 = -agentManager.getAgent(atom).dot(vel);
                        double a2 =  temperature * rdot * dlnpridr;
      //                  System.out.println(a1+", "+a2+", "+(a1+a2));
                        ytemporary[i][j][k] += (pr / q) - agentManager.getAgent(atom).dot(vel) + (temperature * rdot * dlnpridr);
                        if (Double.isNaN(ytemporary[i][j][k])) {
                            Singlet3DmappingDelta0.xyzDot(rivector, deltaVec, sigma);
                            System.out.println("oops");
                        }
                    }
                }
            }
        }
        int n = 0;
        for (int i = 0; i < rnumberofbins; i++) {
            for (int j = 0; j < costhetaphinumberofbins; j++) {
                for (int k = 0; k < costhetaphinumberofbins; k++) {
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
