/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity.mappedDensityfromlatticesite;

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
import etomica.math.SpecialFunctions;
import etomica.math.function.FunctionDifferentiable;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.units.dimensions.*;

public class MeterMappedAvgDebyeWallerxyz implements IDataSource, DataSourceIndependent, AtomLeafAgentManager.AgentSource<Vector> {

    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective id;
    protected final PotentialCalculationForceSum pc;
    protected final AtomLeafAgentManager<Vector> agentManager;
    protected final Box box;
    protected DataSourceUniform xDataSource;
    protected DataFunction data;
    protected IDataInfo dataInfo;
     protected Vector rivector;
     /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    protected int profileDim;
    /**
     * Meter that defines the property being profiled.
     */
    protected final DataTag tag;
    protected double temperature;
protected CoordinateDefinition latticesite;
    protected final FunctionDifferentiable c;
    protected Behavior behavior;
    protected double zidotz;
    protected double msd;
    protected double[] qvector;
protected int numAtoms;
    public enum Behavior {
        NORMAL, P, ZIDOT, DZIDOT
    }

    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterMappedAvgDebyeWallerxyz(int numAtoms, double[] qvector, double msd, Box box, PotentialMaster potentialMaster, double temperature, FunctionDifferentiable c, CoordinateDefinition latticesite) {
        this.box = box;
        this.temperature = temperature;
        this.c = c;
        this.qvector = qvector;
this.numAtoms=numAtoms;
        this.msd = msd;
 this.latticesite=latticesite;
this.rivector =box.getSpace().makeVector();
        xDataSource = new DataSourceUniform("x", Length.DIMENSION);
        tag = new DataTag();
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pc = new PotentialCalculationForceSum();
        agentManager = new AtomLeafAgentManager<Vector>(this, box);
        pc.setAgentManager(agentManager);
    }

    public void setBehavior(Behavior b) {
        behavior = b;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
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

             for (IAtom atom : atoms) {
                 rivector.Ev1Mv2(atom.getPosition(),latticesite.getLatticePosition(atom)) ;
                 box.getBoundary().nearestImage(rivector);
                double ri = Math.sqrt(rivector.squared()) ;
                double fr = agentManager.getAgent(atom).dot(rivector)/ri;
                double beta=1/temperature;
                double erfx = 1- SpecialFunctions.erfc(ri*Math.sqrt(3/(2*msd)));
                double qdotrcap =0.0;

                for (int i=0;i<3;i++){
                    qdotrcap =qdotrcap +(qvector[i]*rivector.getX(i))/ri;         //  double qdotrcap = qvector.dot(rivector)/ri;
                }


double er5 = 1- SpecialFunctions.erfc(5*Math.sqrt(3.0/(2.0*msd)));
double ex=Math.exp(3*ri*ri/(2*msd));

double Drxyz=1;
double Nrxyz=qdotrcap*qdotrcap*( 0.9999999999999999*msd*ri + 1.6825506805955375*Math.pow(10,-16)*msd*ri + 0.9999999999999999*ri*ri*ri + beta*fr*(0.3333333333333333*msd*ri*ri) )/ri;
//double Nrxyznew=( 2*Math.PI*qdotrcap*qdotrcap* (ri*(0.15915494309189532*msd*ri + 2.6778625781941256*Math.pow(10,-17)*msd*ri + 0.15915494309189532*ri*ri*ri + beta*fr*(0.05305164769729844*msd* ri*ri))) )/(ri*ri);

//double Nrxyznewnew=3*qdotrcap*qdotrcap*(ri *(-0.111111*beta*fr*msd*msd + 4.54608*Math.pow(10,-17)*msd*ri) + ex*(0.0804001*beta*fr*Math.pow(msd,2.5)+ 0.2412*Math.pow(msd,1.5)*ri) *erfx)/(ri*ri);

double Drxyzminus5to5=er5*er5*er5*(ri*(-0.333333*beta*fr*Math.pow(msd,2.5) + 1.36382*Math.pow(10,-16)*Math.pow(msd,1.5)*ri) +  ex*erfx*(0.2412*beta*fr*msd*msd*msd + 0.723601* msd*msd *ri)  )/(ri*ri*Math.pow(msd,1.5));

//double Nrxyzminus5to5=;

//System.out.println(Nrwithr2tillinfinity+" "+Nrwithr2till5+" "+Drwithr2tillinfinity+" "+Drwithr2till5+" "+Nrxyz+" "+Drxyzminus5to5+" "+Nrxyznewnew);

  //              y[0]=y[0]+ (Nrxyzminus5to5/Drxyzminus5to5);

                y[0]=y[0]+ (Nrxyz/Drxyz);

                 //            System.out.println(y[0]);

            }
        y[0]=y[0]/numAtoms;
        return data;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray) xDataSource.getData();
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray) xDataSource.getDataInfo();
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataTag getIndependentTag() {
        return xDataSource.getTag();
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    public void reset() {
        if (box == null) return;

         xDataSource.setXMin(0);

        data = new DataFunction(new int[]{xDataSource.getNValues()});
        dataInfo = new DataInfoFunction("Mapped Average Profile", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{1, -1}), this);
        dataInfo.addTag(tag);
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

    @Override
    public Vector makeAgent(IAtom a, Box agentBox) {
        return box.getSpace().makeVector();
    }

    @Override
    public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {

    }
}
