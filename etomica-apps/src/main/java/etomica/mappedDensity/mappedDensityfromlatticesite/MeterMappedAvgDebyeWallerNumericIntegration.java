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

public class MeterMappedAvgDebyeWallerNumericIntegration implements IDataSource, DataSourceIndependent, AtomLeafAgentManager.AgentSource<Vector> {

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
    public MeterMappedAvgDebyeWallerNumericIntegration(int numAtoms, double[] qvector, double msd, Box box, PotentialMaster potentialMaster, double temperature, FunctionDifferentiable c, CoordinateDefinition latticesite) {
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
              double  Nrold=((qdotrcap*qdotrcap/(ri*ri*ri))*((ri*Math.exp(3*ri*ri/(2*msd))*(-(0.0575824*Math.sqrt(msd)*ri)-(0.0191941*beta*Math.pow(msd,1.5)*fr))*erfx)+(ri*( ri* (-0.0795775*ri-0.0265258*beta*msd*fr+0.0795775*ri)+(Math.exp(3*ri*ri/(2*msd))*(0.0575824*Math.pow(msd,0.5)*ri+0.0191941*beta*Math.pow(msd,1.5)*fr)*erfx)+(ri*(0.0795775*ri+0.0265258*beta*msd*fr))  ))));
              double  Drold=  ((1/(Math.pow(msd,1.5)*ri*ri))*((ri*((-0.0795775*beta*Math.pow(msd,1.5)*fr)+(Math.pow(msd,0.5)*ri*(-0.238732+0.238732))))   ));

  double Nr=(qdotrcap*qdotrcap/(ri*ri))*( ((-0.0265258*beta*fr*msd+(1.33893*Math.pow(10,-17))*ri)*ri) + (Math.exp(3*ri*ri/(2*msd))*(0.0191941*beta*fr*Math.pow(msd,1.5)+0.0575824*Math.pow(msd,0.5)*ri)*erfx ) + ( (0.0265258*beta*fr*msd+0.0795775*ri)*ri )+(Math.exp(3*ri*ri/(2*msd))*Math.sqrt(msd)*erfx*(-0.0191941*beta*fr*msd-0.0575824*ri)    ));

  double Dr=(1/(Math.pow(msd,1.5)*ri*ri))*( (ri*(-0.0795775*beta*fr*Math.pow(msd,1.5)+4.01679*Math.pow(10,-17)*Math.pow(msd,0.5)*ri))   );

  //with rsquare dr integration




                 System.out.println(Nr+" "+Nrold+" "+Dr+" "+Drold);

                y[0]=y[0]+ (Nr/Dr);
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
