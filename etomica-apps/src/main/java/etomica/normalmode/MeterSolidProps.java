/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterSolidProps implements IDataSource, AgentSource<MyAgent> {

    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final CoordinateDefinition coordinateDefinition;
    protected final DataSourceScalar meterPE;
    protected final PotentialMaster potentialMaster;
    protected final AtomLeafAgentManager<MyAgent> forceManager;
    protected final IteratorDirective id;
    protected final Vector dr;
    protected double ULat, PLat;
    protected double volume, density;
    protected int nMol;
    protected double dP, ddP, f1, f2;
    protected final double temperature;
    private final PotentialCalculation pcSolidProps;
    protected boolean isLJ = false;
    
    public MeterSolidProps(Space space, DataSourceScalar meterPE, PotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition, double temperature, double dP, double ddP, double ULat, double PLat, boolean isLS) {
        int nData = 2;
//        int nData = 10;
        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.space = space;
        this.coordinateDefinition = coordinateDefinition;
        this.meterPE = meterPE;
        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
//        pcForceSum = new PotentialCalculationForceSum();
//        forceManager = new AtomLeafAgentManager<MyAgent>(this, coordinateDefinition.getBox(), MyAgent.class);
//        pcForceSum.setAgentManager(forceManager);
        dr = space.makeVector();
        this.temperature = temperature;
        this.dP = dP;
        this.ddP = ddP;
        this.ULat = ULat;
        this.PLat = PLat;
        volume = coordinateDefinition.getBox().getBoundary().volume();
        nMol = coordinateDefinition.getBox().getLeafList().getAtomCount();
        density = nMol/volume;
    	f1 = (-1.0/volume + dP/temperature)/3.0/(nMol-1.0);
    	f2 = (1.0/volume/volume + ddP/temperature)/3.0/(nMol-1.0)  +  f1*f1;
    	if(isLJ){
        	double fe = 1.00000;
        	double fee = 1.00000;
    		pcSolidProps = new PotentialCalculationLJSP(space,coordinateDefinition.getBox(),coordinateDefinition,temperature,dP, f1, fe,fee);
    	}else{
    		pcSolidProps = new PotentialCalculationEFSSP(space,coordinateDefinition.getBox(),coordinateDefinition,temperature, f1, isLS);
    	}
        forceManager = new AtomLeafAgentManager<MyAgent>(this, coordinateDefinition.getBox());
//        pcUP.setAgentManager(forceManager);


    }
    

    //U_direct is computed in the main classes using MeterPotentialEnergyFromIntegrator.
    public IData getData() {
        Box box = coordinateDefinition.getBox();
        double[] sum;
        if(isLJ){
            ((PotentialCalculationLJSP)pcSolidProps).reset();
            potentialMaster.calculate(box, id, pcSolidProps);
            sum = ((PotentialCalculationLJSP)pcSolidProps).getSum();
        }else{
            ((PotentialCalculationEFSSP)pcSolidProps).reset();
            potentialMaster.calculate(box, id, pcSolidProps);
            sum = ((PotentialCalculationEFSSP)pcSolidProps).getSum();
        }
        double[] x = data.getData();
        //ALL with "ij" pairs.
        double fdr     = sum[0];
//        double fr      = sum[1];
//        double fR      = sum[2];
//        double T1PhiT1 = sum[3];
//        double T1Phidr = sum[4];
//        double rPhir   = sum[5];
//        double drPhidr = sum[6];
        double dU    = meterPE.getDataAsScalar() - ULat;
        double Ur   = dU + 0.5*fdr;
//        double Pvir = fr/3.0/volume;
//        double Pr   = 1.0/3.0/volume*fR + f1*fdr - PLat;
//        System.out.println(" dUc = " + dU/nMol + "  dUr = " + ((3.0/2.0*(nMol-1.0)*temperature +Ur)/nMol));
//        System.out.println(ElectronVolt.UNIT.fromSim(dU/nMol) + "  " + ElectronVolt.UNIT.fromSim((3.0/2.0*(nMol-1.0)*temperature +Ur)/nMol));

        x[0] =  dU;//U
        x[1] = Ur; //Um   
//        x[2] = Pvir;//P	    
//        x[3] = Pr; //Pm
//        x[4] = dU*dU/temperature/temperature;//Cv
//        x[5] = -1.0/4.0/temperature*(fdr + drPhidr) + Ur*Ur/temperature/temperature;//Cvm
//        x[6] = 2.0/3.0*Pvir + 1.0/9.0/volume*rPhir - volume/temperature*Pvir*Pvir; //B
//        x[7] = -volume*(-2.0/9.0/volume/volume*fR + f2*fdr  - T1PhiT1) - volume/temperature*Pr*Pr; //Bm
//        x[8] = 1.0/temperature/temperature*dU*Pvir; //dPdT
//        x[9] = 1.0/temperature*(f1/2.0*fdr - 1.0/2.0*T1Phidr) + 1.0/temperature/temperature*Ur*Pr; //dPdTm
        return data; //U_New
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public final MyAgent makeAgent(IAtom a, Box oox) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(MyAgent agent, IAtom atom, Box box) {}
}
