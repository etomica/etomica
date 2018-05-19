/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.gaussianwork;

import etomica.data.*;
import etomica.data.histogram.HistogramReweightedData;
import etomica.data.histogram.HistogramSimple;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.integrator.IntegratorBox;
import etomica.math.DoubleRange;
import etomica.potential.PotentialMaster;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Null;

/**
 * Meter used for multiharmonic overlap sampling in the A-sampled system.  \
 * The meter measures the energy difference for the A and B potentials.
 * 
 * the sampling of <beta*U_AB>A
 * gives us the information of
 * 	a. < beta*U_AW >A : notation ---- [ betaUBAf_avg ]
 * 	b. < beta*U_AW >W : notation ---- [ betaUBAr_avg ] 
 *  ---------- U_AW = U_A - U_W
 *  
 * output the probability distribution of the histogram for a. and b.
 * 
 * Key:
 * 	a. eA = e^(-beta*U_A)
 * 	b. eB = e^(-beta*U_B)
 * 	c. eW = e^(-beta*U_W) 
 * 
 * @author Tai Boon Tan
 */
public class MeterWorkABandReweighting implements IDataSource {
    
    public MeterWorkABandReweighting(IntegratorBox integratorA, PotentialMaster potentialMasterB, double ref) {
    	meterA = new MeterPotentialEnergyFromIntegrator(integratorA);
    	meterB = new MeterPotentialEnergy(potentialMasterB, integratorA.getBox());
    	
        this.refPref = ref;
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Scaled Multiharmonic System A and System B Energies", Null.DIMENSION);
        histogramUAWf = new HistogramSimple(4000, new DoubleRange(-200, 200));
        histogramUAWr = new HistogramReweightedData(4000, new DoubleRange(-200, 200));
        nSamples = 0;
        numSum = 0;
        denomSum = 0;
        betaUAWf = 0;        
        
        tag = new DataTag();
    }

    public IData getData() {
    	betaUAB =	(meterB.getDataAsScalar()  - meterA.getDataAsScalar()) / temperature;
    	double eA = Math.exp(-meterA.getDataAsScalar()/temperature);
    	double eB = Math.exp(-meterB.getDataAsScalar()/temperature);
        double eW = eA/(1+refPref*(eA/eB));
    	
    	double eAB = Math.exp(betaUAB);	
    	
        /*
         * Histogram Transformation
         * 
         * Notes: If UB is infinity, there is nothing you can do about it in term 
         * 			of betaUAW expression; you will get an infinity value for the
         * 			term.
         * 		  The numSum will become NaN because [infinity x 0].
         * 			betaUAW * (eW / eA) because eW is zero.
         * 
         */
        
        // < beta*U_AW >A
        betaUAW = Math.log(1 + refPref*eAB); 
        
        betaUAWf += betaUAW;
        nSamples++;
        
    	histogramUAWf.addValue(betaUAW);
    	
    	// < beta*U_AW >W
   		numSum += betaUAW*(eW / eA);    // ( beta*U_AW )*(eW/eA)
   		denomSum += (eW/eA); // eW/eA
    	    	
    	if (!Double.isInfinite(betaUAW)){
    		histogramUAWr.addValue(betaUAW, (eW/eA));
    	}
    	
    	if(Double.isInfinite(betaUAW)){
    		System.err.println(nSamples + " betaUAB: " + betaUAB + " ;betaUAW: " + betaUAW + " ;betaUAW*(eW / eA): " + betaUAW*(eW / eA));
    	}
    	//System.out.println("energy diff: " + betaUAB);
    	data.x = betaUAB;
    	return data;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }
    
    
    public IData getDataHistogramBetaUAWf(){
    	dataHistogram = new DataFunction(new int[]{histogramUAWf.getNBins()}, histogramUAWf.getHistogram());
    	return dataHistogram;
    }
    
    public IDataInfo getDataInfoHistogramBetaUAWf(){
    	DataInfoDoubleArray independentInfo = new DataInfoDoubleArray("Energy", Energy.DIMENSION, new int[]{histogramUAWf.getNBins()});
    	xDataSource = new DataSourceIndependentSimple(histogramUAWf.xValues(), independentInfo);
    	dataInfoHistogram = new DataInfoFunction("Energy Histogram", Null.DIMENSION,xDataSource);
    	return dataInfoHistogram;
    }
    
    public IData getDataHistogramBetaUAWr(){
    	dataHistogram = new DataFunction(new int[]{histogramUAWr.getNBins()}, histogramUAWr.getHistogram());
    	return dataHistogram;
    }
    
    public IDataInfo getDataInfoHistogramBetaUAWr(){
    	DataInfoDoubleArray independentInfo = new DataInfoDoubleArray("Energy", Energy.DIMENSION, new int[]{histogramUAWr.getNBins()});
    	xDataSource = new DataSourceIndependentSimple(histogramUAWr.xValues(), independentInfo);
    	dataInfoHistogram = new DataInfoFunction("Energy Histogram", Null.DIMENSION,xDataSource);
    	return dataInfoHistogram;
    }
    
    public double getBetaUAWf(){
    	return betaUAWf/nSamples;
    	
    }
    
    public double getBetaUAWr(){
    	return numSum/denomSum;
    	
    }
    
    
    protected double temperature;
    protected MeterPotentialEnergyFromIntegrator meterA;
    protected MeterPotentialEnergy meterB;
    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected double latticeEnergy;
    protected double betaUAW, betaUAB; 
    protected double betaUAWf;
    protected final double refPref;
    protected HistogramSimple histogramUAWf;
    protected HistogramReweightedData histogramUAWr;
    protected double nSamples;
    protected double numSum, denomSum;
    
    protected DataSourceIndependentSimple xDataSource;
    private DataFunction dataHistogram;
    private DataInfoFunction dataInfoHistogram;
    
}
