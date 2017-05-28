/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataSourceIndependentSimple;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.integrator.IntegratorBox;
import etomica.units.Energy;
import etomica.units.Null;
import etomica.math.DoubleRange;
import etomica.data.histogram.HistogramReweightedData;
import etomica.data.histogram.HistogramSimple;

/**
 * Meter used for overlap sampling in the target-sampled system.  The meter
 * measures the energy difference for the harmonic and target
 * potentials.
 * 
 * the sampling of <beta*U_BA>B
 * gives us the information of
 * 	a. < beta*U_BW >B : notation ---- [ betaUBWf_avg ]
 * 	b. < beta*U_BW >W : notation ---- [ betaUBWr_avg ] 
 *  ---------- U_BW = U_B - U_W
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
public class MeterWorkTargetHarmonicandReweighting implements IEtomicaDataSource {
    
    public MeterWorkTargetHarmonicandReweighting(IntegratorBox integrator, MeterHarmonicEnergy meterHarmonicEnergy, double ref) {
        meterEnergy = new MeterPotentialEnergyFromIntegrator(integrator);
        this.integrator = integrator;
        this.meterHarmonicEnergy = meterHarmonicEnergy;
        this.refPref = ref;
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Scaled Harmonic and hard sphere Energies", Null.DIMENSION);
        histogramUBWf = new HistogramSimple(4000, new DoubleRange(-200, 200));
        histogramUBWr = new HistogramReweightedData(4000, new DoubleRange(-200, 200));
        nSamples = 0;
        numSum = 0;
        denomSum = 0;
        betaUBWf = 0;
        
        tag = new DataTag();
    }

    public IData getData() {
    	betaUBA =	(meterHarmonicEnergy.getDataAsScalar() -
                	(meterEnergy.getDataAsScalar() - latticeEnergy)) / integrator.getTemperature();
    	double eA = Math.exp(-meterHarmonicEnergy.getDataAsScalar()/integrator.getTemperature());
    	double eB = Math.exp(-(meterEnergy.getDataAsScalar() - latticeEnergy) / integrator.getTemperature());
        double eW = eA/(1+refPref*(eA/eB));
    	
    	double eBA = Math.exp(betaUBA);	
    	
        /*
         * Histogram Transformation
         */
        
        // < beta*U_BW >B
        betaUBW = Math.log(eBA+ refPref);
        
        betaUBWf += betaUBW;
        nSamples++;
        
    	histogramUBWf.addValue(betaUBW);
    	
    	// < beta*U_BW >W
    	numSum += betaUBW*(eA / (eB + refPref*eA));    // ( beta*U_AW )*(eW/eB)
    	denomSum += eA /(eB + refPref*eA); // eW/eB
    	
    	histogramUBWr.addValue(betaUBW, (eW/eB));
    	
    	if(Double.isInfinite(betaUBW)){
    		System.out.println(nSamples + " betaUBA: " + betaUBA + " ;betaUBW: " + betaUBW + " ;(eA /(eB + refPref*eA)): " + (eA /(eB + refPref*eA)));
    	}
    	
    	data.x = betaUBA;
    	return data;
    }

    public void setLatticeEnergy(double newLatticeEnergy) {
        latticeEnergy = newLatticeEnergy;
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }
    
    
    public IData getDataHistogramBetaUBWf(){
    	dataHistogram = new DataFunction(new int[]{histogramUBWf.getNBins()}, histogramUBWf.getHistogram());
    	return dataHistogram;
    }
    
    public IEtomicaDataInfo getDataInfoHistogramBetaUBWf(){
    	DataInfoDoubleArray independentInfo = new DataInfoDoubleArray("Energy", Energy.DIMENSION, new int[]{histogramUBWf.getNBins()});
    	xDataSource = new DataSourceIndependentSimple(histogramUBWf.xValues(), independentInfo);
    	dataInfoHistogram = new DataInfoFunction("Energy Histogram", Null.DIMENSION,xDataSource);
    	return dataInfoHistogram;
    }
    
    public IData getDataHistogramBetaUBWr(){
    	dataHistogram = new DataFunction(new int[]{histogramUBWr.getNBins()}, histogramUBWr.getHistogram());
    	return dataHistogram;
    }
    
    public IEtomicaDataInfo getDataInfoHistogramBetaUBWr(){
    	DataInfoDoubleArray independentInfo = new DataInfoDoubleArray("Energy", Energy.DIMENSION, new int[]{histogramUBWr.getNBins()});
    	xDataSource = new DataSourceIndependentSimple(histogramUBWr.xValues(), independentInfo);
    	dataInfoHistogram = new DataInfoFunction("Energy Histogram", Null.DIMENSION,xDataSource);
    	return dataInfoHistogram;
    }
    
    public double getBetaUBWf(){
    	return betaUBWf/nSamples;
    	
    }
    
    public double getBetaUBWr(){
    	return numSum/denomSum;
    	
    }
    
    
    protected final MeterPotentialEnergyFromIntegrator meterEnergy;
    protected final MeterHarmonicEnergy meterHarmonicEnergy;
    protected final IntegratorBox integrator;
    protected final DataDouble data;
    protected final DataInfoDouble dataInfo;
    protected final DataTag tag;
    protected double latticeEnergy;
    protected double betaUBW, betaUBA; 
    protected double betaUBWf;
    protected final double refPref;
    protected HistogramSimple histogramUBWf;
    protected HistogramReweightedData histogramUBWr;
    protected double nSamples;
    protected double numSum, denomSum;
    
    protected DataSourceIndependentSimple xDataSource;
    private DataFunction dataHistogram;
    private DataInfoFunction dataInfoHistogram;
    
}
