/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.data.DataSourceIndependentSimple;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.potential.PotentialMaster;
import etomica.units.Energy;
import etomica.units.Null;
import etomica.math.DoubleRange;
import etomica.data.histogram.HistogramReweightedData;
import etomica.data.histogram.HistogramSimple;

/**
 * Meter used for overlap sampling in the harmonic-sampled system.  The meter
 * measures the energy difference for the harmonic and target
 * potentials.
 * 
 * the sampling of <beta*U_AB>A
 * gives us the information of
 * 	a. < beta*U_AW >A : notation ---- [ betaUABf_avg ]
 * 	b. < beta*U_AW >W : notation ---- [ betaUABr_avg ] 
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
public class MeterWorkHarmonicTargetandReweighting implements IEtomicaDataSource {
    
    public MeterWorkHarmonicTargetandReweighting(MCMoveHarmonic mcMoveHarmonic, PotentialMaster potentialMaster, double ref) {
    	this.mcMoveHarmonic = mcMoveHarmonic;
        meterEnergy = new MeterPotentialEnergy(potentialMaster);
        meterEnergy.setBox(mcMoveHarmonic.getBox());
        this.refPref = ref;
        data = new DataDouble();
        dataInfo = new DataInfoDouble("Scaled Harmonic and soft sphere Energies", Null.DIMENSION);
        histogramUAWf = new HistogramSimple(4000, new DoubleRange(-200, 200));
        histogramUAWr = new HistogramReweightedData(4000, new DoubleRange(-200, 200));
        nSamples = 0;
        numSum = 0;
        denomSum = 0;
        betaUAWf = 0;        
        
        tag = new DataTag();
    }

    public IData getData() {
    	betaUAB =	((meterEnergy.getDataAsScalar() - latticeEnergy) -
        		mcMoveHarmonic.getLastTotalEnergy()) / temperature;
    	double eA = Math.exp(-mcMoveHarmonic.getLastTotalEnergy()/temperature);
    	double eB = Math.exp(-(meterEnergy.getDataAsScalar() - latticeEnergy) / temperature);
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
    		System.out.println(nSamples + " betaUAB: " + betaUAB + " ;betaUAW: " + betaUAW + " ;betaUAW*(eW / eA): " + betaUAW*(eW / eA));
    	}
    	data.x = betaUAB;
    	return data;
    }

    public void setLatticeEnergy(double newLatticeEnergy) {
        latticeEnergy = newLatticeEnergy;
    }
    
    public void setTemperature(double newTemperature) {
        temperature = newTemperature;
    }
    
    public double getTemperature() {
        return temperature;
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }
    
    
    public IData getDataHistogramBetaUAWf(){
    	dataHistogram = new DataFunction(new int[]{histogramUAWf.getNBins()}, histogramUAWf.getHistogram());
    	return dataHistogram;
    }
    
    public IEtomicaDataInfo getDataInfoHistogramBetaUAWf(){
    	DataInfoDoubleArray independentInfo = new DataInfoDoubleArray("Energy", Energy.DIMENSION, new int[]{histogramUAWf.getNBins()});
    	xDataSource = new DataSourceIndependentSimple(histogramUAWf.xValues(), independentInfo);
    	dataInfoHistogram = new DataInfoFunction("Energy Histogram", Null.DIMENSION,xDataSource);
    	return dataInfoHistogram;
    }
    
    public IData getDataHistogramBetaUAWr(){
    	dataHistogram = new DataFunction(new int[]{histogramUAWr.getNBins()}, histogramUAWr.getHistogram());
    	return dataHistogram;
    }
    
    public IEtomicaDataInfo getDataInfoHistogramBetaUAWr(){
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
    protected final MeterPotentialEnergy meterEnergy;
    protected final MCMoveHarmonic mcMoveHarmonic;
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
