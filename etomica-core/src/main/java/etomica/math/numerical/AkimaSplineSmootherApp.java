/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.data.AccumulatorHistory;
import etomica.data.DataSourceCountSteps;
import etomica.data.DataSourceIndependentSimple;
import etomica.data.IDataSink;
import etomica.data.history.HistoryScrolling;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.graphics.*;
import etomica.graphics.DeviceBox.LabelType;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorListenerAction;
import etomica.modifier.ModifierGeneral;
import etomica.units.dimensions.Null;
import etomica.util.Constants.CompassDirection;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;
import net.miginfocom.swing.MigLayout;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * App to drive AkimaSplineSmoother
 * 
 * @author Andrew Schultz
 */
public class AkimaSplineSmootherApp {
    
    public final JPanel panel;
    protected final AccumulatorHistory historyEy, historyED2, historyED2D, historyED3, historyED3D;
    protected final AkimaSplineSmoother fitter;
    protected final IDataSink sinkDy, sinkDy2, sink0, sinkSmoothed, sinkEy;
    protected DataFunction ySmooth, dataDy, dataDy2, y0, dataEy;
    protected double[] xd;
    protected double[] x, y, dy;
    protected PlotAction plotAction;

    public AkimaSplineSmootherApp(String infile) {
        double d2fac = 0;
        double d2dfac = .01;
        double d3fac = 0;
        double d3dfac = 0;
        if (infile != null) {
            double[][] readx = new double[3][0];
            readFile(infile, readx);
            x = readx[0];
            y = readx[1];
            dy = readx[2];
        }
        else {
            int N = 50;
            x = new double[N];
            y = new double[N];
            dy = new double[N];
            for (int i=0; i<N; i++) {
                x[i] = i;
                dy[i] = i*i;
                y[i] = i*i*i*i;
            }
        }
        fitter = new AkimaSplineSmoother(new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray()));
        fitter.setD2fac(d2fac);
        fitter.setD2dfac(d2dfac);
        fitter.setD3fac(d3fac);
        fitter.setD3dfac(d3dfac);
        fitter.setInputData(x, y, dy);

        panel = new JPanel(new MigLayout());
        
        JPanel plotPanel = new JPanel(new MigLayout("flowy"));

        DisplayPlotXChart plot = new DisplayPlotXChart();
        plotPanel.add(plot.graphic());

        final int nSubPoints = 10;
        
        JTabbedPane morePlotPanel = new JTabbedPane();
        DisplayPlotXChart plotdy = new DisplayPlotXChart();
        morePlotPanel.add(plotdy.graphic(), "dy");
        plotdy.setDoLegend(false);
        DisplayPlotXChart plotdy2 = new DisplayPlotXChart();
        plotdy2.setDoLegend(false);
        morePlotPanel.add(plotdy2.graphic(), "dy2");
        DisplayPlotXChart plotey = new DisplayPlotXChart();
        plotey.setDoLegend(false);
        morePlotPanel.add(plotey.graphic(), "ey");
        
        IntegratorSmoother integrator = new IntegratorSmoother(fitter);
        ActivityIntegrate ai = new ActivityIntegrate(integrator);
        Controller controller = new Controller();
        controller.addAction(ai);
        
        DisplayPlotXChart ePlot = new DisplayPlotXChart();
        ePlot.getPlot().setYLog(true);
        morePlotPanel.add(ePlot.graphic(), "err");
        DataSourceCountSteps counter = new DataSourceCountSteps(integrator);
        historyEy = new AccumulatorHistory(new HistoryScrolling());
        historyEy.setTimeDataSource(counter);
        DataInfoDouble eyInfo = new DataInfoDouble("ey", Null.DIMENSION);
        historyEy.putDataInfo(eyInfo);
        historyED2 = new AccumulatorHistory(new HistoryScrolling());
        historyED2.setTimeDataSource(counter);
        DataInfoDouble ed2Info = new DataInfoDouble("ed2", Null.DIMENSION);
        historyED2.putDataInfo(ed2Info);
        historyED2D = new AccumulatorHistory(new HistoryScrolling());
        historyED2D.setTimeDataSource(counter);
        DataInfoDouble ed2dInfo = new DataInfoDouble("ed2d", Null.DIMENSION);
        historyED2D.putDataInfo(ed2dInfo);
        historyED3 = new AccumulatorHistory(new HistoryScrolling());
        historyED3.setTimeDataSource(counter);
        DataInfoDouble eyd3Info = new DataInfoDouble("ed3", Null.DIMENSION);
        historyED3.putDataInfo(eyd3Info);
        historyED3D = new AccumulatorHistory(new HistoryScrolling());
        historyED3D.setTimeDataSource(counter);
        DataInfoDouble ed3dInfo = new DataInfoDouble("ed3d", Null.DIMENSION);
        historyED3D.putDataInfo(ed3dInfo);
        historyEy.setDataSink(ePlot.getDataSet().makeDataSink());
        historyED2.setDataSink(ePlot.getDataSet().makeDataSink());
        historyED2D.setDataSink(ePlot.getDataSet().makeDataSink());
        historyED3.setDataSink(ePlot.getDataSet().makeDataSink());
        historyED3D.setDataSink(ePlot.getDataSet().makeDataSink());

        plotPanel.add(morePlotPanel);

        sink0 = plot.getDataSet().makeDataSink();
        sinkDy = plotdy.getDataSet().makeDataSink();
        sinkDy2 = plotdy2.getDataSet().makeDataSink();
        sinkEy = plotey.getDataSet().makeDataSink();
        sinkSmoothed = plot.getDataSet().makeDataSink();
        init(nSubPoints);

        plotAction = new PlotAction(nSubPoints);
        integrator.getEventManager().addListener(new IntegratorListenerAction(plotAction));
        plotAction.actionPerformed();
        
        JPanel controlPanel = new JPanel(new MigLayout("flowy"));

        JButton loadButton = new JButton();
        loadButton.setAction(new AbstractAction("Load data") {
            
            public void actionPerformed(ActionEvent e) {
                JFileChooser fileChooser = new JFileChooser();
                int rv = fileChooser.showOpenDialog(panel);
                
                if (rv == JFileChooser.APPROVE_OPTION) {
                    File file = fileChooser.getSelectedFile();
                    double[][] readx = new double[3][0];
                    readFile(file.getAbsolutePath(), readx);
                    x = readx[0];
                    y = readx[1];
                    dy = readx[2];
                    fitter.setInputData(x, y, dy);
                    AkimaSplineSmootherApp.this.init(nSubPoints);
                    plotAction.actionPerformed();
                }
           }
        });
        controlPanel.add(loadButton);

        JButton padButton = new JButton();
        padButton.setAction(new AbstractAction("Double data") {

            public void actionPerformed(ActionEvent e) {
                AkimaSpline spline = new AkimaSpline();
                spline.setInputData(x, fitter.y0);

                double[] newx = new double[x.length*2-1];
                double[] newdy = new double[newx.length];
                for (int i=0; i<x.length-1; i++) {
                    newx[2*i] = x[i];
                    newdy[2*i] = dy[i];
                    newx[2*i+1] = 0.5*(x[i]+x[i+1]);
                    newdy[2*i+1] = Double.MAX_VALUE;
                }
                newx[newx.length-1] = x[x.length-1];
                newdy[newx.length-1] = dy[x.length-1];
                y = spline.doInterpolation(newx);
                x = newx;
                dy = newdy;
                fitter.setInputData(x, y, dy);
                AkimaSplineSmootherApp.this.init(nSubPoints);
                plotAction.actionPerformed();
           }
        });
        controlPanel.add(padButton);

        DeviceControllerButton startButton = new DeviceControllerButton(controller);
        controlPanel.add(startButton.graphic());
        
        JPanel dPanel = new JPanel(new GridLayout(2, 2));
        
        DeviceBox d2Box = new DeviceBox();
        d2Box.setController(controller);
        d2Box.setModifier(new ModifierGeneral(fitter, "d2fac"));
        d2Box.setLabel("d2");
        d2Box.setLabelPosition(CompassDirection.NORTH);
        d2Box.setLabelType(LabelType.BORDER);
        dPanel.add(d2Box.graphic());
        
        DeviceBox d2dBox = new DeviceBox();
        d2dBox.setController(controller);
        d2dBox.setModifier(new ModifierGeneral(fitter, "d2dfac"));
        d2dBox.setLabel("d2d");
        d2dBox.setLabelPosition(CompassDirection.NORTH);
        d2dBox.setLabelType(LabelType.BORDER);
        dPanel.add(d2dBox.graphic());

        DeviceBox d3Box = new DeviceBox();
        d3Box.setController(controller);
        d3Box.setModifier(new ModifierGeneral(fitter, "d3fac"));
        d3Box.setLabel("d3");
        d3Box.setLabelPosition(CompassDirection.NORTH);
        d3Box.setLabelType(LabelType.BORDER);
        dPanel.add(d3Box.graphic());
        
        DeviceBox d3dBox = new DeviceBox();
        d3dBox.setController(controller);
        d3dBox.setModifier(new ModifierGeneral(fitter, "d3dfac"));
        d3dBox.setLabel("d3d");
        d3dBox.setLabelPosition(CompassDirection.NORTH);
        d3dBox.setLabelType(LabelType.BORDER);
        dPanel.add(d3dBox.graphic());
        
        controlPanel.add(dPanel);
        
        DeviceSlider tiltSlider = new DeviceSlider(controller);
        tiltSlider.setModifier(new ModifierGeneral(plotAction, "tilt"));
        tiltSlider.setMinimum(-10);
        tiltSlider.setMaximum(10);
        tiltSlider.setNMajor(5);
        tiltSlider.setPrecision(1);
        tiltSlider.setPostAction(plotAction);
        tiltSlider.setShowBorder(true);
        tiltSlider.setShowValues(true);
        tiltSlider.setEditValues(true);
        controlPanel.add(tiltSlider.graphic());
        
        
        panel.add(plotPanel);
        panel.add(controlPanel);
    }
    
    protected void init(int nSubPoints) {
        DataSourceIndependentSimple dataSourceX0 = new DataSourceIndependentSimple(x, new DataInfoDoubleArray("x0", Null.DIMENSION, new int[]{x.length}));
        DataInfoFunction dataInfo0 = new DataInfoFunction("y0", Null.DIMENSION, dataSourceX0);
        sink0.putDataInfo(dataInfo0);
        y0 = new DataFunction(new int[]{y.length});
        System.arraycopy(y, 0, y0.getData(), 0, y.length);
        sink0.putData(y0);

        xd = new double[(x.length-1)*nSubPoints+1];
        for (int i=0; i<x.length; i++) {
            xd[i*nSubPoints] = x[i];
            if (i<x.length-1) {
                for (int j=0; j<nSubPoints; j++) {
                    xd[i*nSubPoints+j] = x[i] + j*(x[i+1]-x[i])/nSubPoints;
                }
            }
        }

        double[][] dy12 = fitter.getDy12(nSubPoints);
        DataSourceIndependentSimple dataSourceXd = new DataSourceIndependentSimple(xd, new DataInfoDoubleArray("xd", Null.DIMENSION, new int[]{(x.length-1)*nSubPoints+1}));
        DataInfoFunction dataInfoDy = new DataInfoFunction("dy", Null.DIMENSION, dataSourceXd);
        sinkDy.putDataInfo(dataInfoDy);
        dataDy = new DataFunction(new int[]{xd.length}, dy12[0]);
        DataInfoFunction dataInfoDy2 = new DataInfoFunction("dy2", Null.DIMENSION, dataSourceXd);
        sinkDy2.putDataInfo(dataInfoDy2);
        dataDy2 = new DataFunction(new int[]{xd.length}, dy12[1]);
        DataInfoFunction dataInfoEy = new DataInfoFunction("ey", Null.DIMENSION, dataSourceX0);
        sinkEy.putDataInfo(dataInfoEy);
        dataEy = new DataFunction(new int[]{x.length});
        
        DataInfoFunction dataInfoSmooth = new DataInfoFunction("smoothed", Null.DIMENSION, dataSourceX0);
        sinkSmoothed.putDataInfo(dataInfoSmooth);
        ySmooth = new DataFunction(new int[]{y.length});
        System.arraycopy(fitter.y, 0, ySmooth.getData(), 0, fitter.y.length);

        historyEy.reset();
        historyED2.reset();
        historyED2D.reset();
        historyED3.reset();
        historyED3D.reset();
        fitter.calcErr(0, x.length-1);
    }
    
    public void go() {

        JFrame f = new JFrame();
        f.setSize(700,500);
        f.getContentPane().add(panel);
        f.pack();
        f.setTitle("Smoother");
        f.setVisible(true);
        f.addWindowListener(SimulationGraphic.WINDOW_CLOSER);
    }
    
    public static void readFile(String infile, double[][]x) {
        FileReader fileReader;
        try {
            fileReader = new FileReader(infile);
        }catch(IOException e) {
            throw new RuntimeException("Cannot open "+infile+", caught IOException: " + e.getMessage());
        }
        ArrayList<Double>[] xLists = new ArrayList[x.length];
        for (int i=0; i<x.length; i++) {
            xLists[i] = new ArrayList<Double>();
        }
        try {
            BufferedReader bufReader = new BufferedReader(fileReader);
            while (true) {
                String line = bufReader.readLine();
                if (line == null) {
                    break;
                }
                String[] xstr = line.split("[	 ]+");
                for (int i=0; i<x.length; i++) {
                    xLists[i].add(Double.parseDouble(xstr[i]));
                }
            }
            fileReader.close();
        } catch(IOException e) {
            throw new RuntimeException("Problem reading d.dat, caught IOException: " + e.getMessage());
        }

        for (int j=0; j<x.length; j++) {
            x[j] = new double[xLists[j].size()];
            for (int i=0; i<xLists[j].size(); i++) {
                x[j][i] = xLists[j].get(i);
            }
        }
    }

    public final class PlotAction implements IAction {
        protected final DataDouble dd = new DataDouble();
        protected final int nSubPoints;

        public PlotAction(int nSubPoints) {
            this.nSubPoints = nSubPoints;
        }

        public void actionPerformed() {
            double[] smoothy = ySmooth.getData();
            double[] originalY = y0.getData();
            for (int i=0; i<smoothy.length; i++) {
                smoothy[i] = fitter.y[i] + tilt*fitter.x[i];
            }
            sinkSmoothed.putData(ySmooth);
            for (int i=0; i<smoothy.length; i++) {
                originalY[i] = y[i] + tilt*fitter.x[i];
            }
            sink0.putData(y0);
            fitter.getDy12(nSubPoints);
            sinkDy.putData(dataDy);
            sinkDy2.putData(dataDy2);
            
            double[] ey = dataEy.getData();
            for (int i=0; i<ey.length; i++) {
                ey[i] = (fitter.y[i] - fitter.y0[i]) / fitter.ey[i];
            }
            sinkEy.putData(dataEy);

            dd.x = Math.min(fitter.sumSqDy, 4*x.length);
            if (dd.x != 0) historyEy.putData(dd);
            dd.x = Math.min(fitter.d2fac*fitter.sumSqD2, 4*x.length);
            if (dd.x != 0) historyED2.putData(dd);
            dd.x = Math.min(fitter.d2dfac*fitter.sumSqD2D, 4*x.length);
            if (dd.x != 0) historyED2D.putData(dd);
            dd.x = Math.min(fitter.d3fac*fitter.sumSqD3, 4*x.length);
            if (dd.x != 0) historyED3.putData(dd);
            dd.x = Math.min(fitter.d3dfac*fitter.sumSqD3D, 4*x.length);
            if (dd.x != 0) historyED3D.putData(dd);
        }
        
        public double getTilt() {
            return tilt;
        }
        
        public void setTilt(double newTilt) {
            tilt = newTilt;
        }
        
        protected double tilt;
    }

    public static class IntegratorSmoother extends Integrator {
        
        protected final AkimaSplineSmoother fitter;

        public IntegratorSmoother(AkimaSplineSmoother fitter) {
            this.fitter = fitter;
        }

        protected void doStepInternal() {
            fitter.doStep();
        }
    }
    
    public static void main(String[] args) {
        SimulationGraphic.initGraphics();
        String infile = null;
        if (args.length > 0) {
            infile = args[0];
        }
        AkimaSplineSmootherApp smoother = new AkimaSplineSmootherApp(infile);
        smoother.go();
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            getRootPane().putClientProperty(
                            "defeatSystemEventQueueCheck", Boolean.TRUE);
            AkimaSplineSmootherApp smoother = new AkimaSplineSmootherApp(null);
            getContentPane().add(smoother.panel);
            
            smoother.x = new double[]{-1.6666666667e+00,-1.6129032258e+00,-1.5625000000e+00,-1.5151515152e+00,-1.4705882353e+00,-1.4285714286e+00,-1.3888888889e+00,-1.3513513514e+00,-1.3157894737e+00,-1.2820512821e+00,-1.2500000000e+00,-1.2195121951e+00,-1.1904761905e+00,-1.1627906977e+00,-1.1363636364e+00,-1.1111111111e+00,-1.0869565217e+00,-1.0638297872e+00,-1.0416666667e+00,-1.0204081633e+00,-1.0000000000e+00,-9.8039215686e-01,-9.6153846154e-01,-9.4339622642e-01,-9.2592592593e-01,-9.0909090909e-01,-8.9285714286e-01,-8.7719298246e-01,-8.6206896552e-01,-8.4745762712e-01,-8.3333333333e-01,-8.1967213115e-01,-8.0645161290e-01,-7.9365079365e-01,-7.8125000000e-01,-7.6923076923e-01,-7.5757575758e-01,-7.4626865672e-01,-7.3529411765e-01,-7.1428571429e-01,-6.9444444444e-01,-6.7567567568e-01,-6.5789473684e-01,-6.4102564103e-01,-6.2500000000e-01,-6.0975609756e-01,-5.9523809524e-01,-5.8139534884e-01,-5.6818181818e-01,-5.5555555556e-01,-5.4347826087e-01,-5.3191489362e-01,-5.2083333333e-01,-5.1020408163e-01,-5.0000000000e-01,-4.7619047619e-01,-4.5454545455e-01,-4.3478260870e-01,-4.1666666667e-01,-4.0000000000e-01,-3.8461538462e-01,-3.7037037037e-01,-3.5714285714e-01,-3.4482758621e-01,-3.3333333333e-01,-3.2258064516e-01,-3.1250000000e-01,-3.0303030303e-01,-2.9411764706e-01,-2.8571428571e-01,-2.7777777778e-01,-2.7027027027e-01,-2.6315789474e-01,-2.5641025641e-01,-2.5000000000e-01,-2.3809523810e-01,-2.2727272727e-01,-2.0833333333e-01,-1.9230769231e-01,-1.7857142857e-01,-1.6666666667e-01,-1.4285714286e-01,-1.2500000000e-01,-1.1111111111e-01,-1.0000000000e-01,-8.3333333333e-02,-7.1428571429e-02,-5.5555555556e-02,-4.1666666667e-02,-3.3333333333e-02,-2.5000000000e-02};
            smoother.y = new double[]{-6.1422267036e-01,-7.0979820885e-01,-8.1214309243e-01,-9.0264812602e-01,-9.8870003338e-01,-1.0723191225e+00,-1.1468181632e+00,-1.2044325134e+00,-1.2539617400e+00,-1.2928189928e+00,-1.2785271507e+00,-1.3001540490e+00,-1.2920628768e+00,-1.2451229053e+00,-1.2063836031e+00,-1.1491492624e+00,-1.0905976306e+00,-1.0142341575e+00,-9.5546918915e-01,-8.5410991248e-01,-7.6412883785e-01,-6.4901129665e-01,-5.3996084217e-01,-4.6705529164e-01,-3.7523523968e-01,-2.9056816680e-01,-2.0809476735e-01,-1.3830084822e-01,-7.0736258166e-02,5.2635919661e-03,4.0404936358e-02,3.7147126461e-02,1.0602108325e-01,5.3659697495e-02,1.1511666345e-01,7.1022588570e-02,4.5786096319e-02,8.1176064841e-02,-5.2762396897e-03,-1.4710589216e-02,-4.7816452212e-02,-1.1814040674e-01,-1.3006561546e-01,-1.6176291204e-01,-2.0888008981e-01,-2.0285001956e-01,-2.1131123419e-01,-2.0410598972e-01,-2.0556782654e-01,-1.7979625444e-01,-1.6945712322e-01,-1.2359789034e-01,-1.0067239345e-01,-9.0908092690e-02,-6.9973814412e-02,1.1932950580e-02,6.8034845460e-02,1.4913473169e-01,2.2331032865e-01,2.7691581834e-01,3.3694051814e-01,3.8246193355e-01,4.3293595540e-01,4.7179680273e-01,5.1230431321e-01,5.5065624927e-01,5.8458979087e-01,6.1519849731e-01,6.4633851512e-01,6.7001783335e-01,6.9227275603e-01,7.1345377858e-01,7.3676633372e-01,7.5073589005e-01,7.6804197817e-01,8.0030352699e-01,8.2716751223e-01,8.6828788507e-01,9.0183891917e-01,9.2608239449e-01,9.4735929061e-01,9.8024220986e-01,1.0007866061e+00,1.0147891716e+00,1.0223257355e+00,1.0287586554e+00,1.0306076898e+00,1.0262115347e+00,1.0130388579e+00,1.0025311754e+00,9.8343399712e-01};
            smoother.dy = new double[]{1.7360804819e-03,2.0466758601e-03,2.4485061902e-03,2.8814900151e-03,3.4257619928e-03,4.0006538066e-03,4.5151713859e-03,5.1317461628e-03,5.9176354220e-03,6.7556728761e-03,7.2693400606e-03,7.9790033741e-03,8.9730856527e-03,1.0121504787e-02,1.0479268868e-02,1.1412959236e-02,1.2779816970e-02,1.3753763546e-02,1.4056585859e-02,1.5552399776e-02,1.6692384991e-02,1.7256220816e-02,1.8881300049e-02,1.9907941008e-02,2.0930192394e-02,2.2080791161e-02,2.3561208000e-02,2.3558356190e-02,2.5055564352e-02,2.6157097218e-02,2.5928497256e-02,2.6676915480e-02,2.6891526216e-02,2.6487185185e-02,2.7937294540e-02,2.7741421475e-02,2.7024527968e-02,2.5779437406e-02,2.5732147904e-02,2.2856090290e-02,2.0673853423e-02,2.0516550173e-02,1.8134873638e-02,1.6513727366e-02,1.5122436429e-02,1.3799857749e-02,1.2599665037e-02,1.1563311282e-02,1.0441848589e-02,9.5320690205e-03,8.8588556978e-03,8.1429740497e-03,7.6943319833e-03,7.2188886652e-03,6.3976848549e-03,5.6899057085e-03,4.8899620577e-03,4.3644665621e-03,3.7422871973e-03,3.3730655671e-03,3.0882221485e-03,2.7326201723e-03,2.5629492804e-03,2.3761776766e-03,2.1724281321e-03,2.0470209385e-03,1.9103601223e-03,1.8128686801e-03,1.7069891637e-03,1.6458584690e-03,1.5791489363e-03,1.5068887470e-03,1.4291046617e-03,1.4043361866e-03,1.3169259316e-03,1.2508504933e-03,1.2390225059e-03,1.1282162067e-03,1.0607846124e-03,1.0528182225e-03,1.0339660066e-03,9.4012110671e-04,9.3719575428e-04,9.3016362881e-04,9.1342685013e-04,8.9931102371e-04,9.0375627464e-04,9.2459444457e-04,9.2784447732e-04,9.1432213075e-04,9.4713002952e-04};
            smoother.fitter.setInputData(smoother.x, smoother.y, smoother.dy);
            smoother.init(10);
            smoother.plotAction.actionPerformed();
        }
    }
}
