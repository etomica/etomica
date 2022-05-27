/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.controller.Controller;
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
public class AkimaSplineSmootherDyApp {
    
    final JPanel panel;
    protected final AccumulatorHistory historyEy, historyED2, historyED2D, historyED3, historyED3D, historyTrueEy;
    protected final AkimaSplineSmootherDy fitter;
    protected IDataSink sinkDy, sinkDy2, sinkDy3, sink0, sinkSmoothed, sinkEy, sinkTrue, sinkTrueEy;
    protected DataFunction ySmooth, dataDy, dataDy2, dataDy3, y0, dataEy, yTrue, dataTrueEy;
    protected double[] xd;
    protected double[] x, y, dy, trueY;
    protected boolean haveTrueData = false;
    protected int nPretendPoints;
    protected final PlotAction plotAction;

    public AkimaSplineSmootherDyApp(String infile) {
        nPretendPoints = 1;
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
        fitter = new AkimaSplineSmootherDy(new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray()));
        fitter.setD2fac(d2fac);
        fitter.setD2dfac(d2dfac);
        fitter.setD3fac(d3fac);
        fitter.setD3dfac(d3dfac);
        fitter.setInputData(x, y, dy);

        panel = new JPanel();
        
        JPanel plotPanel = new JPanel(new GridBagLayout());
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.gridx = 0;

        final DisplayPlot plot = new DisplayPlot();
        plot.getDataSet().setUpdatingOnAnyChange(true);
        plotPanel.add(plot.graphic(), gbc);

        final int nSubPoints = 10;
        
        JTabbedPane morePlotPanel = new JTabbedPane();
        DisplayPlot plotdy = new DisplayPlot();
        morePlotPanel.add(plotdy.graphic(), "dy");
        plotdy.setDoLegend(false);
        DisplayPlot plotdy2 = new DisplayPlot();
        morePlotPanel.add(plotdy2.graphic(), "dy2");
        plotdy2.setDoLegend(false);
        DisplayPlot plotdy3 = new DisplayPlot();
        plotdy3.setDoLegend(false);
        morePlotPanel.add(plotdy3.graphic(), "dy3");
        DisplayPlot plotey = new DisplayPlot();
        morePlotPanel.add(plotey.graphic(), "ey");
        
        IntegratorSmoother integrator = new IntegratorSmoother(fitter);

        // TODO !!
        Controller controller = new Controller();
        controller.addActivity(new ActivityIntegrate(integrator));
        controller.start();

        DisplayPlot ePlot = new DisplayPlot();
        ePlot.getPlot().setYLog(true);
        ePlot.getDataSet().setUpdatingOnAnyChange(true);
        morePlotPanel.add(ePlot.graphic(), "err");
        DataSourceCountSteps counter = new DataSourceCountSteps(integrator);
        historyEy = new AccumulatorHistory(new HistoryScrolling());
        historyEy.setTimeDataSource(counter);
        DataInfoDouble eyInfo = new DataInfoDouble("ey", Null.DIMENSION);
        historyEy.putDataInfo(eyInfo);
        historyTrueEy = new AccumulatorHistory(new HistoryScrolling());
        historyTrueEy.setTimeDataSource(counter);
        DataInfoDouble eTrueYInfo = new DataInfoDouble("true ey", Null.DIMENSION);
        historyTrueEy.putDataInfo(eTrueYInfo);
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
        historyTrueEy.setDataSink(ePlot.getDataSet().makeDataSink());
        historyED2.setDataSink(ePlot.getDataSet().makeDataSink());
        historyED2D.setDataSink(ePlot.getDataSet().makeDataSink());
        historyED3.setDataSink(ePlot.getDataSet().makeDataSink());
        historyED3D.setDataSink(ePlot.getDataSet().makeDataSink());

        plotPanel.add(morePlotPanel, gbc);

        sink0 = plot.getDataSet().makeDataSink();
        sinkDy = plotdy.getDataSet().makeDataSink();
        sinkDy2 = plotdy2.getDataSet().makeDataSink();
        sinkDy3 = plotdy3.getDataSet().makeDataSink();
        sinkEy = plotey.getDataSet().makeDataSink();
        sinkTrueEy = plotey.getDataSet().makeDataSink();
        sinkSmoothed = plot.getDataSet().makeDataSink();
        sinkTrue = plot.getDataSet().makeDataSink();
        init(nSubPoints);

        plotAction = new PlotAction(nSubPoints);
        integrator.getEventManager().addListener(new IntegratorListenerAction(plotAction));
        plotAction.actionPerformed();
        
        JPanel controlPanel = new JPanel(new GridBagLayout());

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
                    AkimaSplineSmootherDyApp.this.init(nSubPoints);
                    plotAction.actionPerformed();
                }
           }
        });
        controlPanel.add(loadButton, gbc);

        JButton padButton = new JButton();
        padButton.setAction(new AbstractAction("Double data") {
            
            public void actionPerformed(ActionEvent e) {
                nPretendPoints *= 2;
                AkimaSpline spline = new AkimaSpline();
                spline.setInputData(x, fitter.y);
                
                double[] newx = new double[x.length*2-1];
                double[] newdy = new double[newx.length];
                double[] newy = new double[newx.length];
                for (int i=0; i<x.length-1; i++) {
                    newx[2*i] = x[i];
                    newdy[2*i] = dy[i];
                    newy[2*i] = y[i];
                    newx[2*i+1] = 0.5*(x[i]+x[i+1]);
                    newdy[2*i+1] = Double.MAX_VALUE;
                    newy[2*i+1] = 0.5*(y[i]+y[i+1]);
                }
                newx[newx.length-1] = x[x.length-1];
                newdy[newx.length-1] = dy[x.length-1];
                newy[newx.length-1] = y[x.length-1];
                double[] dydx = spline.doInterpolation(newx);
                x = newx;
                y = newy;
                dy = newdy;
                fitter.setInputData(x, y, dydx, dy);
                AkimaSplineSmootherDyApp.this.init(nSubPoints);
                plotAction.actionPerformed();
           }
        });
        controlPanel.add(padButton, gbc);

        JButton loadTrueButton = new JButton();
        loadTrueButton.setAction(new AbstractAction("Load true data") {
            
            public void actionPerformed(ActionEvent e) {
                JFileChooser fileChooser = new JFileChooser();
                int rv = fileChooser.showOpenDialog(panel);
                
                if (rv == JFileChooser.APPROVE_OPTION) {
                    File file = fileChooser.getSelectedFile();
                    double[][] readx = new double[2][0];
                    readFile(file.getAbsolutePath(), readx);
                    trueY = readx[1];

                    System.arraycopy(trueY, 0, yTrue.getData(), 0, trueY.length);
                    sinkTrue.putData(yTrue);
                    haveTrueData = true;
                    plotAction.actionPerformed();
                }
           }
        });
        controlPanel.add(loadTrueButton, gbc);
        
        DeviceRunControls startButton = new DeviceRunControls(controller);
        controlPanel.add(startButton.graphic(), gbc);
        
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
        
        controlPanel.add(dPanel, gbc);
        
        DeviceSlider tiltSlider = new DeviceSlider(controller);
        tiltSlider.setModifier(new ModifierGeneral(plotAction, "tilt"));
        tiltSlider.setNMajor(5);
        tiltSlider.setPrecision(2);
        tiltSlider.setMinimum(-20);
        tiltSlider.setMaximum(20);
//        tiltSlider.setPostAction(plotAction);
        tiltSlider.setShowBorder(true);
        tiltSlider.setShowValues(true);
        tiltSlider.setEditValues(true);
        controlPanel.add(tiltSlider.graphic(), gbc);
        
        
        panel.add(plotPanel);
        panel.add(controlPanel);
    }
    
    protected void init(int nSubPoints) {
        DataSourceIndependentSimple dataSourceX0 = new DataSourceIndependentSimple(x, new DataInfoDoubleArray("x0", Null.DIMENSION, new int[]{x.length}));
        double[] x0real = new double[(x.length+nPretendPoints-1)/nPretendPoints];
        for (int i=0; i<x0real.length; i++) {
            x0real[i] = x[i*nPretendPoints];
        }
        DataSourceIndependentSimple dataSourceX0Real = new DataSourceIndependentSimple(x0real, new DataInfoDoubleArray("x0Real", Null.DIMENSION, new int[]{(x.length+nPretendPoints-1)/nPretendPoints}));
        DataInfoFunction dataInfo0 = new DataInfoFunction("y0", Null.DIMENSION, dataSourceX0Real);
        sink0.putDataInfo(dataInfo0);
        y0 = new DataFunction(new int[]{dataInfo0.getLength()});
        double[] y0d = y0.getData();
        for (int i=0; i<y0d.length; i++) {
            y0d[i] = y[i*nPretendPoints];
        }
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
        DataInfoFunction dataInfoDy = new DataInfoFunction("dy", Null.DIMENSION, dataSourceX0);
        sinkDy.putDataInfo(dataInfoDy);
        dataDy = new DataFunction(new int[]{fitter.y.length}, fitter.y);
        DataInfoFunction dataInfoDy2 = new DataInfoFunction("dy2", Null.DIMENSION, dataSourceXd);
        sinkDy2.putDataInfo(dataInfoDy2);
        dataDy2 = new DataFunction(new int[]{xd.length}, dy12[0]);
        DataInfoFunction dataInfoDy3 = new DataInfoFunction("dy3", Null.DIMENSION, dataSourceXd);
        sinkDy3.putDataInfo(dataInfoDy3);
        dataDy3 = new DataFunction(new int[]{xd.length}, dy12[1]);
        DataInfoFunction dataInfoEy = new DataInfoFunction("ey", Null.DIMENSION, dataSourceX0);
        sinkEy.putDataInfo(dataInfoEy);
        dataEy = new DataFunction(new int[]{(x.length+nPretendPoints-1)/nPretendPoints});
        DataInfoFunction dataInfoTrueEy = new DataInfoFunction("true ey", Null.DIMENSION, dataSourceX0Real);
        sinkTrueEy.putDataInfo(dataInfoTrueEy);
        dataTrueEy = new DataFunction(new int[]{dataInfoTrueEy.getLength()});
        
        DataInfoFunction dataInfoSmooth = new DataInfoFunction("smoothed", Null.DIMENSION, dataSourceX0);
        sinkSmoothed.putDataInfo(dataInfoSmooth);
        ySmooth = new DataFunction(new int[]{y.length});
        DataInfoFunction dataInfoTrue = new DataInfoFunction("true", Null.DIMENSION, dataSourceX0Real);
        sinkTrue.putDataInfo(dataInfoTrue);
        yTrue = new DataFunction(new int[]{dataInfoTrue.getLength()});

        historyEy.reset();
        historyTrueEy.reset();
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
        f.setTitle("Dy Smoother");
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
            double[] iy = fitter.getIy();
            for (int i=0; i<smoothy.length; i++) {
                smoothy[i] = iy[i] + tilt*fitter.x[i];
            }
            sinkSmoothed.putData(ySmooth);
            for (int i=0; i<originalY.length; i++) {
                originalY[i] = y[nPretendPoints*i] + tilt*fitter.x[nPretendPoints*i];
            }
            sink0.putData(y0);
            if (!doJustTilt) {
                fitter.getDy12(nSubPoints);
                sinkDy.putData(dataDy);
                sinkDy2.putData(dataDy2);
                sinkDy3.putData(dataDy3);
                
                double[] ey = dataEy.getData();
                for (int i=0; i<ey.length; i++) {
                    ey[i] = (iy[nPretendPoints*i] - y[nPretendPoints*i]) / fitter.ey[nPretendPoints*i];
                }
                sinkEy.putData(dataEy);
            }

            if (haveTrueData) {
                double sumSqTrueDy = 0;
                double[] trueEy = dataTrueEy.getData();
                double[] yTrueArray = yTrue.getData();
                for (int i=0; i<trueEy.length; i++) {
                    trueEy[i] = (iy[nPretendPoints*i] - trueY[i]) / fitter.ey[nPretendPoints*i];
                    sumSqTrueDy += trueEy[i]*trueEy[i];
                    yTrueArray[i] = trueY[i] + tilt*fitter.x[nPretendPoints*i];
                }
                sinkTrue.putData(yTrue);
                if (!doJustTilt) {
                    sinkTrueEy.putData(dataTrueEy);
                    dd.x = Math.min(sumSqTrueDy, 400*x.length);
                    if (dd.x != 0) historyTrueEy.putData(dd);
                }
            }
            if (doJustTilt) {
                doJustTilt = false;
                return;
            }

            dd.x = Math.min(fitter.sumSqDy, 400*x.length);
            if (dd.x == 0) {
                dd.x = Double.NaN;
            }
            historyEy.putData(dd);
            dd.x = Math.min(fitter.d2fac*fitter.sumSqD2, 400*x.length);
            if (dd.x == 0) {
                dd.x = Double.NaN;
            }
            historyED2.putData(dd);
            dd.x = Math.min(fitter.d2dfac*fitter.sumSqD2D, 400*x.length);
            if (dd.x == 0) {
                dd.x = Double.NaN;
            }
            historyED2D.putData(dd);
            dd.x = Math.min(fitter.d3fac*fitter.sumSqD3, 400*x.length);
            if (dd.x == 0) {
                dd.x = Double.NaN;
            }
            historyED3.putData(dd);
            dd.x = Math.min(fitter.d3dfac*fitter.sumSqD3D, 400*x.length);
            if (dd.x == 0) {
                dd.x = Double.NaN;
            }
            historyED3D.putData(dd);
        }
        
        public double getTilt() {
            return tilt;
        }
        
        public void setTilt(double newTilt) {
            tilt = newTilt;
            doJustTilt = true;
            actionPerformed();
        }
        
        protected double tilt;
        protected boolean doJustTilt = false;
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
        String infile = null;
        if (args.length > 0) {
            infile = args[0];
        }
        AkimaSplineSmootherDyApp smoother = new AkimaSplineSmootherDyApp(infile);
        smoother.go();
    }

}
