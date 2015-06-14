package etomica.normalmode.nptdemo.fluid;
import javax.swing.JFrame;

import etomica.action.IAction;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.atom.DiameterHashByType;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.IData;
import etomica.data.IEtomicaDataSource;
import etomica.graphics.ColorScheme;
import etomica.graphics.ColorSchemeOverlap;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DeviceToggleButton;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayCanvas;
import etomica.graphics.SimulationGraphic;
import etomica.listener.IntegratorListenerAction;
import etomica.modifier.ModifierBoolean;

/**
 * Simple hard-sphere molecular dynamics simulation in 2D.
 *
 * @author David Kofke
 */
 
public class HSNPT2DGraphic extends SimulationGraphic {
    
    public double nominalScale;
    public DisplayBox displayBoxScaled;
    public HS2DFluid sim;
    protected final ColorScheme redColorScheme;
    protected final ColorSchemeOverlap colorSchemeOverlap;
    protected final ColorSchemeFluidOverlap colorSchemeScaled;

    public HSNPT2DGraphic(final HS2DFluid sim) {
        super(sim, TABBED_PANE, "HS Fluid", 1, sim.getSpace(), sim.getController());
        this.sim = sim;
        getDisplayBox(sim.box).setScale(0.8);
        redColorScheme = getDisplayBox(sim.box).getColorScheme();

        displayBoxScaled = new DisplayBox(sim, sim.box, sim.getSpace(), sim.getController());
        displayBoxScaled.setScale(0.9);
        System.out.println("displayBoxScaled created");
        DisplayCanvas oldCanvas = displayBoxScaled.canvas;
        final IStuff stuff = new Stuff2(space, sim.potential, sim.integrator.getTemperature(), sim.box, 40, 1000);
        ((Stuff2)stuff).setRCut(6);
        stuff.setTemperature(1);
        displayBoxScaled.setBoxCanvas(new DisplayBoxCanvas2DFluidScaling(displayBoxScaled, sim.getSpace(), sim.getController(), stuff));
        oldCanvas.dispose();
        System.out.println("displayBoxScaled has a new canvas");
        addAsTab(displayBoxScaled.graphic(), "scaled", true);
        IntegratorListenerAction repaintAction = new IntegratorListenerAction(new IAction() {
            public void actionPerformed() {
                displayBoxScaled.repaint();
            }
        });
        repaintAction.setInterval(1);
        sim.integrator.getEventManager().addListener(repaintAction);
        displayBoxScaled.setDiameterHash(getDisplayBox(sim.box).getDiameterHash());
        
        colorSchemeOverlap = new ColorSchemeOverlap(sim.getSpace(), sim.potentialMaster, sim.box);
        getDisplayBox(sim.box).setColorScheme(colorSchemeOverlap);
        
        colorSchemeScaled = new ColorSchemeFluidOverlap(sim.getSpace(), sim.potentialMaster, sim.box, stuff);
        displayBoxScaled.setColorScheme(colorSchemeScaled);

        getController().getReinitButton().setPostAction(new IAction() {
            public void actionPerformed() {
                getDisplayBox(simulation.getBox(0)).repaint();
                displayBoxScaled.repaint();
            }
        });
        
        nominalScale = getDisplayBox(sim.box).getScale();

        DeviceSlider rhoSlider = new DeviceSlider(sim.getController(), this, "displayDensity");
        rhoSlider.setMinimum(0.1);
        rhoSlider.setMaximum(0.9);
        rhoSlider.setPrecision(3);
        rhoSlider.setMinimum(0.1);
        rhoSlider.setMaximum(0.9);
        rhoSlider.setNMajor(2);
        rhoSlider.setNMinor(3);
        rhoSlider.doUpdate();
        rhoSlider.setShowBorder(true);
        rhoSlider.setShowValues(true);
        rhoSlider.setLabel("density");
        add(rhoSlider);
        setDisplayDensity(0.2);
        
        DeviceToggleButton colorButton = new DeviceToggleButton(sim.getController(), new ModifierBoolean() {
            
            public void setBoolean(boolean b) {
                if (b) {
                    getDisplayBox(simulation.getBox(0)).setColorScheme(redColorScheme);
                    displayBoxScaled.setColorScheme(redColorScheme);
                }
                else {
                    getDisplayBox(simulation.getBox(0)).setColorScheme(colorSchemeOverlap);
                    displayBoxScaled.setColorScheme(colorSchemeScaled);
                }
            }
            
            public boolean getBoolean() {
                return displayBoxScaled.getColorScheme() == redColorScheme;
            }
        }, "Color: overlaps", "Color: red");
        add(colorButton);
        
        final MeterPressureByVolumeChangeStuff meterPBVS = new MeterPressureByVolumeChangeStuff(sim, space, stuff);
        meterPBVS.setX(Math.log(1/1.05), Math.log(1/0.95), 2);
        meterPBVS.setIntegrator(sim.integrator);
        final AccumulatorAverageFixed accPBVS = new AccumulatorAverageFixed(100);
        DataPumpListener pumpPBVS = new DataPumpListener(meterPBVS, accPBVS, 1);
        sim.integrator.getEventManager().addListener(pumpPBVS);

        sim.integrator.getEventManager().addListener(new IIntegratorListener() {
            int countdown = 100;
            public void integratorStepStarted(IIntegratorEvent e) {}
            
            public void integratorStepFinished(IIntegratorEvent e) {
                if (--countdown>0) return;
                countdown = 100;
                IEtomicaDataSource xDataSource = meterPBVS.getScalingDataSource();
                IData xData = xDataSource.getData();
                int last = xData.getLength()-1;

                double avg0Stuff= accPBVS.getData(accPBVS.AVERAGE).getValue(0);
                double err0Stuff = accPBVS.getData(accPBVS.ERROR).getValue(0);
                double avg10Stuff = accPBVS.getData(accPBVS.AVERAGE).getValue(last);
                double err10Stuff = accPBVS.getData(accPBVS.ERROR).getValue(last);
                double dv0 = sim.box.getBoundary().volume()*(xData.getValue(0)-1);
                double dvN = sim.box.getBoundary().volume()*(xData.getValue(last)-1);

                System.out.print(String.format("stuff   % 7.4f  %8.2e    % 7.4f  %8.2e\n", Math.log(avg0Stuff)/dv0, err0Stuff/avg0Stuff/Math.abs(dv0), Math.log(avg10Stuff)/dvN, err10Stuff/avg10Stuff/Math.abs(dvN)));
            }
            
            public void integratorInitialized(IIntegratorEvent e) {}
        });
    }
    
    public void setDisplayDensity(double newDisplayDensity) {
        double density = sim.box.getMoleculeList().getMoleculeCount() / sim.box.getBoundary().volume();
        double newScale = Math.sqrt(density/newDisplayDensity) * nominalScale;
        DisplayBox unscaledDisplayBox = getDisplayBox(sim.box);
        unscaledDisplayBox.setScale(newScale);
        displayBoxScaled.setScale(newScale);
        ((DisplayBoxCanvas2DFluidScaling)displayBoxScaled.canvas).setDisplayDensity(newDisplayDensity);
        colorSchemeScaled.setDisplayDensity(newDisplayDensity);
        double displayDiameter = nominalScale / newScale;
        ((DiameterHashByType)unscaledDisplayBox.getDiameterHash()).setDiameter(sim.getSpecies(0).getAtomType(0), displayDiameter);
        colorSchemeOverlap.setSigma(displayDiameter);
        unscaledDisplayBox.repaint();
        displayBoxScaled.repaint();
    }
    
    public double getDisplayDensity() {
        double scale = getDisplayBox(sim.box).getScale();
        double density = sim.box.getMoleculeList().getMoleculeCount() / sim.box.getBoundary().volume();
        double s = scale/nominalScale;
        return density /(s*s); 
    }

    public static class Applet extends javax.swing.JApplet {
        public void init() {
            HS2DFluid sim = new HS2DFluid();
            sim.ai.setSleepPeriod(1);
            
            HSNPT2DGraphic graphic = new HSNPT2DGraphic(sim);
            getContentPane().add(graphic.getPanel());
        }

        private static final long serialVersionUID = 1L;
    }

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        HS2DFluid sim = new HS2DFluid();
        sim.ai.setSleepPeriod(1);
        HSNPT2DGraphic graphic = new HSNPT2DGraphic(sim);
        JFrame f = graphic.makeAndDisplayFrame();
        f.setSize(700, 500);

//        MeterPressureHard meterP = new MeterPressureHard(sim.space);
//        meterP.setIntegrator(sim.integrator);
//        sim.integrator.addCollisionListener(meterP);
//        sim.ai.setMaxSteps(10000);
//        sim.getController().actionPerformed();
//        System.out.println("P: "+meterP.getDataAsScalar());
//        System.exit(1);

   }
}