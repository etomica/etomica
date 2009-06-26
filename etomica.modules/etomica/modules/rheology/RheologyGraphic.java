package etomica.modules.rheology;

import java.util.ArrayList;

import etomica.action.IAction;
import etomica.atom.AtomPair;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPumpListener;
import etomica.graphics.DeviceDelaySlider;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.modifier.ModifierGeneral;
import etomica.space3d.Space3D;
import g3dsys.images.Figure;

/**
 * Module for "Macromolecular dynamics related to rheological properties
 * through Brownian dynamics simulations", designed by Lewis Wedgewood.
 * 
 * @author Andrew Schultz
 */
public class RheologyGraphic extends SimulationGraphic {

    public RheologyGraphic(final SimulationRheology sim) {
        super(sim, SimulationGraphic.TABBED_PANE, "Polymer Rheology", 1, sim.getSpace(), sim.getController());
        sim.setChainLength(10);
        final Object bondObject = new Object();
        final ArrayList<Figure> bondList = new ArrayList<Figure>();
        for (int i=0; i<9; i++) {
            bondList.add((Figure)((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).makeBond(new AtomPair(sim.box.getLeafList().getAtom(i), sim.box.getLeafList().getAtom(i+1)), bondObject));
        }
        
        final DeviceSlider sliderA = new DeviceSlider(sim.getController());
        sliderA.setLabel("a");
        sliderA.setShowBorder(true);
        ModifierGeneral modifierA = new ModifierGeneral(sim.integrator, "a");
        sliderA.setModifier(modifierA);
        sliderA.setPrecision(1);
        sliderA.setMinimum(-1);
        sliderA.setMaximum(1);
        sliderA.setNMajor(4);
        sliderA.setShowValues(true);
        add(sliderA);

        final DeviceSlider sliderShear = new DeviceSlider(sim.getController());
        sliderShear.setShowBorder(true);
        ModifierGeneral modifierShear = new ModifierGeneral(sim.integrator, "shearRateNumber");
        sliderShear.setModifier(modifierShear);
        sliderShear.setLabel("Shear rate");
        sliderShear.setPrecision(2);
        sliderShear.setMinimum(0);
        sliderShear.setMaximum(10);
        sliderShear.setNMajor(4);
        sliderShear.setShowValues(true);
        sliderShear.setEditValues(true);
        add(sliderShear);
        
        DeviceSlider sliderPolymerLength= new DeviceSlider(sim.getController());
        sliderPolymerLength.setShowBorder(true);
        ModifierGeneral modifierLength = new ModifierGeneral(sim, "chainLength");
        sliderPolymerLength.setModifier(modifierLength);
        sliderPolymerLength.setLabel("Chain length");
        sliderPolymerLength.setMinimum(0);
        sliderPolymerLength.setMaximum(200);
        sliderPolymerLength.setNMajor(4);
        sliderPolymerLength.setShowValues(true);
        add(sliderPolymerLength);
        sliderPolymerLength.setPreAction(new IAction() {
            public void actionPerformed() {
                for (int i=0; i<bondList.size(); i++) {
                    ((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).releaseBond(bondList.get(i));
                }
                bondList.clear();
            }
        });
        sliderPolymerLength.setPostAction(new IAction() {
            public void actionPerformed() {
                for (int i=0; i<sim.box.getLeafList().getAtomCount()-1; i++) {
                    bondList.add((Figure)((DisplayBoxCanvasG3DSys)getDisplayBox(sim.box).canvas).makeBond(new AtomPair(sim.box.getLeafList().getAtom(i), sim.box.getLeafList().getAtom(i+1)), bondObject));
                }

                getPaintAction(sim.box).actionPerformed();
            }
        });

        final MeterViscosity meterViscosity = new MeterViscosity(sim.getSpace());
        meterViscosity.setIntegrator(sim.integrator);
        meterViscosity.setBox(sim.box);
        AccumulatorAverageCollapsing avgViscosity = new AccumulatorAverageCollapsing();
        DataPumpListener viscosityPump = new DataPumpListener(meterViscosity, avgViscosity, 10);
        sim.integrator.getEventManager().addListener(viscosityPump);
        DisplayTextBoxesCAE viscosityDisplay = new DisplayTextBoxesCAE();
        viscosityDisplay.setAccumulator(avgViscosity);
        add(viscosityDisplay);
        getController().getDataStreamPumps().add(viscosityPump);

        final MeterNormalStress meterNormalStress1 = new MeterNormalStress(sim.getSpace());
        meterNormalStress1.setIntegrator(sim.integrator);
        meterNormalStress1.setBox(sim.box);
        meterNormalStress1.setDims(new int[]{0,1});
        AccumulatorAverageCollapsing avgNormalStress1 = new AccumulatorAverageCollapsing();
        DataPumpListener normalStress1Pump = new DataPumpListener(meterNormalStress1, avgNormalStress1, 10);
        sim.integrator.getEventManager().addListener(normalStress1Pump);
        DisplayTextBoxesCAE normalStress1Display = new DisplayTextBoxesCAE();
        normalStress1Display.setLabel("first normal stress coefficient");
        normalStress1Display.setAccumulator(avgNormalStress1);
        add(normalStress1Display);
        getController().getDataStreamPumps().add(normalStress1Pump);
        sliderA.setPostAction(new IAction() {
            public void actionPerformed() {
                meterNormalStress1.setDoDouble(sliderA.getValue() < 0);
            }
        });


        final MeterNormalStress meterNormalStress2 = new MeterNormalStress(sim.getSpace());
        meterNormalStress2.setIntegrator(sim.integrator);
        meterNormalStress2.setBox(sim.box);
        meterNormalStress2.setDims(new int[]{1,2});
        AccumulatorAverageCollapsing avgNormalStress2 = new AccumulatorAverageCollapsing();
        DataPumpListener normalStress2Pump = new DataPumpListener(meterNormalStress2, avgNormalStress2, 10);
        sim.integrator.getEventManager().addListener(normalStress2Pump);
        DisplayTextBoxesCAE normalStress2Display = new DisplayTextBoxesCAE();
        normalStress2Display.setLabel("second normal stress coefficient");
        normalStress2Display.setAccumulator(avgNormalStress2);
        add(normalStress2Display);
        getController().getDataStreamPumps().add(normalStress2Pump);
        
        DeviceDelaySlider delaySlider = new DeviceDelaySlider(sim.getController(), sim.activityIntegrate);
        getPanel().controlPanel.add(delaySlider.graphic(), SimulationPanel.getVertGBC());
    }
    
    public static void main(String[] args) {
        SimulationRheology sim = new SimulationRheology(Space3D.getInstance());
        RheologyGraphic graphic = new RheologyGraphic(sim);
        graphic.makeAndDisplayFrame();
    }
}
