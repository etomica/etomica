package etomica.modules.materialfracture;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.util.HistoryScrolling;

public class MaterialFractureGraphic extends SimulationGraphic {

    
//********************start of Composite Manager Panel*************************************************
    public static JPanel compositePanel = new JPanel(new java.awt.GridBagLayout()); //default shape of panel
    public static GridBagLayout gbLayout = new GridBagLayout();
    public static GridBagConstraints gbConst = new GridBagConstraints();    
    public static JPanel temperaturePanel = new JPanel(new java.awt.GridLayout(1, 2)); //default shape of panel
    public static JPanel pdPanel = new JPanel(new java.awt.GridLayout()); //default shape of panel
    public static JTabbedPane sliderTabbedPane = new JTabbedPane();
//********************end of Composite Manager Panel***************************************************

/**
/* constructor
 */
  public MaterialFractureGraphic(MaterialFracture sim) {
    super(sim, SimulationGraphic.TABBED_PANE, "Material Fracture", 1, sim.getSpace());
    gbLayout = (GridBagLayout)compositePanel.getLayout();
//    SimulationResetForSSF srssf = new SimulationResetForSSF(); // this is different from SimulationRestart
//        srssf.setSpecies(sim.species);
//        srssf.setDisplay(displayPhase0);
//        srssf.setPhase(phase0);
//    etomica.graphics.DeviceTrioControllerButton dscb = new etomica.graphics.DeviceTrioControllerButton();//start, reset, reset average
//        dscb.setShape("VERTICAL");
//        dscb.resetButton().setAction(new ActionGraphic(srssf));
             
    final StrainColorScheme strainColor = new StrainColorScheme(sim);    
    getDisplayBox(sim.box).setColorScheme(strainColor);
    strainColor.setNumber(37);

//********************************* Simulation Configuration Data*********************************************************
//    JPanel textPanel = new JPanel();
//    textPanel.setBorder(new javax.swing.border.TitledBorder(new javax.swing.border.EtchedBorder(
//                                javax.swing.border.EtchedBorder.RAISED, java.awt.Color.gray, java.awt.Color.black)));
//    final JTextField textField1[] = new JTextField[4];
//    final JTextField textField2[] = new JTextField[4];
//        textPanel.setLayout(new java.awt.GridLayout(1, (textField1.length+textField2.length) ));
//        for(int i=0; i<textField1.length;i++){
//                textField1[i] = new JTextField();
//                textField1[i].setHorizontalAlignment(SwingConstants.CENTER);textField1[i].setEditable(false); 
//                textField1[i].setFont(new java.awt.Font(null,java.awt.Font.BOLD,12));
//                textField2[i] = new JTextField();
//                textField2[i].setHorizontalAlignment(SwingConstants.CENTER);textField2[i].setEditable(true);                
//                textPanel.add(textField1[i]);textPanel.add(textField2[i]);
//            }
//     textField1[0].setText(" Number of Atom ");   
//     textField1[1].setText(" Dimension X ");   
//     textField1[2].setText(" Dimension Y ");
//     textField1[3].setText(" Gage Length ");
//     srssf.setTextFields(textField2);     
// *********************************Simulation Configuraiton Data*********************************************************
//************************ start of Display Plot related ********************************************            

    final MeterStrain meterStrain = new MeterStrain();
    meterStrain.setBox(sim.box);
    meterStrain.setAtomNumber(37);
    DataFork strainFork = new DataFork();
    DataPump strainPump = new DataPump(meterStrain, strainFork);
    getController().getDataStreamPumps().add(strainPump);
    sim.integrator.addIntervalAction(strainPump);
    AccumulatorAverageCollapsing strainAverage = new AccumulatorAverageCollapsing();
    strainFork.addDataSink(strainAverage);
    
    final MeterStress meterStress = new MeterStress(sim.pc);
    meterStress.setBox(sim.box);
    DataFork stressFork = new DataFork();
    DataPump stressPump = new DataPump(meterStress, stressFork);
    getController().getDataStreamPumps().add(stressPump);
    sim.integrator.addIntervalAction(stressPump);
    AccumulatorHistory stressHistory = new AccumulatorHistory(new HistoryScrolling(400));
    stressHistory.setTimeDataSource(meterStrain);
    stressFork.addDataSink(stressHistory);
    AccumulatorAverageCollapsing stressAverage = new AccumulatorAverageCollapsing();
    stressAverage.setPushInterval(10);
    stressFork.addDataSink(stressAverage);
    
    DisplayPlot stressHistoryPlot = new DisplayPlot();
    stressHistory.setDataSink(stressHistoryPlot.getDataSet().makeDataSink());
    stressHistoryPlot.setLabel("Stress");
    stressHistoryPlot.setDoLegend(false);
    
    add(stressHistoryPlot);
    
    final MeterElongation meterElongation = new MeterElongation();
    meterElongation.setBox(sim.box);
    meterElongation.setAtomNumber(37);
    DataFork elongationFork = new DataFork();
    DataPump elongationPump = new DataPump(meterElongation, elongationFork);
    getController().getDataStreamPumps().add(elongationPump);
    sim.integrator.addIntervalAction(elongationPump);
    
    final MeterLoad meterLoad = new MeterLoad(sim.pc);
    DataFork loadFork = new DataFork();
    DataPump loadPump = new DataPump(meterLoad, loadFork);
    getController().getDataStreamPumps().add(loadPump);
    sim.integrator.addIntervalAction(loadPump);
    AccumulatorHistory loadHistory = new AccumulatorHistory(new HistoryScrolling(400));
    loadHistory.setTimeDataSource(meterElongation);
    loadFork.addDataSink(loadHistory);
    AccumulatorAverageCollapsing loadAverage = new AccumulatorAverageCollapsing();
    loadAverage.setPushInterval(10);
    loadFork.addDataSink(loadAverage);
    
    DisplayPlot loadHistoryPlot = new DisplayPlot();
    loadHistory.setDataSink(loadHistoryPlot.getDataSet().makeDataSink());
    loadHistoryPlot.setLabel("Load");
    stressHistoryPlot.setDoLegend(false);
    
    add(loadHistoryPlot);
    
    DisplayTextBoxesCAE stressDisplay = new DisplayTextBoxesCAE();
    stressDisplay.setAccumulator(stressAverage);
    add(stressDisplay);
    
    DisplayTextBoxesCAE strainDisplay = new DisplayTextBoxesCAE();
    strainDisplay.setAccumulator(strainAverage);
    add(strainDisplay);

     
//Potential Sliders
//    etomica.ModulatorAbstract[] modul1 = new etomica.ModulatorAbstract[4];    
//    DiameterModulator diaModulator1 = this.new DiameterModulator(P2LennardJones1,speciesSpheres0 );
//        diaModulator1.setDisplay(displayPhase0);
//         modul1[0] = diaModulator1;
//    EpsilonModulator epsilonModulator1 = this.new EpsilonModulator(P2LennardJones1);
//         modul1[1] = epsilonModulator1;
//    SpringConstantModulator scModulator1 = this.new SpringConstantModulator(p1, meterSS);
//        scModulator1.setDisplay(displayPhase0);
//         modul1[2] = scModulator1;
//    PotentialTruncationModulator ptModulator1 = this.new PotentialTruncationModulator(pt);
//        ptModulator1.setDisplay(displayPhase0);
//         modul1[3] = ptModulator1;
//    DevicePotentialSlider potSlider = this.new DevicePotentialSlider(this);  
//        potSlider.setModulator(modul1);
        
//Integrator Sliders        
//    etomica.ModulatorAbstract[] modul2 = new etomica.ModulatorAbstract[2];    
//    TimeStepModulator timeModulator1 = this.new TimeStepModulator(integratorVelocityVerlet1);
//    AndersonNuModulator annuModulator1 = this.new AndersonNuModulator(integratorVelocityVerlet1);
//         modul2[0] = timeModulator1;
//         modul2[1] = annuModulator1;
//     DeviceIntegratorSlider intSlider = this.new DeviceIntegratorSlider(this);  
//         intSlider.setModulator(modul2);
//Auxiliary Sliders and Buttons
//    etomica.ModulatorAbstract[] modul3 = new etomica.ModulatorAbstract[2];    
//    final GageLengthModulator glModulator1 = this.new GageLengthModulator(strainColor, meterSS);
//         glModulator1.setDisplay(displayPhase0);
//         glModulator1.setTextField(textField2[3]);
//         modul3[0] = glModulator1;
//    DisplayModulator dModulator1 = this.new DisplayModulator(displayPhase0);
//         modul3[1] = dModulator1;
//    final DeviceAuxiliary auxDevices = this.new DeviceAuxiliary(this);     
//         auxDevices.setModulator(modul3);
         
//************************start of Atom Configuration related************************************
//     etomica.ModulatorBoolean modulator1 = new etomica.ModulatorBoolean() {
//           public void setBoolean(boolean b) {
//               if(b) {
//                      ((etomica.ConfigurationSequential)phase0.getConfiguration()).setSquareConfig(true);
//                      speciesSpheres0.setNMolecules(200);
//                      glModulator1.setHexagonal(false);strainColor.setHexagonal(false); meterSS.setHexagonal(false);
//                      strainColor.setNumber(40); glModulator1.setValue(50);auxDevices.getSlider().setValue(50);
//                      textField2[0].setText(String.valueOf(speciesSpheres0.getNMolecules()));   
//                      display.repaint();
//                } else {
//                      ((etomica.ConfigurationSequential)phase0.getConfiguration()).setSquareConfig(false); 
//                      speciesSpheres0.setNMolecules(198);
//                      glModulator1.setHexagonal(true);strainColor.setHexagonal(true); meterSS.setHexagonal(true);
//                      strainColor.setNumber(45); glModulator1.setValue(50);auxDevices.getSlider().setValue(50);
//                      textField2[0].setText(String.valueOf(speciesSpheres0.getNMolecules()));   
//                      display.repaint();
//                }
//           }
//                public boolean getBoolean() {return false;}
//            };
//     etomica.graphics.DeviceToggleRadioButtons configButton = new etomica.graphics.DeviceToggleRadioButtons(this, modulator1,"", "Square", "Triangular");
//               configButton.setName("Config Shape");  
//               auxDevices.setGraphicComponents(configButton);
//*************************end of Atom Configuration related*****************************   
               
//************************start of Display Tension Plot related************************************
//     final SSF2DPlot ssf2dPlot = new SSF2DPlot(displayPhase0);
//     final SSF2DPlotPotentialCutOff ssf2dPlotPT = new SSF2DPlotPotentialCutOff(displayPhase0, phase0);
//     final etomica.ModulatorBoolean modulator2;
//     final etomica.ModulatorBoolean modulator3;
//     
//        ssf2dPlot.setPotentialTension(p1);
//        modulator2 = new etomica.ModulatorBoolean() {
//           public void setBoolean(boolean b) {
//               if(b) {
//                     display.addDrawable(ssf2dPlot); //ssf2dPlotPT.setShowLJPotential(false); 
//                     if(display.getDrawables().size()!=1){ssf2dPlotPT.setShowLJPotential(false);}
//                     display.repaint();
//                } else {
//                     display.removeDrawable(ssf2dPlot); ssf2dPlotPT.setShowLJPotential(true); 
//                     display.repaint();
//                }
//           }
//                public boolean getBoolean() {return false;}
//            };
//     final etomica.graphics.DeviceToggleRadioButtons plotButton = new etomica.graphics.DeviceToggleRadioButtons(this, modulator2,"", "On", "Off");
//               plotButton.setName("Potential Field");  
//               auxDevices.setGraphicComponents(plotButton);
//*************************end of Display Tension plot related*****************************     

//************************start of Display Tension Plot related************************************
//     ssf2dPlotPT.setPotentialCutoff(pt); //passing PotentialTruncationSimple instance
//     ssf2dPlotPT.setPotential(P2LennardJones1);//passin LJ Potential
//        modulator3= new etomica.ModulatorBoolean() {
//           public void setBoolean(boolean b) {
//               if(b) { 
//                     display.addDrawable(ssf2dPlotPT); //System.out.println(" here "+ssf2dPlotPT.getShowPotential());//ssf2dPlotPT.setShowLJPotential(true); 
//                     if(display.getDrawables().size()!=1){ssf2dPlotPT.setShowLJPotential(false);}
//                     display.repaint();
//               } else {
//                     display.removeDrawable(ssf2dPlotPT);  ssf2dPlotPT.setShowLJPotential(true);
//                     display.repaint();  
//                }
//           }
//                public boolean getBoolean() {return false;}
//            };
//     etomica.graphics.DeviceToggleRadioButtons plotPtButton = new etomica.graphics.DeviceToggleRadioButtons(this, modulator3,"", "On", "Off");
//               plotPtButton.setName("Potential Cutoff");  
//               auxDevices.setGraphicComponents(plotPtButton);
//*************************end of Display Tension plot related*****************************     
               
                  
//********************start of temperature related**************************************************
//    etomica.graphics.DeviceThermoSelector dtSelector = new etomica.graphics.DeviceThermoSelector();
//            dtSelector.setTemperatures(new double [] {50, 100, 200, 300, 400, 500, 600 });
//            dtSelector.setUnit(new etomica.units.Unit(Kelvin.UNIT) ); 
//            dtSelector.setSelected(4);
//            dtSelector.getLabel().setText("Set Value(K)");
//
//            temperaturePanel.setBorder(new javax.swing.border.TitledBorder(
//                            new javax.swing.border.EtchedBorder(
//                                javax.swing.border.EtchedBorder.RAISED, java.awt.Color.red, java.awt.Color.blue) 
//                                ,"Temperature"
//                                ,javax.swing.border.TitledBorder.LEFT
//                                ,javax.swing.border.TitledBorder.TOP
//                                ,new java.awt.Font(null,java.awt.Font.BOLD,15)
//                                ,java.awt.Color.black));
//            temperaturePanel.add(dtSelector.graphic());
//*********************** end of temperature related **********************************************
                        
//*********************** start of Main Composite Board *******************************************
//            compositePanel.setLayout(gbLayout);
//            compositePanel.setBorder(new javax.swing.border.TitledBorder(
//                            new javax.swing.border.EtchedBorder(
//                                javax.swing.border.EtchedBorder.RAISED, java.awt.Color.gray, java.awt.Color.black) 
//                                ,"Simulation Manager"
//                                ,javax.swing.border.TitledBorder.CENTER 
//                                ,javax.swing.border.TitledBorder.TOP
//                                ,new java.awt.Font(null,java.awt.Font.BOLD,18)
//                                ,java.awt.Color.black));
//                                
//            sliderTabbedPane.add("Potential",potSlider.graphic());                               
//            sliderTabbedPane.add("Integrator",intSlider.graphic());                               
//            sliderTabbedPane.add("Aux", auxDevices.graphic());                               
//            gbConst.fill = GridBagConstraints.BOTH;
//            gbConst.gridx = 1;gbConst.gridy = 0; gbConst.gridwidth = 1;gbConst.gridheight = 1;
//                gbLayout.setConstraints(dscb.graphic(), gbConst);
//                  compositePanel.add(dscb.graphic()); 
//            gbConst.fill = GridBagConstraints.HORIZONTAL;
//            gbConst.gridx = 1;gbConst.gridy = 1; gbConst.gridwidth = 1;gbConst.gridheight = 1;
//                gbLayout.setConstraints(pdPanel, gbConst);
//                  compositePanel.add(pdPanel);
//            gbConst.fill = GridBagConstraints.HORIZONTAL;
//            gbConst.gridx = 1;gbConst.gridy = 2; gbConst.gridwidth = 1;gbConst.gridheight = 1;
//                gbLayout.setConstraints(temperaturePanel, gbConst);
//                  compositePanel.add(temperaturePanel);
//            gbConst.fill = GridBagConstraints.HORIZONTAL;
//            gbConst.gridx = 1;gbConst.gridy = 3; gbConst.gridwidth = 1;gbConst.gridheight = 1;
//                gbLayout.setConstraints(sliderTabbedPane, gbConst);
//                  compositePanel.add(sliderTabbedPane);             
//************************ end of Main Composite Board ********************************************

//    this.mediator().addMediatorPair(new MediatorGraphic.DisplayNull.NoAction(this.mediator()));
//    this.mediator().addMediatorPair(new MediatorGraphic.DeviceNull.NoAction(this.mediator()));    
//
//    this.mediator().go();
//                       
//    ((SimulationGraphic)this).panel().setLayout(new java.awt.BorderLayout());
//     javax.swing.JTabbedPane displayPanel = new javax.swing.JTabbedPane();
//     displayPanel.add("Simulation Configuration", displayPhase0.graphic());
//     displayPanel.add("Plot-LE",plot.getPlot());
//     displayPanel.add("Plot-SS",plot1.getPlot());
//     displayPanel.add("Simulation Data", displayTable.graphic());
//     
//     textField2[0].setText(String.valueOf(phase0.getAgent(speciesSpheres0).getNMolecules()));   
//     textField2[1].setText(String.valueOf(phase0.dimensions().x(0)));
//     textField2[2].setText(String.valueOf(phase0.dimensions().x(1)));
//     textField2[3].setText(String.valueOf(glModulator1.getValue()));     
//     
//     ((SimulationGraphic)this).panel().add(displayPanel, BorderLayout.CENTER);
//     ((SimulationGraphic)this).panel().add(compositePanel, BorderLayout.WEST); 
//     ((SimulationGraphic)this).panel().add(textPanel, BorderLayout.SOUTH);
         
  } //end of constructor
  
/**
/* main method
 */
  public static void main(String[] args) {
     MaterialFracture sim = new MaterialFracture();
     MaterialFractureGraphic simGraphic = new MaterialFractureGraphic(sim);

     simGraphic.makeAndDisplayFrame();

  }//end of main

  public static class Applet extends javax.swing.JApplet{
        public void init(){ 
            getContentPane().add(new MaterialFractureGraphic(new MaterialFracture()).getPanel());
        }
  }
    
}//end of class




