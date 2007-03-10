package etomica.graphics;
import java.awt.Color;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.atom.AtomLeaf;
import etomica.integrator.IntegratorHard;
import etomica.modifier.ModifierGeneral;
import etomica.nbr.list.PotentialMasterList;
import etomica.util.Default;

/**
 * Class that defines the algorithm used to determine atoms colors when drawn to DisplayPhase.
 * The atomColor method is called just before the atom is drawn to set the graphics color.
 *
 * @author David Kofke
 */
 
public abstract class ColorScheme implements java.io.Serializable {

    protected Color defaultColor;
    
    public ColorScheme() {
        this(DEFAULT_ATOM_COLOR);
    }
    public ColorScheme(Color color) {
        defaultColor = color;
    }

    public abstract Color getAtomColor(AtomLeaf a);
    
    public final void setDefaultColor(Color c) {defaultColor = c;}
    public final Color getDefaultColor() {return defaultColor;}

    public static Color DEFAULT_ATOM_COLOR = Color.red;
    
    /**
     * Colors all atoms with baseColor.
     */
    public static class Simple extends ColorScheme {
        private static final long serialVersionUID = 1L;
        public Simple() {super();}
        public Simple(java.awt.Color color) {super(color);}
        public Color getAtomColor(AtomLeaf a) {return defaultColor;}
    }//end of Simple
    
    //for color scheme demo
    public static void main(String args[]) {
      Default defaults = new Default();
      defaults.doSleep = false;
      defaults.ignoreOverlap = true;
      final etomica.simulation.prototypes.HSMD3D sim = new etomica.simulation.prototypes.HSMD3D();
      final SimulationGraphic simGraphic = new SimulationGraphic(sim);
      DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
      nSelector.setResetAction(new SimulationRestart(sim));
      nSelector.setSpeciesAgent(sim.phase.getAgent(sim.species));
      simGraphic.add(nSelector);
      
      final ColorSchemeByType ct = new ColorSchemeByType();
      final ColorSchemeTemperature ctemp = new ColorSchemeTemperature(0,5);
      final ColorSchemeColliders ccld =
        new ColorSchemeColliders((IntegratorHard)sim.getIntegratorList().getFirst());
      final ColorSchemeNeighbor nghb = new ColorSchemeNeighbor(sim,sim.phase);
      nghb.setAtom(sim.phase.firstAtom());
      final ColorSchemeRandom rand = new ColorSchemeRandom(sim.phase, sim.getRandom());
      final ColorSchemeCell cell = new ColorSchemeCell(sim,sim.phase);
      cell.setLattice(((PotentialMasterList)sim.getPotentialMaster()).getNbrCellManager(sim.phase).getLattice());
      
      Action act = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          ct.setColor(dp.getPhase().firstAtom().getType(), 
              new java.awt.Color(sim.getRandom().nextInt(256),
                  sim.getRandom().nextInt(256),
                  sim.getRandom().nextInt(256))
              );
          dp.setColorScheme(ct);
        }

        public String getLabel() {
          return "Global Random";
        }        
      };
      Action act2 = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(ctemp);
        }
        public String getLabel() {
          return "By Temperature (0-5)";
        }
      };
      Action act3 = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(ccld);
        }
        public String getLabel() {
          return "Colliders";
        }
      };
      Action act4 = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(nghb);
        }
        public String getLabel() {
          return "Neighbors";
        }
      };
      Action act5 = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(rand);
        }
        public String getLabel() {
          return "Unique Random";
        }
      };
      Action act6 = new Action() {
        public void actionPerformed() {
          ct.setColor(sim.phase.firstAtom().getType(), java.awt.Color.red);
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(ct);
        }
        public String getLabel() {
          return "Default Red";
        }
      };
      Action act7 = new Action() {
        public void actionPerformed() {
          DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
          dp.setColorScheme(cell);
        }
        public String getLabel() {
          return "Cell";
        }
      };

      DeviceButton colorer = new DeviceButton(sim.getController(),act);
      DeviceButton tempcolorer = new DeviceButton(sim.getController(),act2);
      DeviceButton colliders = new DeviceButton(sim.getController(),act3);
      DeviceButton neighbors = new DeviceButton(sim.getController(),act4);
      DeviceButton randomcol = new DeviceButton(sim.getController(),act5);
      DeviceButton def = new DeviceButton(sim.getController(),act6);
      DeviceButton cellbtn = new DeviceButton(sim.getController(),act7);
      
      DeviceSlider slabslide = new DeviceSlider(sim.getController());
      slabslide.setMinimum(0);
      slabslide.setMaximum(100);
      final DisplayPhaseCanvasG3DSys dpg3d = (DisplayPhaseCanvasG3DSys)simGraphic.getDisplayPhase(sim.phase).canvas;
      ModifierGeneral mg = new ModifierGeneral(dpg3d,"slab");
      slabslide.setModifier(mg);
      slabslide.setName("Slab");
      
      DeviceSlider depthslide = new DeviceSlider(sim.getController());
      depthslide.setMinimum(0);
      depthslide.setMaximum(100);
      ModifierGeneral mg2 = new ModifierGeneral(dpg3d,"depth");
      depthslide.setModifier(mg2);
      depthslide.setName("Depth");
      
      simGraphic.add(def);
      simGraphic.add(colorer);
      simGraphic.add(tempcolorer);
      simGraphic.add(colliders);
      simGraphic.add(neighbors);
      simGraphic.add(randomcol);
      simGraphic.add(cellbtn);
      simGraphic.add(slabslide);
      simGraphic.add(depthslide);
      
      
//      d.get
//      Action setslab = new Action() {
//
//        public void actionPerformed() {
//          //d.setS
//        }
//
//        public String getLabel() {
//          // TODO Auto-generated method stub
//          return null;
//        }
//        
//      };
      
      simGraphic.makeAndDisplayFrame();
      
      //ColorSchemeTemperature colorScheme = new ColorSchemeTemperature(0.0f,5.0f);

      ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simGraphic.displayList().getFirst()).getColorScheme());
      colorScheme.setColor(sim.species.getMoleculeType(), java.awt.Color.red);

      DisplayPhase dp = (DisplayPhase)simGraphic.displayList().getFirst();
      dp.setColorScheme(colorScheme);
      simGraphic.panel().setBackground(java.awt.Color.yellow);        
    }
}
