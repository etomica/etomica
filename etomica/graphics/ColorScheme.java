package etomica.graphics;
import java.awt.Color;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.atom.IAtom;
import etomica.modifier.ModifierGeneral;
import etomica.nbr.list.PotentialMasterList;

/**
 * Class that defines the algorithm used to determine atoms colors when drawn to DisplayBox.
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

    public abstract Color getAtomColor(IAtom a);
    
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
        public Color getAtomColor(IAtom a) {return defaultColor;}
    }
    
    //for color scheme demo
    public static void main(String args[]) {
      final String APP_NAME = "Color Scheme";

      final etomica.simulation.prototypes.HSMD3D sim = new etomica.simulation.prototypes.HSMD3D();
      final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME);

      Action repaintAction = simGraphic.getPaintAction(sim.box);

      DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
      nSelector.setResetAction(new SimulationRestart(sim));
      nSelector.setSpecies(sim.species);
      nSelector.setBox(sim.box);
      nSelector.setPostAction(repaintAction);
      simGraphic.add(nSelector);


      simGraphic.getController().getReinitButton().setPostAction(repaintAction);

      final ColorSchemeByType ct = new ColorSchemeByType();
      final ColorSchemeTemperature ctemp = new ColorSchemeTemperature(0,5);
      final ColorSchemeColliders ccld = new ColorSchemeColliders(sim.integrator);
      final ColorSchemeNeighbor nghb = new ColorSchemeNeighbor(sim, (PotentialMasterList)sim.potentialMaster, sim.box);
      nghb.setAtom(sim.box.getLeafList().getAtom(0));
      final ColorSchemeRandom rand = new ColorSchemeRandom(sim.box, sim.getRandom());
      final ColorSchemeCell cell = new ColorSchemeCell((PotentialMasterList)sim.potentialMaster,sim.getRandom(),sim.box);
      cell.setLattice(((PotentialMasterList)sim.potentialMaster).getNbrCellManager(sim.box).getLattice());
      
      Action act = new Action() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          ct.setColor(sim.species.getMoleculeType(), 
              new java.awt.Color(sim.getRandom().nextInt(256),
                  sim.getRandom().nextInt(256),
                  sim.getRandom().nextInt(256))
              );
          dp.setColorScheme(ct);
        }
      };
      Action act2 = new Action() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(ctemp);
        }
      };
      Action act3 = new Action() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(ccld);
        }
      };
      Action act4 = new Action() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(nghb);
        }
      };
      Action act5 = new Action() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(rand);
        }
      };
      Action act6 = new Action() {
        public void actionPerformed() {
          ct.setColor(sim.species.getMoleculeType(), java.awt.Color.red);
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(ct);
        }
      };
      Action act7 = new Action() {
        public void actionPerformed() {
          DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
          dp.setColorScheme(cell);
        }
      };

      DeviceButton colorer = new DeviceButton(sim.getController(),act);
      colorer.setLabel("Global Random");
      colorer.setPostAction(repaintAction);
      DeviceButton tempcolorer = new DeviceButton(sim.getController(),act2);
      tempcolorer.setLabel("By Temperature (0-5)");
      DeviceButton colliders = new DeviceButton(sim.getController(),act3);
      colliders.setLabel("Colliders");
      DeviceButton neighbors = new DeviceButton(sim.getController(),act4);
      neighbors.setLabel("Neighbors");
      neighbors.setPostAction(repaintAction);
      DeviceButton randomcol = new DeviceButton(sim.getController(),act5);
      randomcol.setLabel("Unique Random");
      randomcol.setPostAction(repaintAction);
      DeviceButton def = new DeviceButton(sim.getController(),act6);
      def.setLabel("Default Red");
      def.setPostAction(repaintAction);
      DeviceButton cellbtn = new DeviceButton(sim.getController(),act7);
      cellbtn.setLabel("Cell");
      
      DeviceSlider slabslide = new DeviceSlider(sim.getController());
      slabslide.setMinimum(0);
      slabslide.setMaximum(100);
      final DisplayBoxCanvasG3DSys dpg3d = (DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas;
      ModifierGeneral mg = new ModifierGeneral(dpg3d,"slab");
      slabslide.setModifier(mg);
      slabslide.setLabel("Slab");
      
      DeviceSlider depthslide = new DeviceSlider(sim.getController());
      depthslide.setMinimum(0);
      depthslide.setMaximum(100);
      ModifierGeneral mg2 = new ModifierGeneral(dpg3d,"depth");
      depthslide.setModifier(mg2);
      depthslide.setLabel("Depth");
      
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
//      };
      
      simGraphic.makeAndDisplayFrame(APP_NAME);
      
      //ColorSchemeTemperature colorScheme = new ColorSchemeTemperature(0.0f,5.0f);

      ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
      colorScheme.setColor(sim.species.getMoleculeType(), java.awt.Color.red);

      DisplayBox dp = (DisplayBox)simGraphic.displayList().getFirst();
      dp.setColorScheme(colorScheme);
    }
}
