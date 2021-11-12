package etomica.modules;

import etomica.modules.adsorption.AdsorptionGraphic;
import etomica.modules.catalysis.CatalysisGraphic;
import etomica.modules.chainequilibrium.ChainEquilibriumGraphic;
import etomica.modules.colloid.ColloidGraphic;
import etomica.modules.dcvgcmd.DCVGCMDGraphic;
import etomica.modules.droplet.DropletAtomicGraphic;
import etomica.modules.droplet.DropletGraphic;
import etomica.modules.ensembles.LJMCGraphic;
import etomica.modules.glass.GlassGraphic;
import etomica.modules.insertion.InsertionGraphic;
import etomica.modules.interfacial.InterfacialSWGraphic;
import etomica.modules.ljmd.LjmdGraphic;
import etomica.modules.materialfracture.MaterialFractureGraphic;
import etomica.modules.mu.MuGraphic;
import etomica.modules.multiharmonic.MultiharmonicGraphicMC;
import etomica.modules.osmosis.Osmosis;
import etomica.modules.pistoncylinder.PistonCylinderGraphic;
import etomica.modules.reactionequilibrium.ReactionEquilibriumGraphic;
import etomica.modules.rosmosis.ReverseOsmosisGraphic;
import etomica.modules.rosmosis.ReverseOsmosisWaterGraphic;
import etomica.modules.sam.SamGraphic;
import etomica.modules.statistics.StatisticsMCGraphic;
import etomica.modules.swmd.SwmdGraphic;
import etomica.modules.vle.VLE;

public final class EtomicaModuleInfo {
    public final String moduleName;
    public final String description;
    public final Class<?> moduleClass;
    public final String args;

    public EtomicaModuleInfo(String moduleName, String description, Class<?> moduleClass, String args) {
        this.moduleName = moduleName;
        this.description = description;
        this.moduleClass = moduleClass;
        this.args = args;
    }

    public static EtomicaModuleInfo mod(String moduleName, String description, Class<?> moduleClass) {
        return new EtomicaModuleInfo(moduleName, description, moduleClass, "");
    }

    public static EtomicaModuleInfo mod(String moduleName, String description, Class<?> moduleClass, String args) {
        return new EtomicaModuleInfo(moduleName, description, moduleClass, args);
    }

    public String toString() {
        return this.moduleName;
    }

    public static final EtomicaModuleInfo[] ETOMICA_MODULES = {
            mod("Adsorption", "No description yet", AdsorptionGraphic.class),
            mod("Catalysis", "No description yet", CatalysisGraphic.class),
            mod("Chain Equilibrium", "No description yet", ChainEquilibriumGraphic.class),
            mod("Colloid", "No description yet", ColloidGraphic.class),
            mod("Crystal Viewer", "Permits viewing of static 3-dimensional lattices and the planes they define. Various elementary lattices can be selected, and any plane or surface can be viewed by specifying it via its Miller indices. Image can be rotated to permit viewing from any angle.", etomica.modules.crystalviewer.CrystalViewer.class),
            mod("Discontinuous MD 2D", "No description yet", SwmdGraphic.class, "2"),
            mod("Discontinuous MD 3D", "No description yet", SwmdGraphic.class, "3"),
            mod("DCVGCMD", "No description yet", DCVGCMDGraphic.class),
            mod("Droplet", "No description yet", DropletGraphic.class),
            mod("Droplet (atomic)", "No description yet", DropletAtomicGraphic.class),
            mod("Ensembles", "No description yet", LJMCGraphic.class),
            mod("Glass", "No description yet", GlassGraphic.class),
            mod("Insertion", "No description yet", InsertionGraphic.class),
            mod("Interfacial", "No description yet", InterfacialSWGraphic.class),
            mod("LJMD", "No description yet", LjmdGraphic.class),
            mod("Material Fracture", "No description yet", MaterialFractureGraphic.class),
            mod("Chemical Potential", "No description yet", MuGraphic.class),
            mod("Multiharmonic", "No description yet", MultiharmonicGraphicMC.class),
            mod("Nucleation", "No description yet", etomica.modules.nucleation.NucleationGraphic.class),
            mod("Osmosis", "No description yet", Osmosis.class),
            mod("Piston-Cylinder 2D", "The piston-cylinder apparatus is a standard tool used to illustrate thermodynamic concepts involving heat, work, and internal energy. This module simulates a collection of atoms in a chamber with a movable wall under external pressure.", PistonCylinderGraphic.class, "-dim 2"),
            mod("Piston-Cylinder 3D", "The piston-cylinder apparatus is a standard tool used to illustrate thermodynamic concepts involving heat, work, and internal energy. This module simulates a collection of atoms in a 3D chamber with a movable wall under external pressure.", PistonCylinderGraphic.class, "-dim 3"),
            mod("Piston-Cylinder-SWMD", "The piston-cylinder apparatus is a standard tool used to illustrate thermodynamic concepts involving heat, work, and internal energy. This module simulates a collection of atoms in a chamber with a movable wall under external pressure.  This version of the simulation has extra features used for the SWMD module", PistonCylinderGraphic.class, "-dim 2 -configButton -densityInput -nMoleculeSlider -rdf -fastButton"),
            mod("Reaction Equilibrium", "Simple reaction equlibrium involving two atomic species and the three dimeric molecules they can form. Atoms move about via 2-D molecular dynamics, and can \"react\" to form dimers. Dynamic equilibria is demonstrated through constant recombining of atoms, and equilibria can be quantified and analyzed with thermodynamic reaction equilibria models.", ReactionEquilibriumGraphic.class),
            mod("Rheology", "No description yet", etomica.modules.rheology.RheologyGraphic.class),
            mod("Reverse Osmosis", "No description yet", ReverseOsmosisGraphic.class),
            mod("Reverse Osmosis (water)", "No description yet", ReverseOsmosisWaterGraphic.class),
            mod("Self-Assembled Monolayer", "No description yet", SamGraphic.class),
            mod("Statistics", "No description yet", StatisticsMCGraphic.class, "-moduleNum 1"),
            mod("Statistics Part 2", "No description yet", StatisticsMCGraphic.class, "-moduleNum 2"),
            mod("VLE", "No description yet", VLE.class),
            mod("B2 fit (for VLE)", "No description yet", etomica.modules.vle.B2Fit.class),
    };
}

