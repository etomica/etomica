package etomica.units;

public abstract class BaseUnitPseudo3D extends BaseUnit implements BaseUnit.D3 {
    
    public static final double FALSE_DEPTH = 5.0;  //Angstroms
    
    public BaseUnitPseudo3D(BaseUnit unit3D, double mult) {
    	super(unit3D.toSim(1.0)*mult, unit3D.toString(), unit3D.symbol(), unit3D.dimension(), unit3D.prefixAllowed());
//        name = unit3D.toString();
//        symbol = unit3D.symbol();
//        prefixAllowed = unit3D.prefixAllowed();
    }
      
    public static class Pressure extends BaseUnitPseudo3D {
        public Pressure(BaseUnit.Pressure baseUnit) {
            super(baseUnit, FALSE_DEPTH);
        }
    }//end of Pressure
    
    public static class Volume extends BaseUnitPseudo3D {
        public Volume(BaseUnit.Volume baseUnit) {
            super(baseUnit, 1.0/FALSE_DEPTH);
        }
    }//end of Volume
    
}//end of BaseUnitPseudo3D