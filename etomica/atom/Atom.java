package etomica.atom;

import etomica.api.IAtom;

 /**
  * Object corresponding to one physical atom or group of atoms. Each atom holds
  * the following publicly accessible fields:
  * <ul>
  * <li>an AtomType instance (fieldname: type) that holds information this atom
  * has in common with other atoms made by the same factory
  * </ul>
  * <p>
  * @author David Kofke, Andrew Schultz, and C. Daniel Barnes
  * 
  */
public abstract class Atom implements IAtom, java.io.Serializable {

    public Atom() {
    }
    
    public abstract String signature();
    
    private static final long serialVersionUID = 1L;
}
