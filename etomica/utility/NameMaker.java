/*
 * History
 * Created on Nov 15, 2004 by kofke
 */
package etomica.utility;

import java.util.HashMap;

import etomica.Integrator;
import etomica.Phase;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class NameMaker {

    /**
     * 
     */
    private NameMaker() { }

    public static String makeName(Class className) {
        ClassCount thisCount = (ClassCount)classList.get(className);
        if (thisCount == null) {
            thisCount = new ClassCount();
            classList.put(className,thisCount);
        }
        return className.getName().substring(8)+thisCount.count++;
    }
    
    private static HashMap classList = new HashMap();
    
    private static class ClassCount {
        int count;
    }
    
    public static void main(String[] argv) {
        System.out.println(NameMaker.makeName(Phase.class));
        System.out.println(NameMaker.makeName(Integrator.class));
        System.out.println(NameMaker.makeName(Phase.class));
    }
}
