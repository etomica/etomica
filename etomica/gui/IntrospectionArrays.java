package etomica.gui;

public class IntrospectionArrays {
    static final Class[] spaceClasses = introspect("Space", true);
    static final Class[] speciesClasses = introspect("Species", true);
    static final Class[] potentialClasses = introspect("P2", true);
    static final Class[] p1Classes = introspect("P1", true);
    static final Class[] integratorClasses = introspect("Integrator", true);
    static final Class[] phaseClasses = introspect("Phase", true);
    static final Class[] controllerClasses = introspect("Controller", true);
    static final Class[] displayClasses = introspect("Display", true);
    static final Class[] meterClasses = introspect("Meter", true);
    static final Class[] deviceClasses = introspect("Device", true);
    static final Class[] actionClasses = introspect(etomica.Default.CLASS_DIRECTORY+"/action","",false);
    static Class[] workingClass = null;
    static Class[] validClasses = null;
    static int validCount = 0;

    public static Class[] introspect(String t, boolean b){
        return introspect(etomica.Default.CLASS_DIRECTORY, t, b);
    }
    
    public static Class[] introspect(String path, String t, boolean b){
        boolean saveGL = etomica.Default.DISPLAY_USE_OPENGL;
        etomica.Default.DISPLAY_USE_OPENGL = false;
        final String title = t;
        final boolean checkInterfaces = b;
        validCount = 0;
        workingClass = null;
        validClasses = null;
	    java.io.File dir = new java.io.File(path);
	    String[] files = dir.list(new java.io.FilenameFilter() {
	        public boolean accept(java.io.File d, String name) {
	                return name.startsWith(title)
	                && name.endsWith("class")
	                && name.indexOf("$") == -1
	                && !name.endsWith("BeanInfo.class")
	                && !name.endsWith("Editor.class")
	                && !name.endsWith("Customizer.class");}
	        });
	    workingClass = new Class[files.length];
	    for(int i=0; i<files.length; i++) {
	        int idx = files[i].lastIndexOf(".");
	        files[i] = files[i].substring(0,idx);
	        workingClass[validCount] = null;
	        String packageName = "";
	        try{
	            packageName = path.substring(etomica.Default.WORKING_DIRECTORY.length());
	            packageName = packageName.replace('/', '.');
	            workingClass[validCount] = Class.forName(packageName + "." + files[i]);
	            if (checkInterfaces){
                    checkIfGUI(workingClass[validCount]);
                }
                else { validCount++; }
	        } catch(ClassNotFoundException e) {System.out.println("Failed for "+files[i]);}
	        catch(java.lang.SecurityException se) {System.out.println("security exc");}
	        catch(java.lang.NoClassDefFoundError e) {
	            System.out.println(packageName);
	            System.out.println(files[i]);
	            System.out.println("no class def found error");
	            e.printStackTrace();
	        }
	    }// End of initialization of Classes array
	    if (validCount !=0) {
	        validClasses = new Class[validCount];
            for (int i = 0; i < validCount; i++){
                validClasses[i] = workingClass[i];
            }
        }
        etomica.Default.DISPLAY_USE_OPENGL = saveGL;
	    return validClasses;
    }// end of introspect method
    
    private static void checkIfGUI(Class wClass){
        java.lang.Class[] iFaces = wClass.getInterfaces();
//        if (checkForEtomicaInfo(wClass)) { validCount++; }
//        else {
            for (int j = 0; j < iFaces.length; j++){
                if (iFaces[j].toString().equals(String.valueOf("interface etomica.EtomicaElement"))){
                    validCount++;
                }
            }
//        }
    }// end of checkIfGUI
    
    private static boolean checkForEtomicaInfo(Class wClass){
        try {
            java.lang.reflect.Method method = wClass.getMethod("getEtomicaInfo",null);
            etomica.EtomicaInfo info = (etomica.EtomicaInfo)method.invoke(wClass,null);
  //          java.util.Vector compatElements = info.getCompatibility();
            return true;
        }
        catch(java.lang.SecurityException se){return false;}
        catch(java.lang.IllegalAccessException iae){System.out.println("illegal access exception");return false;}
        catch(java.lang.reflect.InvocationTargetException ite){System.out.println("invocation target exception");return false;}
        catch(java.lang.NoSuchMethodException nsme){return false;}
    }
}// end of IntrospectionArrays class