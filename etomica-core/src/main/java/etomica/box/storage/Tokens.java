package etomica.box.storage;

public class Tokens {
    public static final Object POSITION = new Object();
    public static final Object VELOCITY = new Object();
    public static final Object ANGULAR_VELOCITY = new Object();
    public static final Object ANGULAR_MOMENTUM = new Object();

    public static final OrientationsToken ORIENTATION_FULL = () -> false;
    public static final OrientationsToken ORIENTATION = () -> true;
    public static final Object BOND_LENGTH = new Object();

    public interface OrientationsToken {
        boolean isAxisSymmetric();
    }
}
