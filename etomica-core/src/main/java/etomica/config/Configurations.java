package etomica.config;

public final class Configurations {
    private Configurations() {}

    public static Configuration fromResourceFile(String filename, Class<?> callingClass) {
        return new ConfigurationResourceFile(filename, callingClass);
    }
}
