package etomica.util;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public final class BuildProperties {
    public static final Properties PROPERTIES = loadProperties();

    private static Properties loadProperties() {
        try (InputStream input = BuildProperties.class.getResourceAsStream("/build.properties")) {
            Properties props = new Properties();
            if (input == null) {
                return props;
            }
            props.load(input);
            return props;
        } catch (IOException e) {
            return new Properties();
        }
    }

    public static void main(String[] args) {
        System.out.println(PROPERTIES);
    }
}
