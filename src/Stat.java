/**
 * A static class that performs rudimentary statistics operations.
 * <p>
 *     This class and all functions within were added after I submitted my application to your institution.
 * </p>
 *
 * @author Keeler Spear
 * @version %I%, %G%
 * @since 1.0
 */
public class Stat {

    /**
     * Calculates the mean (average) of the provided data set.
     *
     * @param data The data set to be used.
     * @return The mean of the provided data set.
     */
    public static double mean(double[] data) {
        double mean = 0.0;

        for (int i = 0; i < data.length; i++) {
            mean += data[i];
        }

        return mean / data.length;
    }

    /**
     * Calculates the sample standard deviation of the provided data set.
     *
     * @param data The data set to be used.
     * @return The sample standard deviation of the provided data set.
     */
    public static double stDev(double[] data) {
        double stDev = 0.0;
        double mean = mean(data);

        for (int i = 0; i < data.length; i++) {
            stDev += Math.pow(data[i] - mean, 2);
        }

        return Math.sqrt(stDev / (data.length - 1));
    }

    /**
     * Calculates the variance of the provided data set.
     *
     * @param data The data set to be used.
     * @return The variance of the provided data set.
     * @see #stDev(double[])
     */
    public static double var(double[] data) {
        return Math.pow(stDev(data), 2);
    }

    /**
     * Standardizes the data set.
     *
     * @param data The data set to be used.
     * @return The standardized data set.
     */
    public static double[] standardize(double[] data) {
        double[] sData = new double[data.length];
        double mean = mean(data);
        double stDev = stDev(data);

        for (int i = 0; i < data.length; i++) {
            sData[i] = (data[i] - mean) / stDev;
        }

        return sData;
    }
}
