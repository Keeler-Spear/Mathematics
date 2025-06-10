package Mathematics;

@FunctionalInterface
public interface NFunction<R> {
    R apply(Double... args);
}
