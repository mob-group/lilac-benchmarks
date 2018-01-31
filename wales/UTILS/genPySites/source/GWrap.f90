DOUBLE PRECISION FUNCTION WRAP(value, ValRange, ZeroPhase)

        DOUBLE PRECISION :: value, valRange, zeroPhase
        INTEGER :: n

        n=FLOOR((value-ZeroPhase)/ValRange)

        WRAP=value - n * ValRange

        RETURN

        END
