NAME
    MIS_Project_2 - Performs various Compression Types to csv data sets

SYNOPSIS
    java Main [options] <infiles' root>

DESCRIPTION
    Reads the data sets from 20 line csv files and performs the requested actions such as the following: Quantization,
    Predictive Coding, Noise Visualization, Encoding, Decoding, and finding the most similar sensors by different measures.

OPTIONS
    -h Shows the help dialogue.

    -q <i> Sets the Quantizer Scheme, defaults to uniform (1), input range of [0,3]

    -b <i> Sets the number of bits per symbol to the integer i (default 8)

    -p <i> Sets the prediction scheme, defaults to no predictor (1), input range of [0,6]

    -e <i> Sets the encoder/decoder scheme, defaults to no encoding (1), input range of [1,4]

    -c <i> Sets the similarity measure, defaults to L1 (1), input range of [0,2]


DEPENDENCIES
    JVM 7.0+

EXAMPLES
    java Main -q 0 -e 0 -c 0 -p 0 Data/sampledata
        Reads from the root folder ./Data/sampledata, uses no quantization, uses no encoder, uses no 
            predictor, uses no similarity measure (calculates the three regardless), and uses 
            8 bits per symbol.

    java Main -q 3 -e 2 -c 2 -p 6 -b 12 Data/sampledata
        Reads from the root folder ./Data/sampledata, uses reverse gaussian quantization, uses run length encoding,
            uses predictor 6, uses correlation similarity (calculates the three regardless), and uses 12
            bits per symbol.

COMPILATION
    Compiles using the following command:
        javac Main.java

AUTHORS
    Steven Brown
    Kyle St. Leger-Barter
    
