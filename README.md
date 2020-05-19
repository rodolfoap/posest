# POSEST v1.2

* Original is http://users.ics.forth.gr/~lourakis/posest/. Will try to keep this updated.
* Included **levmar** (http://users.ics.forth.gr/~lourakis/levmar/) dependency.
* Included history taken from https://github.com/imbinwang/posest
* **CMakeLists.txt** slightly refactored.
* Still warnings from <math.h> _sincos()_ function (see `compiler.h`)
* See the original [README.txt](README.txt)

## Requirements

```
apt -y install liblapack-dev libblas-dev libf2c2-dev
```
## Example

```
$ head 32D.txt

   8.1878000e-02  -1.1523000e-01   2.1624960e+00   7.3585876e+02   3.5842645e+02
   4.5190000e-02   8.9986000e-02   2.1821140e+00   5.5576343e+02   5.5143250e+02
  -1.1566500e-01  -1.5556100e-01   2.0953220e+00   6.0760284e+02   1.2843721e+02
  -4.3963000e-02  -2.1449400e-01   2.1593870e+00   7.5831317e+02   1.5376729e+02
   6.0513000e-02   1.4315800e-01   2.2882420e+00   6.2063843e+02   6.6775958e+02
   1.8366000e-02  -1.4135400e-01   2.0753480e+00   6.3892877e+02   2.2926468e+02
   6.7479000e-02  -1.1182500e-01   2.1084750e+00   6.7085785e+02   3.2204453e+02
  -1.0169000e-02  -1.8204900e-01   2.1258330e+00   7.1237640e+02   1.9444815e+02
   6.2386000e-02  -1.5903000e-02   2.1521210e+00   6.2569055e+02   4.4308145e+02
   4.0083000e-02  -2.6163000e-02   2.0726150e+00   5.4278918e+02   3.7011740e+02
   ... (891 lines)

$ cat K.txt

	2980	0 	600
	0.0 	3000 	450
	0.0 	0.0 	1.0

$ ./posest_demo K.txt 32D.txt

Camera pose estimation using 891 image matches
PnP pose: -0.0534043 0.954171 0.785192 -1.54126 -0.825812 0.852895

Refinement using 1432 measurements, 6 variables
LM returned 8 in 8 iter, reason 2, error 0.684658 [initial 0.822646], 8/8 func/fjac evals

Estimated motion ([rv t]) [175 outliers, 19.64%]
-0.05213811 0.9536027 0.7843842 -1.54212 -0.8226161 0.8482345

Elapsed time: 0.01 seconds, 6.84 msecs

Camera projection center: 1.62 -0.60 0.89

Reprojection RMS and RMedS errors for input points: 110.165 0.559429
	25%, 50% and 75% quartiles: 0.308074 0.559429 1.59854

```

### Credits

All the credit to the author, Manolis Lourakis.
