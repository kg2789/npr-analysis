# npr-analysis
These are the scripts to run the updated nor-analysis which included Masaaaki's qslashed scheme

All the files are borrowed from Ryan Abbott. I have only made some slight modifications to his code to run the Masaaki-qslashed scheme.

To run these files, please use CMake.

Download these files from the repository. 

Then build a directory called build and run the following : 
```
cd build
cmake ..
make convert
make analysis
```
This should do it! 

To run the code, we use two parts of Ryan's code : 

Converting the data :
```
./convert -i ./run_1/ -o ./output 
```
Analysing the data to finally get the Z-factors :
```
./analysis ./output masaaki qslash --known-ZV 0.71408 --known-ZA 0.71408
```
