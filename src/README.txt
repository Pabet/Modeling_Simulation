/////////////////////////////////////////////////////////
Important notes to Input/Output
/////////////////////////////////////////////////////////
The program can be called with 3 parameters input file(string,  input file must be in MolSim-master),
end time(double), delta time(double).If the number of parameters is wrong or any of the inputs are 
invalid(for ex. input file can not be opened), the program will take the hard coded default values
for each of the invalid files, or all of them for the wrong number of parameters.
Default input file is eingabe-sonne.txt
//////////////////////////////////////////////////////////
Output is done in 3 forms. A group of XYZ files, a group of VTU Files and a pvd file which is a collection
of the vtu files and is usefull for linking the data in paraview. They are all automatically generated in the 	 
program and no implementation to swich the generation on and off was undertaken. 
//////////////////////////////////////////////////////////