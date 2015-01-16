Tool for writing DOS from Gaussian output.

prerequisites:

	input_energy_mesh.txt : A file including the information about energy mesh.
			       Same format as is required for elses program to write DOS.
	levels.txt            : Truncated Gaussian output file.
			      (This (elses-dos) program searches the first match
			      to the string "Orbital energies and kinetic energies"
			      and begin to read eigen levels.
			      If you let it read other block (the 2nd match and after) of Gaussian output file,
			      truncate all the lines before the target block beginning with
			      "Orbital energies and kinetic energies".)

	* You must know how many levels there are in the Gaussian output file.
	This program ask you to input the number of levels from the standard input.
