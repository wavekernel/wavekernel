all: elses_optimize_param elses_generate_molfile \
	elses_obtain_data elses_compare_data elses_convert_dzparam \
	elses_read_config

elses_optimize_param: elses_optimize_param.cc
	g++ elses_optimize_param.cc -o elses_optimize_param

elses_generate_molfile: elses_generate_molfile.cc
	g++ elses_generate_molfile.cc -o elses_generate_molfile

elses_obtain_data: elses_obtain_data.cc
	g++ elses_obtain_data.cc -o elses_obtain_data

elses_compare_data: elses_compare_data.cc
	g++ elses_compare_data.cc -o elses_compare_data

elses_convert_dzparam: elses_convert_dzparam.cc
	g++ elses_convert_dzparam.cc -o elses_convert_dzparam

elses_read_config: elses_read_config.cc
	g++ elses_read_config.cc -o elses_read_config

clean:
	rm -f elses_optimize_param elses_generate_molfile \
	elses_obtain_data elses_compare_data elses_convert_dzparam \
	elses_read_config

