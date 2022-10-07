import os, sys
working_dir = os.getcwd()
with open("PIQLE.py", "rt") as input_file:
	with open("PIQLE_tmp.py", "wt") as output_file:
		for line in input_file:
			output_file.write(line.replace('change/to/your/current/directory', working_dir + '/'))
os.system('mv PIQLE_tmp.py PIQLE.py')

with open("PIQLE.py", "rt") as input_file:
        with open("PIQLE_tmp.py", "wt") as output_file:
                for line in input_file:
                        output_file.write(line.replace('configured = 0', 'configured = 1'))
os.system('mv PIQLE_tmp.py PIQLE.py')
os.system('chmod -R a+x ' + working_dir + '/apps/*')

#configure glinter
os.system('chmod -R a+x ' + working_dir + '/apps/glinter/*')
os.chdir(working_dir + '/scripts/')
with open("msa_concat.py", "rt") as input_file:
	with open("msa_concat_tmp.py", "wt") as output_file:
		for line in input_file:
			output_file.write(line.replace('change/to/your/current/directory', working_dir + '/'))
os.system('mv msa_concat_tmp.py msa_concat.py')
os.chdir(working_dir)


print('\nConfigured successfully!\n')
